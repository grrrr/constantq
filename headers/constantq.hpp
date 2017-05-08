#ifndef __CONSTANTQ_HPP
#define __CONSTANTQ_HPP

#include <complex>
#include <ostream>
#include <vector>
#include <cstdlib>
#include <blitz/array.h>
#include "fft.hpp"
#include "window.hpp"

using namespace std;
using namespace blitz;


template<typename T,bool banded> class Kernel;

// template specialization for individual sparse nodes
template<typename T> 
class Kernel<T,false>
{
public:
    void init(Array<complex<T>,1> const &ckr,T thr)
    {
        for(int c = 0; c < ckr.size(); ++c)
            if(abs(ckr(c)) >= thr)
                sp.push_back(SparseNode(c,ckr(c)));
    }
    
    complex<T> operator()(Array<complex<T>,1> const &,Array<complex<T>,1> const &hft) const
    {
        complex<T> s(0);
        for(typename vector<SparseNode>::const_iterator it1 = sp.begin(); it1 != sp.end(); ++it1)
            s += hft(it1->p)*it1->v;
        return s;
    }

protected:
    struct SparseNode 
    {
        SparseNode(int _p,complex<T> _v): p(_p),v(_v) {}
        int p; complex<T> v; 
    };
    
    vector<SparseNode> sp;
};

// template specialization for banded sparsity
template<typename T> 
class Kernel<T,true>
{
public:
    void init(Array<complex<T>,1> const &ckr,T thr)
    {
        int c;
        for(c = 0; c < ckr.size() && abs(ckr(c)) < thr; ++c);
        if(c < ckr.size()) {
            lo = c;
            for(c = ckr.size()-1; c >= 0 && abs(ckr(c)) < thr; --c);
            hi = c;
        }
        else
            lo = -1;
    }

    complex<T> operator()(Array<complex<T>,1> const &ckernel,Array<complex<T>,1> const &hft) const
    {
        if(lo >= 0) {
            Range r(lo,hi);
            return sum(ckernel(r)*hft(r));
        }
        else
            return 0;
    }

protected:
    int lo,hi;
};



template<typename T,bool banded = false>
class ConstantQ
{
public:
    ConstantQ(float sr,Array<T,1> const &fr,Array<T,1> const &q,Window::Base<T> const &window,T thr=0.0001,int wndalign=0,int minhop=0)
        : srate(sr),freqs(fr)
        , threshold(thr)
        , fft(NULL),sparse(NULL)
    {    
        int nfreqs = freqs.size();
        Array<T,1> digfreq(freqs*T(2.*M_PI/srate));

        int i;
        Array<int,1> wndsz(q/freqs*srate+0.5);
        wndsz = max(wndsz,minhop*2);
        
//        cerr << "Windows " << wndsz << endl;
        
        wndszmax = max(wndsz);
        wndszmin = min(wndsz);
        datasize = nextpow2(wndszmax);
        
        if(fft && fft->length() != datasize) { 
            delete fft; fft = NULL;
        }
        if(!fft)
            fft = new RFFT<T>(datasize,true);

        ckernel[0].resize(nfreqs,datasize/2);
        ckernel[1].resize(nfreqs,datasize/2);
        
//        cerr << "ckernel " << ckernel[0].shape() << endl;

        for(i = 0; i < nfreqs; ++i) {
            int sz = wndsz(i);
            
            firstIndex ix;
            Array<T,1> fr(sz);
            fr = (ix-sz/2.)*digfreq(i);
            
//            cerr << i << " WAVE " << fr << endl;
//            cerr << i << " -> " << wndszmax << "," << sz << endl;

            int numz;
            if(!wndalign)
                numz = int(ceil((wndszmax-sz)/2.));  // center-align
            else if(wndalign < 0)
                numz = 0;   // front-align
            else
                numz = wndszmax-sz;  // rear-align
                
            Array<complex<T>,1> kernel(datasize);
            kernel = 0;
            // fill with windowed complex waves
            Array<complex<T>,1> focus(kernel,Range(numz,numz+sz-1));
//            cerr << "FOCUS" << focus.shape() << '/' << numz << ',' << (numz+sz) << endl;
            
            Array<T,1> wnd(window(sz)/sz);
            
            real(focus) = cos(fr);
            imag(focus) = sin(fr);
            
//            cerr << "Focus" << i << " ... " << focus << endl;
            
            focus *= wnd;
            
            CFFT<T> fft(datasize,true);
            Array<complex<T>,1> ft(fft(kernel)/T(datasize)); 
            
//            cerr << "-------------------------- " << endl << i << " " << ft << endl;
            
            ckernel[0](i,Range::all()) = conj(ft(Range(fromStart,datasize/2-1)));
            ckernel[1](i,Range::all()) = ft(Range(datasize/2,toEnd));
        }
        
//        cerr << "KERNEL1" << endl << ckernel[0] << endl;
//        cerr << "KERNEL2" << endl << ckernel[1] << endl;

//        cerr << "Ready" << endl;
        
        if(threshold) {
            Array<T,2> ka1(abs(ckernel[0]));
            Array<T,2> ka2(abs(ckernel[1]));
            T const thr = max(max(ka1),max(ka2))*threshold; // compute relative threshold

//            cerr << "Threshold = " << thr << endl;
            
            if(sparse) delete[] sparse;
            sparse = new SparseC[nfreqs];
            
            for(int r = 0; r < nfreqs; ++r)
                for(int k = 0; k < 2; ++k)
                    sparse[r][k].init(ckernel[k](r,Range::all()),thr);
        }

    }
    
    int hopsize() const { return wndszmin/2; }

    int length() const { return datasize; }

    int windowmin() const { return wndszmin; }
    int windowmax() const { return wndszmax; }
    
    Array<complex<T>,1> operator()(const Array<T,1> &data)
    {
        // transform
        Array<complex<T>,1> ft((*fft)(data));
        // first and second half
        Array<complex<T>,1> hft1(ft,Range(0,datasize/2-1));
        Array<complex<T>,1> hft2(ft,Range(datasize/2,1,-1));
        
//        cerr << endl << ft.shape() << " -> " << hft1.shape() << " " << hft2.shape() << endl;
//        cerr << endl << "F1" << endl << hft1 << endl << "F2" << endl << hft2 << endl;
        
        Array<complex<T>,1> ret(freqs.size());
        
        if(threshold) {
            // sparse multiplication of kernel and spectrum

            for(int f = 0; f < freqs.size(); ++f) {
                Array<complex<T>,1> hkernel1 = ckernel[0](f,Range::all());
                Array<complex<T>,1> hkernel2 = ckernel[1](f,Range::all());
                ret(f) = sparse[f][0](hkernel1,hft1)+conj(sparse[f][1](hkernel2,hft2));
            }
        }
        else {
            // dot product of kernel and spectrum
            for(int f = 0; f < freqs.size(); ++f) {
                Array<complex<T>,1> hkernel1 = ckernel[0](f,Range::all());
                Array<complex<T>,1> hkernel2 = ckernel[1](f,Range::all());
                ret(f) = sum(hkernel1*hft1)+conj(sum(hkernel2*hft2));
            }
        }
        return ret;
    }
    
    friend ostream &operator <<(ostream &s,const ConstantQ &cq)
    {
        s << "ConstantQ< [";
        for(typename ConstantQ::Sparse::const_iterator it1 = cq.sparse1.begin(); it1 != cq.sparse1.end(); ++it1)
            s << "(" << it1->r << "," << it1->c << ":" << it1->v << ") ";
        s << "] [";
        for(typename ConstantQ::Sparse::const_iterator it2 = cq.sparse2.begin(); it2 != cq.sparse2.end(); ++it2)
            s << "(" << it2->r << "," << it2->c << ":" << it2->v << ") ";
        s << ">";
        return s;
    }
    
protected:

    static int nextpow2(int n)
    {
	   --n;
	   n |= n>>1;
	   n |= n>>2;
	   n |= n>>4;
	   n |= n>>8;
	   n |= n>>16;
	   return n+1;
	}

    float srate;
    Array<T,1> freqs;
    Array<complex<T>,2> ckernel[2];
    int wndszmin,wndszmax,datasize;
    T threshold;

    RFFT<T> *fft;

    typedef Kernel<T,banded> Sparse;
    typedef Sparse SparseC[2];   
    SparseC *sparse;
};

#endif
