/*
constantq~ - constant-q analysis

Copyright (c) 2008-2009 Thomas Grill (gr@grrrr.org)
For information on usage and redistribution, and for a DISCLAIMER OF ALL
WARRANTIES, see the file, "license.txt," in this distribution.  

$LastChangedRevision$
$LastChangedDate$
$LastChangedBy$
*/

#define FLEXT_ATTRIBUTES 1

#include <flext.h>

#if FLEXT_VERSION < 502
#error flext version too old.
#endif

#ifdef _REENTRANT
#undef _REENTRANT
#endif

#include "constantq.hpp"

namespace { // anonymous namespace


template<typename T>
static T hz2bark(T f)
{
    // return N.arctan(f*0.00076)*13.+N.arctan(N.square(f/7500.))*3.5
    // approximation (invertible!)
    T f1 = f/600.;
    return log(f1+sqrt(f1*f1+1))*6.;
}

template<typename T> 
static T bark2hz(T b) { return (exp(b/6.)-exp(b/-6.))*300.; }

template<typename T>
static T hz2mel(T f) { return hz2bark(f)*100.; }

template<typename T>
static T mel2hz(T m) { return bark2hz(m/100.); }


class constantq
    : public flext_dsp
{
    FLEXT_HEADER_S(constantq,flext_dsp,setup);
    
public:
    constantq()
        : cq(NULL),offset(0)
        , threshold(0.001),wndalign(1),minhop(0)
        , srate(0),rate(-1)
        , buf(NULL),bufupd(true)
    {
        AddInSignal();
        AddOutList();

	    FLEXT_ADDTIMER(tmr,m_bang);
    }
    
    ~constantq()
    {
        delete cq;
        delete buf;
    }

protected:

    static const t_symbol *sym_reset;
    static const t_symbol *sym_triangular,*sym_hamming,*sym_hann,*sym_cosine,*sym_bartlett,*sym_blackman;

    static void setup(t_classid c)
    {
        sym_reset = MakeSymbol("reset");

        sym_triangular = MakeSymbol("triangular");
        sym_hamming = MakeSymbol("hamming");
        sym_hann = MakeSymbol("hann");
        sym_cosine = MakeSymbol("cosine");
        sym_bartlett = MakeSymbol("bartlett");
        sym_blackman = MakeSymbol("blackman");

        FLEXT_CADDMETHOD(c,0,m_bang);
        FLEXT_CADDMETHOD(c,0,m_frqs);
        FLEXT_CADDMETHOD_(c,0,sym_reset,m_reset); // recalculate kernel
        FLEXT_CADDATTR_VAR(c,"buffer",mg_buffer,ms_buffer);
        FLEXT_CADDATTR_VAR1(c,"bufupd",bufupd);
        FLEXT_CADDATTR_VAR(c,"frqs",mg_frqs,ms_frqs);
        FLEXT_CADDMETHOD_(c,0,"logfrqs",m_logfrqs);
        FLEXT_CADDMETHOD_(c,0,"melfrqs",m_melfrqs);
        FLEXT_CADDATTR_VAR(c,"q",mg_qfactors,ms_qfactors);
        FLEXT_CADDATTR_VAR1(c,"threshold",threshold);
        FLEXT_CADDATTR_GET(c,"granularity",mg_granularity);
        FLEXT_CADDATTR_GET(c,"length",mg_length);
        FLEXT_CADDATTR_VAR(c,"window",wndname,ms_window);
        FLEXT_CADDATTR_VAR1(c,"wndalign",wndalign);
        FLEXT_CADDATTR_VAR1(c,"minhop",minhop);
        FLEXT_CADDATTR_VAR(c,"rate",rate,ms_rate);       
    }

    virtual bool CbDsp()
    {
        if(srate != Samplerate()) {
            srate = Samplerate();
            refreshkernel();
        }
        else
            resetbuffer();
        return true;
    }
    
    virtual void CbSignal()
    {
        if(!cq) return;
        
        int const n = Blocksize();
        Array<t_sample,1> vec(InSig(0),n,neverDeleteData);
        size_t const len = cq->length();
        data(Range(offset,offset+n-1)) = vec;
        if(offset < len-n) {
            data(Range(offset+len,offset+len+n-1)) = vec;
            offset += n;
        }
        else {
            // wrap around
            int r = offset-len; // < 0
            data(Range(offset+len,toEnd)) = vec(Range(0,-r-1));
            data(Range(0,n-1+r)) = vec(Range(-r,n-1));
            offset = r+n;
        }        

        if(rate < 0)
            tmr.Now(); // trigger analysis
    }
    
    void m_bang(void * = NULL)
    {
        if(!cq) return;
    
        Array<t_sample,1> ref(data,Range(offset,offset+cq->length()-1));
        
        // analyze
        Array<t_sample,1> cqspec(abs((*cq)(ref)));
        
        if(buf && buf->Ok() && buf->Valid()) {
            {
                buffer::Locker buflock(*buf);
                buf->Update();
                
                Array<t_sample,2> barr(buf->Data(),TinyVector<int,2>(buf->Frames(),buf->Channels()),neverDeleteData);
                Range rng(0,min(cqspec.size(),buf->Frames())-1);
                barr(rng,0) = cqspec(rng);
                
                buf->Dirty(bufupd);
            }
            // unlock before bang
            ToOutBang(0);
        }
        else {
            AtomListStatic<256> ret(cqspec.size());
            for(int i = 0; i < ret.Count(); ++i)
                SetFloat(ret[i],cqspec(i));
            ToOutList(0,ret);
        }
    }
    
    void m_reset() { refreshkernel(); }
    
    void m_frqs(int argc,const t_atom *argv)
    {
        freqs.resize(argc);
        for(int i = 0; i < argc; ++i)
            freqs(i) = GetAFloat(argv[i]);
    }

    void ms_frqs(const AtomList &l) { m_frqs(l.Count(),l.Atoms()); }
    
    void mg_frqs(AtomList &l) const
    {
        l(freqs.size());
        for(int i = 0; i < l.Count(); ++i)
            SetFloat(l[i],freqs(i));
    }
    
    void m_logfrqs(int argc,const t_atom *argv)
    {
        if(argc < 2) return;
        float fmin = GetAFloat(argv[0]);
        float fmax = GetAFloat(argv[1]);
        int odiv = argc >= 3?GetAInt(argv[2]):12;
        if(fmin >= fmax || odiv <= 0) return;

        float lmin = log(fmin)/M_LN2;
        float lmax = log(fmax)/M_LN2;
        int bnds = int(floor((lmax-lmin)*odiv)+0.5);
        float pow2n = pow(2,1./odiv);

        freqs.resize(bnds+1);
        qfactors.resize(bnds+1);
        for(int b = 0; b <= bnds; ++b)
            freqs(b) = pow(pow2n,b)*fmin;
        qfactors = sqrt(pow2n)/(pow2n-1.);
    }
    
    void m_melfrqs(float res)
    {
        if(!res) res = 1;
    
        int bnds = int(25./res+0.5);
        freqs.resize(bnds);
        qfactors.resize(bnds);

        for(int b = 0; b < bnds; ++b) {
            float bark = res*(b+0.5);
            freqs(b) = bark2hz(bark);
            float odiv = tanh(bark/6.)*(log(64.)/res);
            float pow2n = pow(2.,1./odiv);
            qfactors(b) = sqrt(pow2n)/(pow2n-1.);
        }
    }
    
    void m_qfactors(int argc,const t_atom *argv)
    {
        qfactors.resize(argc);
        for(int i = 0; i < argc; ++i)
            qfactors(i) = GetAFloat(argv[i]);
    }

    void ms_qfactors(const AtomList &l) { m_qfactors(l.Count(),l.Atoms()); }
    
    void mg_qfactors(AtomList &l) const
    {
        l(qfactors.size());
        for(int i = 0; i < l.Count(); ++i)
            SetFloat(l[i],qfactors(i));
    }
    
    void mg_length(int &l) const { l = cq?cq->length():0; }

    void mg_granularity(int &h) const { h = cq?cq->hopsize():0; }
    
    void ms_rate(float r)
    {
        if(rate != r) {
            rate = r;
            if(rate > 0) 
                tmr.Periodic(rate*0.001);
            else
                tmr.Reset();
        }
    }
    
    void mg_buffer(AtomList &l) const
    {
        if(buf) {
            l(1);
            SetSymbol(l[0],buf->Symbol());
        }
    }
    
    void ms_buffer(const AtomList &l)
    {
        if(l.Count() == 1 && IsSymbol(l[0])) {
            const t_symbol *bn = GetSymbol(l[0]);
            if(!buf || buf->Symbol() != bn) {
                delete buf;
                buf = new buffer(bn);
            }
        }
        else {
            delete buf;
            buf = NULL;
        }
    }
    
    void ms_window(const AtomList &l)
    {
        delete window;
        if(!l.Count())
            wndname(0),window = NULL; // boxcar window
        else if(IsSymbol(l[0])) {
            wndname = l;
            const t_symbol *wn = GetSymbol(l[0]);
            if(wn == sym_triangular)
                window = new Window::Window<Window::Triangular<t_sample> >();
            else if(wn == sym_hamming)
                window = new Window::Window<Window::Hamming<t_sample> >();
            else if(wn == sym_hann)
                window = new Window::Window<Window::Hann<t_sample> >();
            else if(wn == sym_cosine)
                window = new Window::Window<Window::Cosine<t_sample> >();
            else if(wn == sym_bartlett)
                window = new Window::Window<Window::Bartlett<t_sample> >();
            else if(wn == sym_blackman)
                window = new Window::Window<Window::Blackman<t_sample> >();
            else
                window = NULL;
        }
    }
    
    
    Array<t_sample,1> data;
    int offset;
    Array<t_sample,1> freqs,qfactors;
    float threshold;
    int wndalign,minhop;
    float srate;
    ConstantQ<t_sample,false> *cq;
    Timer tmr;
    double rate;
    buffer *buf;
    bool bufupd;
    AtomList wndname;
    Window::Base<t_sample> *window;
    
private:

    void resetbuffer()
    {
        if(!cq) return;
        
        int l = cq->length()*2;
        if(l != data.size()) {
            data.resize(l);
            data = 0; 
            offset = 0;
        }
    }

    void refreshkernel()
    {
        int len = cq?cq->length():0;
        
        delete cq;
        
        const int nfreqs = freqs.length(0);
        if(nfreqs) {
            Array<t_sample,1> qf(nfreqs);

            if(qfactors.length(0) && qfactors(0)) {
                if(nfreqs == qfactors.length(0))
                    qf = qfactors;
                else
                    qf = qfactors(0);            
            }
            else {
                // determination of q from frequencies 
                // \TODO re-check this!!
                qf(0) = 1./(freqs(1)/freqs(0)-1);
                int i;
                for(i = 1; i < nfreqs-1; ++i)
                    qf(i) = 1./(sqrt(freqs(i+1)/freqs(i-1))-1);
                qf(i) = 1./(freqs(i)/freqs(i-1)-1.);
            }
            
            Window::Window<Window::Boxcar<t_sample> > boxcar;
            cq = new ConstantQ<t_sample,false>(srate,freqs,qf,window?*window:boxcar,threshold,wndalign,max(minhop,Blocksize()));
        }
        else
            cq = NULL;
            
        if(cq && cq->length() != len)
            resetbuffer();
            
        ToOutAnything(GetOutAttr(),sym_reset,0,NULL);
    }

    FLEXT_CALLBACK(m_bang)
    FLEXT_CALLBACK(m_reset)
    FLEXT_CALLBACK_V(m_frqs)
    FLEXT_CALLVAR_V(mg_frqs,ms_frqs)
    FLEXT_CALLBACK_V(m_logfrqs)
    FLEXT_CALLBACK_F(m_melfrqs)
    FLEXT_CALLBACK_V(m_qfactors)
    FLEXT_CALLVAR_V(mg_qfactors,ms_qfactors)
    FLEXT_CALLVAR_V(mg_buffer,ms_buffer)
    FLEXT_ATTRVAR_B(bufupd)
    FLEXT_ATTRVAR_F(threshold)
    FLEXT_CALLGET_I(mg_granularity)
    FLEXT_CALLSET_F(ms_rate)
    FLEXT_ATTRGET_F(rate)
    FLEXT_CALLGET_I(mg_length)
    FLEXT_ATTRVAR_I(wndalign)
    FLEXT_CALLSET_V(ms_window)
    FLEXT_ATTRGET_V(wndname)
    FLEXT_ATTRVAR_I(minhop)
    FLEXT_CALLBACK_T(m_bang)
};

const t_symbol *constantq::sym_reset;
const t_symbol *constantq::sym_triangular,*constantq::sym_hamming,*constantq::sym_hann,*constantq::sym_cosine,*constantq::sym_bartlett,*constantq::sym_blackman;

FLEXT_NEW_DSP("constantq~",constantq)

} // namespace
