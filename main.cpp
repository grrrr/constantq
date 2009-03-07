#define FLEXT_ATTRIBUTES 1

#include <flext.h>

#ifdef _REENTRANT
#undef _REENTRANT
#endif

#include "constantq.hpp"


class constantq
    : public flext_dsp
{
    FLEXT_HEADER_S(constantq,flext_dsp,setup);
    
public:
    constantq()
        : cq(NULL),offset(0)
        , threshold(0.001),wndalign(1)
        , srate(0),rate(-1)
    {
        AddInSignal();
        AddOutList();

	    FLEXT_ADDTIMER(tmr,m_bang);
    }
    
    ~constantq()
    {
        delete cq;
    }

protected:
    static void setup(t_classid c)
    {
        FLEXT_CADDMETHOD(c,0,m_bang);
        FLEXT_CADDMETHOD(c,0,m_freqs);
        FLEXT_CADDMETHOD_(c,0,"reset",m_reset); // recalculate kernel
        FLEXT_CADDATTR_VAR(c,"freqs",mg_freqs,ms_freqs);
        FLEXT_CADDATTR_VAR(c,"q",mg_qfactors,ms_qfactors);
        FLEXT_CADDATTR_VAR1(c,"threshold",threshold);
        FLEXT_CADDATTR_GET(c,"hopsize",mg_hopsize);
        FLEXT_CADDATTR_GET(c,"length",mg_length);
        FLEXT_CADDATTR_VAR1(c,"wndalign",wndalign);
        FLEXT_CADDATTR_VAR(c,"rate",rate,ms_rate);
    }

#if 0    
    virtual bool CbFinalize()
    {
        if(flext_base::CbFinalize()) {
            refreshkernel();
            return true;
        }
        else
            return false;
    }
#endif
    
    virtual bool CbDsp()
    {
#ifdef FLEXT_DEBUG
        post("DSP");
#endif
    
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
        if(offset < len-n)
            data(Range(offset+len,offset+len+n-1)) = vec;
        else {
            int r = offset-len; // < 0
            data(Range(offset+len,toEnd)) = vec(Range(0,-r-1));
            data(Range(0,n-1+r)) = vec(Range(-r,n-1));
            offset = r;
        }        
        offset += n;

        if(rate < 0)
            tmr.Now();
    }
    
    void m_bang(void * = NULL)
    {
        if(!cq) return;
    
        Array<t_sample,1> ref(data,Range(offset,offset+cq->length()-1));
        Array<complex<t_sample>,1> tr = (*cq)(ref);
        
        AtomListStatic<256> ret(tr.length(0));
        for(int i = 0; i < ret.Count(); ++i)
            SetFloat(ret[i],abs(tr(i)));
        ToOutList(0,ret);
    }
    
    void m_reset() { refreshkernel(); }
    
    void m_freqs(int argc,const t_atom *argv)
    {
        freqs.resize(argc);
        for(int i = 0; i < argc; ++i)
            freqs(i) = GetAFloat(argv[i]);
    }

    void ms_freqs(const AtomList &l) { m_freqs(l.Count(),l.Atoms()); }
    
    void m_qfactors(int argc,const t_atom *argv)
    {
        qfactors.resize(argc);
        for(int i = 0; i < argc; ++i)
            qfactors(i) = GetAFloat(argv[i]);
    }

    void ms_qfactors(const AtomList &l) { m_qfactors(l.Count(),l.Atoms()); }
    
    void mg_freqs(AtomList &l) const
    {
        l(freqs.length(0));
        for(int i = 0; i < l.Count(); ++i)
            SetFloat(l[i],freqs(i));
    }
    
    void mg_qfactors(AtomList &l) const
    {
        l(qfactors.length(0));
        for(int i = 0; i < l.Count(); ++i)
            SetFloat(l[i],qfactors(i));
    }
    
    void mg_length(int &l) const { l = cq?cq->length():0; }

    void mg_hopsize(int &h) const { h = cq?cq->hopsize():0; }
    
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
    
    Array<t_sample,1> data;
    int offset;
    Array<t_sample,1> freqs,qfactors;
    float threshold;
    int wndalign;
    float srate;
    ConstantQ<t_sample,false> *cq;
    Timer tmr;
    double rate;
    
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
#ifdef FLEXT_DEBUG
        post("Refreshing kernel");
#endif

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
            
/*            
            cerr << "F ";
            for(int i = 0; i < freqs.length(0); ++i)
                cerr << freqs(i) << " ";
            cerr << endl;
            cerr << "Q ";
            for(int i = 0; i < qf.length(0); ++i)
                cerr << qf(i) << " ";
            cerr << endl;
*/            
            cq = new ConstantQ<t_sample,false>(srate,freqs,qf,Window::Window<Window::Hamming<t_sample> >(),threshold,wndalign,Blocksize());
        }
        else
            cq = NULL;
            
        if(cq && cq->length() != len)
            resetbuffer();
    }

    FLEXT_CALLBACK(m_bang)
    FLEXT_CALLBACK(m_reset)
    FLEXT_CALLBACK_V(m_freqs)
    FLEXT_CALLVAR_V(mg_freqs,ms_freqs)
    FLEXT_CALLBACK_V(m_qfactors)
    FLEXT_CALLVAR_V(mg_qfactors,ms_qfactors)
    FLEXT_ATTRVAR_F(threshold)
    FLEXT_CALLGET_I(mg_hopsize)
    FLEXT_CALLSET_F(ms_rate)
    FLEXT_ATTRGET_F(rate)
    FLEXT_CALLGET_I(mg_length)
    FLEXT_ATTRVAR_I(wndalign)
    FLEXT_CALLBACK_T(m_bang)
};

FLEXT_NEW_DSP("constantq~",constantq)
