#ifdef _REENTRANT
#undef _REENTRANT
#endif

#include <flext.h>
#include "constantq.hpp"

class constantq
    : public flext_base
{
    FLEXT_HEADER_S(constantq,flext_base,setup);
    
public:
    constantq(const t_symbol *bname)
        : cq(NULL)
        , threshold(0.001),minhop(64),wndalign(1),srate(44100)
    {
        buf.Set(bname);
        
        AddInAnything();
        AddOutList();
    }
    
    ~constantq()
    {
        delete cq;
    }

protected:
    static void setup(t_classid c)
    {
        FLEXT_CADDBANG(c,0,m_bang);
        FLEXT_CADDMETHOD(c,0,m_freqs);
        FLEXT_CADDATTR_VAR(c,"freqs",mg_freqs,ms_freqs);
        FLEXT_CADDATTR_VAR1(c,"threshold",threshold);
        FLEXT_CADDATTR_VAR1(c,"minhop",minhop);
        FLEXT_CADDATTR_GET(c,"length",mg_length);
        FLEXT_CADDATTR_VAR(c,"wndalign",wndalign,ms_wndalign);
        FLEXT_CADDATTR_VAR(c,"srate",srate,ms_srate);
    }
    
    // transform buffer to list
    void m_bang()
    {
        if(cq && buf.Valid() && buf.Frames() >= cq->length()) {
            blitz::Array<t_sample,1> data(buf.Data(),shape(buf.Frames()),neverDeleteData);
            blitz::Array<t_sample,1> ref(data,Range(0,cq->length()-1));
            blitz::Array<complex<t_sample>,1> tr = (*cq)(ref);
            
            AtomListStatic<64> ret(tr.length(0));
            for(int i = 0; i < ret.Count(); ++i)
                SetFloat(ret[i],abs(tr(i)));
            ToOutList(0,ret);
        }
            
    }
    
    void m_freqs(int argc,const t_atom *argv)
    {
        freqs.resize(argc);
        for(int i = 0; i < argc; ++i)
            freqs(i) = GetAFloat(argv[i]);
            
        refreshkernel();
    }

    void ms_freqs(const AtomList &l) { m_freqs(l.Count(),l.Atoms()); }
    
    void mg_freqs(AtomList &l)
    {
        l(freqs.length(0));
        for(int i = 0; i < l.Count(); ++i)
            SetFloat(l[i],freqs(i));
    }
    
    void mg_length(int &l)
    {
        l = cq?cq->length():0;
    }
    
    void ms_wndalign(int wa)
    {
        if(wndalign != wa) {
            wndalign = wa;
            refreshkernel();
        }
    }
    
    void ms_srate(float sr)
    {
        if(srate != sr) {
            srate = sr;
            refreshkernel();
        }
    }
    
    buffer buf;
    blitz::Array<t_sample,1> freqs;
    float threshold;
    int minhop,wndalign;
    float srate;
    ConstantQ<t_sample> *cq;
    
private:

    void refreshkernel()
    {
        delete cq;
        if(freqs.length(0))
            cq = new ConstantQ<t_sample>(srate,freqs,Window::Window<Window::Hamming<t_sample> >(),threshold,wndalign,minhop);
        else
            cq = NULL;
    }

    FLEXT_CALLBACK(m_bang)
    FLEXT_CALLBACK_V(m_freqs)
    FLEXT_CALLVAR_V(mg_freqs,ms_freqs)
    FLEXT_ATTRVAR_F(threshold)
    FLEXT_ATTRVAR_F(minhop)
    FLEXT_CALLGET_I(mg_length)
    FLEXT_ATTRGET_I(wndalign)
    FLEXT_CALLSET_I(ms_wndalign)
    FLEXT_ATTRGET_F(srate)
    FLEXT_CALLSET_F(ms_srate)
};

FLEXT_NEW_1("constantq",constantq,t_symptr)
