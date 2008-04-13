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
        , threshold(0.01),minhop(64)
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
    }
    
    // transform buffer to list
    void m_bang()
    {
        if(cq && buf.Valid() && buf.Frames() >= cq->length()) {
            blitz::Array<t_sample,1> data(buf.Data(),shape(buf.Frames()),neverDeleteData);
            blitz::Array<t_sample,1> ref(data,Range(0,cq->length()-1));
            blitz::Array<complex<t_sample>,1> tr = (*cq)(ref);
            AtomList ret(tr.length(0));
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
            
        delete cq;
        cq = new ConstantQ<t_sample>(44100,freqs,Window::Window<Window::Hamming<t_sample> >(),threshold,0,minhop);
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
    
    buffer buf;
    blitz::Array<t_sample,1> freqs;
    float threshold;
    int minhop;
    ConstantQ<t_sample> *cq;
    
private:
    FLEXT_CALLBACK(m_bang)
    FLEXT_CALLBACK_V(m_freqs)
    FLEXT_CALLVAR_V(mg_freqs,ms_freqs)
    FLEXT_ATTRVAR_F(threshold)
    FLEXT_ATTRVAR_F(minhop)
    FLEXT_CALLGET_I(mg_length)
};

FLEXT_NEW_1("constantq",constantq,t_symptr)
