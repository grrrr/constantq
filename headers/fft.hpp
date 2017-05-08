#ifndef __FFT_HPP
#define __FFT_HPP

#include <fftw3.h>
#include <complex>
#include <cstdlib>
#include <blitz/array.h>

template<typename T> 
class RFFT
{
public:
    RFFT(int sz,bool optimize = false)
        : size(sz),rdata(sz),cdata(sz/2+1)
    {
        if(sizeof(T) == sizeof(float))
            plan_f = fftwf_plan_dft_r2c_1d(size,(float *)rdata.data(),(fftwf_complex *)cdata.data(),optimize?FFTW_MEASURE:FFTW_ESTIMATE);
        else
            plan_d = fftw_plan_dft_r2c_1d(size,(double *)rdata.data(),(fftw_complex *)cdata.data(),optimize?FFTW_MEASURE:FFTW_ESTIMATE);
    }
    
    ~RFFT()
    {
        if(sizeof(T) == sizeof(float))
            fftwf_destroy_plan(plan_f);
        else
            fftw_destroy_plan(plan_d);
    }
    
    const blitz::Array<std::complex<T>,1> &operator()(const blitz::Array<T,1> &data)
    {
        rdata = data;
        if(sizeof(T) == sizeof(float))
            fftwf_execute(plan_f);
        else
            fftw_execute(plan_d);
        return cdata;
    }
    
    int length() const { return size; }
    
private:
    int size;
    blitz::Array<T,1> rdata;
    blitz::Array<std::complex<T>,1> cdata;
    
    union {
        fftw_plan plan_d;
        fftwf_plan plan_f;
    };
};

template<typename T> 
class RFFTH
{
public:
    RFFTH(int sz,bool optimize = false)
        : size(sz),rdata(sz),cdata(sz)
    {
        if(sizeof(T) == sizeof(float))
            plan_f = fftwf_plan_r2r_1d(size,(float *)rdata.data(),(float *)cdata.data(),FFTW_R2HC,optimize?FFTW_MEASURE:FFTW_ESTIMATE);
        else
            plan_d = fftw_plan_r2r_1d(size,(double *)rdata.data(),(double *)cdata.data(),FFTW_R2HC,optimize?FFTW_MEASURE:FFTW_ESTIMATE);
    }
    
    ~RFFTH()
    {
        if(sizeof(T) == sizeof(float))
            fftwf_destroy_plan(plan_f);
        else
            fftw_destroy_plan(plan_d);
    }
    
    const blitz::Array<T,1> &operator()(const blitz::Array<T,1> &data)
    {
        rdata = data;
        if(sizeof(T) == sizeof(float))
            fftwf_execute(plan_f);
        else
            fftw_execute(plan_d);
        return cdata;
    }
    
    int length() const { return size; }
    
private:
    int size;
    blitz::Array<T,1> rdata,cdata;
    
    union {
        fftw_plan plan_d;
        fftwf_plan plan_f;
    };
};

template<typename T> 
class CFFT
{
public:
    CFFT(int sz,bool optimize = false)
        : size(sz),idata(sz),odata(sz)
    {
        if(sizeof(T) == sizeof(float))
            plan_f = fftwf_plan_dft_1d(size,(fftwf_complex *)idata.data(),(fftwf_complex *)odata.data(),FFTW_FORWARD,optimize?FFTW_MEASURE:FFTW_ESTIMATE); 
        else
            plan_d = fftw_plan_dft_1d(size,(fftw_complex *)idata.data(),(fftw_complex *)odata.data(),FFTW_FORWARD,optimize?FFTW_MEASURE:FFTW_ESTIMATE); 
    }
    
    ~CFFT()
    {
        if(sizeof(T) == sizeof(float))
            fftwf_destroy_plan(plan_f);
        else
            fftw_destroy_plan(plan_d);
    }
    
    const blitz::Array<std::complex<T>,1> &operator()(const blitz::Array<std::complex<T>,1> &data)
    {
        idata = data;
        if(sizeof(T) == sizeof(float))
            fftwf_execute(plan_f);
        else
            fftw_execute(plan_d);
        return odata;
    }
    
    int length() const { return size; }
    
private:
    int size;
    blitz::Array<std::complex<T>,1> idata;
    blitz::Array<std::complex<T>,1> odata;
    
    union {
        fftw_plan plan_d;
        fftwf_plan plan_f;
    };
};


template<typename T> class FFT;

#endif
