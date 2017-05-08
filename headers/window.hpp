#ifndef __WINDOW_HPP
#define __WINDOW_HPP

#include <cstdlib>
#include <cmath>
#include <blitz/array.h>

using namespace blitz;

namespace Window {
// http://en.wikipedia.org/wiki/Window_function

template<typename T> struct Function { typedef T Type; }; 

template<typename T> struct Boxcar: public Function<T> { static T calc(int i,int n) { return 1; } };
template<typename T> struct Triangular: public Function<T> { static T calc(int i,int n) { return 1.-abs(2*i-(T(n)-1)/2)/T(n); } }; // FALSCH
//    template<typename T> struct Hamming: public Function<T> { static T calc(int i,int n) { return 0.53836-0.46164*cos(2*M_PI*i/(n-1)); } };
template<typename T> struct Hamming: public Function<T> { static T calc(int i,int n) { return 0.54-0.46*cos(2.*M_PI*T(i)/T(n-1)); } };
template<typename T> struct Hann: public Function<T> { static T calc(int i,int n) { return 0.5*(1.-cos(2*M_PI*T(i)/T(n-1))); } };
template<typename T> struct Cosine: public Function<T> { static T calc(int i,int n) { return sin(M_PI*T(i)/T(n-1)); } };
template<typename T> struct Bartlett: public Function<T> { static T calc(int i,int n) { return 1.-abs(2.*T(i)/T(n-1)-1); } };

template<typename T> 
struct Blackman: public Function<T>
{ 
    static T calc(int i,int n) 
    {
        static const float a=0.16;
        return (1-a)/2.-0.5*cos(2*M_PI*i/(n-1))+a/2*cos(4*M_PI*i/(n-1));
    } 
};


template<typename T>
class Base
{
public:
    virtual ~Base() {}
    virtual Array<T,1> operator()(int n) const = 0;
};

template<typename F>
class Window
    : public Base<typename F::Type>
{
public:
    virtual Array<typename F::Type,1> operator()(int n) const
    {
        Array<typename F::Type,1> ret(n);
        for(int i = 0; i < n; ++i)
            ret(i) = F::calc(i,n);
        return ret;
    }
};

} // namespace

#endif