//
// Created by alejandro on 10/07/19.
//
#ifndef MONTECARLOENGINE_OPTIONGEN_H
#define MONTECARLOENGINE_OPTIONGEN_H

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <memory>
#include <numeric>
#include <adept.h>
#include <map>
#include <algorithm>

using namespace std;
enum optionType { call, put };
enum asianType {max_, min_, avg_};

template <typename T>
T normalCDF( T x){
    T result = erfc(-x / sqrt(2)) * 0.5;
    return result;
}

template <typename T>
class OptionBS{
public:
    optionType tipo;
    T interesAnual;
    T strike;
    T spot;
    T sigma;
    T tau;

    T forward;
    T factorDescuento;

    OptionBS(optionType tipo,T interesAnual, T strike, T spot, T sigma, T tau):
            tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    }
    T price(){
        T dplus{};
        T dminus{};
        T sigma_sqr = sigma*sigma;
        T inverseSigma = 1 / (sigma * sqrt(tau));
        T logFK = log(forward/strike);
        T medioSigma = 0.5 * sigma_sqr*tau;

        dplus = inverseSigma*(logFK + medioSigma);
        dminus = inverseSigma*(logFK - medioSigma);

        T distribucionPlus{};
        T distribucionMinus{};
        T negDistribucionPlus{};
        T negDistribucionMinus{};

        T result{};
        switch (tipo){
            case call:
                 distribucionPlus = normalCDF(dplus);
                 distribucionMinus = normalCDF(dminus);
                result = factorDescuento * ((distribucionPlus * forward) - (distribucionMinus * strike));
                break;
            case put:
                dplus = -dplus; //correcion para que funcion con adept
                dminus = -dminus;
                 negDistribucionMinus = normalCDF(dminus);
                 negDistribucionPlus = normalCDF(dplus);
                result = factorDescuento * ((negDistribucionMinus * strike) - (negDistribucionPlus * forward));
                break;
            default:
                cout<<"Error"<<endl;
                result = -11111111;
                break;
        }
        return result;
    };
};

template <typename T>
class OptionGen {
public:
    virtual OptionGen<T> *getOptionGen( int )
    {
        return 0;
    }

    int signo;

    virtual map<T,T> evaluate(vector<T> values, double criteria) = 0;
    virtual double getExpiry() = 0;
};

template <typename T>
class Decorator: public OptionGen<T>{
public:
    explicit Decorator (OptionGen<T> * optionAbstract){
        innerOption = optionAbstract;
    }
    OptionGen<T>* innerOption;
    map<T,T> evaluate(vector<T> values, double criteria){
        innerOption->evaluate(values,criteria);
    };
    double getExpiry(){
        innerOption->getExpiry();
    }
};

template <typename T>
class Composite: public OptionGen<T>{
public:
    vector<std::pair<int, OptionGen<T> *>> vectorInst;

    void add (OptionGen<T> *elem){
        vectorInst.push_back(make_pair(elem->signo,elem));
    }


    void remove( const unsigned int index )
    {
        OptionGen<T> *child = vectorInst[ index ];
        vectorInst.erase( vectorInst.begin() + index );
        delete child;
    }

    unsigned long compSize(){
        return vectorInst.size();
    }


    OptionGen<T> *getOptionGen (int n){
        auto elem = get<1>(vectorInst[n]);
        return elem;
    }

    //mapa <U, T>
    map<T,T> evaluate(vector<T> values, double criteria){
        map<T,T> result{};
        map<T,T> result2{};

        for(auto element : vectorInst){

            result2 = get<1>(element)->evaluate(values, criteria);
            for (auto it=result2.begin(); it!=result2.end(); ++it){
                    auto f = result.find(it->first);
                    if(f != result.end()){
                        //multiplicar por signo
                        T suma{};
                        suma = f->second + it->second;
                        result.erase(f);
                        result.insert (std::make_pair(it->first,suma));
                    }else {
                        result.insert(std::make_pair(it->first,it->second));
                    }
            }
        }
        return result;
    }

    double  getExpiry(){
        double  maxExpiry = 0.0;
        for (const auto element : vectorInst){
            maxExpiry = max(maxExpiry,get<1>(element)->getExpiry());
        }
        return maxExpiry;
    }

};

template<typename T>
class Greeks{
public:

    optionType tipo;
    T interesAnual;
    T strike;
    T spot;
    T sigma;
    double tau;

    T forward;
    T factorDescuento;

    Greeks();

    Greeks(optionType tipo,T interesAnual, T strike, T spot, T sigma, double tau):
            tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    };

    T vega(){
        T vega{};
        T dplus{};
        T logFK = log(forward/strike);
        T inverseSigma = 1 / (sigma * sqrt(tau));
        T medioSigma = 0.5 * sigma*sigma*tau;
        dplus = inverseSigma*(logFK + medioSigma);

        T spotSqrTau = spot*sqrt(tau);
        T distribucionPrima = (1/(sqrt(2*M_PI)))*exp(-(dplus*dplus)/2);
        vega = spotSqrTau*distribucionPrima;
        return vega;
    };
    T theta(){
        T theta{};
        T dplus{};
        T dminus{};
        T logFK = log(forward/strike);
        T inverseSigma = 1 / (sigma * sqrt(tau));
        T medioSigma = 0.5 * sigma*sigma*tau;

        dplus = inverseSigma*(logFK + medioSigma);
        T distribucionPrima = (1/(sqrt(2*M_PI)))*exp(-(dplus*dplus)/2);

        dminus = inverseSigma*(logFK - medioSigma);

        boost::math::normal normalDistribution;
        T distribucionMinus{};
        switch (tipo){
            case call:
                distribucionMinus = normalCDF(dminus);
                theta = -(spot*distribucionPrima*sigma)/(2*sqrt(tau)) - (interesAnual*strike*exp(-interesAnual*tau)*(distribucionMinus));
                break;
            case put:
                dminus = -dminus; ///correcion para adouble no permite pasar parametro negativo a la funcion
                distribucionMinus = normalCDF(dminus);
                theta = -((spot*distribucionPrima*sigma)/(2*sqrt(tau))) + (interesAnual*strike*exp(-interesAnual*tau)*(distribucionMinus));
                break;
            default:
                // cout<<"Error"<<endl;
                theta = -11111111;
                break;
        }
        return theta;
    };
    T delta(){
        T delta{};
        T dplus{};
        T logFK = log(forward/strike);
        T inverseSigma = 1 / (sigma * sqrt(tau));
        T medioSigma = 0.5 * sigma*sigma*tau;
        dplus = inverseSigma*(logFK + medioSigma);
        boost::math::normal normalDistribution;
        delta = normalCDF(dplus);
        switch (tipo){
            case call:
                break;
            case put:
                delta = delta -1;
                break;
            default:
                //cout<<"Error"<<endl;
                break;
        }
        return delta;
    };
    T rho(){
        T delta{};
        T dplus{};
        T dminus{};
        T rho{};
        T logFK = log(forward/strike);
        T inverseSigma = 1 / (sigma * sqrt(tau));
        T medioSigma = 0.5 * sigma*sigma*tau;
        dplus = inverseSigma*(logFK + medioSigma);
        dminus = inverseSigma*(logFK - medioSigma);
        boost::math::normal normalDistribution;

        switch (tipo){
            case call:
                rho = strike * tau * exp(-interesAnual*tau) * normalCDF(dminus);
                break;
            case put:
                dminus = -dminus;
                rho = -strike * tau * exp(-interesAnual*tau) * normalCDF(dminus);
                break;
            default:
                //cout<<"Error"<<endl;
                break;
        }
        return rho;
    };

};

template<typename T>
class Call: public OptionGen<T> {
public:
    Greeks<T> griegas;
    int signo;
    T interesAnual;
    T strike;
    T spot;
    T sigma;
    double tau;

    T forward;
    T factorDescuento;

    Call(int signo,T interesAnual, T strike, T spot, T sigma, double tau):
            signo(signo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau), griegas(call,interesAnual,strike,spot,sigma,tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    }
    Call();

    map<T,T> evaluate(vector<T> values, double criteria){
        map<T,T> price{};
        unsigned long  posicion = getExpiry()*criteria;
        price.insert(std::make_pair(tau,signo*max(values.at(posicion) - strike,0.0)));

        return price;
    };
    double getExpiry(){
        return this->tau;
    };
};


template<typename T>
class Put: public OptionGen<T> {
public:
    Greeks<T> griegas;
    int signo;
    T interesAnual;
    T strike;
    T spot;
    T sigma;
    double tau;

    T forward;
    T factorDescuento;

    Put(int signo,T interesAnual, T strike, T spot, T sigma, double tau):
            signo(signo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau), griegas(put,interesAnual,strike,spot,sigma,tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    }
    Put();

    map<T,T> evaluate(vector<T> values, double criteria){
        map<T,T> price{};
        unsigned long posicion = getExpiry()*criteria;
        price.insert(make_pair(tau,signo*max(strike - values.at(posicion),0.0)));

        return price;
    }

    double getExpiry(){
        return this->tau;
    }
};
template<typename T>
class Asian: public Decorator<T>{
public:
    asianType tipo;


    Asian(asianType tipo,Call<T> *optionCall):tipo(tipo),Decorator<T>(optionCall){};

    Asian(asianType tipo,Put<T> *optionPut):tipo(tipo),Decorator<T>(optionPut){};

    map<T,T> evaluate(vector<T> values, double criteria){
        map<T,T> price{};
        double tau = getExpiry();
        unsigned long posicion = (tau*criteria);
        vector<T> result(values.size());

        for (unsigned long i = 0;i<=values.size();++i){
            result[i] = values[i];
        }
        if(!values.empty()) {
            if (tipo == max_) {
                auto it = max_element(values.begin(),values.end());
                result[posicion] = *it;
            } else if (tipo == min_) {
                 auto it = min_element(values.begin(), values.end());
                result[posicion] = *it;
            } else if (tipo == avg_) {
                T init = 0.0;
                auto it = std::accumulate(values.begin(), values.end(), init) / values.size();
                result[posicion] = it;
            } else {
                // cout << "ERROR" << endl;
            }
        }else{
            //cout << "VectorValues vacio" << endl;
        }

        price = Decorator<T>::innerOption->evaluate(result,criteria);
        return price;
    };
    double getExpiry(){
        return Decorator<T>::innerOption->getExpiry();
    };
};

#endif //MONTECARLOENGINE_OPTIONGEN_H
