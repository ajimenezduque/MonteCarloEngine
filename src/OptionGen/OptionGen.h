//
// Created by alejandro on 10/07/19.
//

#ifndef MONTECARLOENGINE_OPTIONGEN_H
#define MONTECARLOENGINE_OPTIONGEN_H

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <memory>
#include <numeric>

using namespace std;
enum optionType { call, put };
enum asianType {max_, min_, avg_};

template <typename T>
class OptionGen {
public:
    virtual OptionGen<T> *getOptionGen( int )
    {
        return 0;
    }

    virtual T evaluate(vector<T> values) = 0;
    virtual T getExpiry() = 0;
};

template <typename T>
class Decorator: public OptionGen<T>{
public:
    Decorator(OptionGen<T> * optionAbstract){
        innerOption = optionAbstract;
    }
    OptionGen<T>* innerOption;
    T evaluate(vector<T> values){
        innerOption->evaluate(values);
    };
    T getExpiry(){
        innerOption->getExpiry();
    }
};

template <typename T>
class Composite: public OptionGen<T>{
public:
    vector<OptionGen<T> *> vectorInst;

    /*  void add (shared_ptr<Option> elem){
          vectorInst.push_back(elem);
      }*/

    void add (OptionGen<T> *elem){
        vectorInst.push_back(elem);
    }

    void remove( const unsigned int index )
    {
        OptionGen<T> *child = vectorInst[ index ];
        vectorInst.erase( vectorInst.begin() + index );
        delete child;
    }

    int compSize(){
        return vectorInst.size();
    }

    OptionGen<T> *getOptionGen (int n){
        auto elem = vectorInst[n];
        return elem;
    }
    //vector de pares<fecha,valor>
    T evaluate(vector<T> values){
        for(auto element : vectorInst){
            element->evaluate(values);
        }
    }
    T  getExpiry(){
        double  maxExpiry = 0.0;
        for (const auto element : vectorInst){
            maxExpiry = max(maxExpiry,element->getExpiry());
        }
        return maxExpiry;
    }

};

template<typename T>
class Greeks{
public:
    //variable de entrada
    optionType tipo;
    double interesAnual;
    double strike;
    double spot;
    double sigma;
    double tau;

//variables intermedias
    double forward;
    double factorDescuento;

    Greeks();
    //Greeks(Greeks const &g);
    Greeks(optionType tipo,double interesAnual, double strike, double spot, double sigma, double tau):
            tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    };

    T vega(){
        double vega{};
        double dplus{};
        double logFK = log(forward/strike);
        double inverseSigma = 1 / (sigma * sqrt(tau));
        double medioSigma = 0.5 * sigma*sigma*tau;
        dplus = inverseSigma*(logFK + medioSigma);

        double spotSqrTau = spot*sqrt(tau);
        double distribucionPrima = (1/(sqrt(2*M_PI)))*exp(-(dplus*dplus)/2);
        vega = spotSqrTau*distribucionPrima;
        return vega;
    };
    T theta(){
        double theta{};
        double dplus{};
        double dminus{};
        double logFK = log(forward/strike);
        double inverseSigma = 1 / (sigma * sqrt(tau));
        double medioSigma = 0.5 * sigma*sigma*tau;

        dplus = inverseSigma*(logFK + medioSigma);
        double distribucionPrima = (1/(sqrt(2*M_PI)))*exp(-(dplus*dplus)/2);


        dminus = inverseSigma*(logFK - medioSigma);

        boost::math::normal normalDistribution;
        double distribucionMinus{};
        switch (tipo){
            case call:
                distribucionMinus =  boost::math::cdf(normalDistribution,dminus);
                //cout<<"Call Distribucion Minus: "<< distribucionMinus<<endl;
                //normalCDF(dminus);
                theta = -(spot*distribucionPrima*sigma)/(2*sqrt(tau)) - (interesAnual*strike*exp(-interesAnual*tau)*(distribucionMinus));
                break;
            case put:
                distribucionMinus =  boost::math::cdf(normalDistribution,-dminus);
                //cout<<" Put Distribucion Minus: "<< distribucionMinus<<endl;
                //normalCDF(-dminus);
                //CFD(-dminus);
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
        double delta{};
        double dplus{};
        double logFK = log(forward/strike);
        double inverseSigma = 1 / (sigma * sqrt(tau));
        double medioSigma = 0.5 * sigma*sigma*tau;
        dplus = inverseSigma*(logFK + medioSigma);
        boost::math::normal normalDistribution;
        delta = boost::math::cdf(normalDistribution,dplus);
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

};

template<typename T>
class Call: public OptionGen<T> {
public:
    Greeks<T> griegas;
//variable de entrada

    T interesAnual;
    T strike;
    T spot;
    T sigma;
    T tau;

//variables intermedias
    T forward;
    T factorDescuento;

    Call(T interesAnual, T strike, T spot, T sigma, T tau):
            interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau), griegas(call,interesAnual,strike,spot,sigma,tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    }
    Call();
    // double price(double interesAnual,double strike, double spot, double sigma, double tau);

    //añadir frecuencia de muestreo
    T evaluate(vector<T> values){
        T price{};
        //sacar valor coincidente en el tiempo
        price = max(values.back() - strike,0.0);
        return price;
    };
    T getExpiry(){
        return this->tau;
    };
};

template<typename T>
class Put: public OptionGen<T> {
public:
    Greeks<T> griegas;
//variable de entrada
    T interesAnual;
    T strike;
    T spot;
    T sigma;
    T tau;

//variables intermedias
    T forward;
    T factorDescuento;

    Put(T interesAnual, T strike, T spot, T sigma, T tau):
            interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau), griegas(put,interesAnual,strike,spot,sigma,tau)
    {
        factorDescuento = exp(-interesAnual*tau);
        forward = spot/factorDescuento;
    }
    Put();
    //double price();
    T evaluate(vector<T> values){
        double price{};
        price =  max(strike - values.back(),0.0);
        return price;
    };
    T getExpiry(){
        return this->tau;
    }
};
template<typename T>
class Asian: public Decorator<T>{
public:
    asianType tipo;

    Asian(T interesAnual,T strike, T spot, T sigma, T tau);

    Asian(asianType tipo,Call<T> *optionCall):tipo(tipo),Decorator<T>(optionCall){}

    Asian(asianType tipo,Put<T> *optionPut):tipo(tipo),Decorator<T>(optionPut){}

    T evaluate(vector<double> values){
        double v{};
        vector<double> result;
        if(!values.empty()) {
            if (tipo == max_) {
                auto it = max_element(values.begin(), values.end());
                result.push_back(*it);
            } else if (tipo == min_) {
                auto it = min_element(values.begin(), values.end());
                result.push_back(*it);
            } else if (tipo == avg_) {
                v = accumulate(values.begin(), values.end(), 0.0) / values.size();
                result.push_back(v);
            } else {
                // cout << "ERROR" << endl;
            }
        }else{
            //cout << "VectorValues vacio" << endl;
        }

        v=Decorator<T>::innerOption->evaluate(result);
        //v=Decorator<T>::evaluate(values);
        return v;
    };
    T getExpiry(){
        return Decorator<T>::innerOption->getExpiry();
    };
};

#endif //MONTECARLOENGINE_OPTIONGEN_H
