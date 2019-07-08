//
// Created by alejandro on 17/06/19.
//

#include "Option.h"
#include <cmath>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <numeric>

using namespace std;
OptionBS::OptionBS (optionType tipo,double interesAnual, double strike, double spot, double sigma, double tau):
        tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau)
        {
            factorDescuento = exp(-interesAnual*tau);
            forward = spot/factorDescuento;
        }

        Call::Call(double interesAnual, double strike, double spot, double sigma, double tau):
        interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau), griegas(call,interesAnual,strike,spot,sigma,tau)
        {
            factorDescuento = exp(-interesAnual*tau);
            forward = spot/factorDescuento;
        }
Greeks::Greeks(const Greeks &p) : tipo(call),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau),factorDescuento(factorDescuento),forward(forward) {}

Put::Put(double interesAnual, double strike, double spot, double sigma, double tau):
        interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau),griegas(put,interesAnual,strike,spot,sigma,tau)
            {
            factorDescuento = exp(-interesAnual*tau);
            forward = spot/factorDescuento;
            }

Greeks::Greeks(optionType tipo,double interesAnual, double strike, double spot, double sigma, double tau):
        tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau)
        {
            factorDescuento = exp(-interesAnual*tau);
            forward = spot/factorDescuento;
        };

double Call::evaluate(vector<double> values){
    double price{};
    price = max(values.back() - strike,0.0);
   // cout<<"evaluate Call"<<endl;
   // cout<<price<<endl;
    return price;
};

double Put::evaluate(vector<double> values){
    double price{};
    price =  max(strike - values.back(),0.0);
    //cout<<"evaluate Put"<<endl;
    //cout<<price<<endl;
    return price;
};

/*vector<double> Asian::valuesAsian(vector<double> values){
    auto it = max_element(values.begin(),values.end());
    vector<double> result;
    result.push_back(*it);
    return result;
};*/
double Put::getExpiry(){
    return this->tau;
}
double Call::getExpiry() {
    return this->tau;
}
double Asian::getExpiry(){
    return innerOption->getExpiry();
}
double Asian::evaluate(vector<double> values){
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
            cout << "ERROR" << endl;
        }
    }else{
        cout << "VectorValues vacio" << endl;
    }

    v=innerOption->evaluate(result);
    //Decorator::evaluate(values);
    return v;
};


double Greeks::delta(){
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
            cout<<"Error"<<endl;
            break;
    }
    return delta;

}

double Greeks::vega(){
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
}



double Greeks::theta(){
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
            cout<<"Error"<<endl;
            theta = -11111111;
            break;
    }
    return theta;
}


double OptionBS::price(){

    double dplus{};
    double dminus{};
    double sigma_sqr = sigma*sigma;
    double inverseSigma = 1 / (sigma * sqrt(tau));
    double logFK = log(forward/strike);
    double medioSigma = 0.5 * sigma_sqr*tau;

    boost::math::normal normalDistribution;

    dplus = inverseSigma*(logFK + medioSigma);
    dminus = inverseSigma*(logFK - medioSigma);

    auto distribucionPlus = boost::math::cdf(normalDistribution,dplus);
    auto distribucionMinus = boost::math::cdf(normalDistribution,dminus);
    auto negDistribucionPlus = boost::math::cdf(normalDistribution,-dplus);
    auto negDistribucionMinus = boost::math::cdf(normalDistribution,-dminus);


double result{};

switch (tipo){
    case call:
       result = factorDescuento * ((distribucionPlus * forward) - (distribucionMinus * strike));
        break;
    case put:
        result = factorDescuento * ((negDistribucionMinus * strike) - (negDistribucionPlus * forward));
        break;
    default:
        cout<<"Error"<<endl;
        result = -11111111;
        break;
}
return result;
}

/*Option::Option(optionType tipo, double interesAnual, double strike, double spot, double sigma, double tau):
    tipo(tipo),interesAnual(interesAnual),strike(strike),spot(spot),sigma(sigma),tau(tau){

    factorDescuento = exp(-interesAnual*tau);
    forward = spot/factorDescuento;

}
Option::Option(){}

Asian::Asian(){}
double normalCDF(double x){
cout<<"normalCDF: "<<erfc(-x / sqrt(2))*0.5<<endl;
return erfc(-x / sqrt(2))*0.5;;
}
double CFD(double x){
const double root = sqrt(0.5);
cout<<"CFD: "<<  0.5*(1.0 + erf(x*root))<<endl;
return 0.5*(1.0 + erf(x*root));
}

double Option::evaluate() {
cout<<"Option"<<endl;
return 55;
}
double Option::delta(){
double delta{};
double dplus{};
double logFK = log(forward/strike);
double inverseSigma = 1 / (sigma * sqrt(tau));
double medioSigma = 0.5 * sigma*sigma*tau;
dplus = inverseSigma*(logFK + medioSigma);
boost::math::normal normalDistribution;
delta = boost::math::cdf(normalDistribution,dplus);
// normalCDF(dplus);
//CFD(dplus);
switch (tipo){
    case call:
        break;
    case put:
        delta = delta -1;
        break;
    default:
        cout<<"Error"<<endl;
        break;
}
return delta;

}

double Option::vega(){
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
}



double Option::theta(){
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
        cout<<"Error"<<endl;
        theta = -11111111;
        break;
}

return theta;
}

double Option::price(){

double dplus{};
double dminus{};
double sigma_sqr = sigma*sigma;
double inverseSigma = 1 / (sigma * sqrt(tau));
double logFK = log(forward/strike);
double medioSigma = 0.5 * sigma_sqr*tau;

boost::math::normal normalDistribution;

dplus = inverseSigma*(logFK + medioSigma);
dminus = inverseSigma*(logFK - medioSigma);

auto distribucionPlus = boost::math::cdf(normalDistribution,dplus);
auto distribucionMinus = boost::math::cdf(normalDistribution,dminus);
auto negDistribucionPlus = boost::math::cdf(normalDistribution,-dplus);
auto negDistribucionMinus = boost::math::cdf(normalDistribution,-dminus);

/*cout<<"factorDescuento"<< factorDescuento<<endl;
cout<<"inverseSigma"<< inverseSigma<<endl;
cout<<"logFK"<< logFK<<endl;
cout<<"medioSigma"<< medioSigma<<endl;
cout<<"dplus"<< dplus<<endl;
cout<<"dminus"<< dminus<<endl;
cout<<"distribucionPlus"<< distribucionPlus<<endl;
cout<<"distribucionMinus"<< distribucionMinus<<endl;
cout<<"negDistribucionPlus"<< negDistribucionPlus<<endl;
cout<<"negDistribucionMinus"<< negDistribucionMinus<<endl;

double result{};

switch (tipo){
    case call:
       result = factorDescuento * ((distribucionPlus * forward) - (distribucionMinus * strike));
        break;
    case put:
        result = factorDescuento * ((negDistribucionMinus * strike) - (negDistribucionPlus * forward));
        break;
    default:
        cout<<"Error"<<endl;
        result = -11111111;
        break;
}
return result;
}
double Asian::evaluate() {
cout<<"Asian"<<endl;
return 5.0;
}

*/