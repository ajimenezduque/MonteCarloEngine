//
// Created by alejandro on 17/06/19.
//

#ifndef MONTECARLOENGINE_OPTION_H
#define MONTECARLOENGINE_OPTION_H

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include<memory>

 using namespace std;
enum optionType { call, put };
enum asianType {max_, min_, avg_};
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

    double delta();
    double theta();
    double vega();

    Greeks();
    Greeks(Greeks const &g);
    Greeks(optionType tipo,double interesAnual, double strike, double spot, double sigma, double tau) ;
};
class Option{
public:
    //Greeks griegas;
    virtual double evaluate(vector<double> values) = 0;
};

class Decorator: public Option{
public:
    Decorator(Option * optionAbstract){
        innerOption = optionAbstract;
    }
    Option* innerOption;
    double evaluate(vector<double> values){
        innerOption->evaluate(values);
    };
};


class Composite: public Option{
public:
    vector<Option *> vectorInst;

  /*  void add (shared_ptr<Option> elem){
        vectorInst.push_back(elem);
    }*/
  void add (Option *elem){
      vectorInst.push_back(elem);
  }

    double evaluate(vector<double> values){
        for(auto element : vectorInst){
        element->evaluate(values);
        }
    }
    //implementar getMaturity

};
class Call: public Option {
public:
    Greeks griegas;
//variable de entrada

    double interesAnual;
    double strike;
    double spot;
    double sigma;
    double tau;

//variables intermedias
    double forward;
    double factorDescuento;

    Call (double interesAnual,double strike, double spot, double sigma, double tau);
    Call();
   // double price(double interesAnual,double strike, double spot, double sigma, double tau);
    double vega();
    double theta();
    double delta();
    double evaluate(vector<double> values);
};


class Put: public Option {
public:
    Greeks griegas;
//variable de entrada
    double interesAnual;
    double strike;
    double spot;
    double sigma;
    double tau;

//variables intermedias
    double forward;
    double factorDescuento;

    Put (double interesAnual,double strike, double spot, double sigma, double tau);
    Put();
    //double price();
    double vega();
    double theta();
    double delta();
    double evaluate(vector<double> values);
};
class Asian: public Decorator{
public:
    asianType tipo;
    Asian(double interesAnual,double strike, double spot, double sigma, double tau);

    Asian(asianType tipo,Call *optionCall):tipo(tipo),Decorator(optionCall){}

    Asian(asianType tipo,Put *optionPut):tipo(tipo),Decorator(optionPut){}

    double evaluate(vector<double> values);
};


/*class Instrumento{
public:

    virtual double evaluate() = 0;
};*/
/*class Option: public Instrumento {
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
    //boost::math::normal distribucion;
    double factorDescuento;

    Option (optionType tipo, double interesAnual,double strike, double spot, double sigma, double tau);
    Option();
    double price();
    double vega();
    double theta();
    double delta();
    double evaluate();
};*/



/*class Asian: public Option{
public:
    Asian ();
    asianType tipo;
    double evaluate();
};*/



#endif //MONTECARLOENGINE_OPTION_H
