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
enum asianType {max, min, avg};

class Instrumento{
public:

    virtual double evaluate() = 0;
};

class Composite: public Instrumento{
public:

    vector<shared_ptr<Instrumento>> vectorInst;

    void add (shared_ptr<Instrumento> elem){
        vectorInst.push_back(elem);
    }

    double evaluate(){

        for(auto element : vectorInst){
        element->evaluate();
        }

    }

};
class Option: public Instrumento {
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
};

class Asian: public Option{
public:
    Asian ();
    asianType tipo;
    double evaluate();
};



#endif //MONTECARLOENGINE_OPTION_H
