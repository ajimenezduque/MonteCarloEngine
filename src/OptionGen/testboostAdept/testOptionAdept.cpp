//
// Created by alejandro on 22/07/19.
//

//
// Created by alejandro on 10/07/19.
//
#define BOOST_TEST_MODULE pruebaOptionAdept

#include <boost/test/included/unit_test.hpp>
//#include <boost/test/included/unit_test_framework.hpp>
#include "OptionGen/OptionGen.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <chrono>
#include <adept.h>

using namespace std;
using namespace adept;

//double
adouble SimpleMonteCarlo2(
        OptionGen<adouble> &option,
        const adouble Spot,
        const map<adouble,adouble> &Vol,
        const adouble r,
        unsigned long NumberOfPaths,
        unsigned long NumberOfSamples)
{

    //use this construction for a random generator with a fixed random sequence
    std::mt19937 gen(1234);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> normalDistribution(0,1);

    double maxExpiry = option.getExpiry();

    adouble dt = maxExpiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<adouble> underlying_values(NumberOfSamples + 1,Spot); //x[0]);

   // adouble variance = Vol * Vol;//x[1] * x[1] ;
    //adouble rootVariance = sqrt(variance * dt);

    //adouble movedSpotFactor = exp((r - (0.5 * variance))*dt);

    map<adouble,adouble> runningSum{};

    adouble thisGaussian {};
    adouble diffusion{};
    adouble drift{};
    map<adouble, adouble> thisPayoff{};
    /*Otra variante
     auto it = Vol.begin();
     posIni = it->first;
     posIni*/

    for (unsigned long i = 0; i < NumberOfPaths; i++) {
        adouble posIni = 0;
        for (unsigned long j = 1; j < underlying_values.size(); j++) {
             //double pos = (j * maxExpiry)/NumberOfSamples;
             posIni = posIni + dt*j;
             if(posIni > (Vol.rbegin())->first){
                 posIni = Vol.rbegin()->first;
             }
             auto it = Vol.lower_bound(posIni);
             cout<<"IT: "<<posIni<<" Frist: "<<it->first <<" Second: "<<it->second<<endl;
             adouble variance = it->second * it->second;
             adouble rootVariance = sqrt(variance * dt);
             adouble movedSpotFactor = exp((r - (0.5 * variance))*dt);
             thisGaussian = normalDistribution(gen);
             diffusion = exp(rootVariance * thisGaussian);
             drift = underlying_values[j - 1] * movedSpotFactor;
             underlying_values[j] = drift * diffusion;
        }

        thisPayoff = option.evaluate(underlying_values, NumberOfSamples/maxExpiry);

        for (auto it = thisPayoff.begin(); it != thisPayoff.end(); ++it) {
           runningSum[it->first] = runningSum[it->first] + (exp(-r * it->first) * (it->second));
        }
    }
    adouble value{};
    for (auto it = runningSum.begin(); it != runningSum.end(); ++it) {
        // cout<< it->first << " => " << it->second << '\n';
        //value += (exp(-r * it->first) * (it->second / NumberOfPaths));
        value+= it->second;
    }

    //value = exp(-r * maxExpiry) * running;

    adouble NumberOfPaths1 = NumberOfPaths;
    return value/NumberOfPaths1;
}

BOOST_AUTO_TEST_CASE(Test_OptionGenAdept){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");

    Stack stack;
    adouble x[3];
    x[0].set_value(17);
    x[1].set_value(0.25);
    x[2].set_value(0.03);
  //  Call<adouble> opcionCallDelta1 (1,0.08,100,305,0.25,4.0);
    //Put<adouble> opcionPutTheta1 (1,0.08,300.0,305.0,0.25,4.0/12.0);
    //Call<adouble> opcionCallVega1 (1,0.08,300,305,0.25,4.0/12.0);
    Composite <adouble> myOptions1;
   /* Composite<adouble> myOptions;
    myOptions.add(&opcionCallDelta1);
    myOptions.add(&opcionPutTheta1);
    myOptions.add(&opcionCallVega1);*/

   /* Call<adouble> opcionCall1 (1,0.08,275,305,0.25,4.0);
    Call<adouble> opcionCall2 (1,0.08,350,305,0.25,4.0);
    Call<adouble> opcionCall3 (1,0.08,325,305,0.25,3.0);
    Call<adouble> opcionCall4 (1,0.08,275,305,0.25,2.0);
    Call<adouble> opcionCall5 (1,0.08,350,305,0.25,2.0);
    Call<adouble> opcionCall6 (1,0.08,325,305,0.25,4.0/12.0);

    Put<adouble> opcionPut1 (1,0.08,325.0,305.0,0.25,4.0);
    Put<adouble> opcionPut2 (1,0.08,275.0,305.0,0.25,3.0);
    Put<adouble> opcionPut3 (1,0.08,350.0,305.0,0.25,3.0);
    Put<adouble> opcionPut4 (1,0.08,325.0,305.0,0.25,2);
    Put<adouble> opcionPut5 (1,0.08,275,305.0,0.25,4.0/12.0);
    Put<adouble> opcionPut6 (1,0.08,350.0,305.0,0.25,4.0/12.0);

    myOptions1.add(&opcionCall1);
    myOptions1.add(&opcionPut1);
    myOptions1.add(&opcionCall2);
    myOptions1.add(&opcionPut2);
    myOptions1.add(&opcionCall3);
    myOptions1.add(&opcionPut3);
    myOptions1.add(&opcionCall4);
    myOptions1.add(&opcionPut4);
    myOptions1.add(&opcionCall5);
    myOptions1.add(&opcionPut5);
    myOptions1.add(&opcionCall6);
    myOptions1.add(&opcionPut6);*/


    adouble interes = 0.08;
    adouble spot = 305.0;
    //adouble sigma = 0.25;
    map<adouble,adouble> sigma;
    sigma[0.5] = 2.8;
    sigma[1.0] = 3.0;
    sigma[2.0] = 3.2;
    sigma [2.5] = 3.4;
    unsigned long paths = 1000;
    unsigned long samples = 365;
    vector<adouble> strikes {275,325,350};
    vector<double> vencimientos {4,3,2,4.0/12.0};


    int par = 0;
    for(int i = 0;i < vencimientos.size();++i){
        for(int j = 0; j < strikes.size();++j){
            if(par%2 == 0){
                Call<adouble> *option = new Call<adouble>(1,0.08,strikes[j],305,0.25,vencimientos[i]);
                myOptions1.add(option);
            }else{
                Put<adouble> *option = new Put<adouble>(1,0.08,strikes[j],305,0.25,vencimientos[i]);
                myOptions1.add(option);
            }
            par++;
        }
    }

    /*Call<adouble> ejercicio (1,0.03,15,17,0.25,4.0/12.0);
    cout<<"Delta: "<<ejercicio.griegas.delta()<<endl;
    cout<<"Vega: "<<ejercicio.griegas.vega()<<endl;*/
    Call<adouble> opcionCallVega(1,0.08, 300, 305, 0.25, 4.0 / 12.0);
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;
    stack.new_recording();

    auto start = std::chrono::high_resolution_clock::now();
   //adouble valoracionBase = SimpleMonteCarlo2(myOptions1,spot,sigma,interes,paths,samples);
    adouble valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma,interes,paths,samples);
    valoracionBase.set_gradient(1.0);
    stack.compute_adjoint();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    cout << "Tiempo de cálculo: "<<diff.count()<<" segundos"<< endl;

    cout<<"Delta: "<<spot.get_gradient()<<endl;
    cout<<"Vega: "<<sigma[0.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[1.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[2.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[2.5].get_gradient()<<endl;
    cout<<"Rho: "<<interes.get_gradient()<<endl;

    cout<<"Valoracion Portfolio: "<<valoracionBase<<endl;



}
