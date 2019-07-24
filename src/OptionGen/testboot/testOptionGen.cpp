//
// Created by alejandro on 10/07/19.
//
#define BOOST_TEST_MODULE pruebaOptionGen

#include <boost/test/included/unit_test.hpp>
//#include <boost/test/included/unit_test_framework.hpp>
#include "../OptionGen.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <chrono>


using namespace std;
//Spot, volatilidad y interes dependen del momento actual por lo que son inputs montcarlo
//double
double SimpleMonteCarlo2(
                                 OptionGen<double> &option,
                                 const double Spot,//dato de mercado dejar como input
                                 const double Vol,//input
                                 const double r,//input
                                 const unsigned long NumberOfPaths,
                                 const unsigned long NumberOfSamples)
{
    //use this construction for a random generator with a fixed random sequence
    std::mt19937 gen(1234);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> nomalDistribution(0,1);

    double maxExpiry = option.getExpiry();

    double dt = maxExpiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<double> underlying_values(NumberOfSamples + 1, Spot);

    auto variance = Vol * Vol ;
    auto rootVariance = std::sqrt(variance * dt);

    auto movedSpotFactor = std::exp((r - (0.5 * variance))*dt);

    map<double,double> runningSum{};
    
    for (unsigned long i = 0; i < NumberOfPaths; i++) {
        for (unsigned long j = 1; j < underlying_values.size(); j++) {
            auto thisGaussian = nomalDistribution(gen);
            auto diffusion = std::exp(rootVariance * thisGaussian);
            auto drift = underlying_values[j - 1] * movedSpotFactor;
            underlying_values[j] = drift * diffusion;
        }

        map<double, double> thisPayoff{};

        thisPayoff = option.evaluate(underlying_values, NumberOfSamples / maxExpiry);
    
    
        for (auto it = thisPayoff.begin(); it != thisPayoff.end(); ++it) {
            runningSum[it->first] = runningSum[it->first] + (exp(-r * it->first) * (it->second / NumberOfPaths));
        }
    }

    double value{};
    for (auto it = runningSum.begin(); it != runningSum.end(); ++it) {
       // cout<< it->first << " => " << it->second << '\n';
        //value += (exp(-r * it->first) * (it->second / NumberOfPaths));
        value += it->second;
    }
    return value;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");
    double interes = 0.08;
    double spot = 305;
    double sigma = 0.25;
    unsigned long paths = 10000;
    unsigned long samples = 365;
    vector<double> strikes {275,325,350};
    vector<double> vencimientos {4,3,2,4.0/12.0};


    Composite<double> myOptions1;
   /* Call<double> opcionCall1 (1,0.08,275,305,0.25,4.0);
    Call<double> opcionCall2 (1,0.08,350,305,0.25,4.0);
    Call<double> opcionCall3 (1,0.08,325,305,0.25,3.0);
    Call<double> opcionCall4 (1,0.08,275,305,0.25,2.0);
    Call<double> opcionCall5 (1,0.08,350,305,0.25,2.0);
    Call<double> opcionCall6 (1,0.08,325,305,0.25,4.0/12.0);

    Put<double> opcionPut1 (1,0.08,325.0,305.0,0.25,4.0);
    Put<double> opcionPut2 (1,0.08,275.0,305.0,0.25,3.0);
    Put<double> opcionPut3 (1,0.08,350.0,305.0,0.25,3.0);
    Put<double> opcionPut4 (1,0.08,325.0,305.0,0.25,2);
    Put<double> opcionPut5 (1,0.08,275,305.0,0.25,4.0/12.0);
    Put<double> opcionPut6 (1,0.08,350.0,305.0,0.25,4.0/12.0);

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



     int par = 0;
     for(int i = 0;i < vencimientos.size();++i){
         for(int j = 0; j < strikes.size();++j){
             if(par%2 == 0){
                 Call<double> *option = new Call<double>(1,interes,strikes[j],spot,sigma,vencimientos[i]);
                 myOptions1.add(option);
             }else{
                 Put<double> *option = new Put<double>(1,interes,strikes[j],spot,sigma,vencimientos[i]);
                 myOptions1.add(option);
             }
                 par++;
         }
     }
     ///Prueba calculo de sensibilidades con unica opcion
    /*Call<double> ejercicio (1,0.03,15,17,0.25,4.0/12.0);
     cout<<"Delta: "<<ejercicio.griegas.delta()<<endl;
    cout<<"Vega: "<<ejercicio.griegas.vega()<<endl;
    */
    Call<double> opcionCallVega(1,0.08, 300, 305, 0.25, 4.0 / 12.0);
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;

    auto start = std::chrono::high_resolution_clock::now();
    double valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma,interes,paths,samples);

    double deltaPrice = 0.00001;

    double bumpedSpot = SimpleMonteCarlo2(opcionCallVega,spot + deltaPrice,sigma,interes,paths,samples);


    double bumpedVol = SimpleMonteCarlo2(opcionCallVega,spot,sigma + deltaPrice,interes,paths,samples);


    double bumpedInt = SimpleMonteCarlo2(opcionCallVega,spot,sigma,interes + deltaPrice,paths,samples);


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff=end-start;
    cout<<"Delta Bumped: "<< ((bumpedSpot - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped: "<< ((bumpedVol - valoracionBase)/deltaPrice)<<endl;
    cout<<"Rho Bumped: "<< ((bumpedInt - valoracionBase)/deltaPrice)<<endl;
    cout<<"Valoracion Portfolio: "<<valoracionBase<<endl;
    cout << "Tiempo de cálculo: "<<diff.count()<<" segundos"<< endl;

    /*Put opcionPut (0.08,300.0,305.0,0.25,4);
    OptionBS putBS (put,0.08,300.0,305.0,0.25,4);

    Asian *optionAsian1 = new Asian(avg_,new Call(0.08,300.0,305.0,0.25,4));
    Asian *optionAsian2 = new Asian(max_,new Call(0.08,300.0,305.0,0.25,4));
    Asian *optionAsian3 = new Asian(min_,new Call(0.08,300.0,305.0,0.25,4));*/

}
