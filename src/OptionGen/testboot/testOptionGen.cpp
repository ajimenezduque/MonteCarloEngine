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

double SimpleMonteCarlo2(
                                 OptionGen<double> &option,
                                 const double Spot,//dato de mercado dejar como input
                                 const map<double,double> &Vol,//input
                                 const double r,//input
                                 const unsigned long NumberOfPaths,
                                 const unsigned long NumberOfSamples)
{

    std::mt19937 gen(1234);

    std::normal_distribution<> normalDistribution(0,1);

    double maxExpiry = option.getExpiry();

    double dt = maxExpiry/NumberOfSamples;

    std::vector<double> underlying_values(NumberOfSamples + 1, Spot);

    map<double,double> runningSum{};

    double thisGaussian {};
    double diffusion{};
    double drift{};
    map<double, double> thisPayoff{};

    for (unsigned long i = 0; i < NumberOfPaths; i++) {
        double posIni = 0;
        for (unsigned long j = 1; j < underlying_values.size(); j++) {
            posIni = posIni + dt*j;
            if(posIni > (Vol.rbegin())->first){
                posIni = Vol.rbegin()->first;
            }
            auto it = Vol.lower_bound(posIni);
            //  cout<<"IT: "<<posIni<<" Frist: "<<it->first <<" Second: "<<it->second<<endl;
            double variance = it->second * it->second;
            double rootVariance = sqrt(variance * dt);
            double movedSpotFactor = exp((r - (0.5 * variance))*dt);
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
    double value{};

    for (auto it = runningSum.begin(); it != runningSum.end(); ++it) {
        value += it->second;
    }

    return value/NumberOfPaths;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen_Portfolio){
    BOOST_TEST_MESSAGE("Se ejecuta test de valoración para un portfolio de opciones.");
    cout<<"///************TEST****************////"<<endl;

    double interes = 0.08;
    double spot = 305;
    double sigma = 0.25;
    map<double,double> sigma_map;
    sigma_map[0.5] = 2.8;
    sigma_map[1.0] = 3.0;
    sigma_map[2.0] = 3.2;
    sigma_map [2.5] = 3.4;
    unsigned long paths = 10000;
    unsigned long samples = 365;
    vector<double> strikes {275,325,350};
    vector<double> vencimientos {4,3,2,4.0/12.0};


    Composite<double> myOptions1;

     int par = 0;
     for(int i = 0;i < vencimientos.size();++i){
         for(int j = 0; j < strikes.size();++j){
             if(par%2 == 0){
                 Call<double> *option = new Call<double>(1,interes,strikes[j],spot,interes,vencimientos[i]);
                 myOptions1.add(option);
             }else{
                 Put<double> *option = new Put<double>(1,interes,strikes[j],spot,interes,vencimientos[i]);
                 myOptions1.add(option);
             }
                 par++;
         }
     }

    auto start = std::chrono::high_resolution_clock::now();
    double valoracionBase = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);

    double deltaPrice = 0.00001;
    double bumpedSpot = SimpleMonteCarlo2(myOptions1,spot + deltaPrice,sigma_map,interes,paths,samples);

    sigma_map[0.5] = sigma_map[0.5] + deltaPrice;
    double bumpedVol = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    sigma_map[0.5] = 2.8;
    sigma_map[1.0] = sigma_map[1.0] + deltaPrice;
    double bumpedVol2 = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    sigma_map[1.0] = 3.0;
    sigma_map[2.0] = sigma_map[2.0] + deltaPrice;
    double bumpedVol3 = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    sigma_map[2.0] = 3.2;
    sigma_map[2.5] = sigma_map[2.5] + deltaPrice;
    double bumpedVol4 = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    sigma_map [2.5] = 3.4;
    double bumpedInt = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes + deltaPrice,paths,samples);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds=end-start;
    cout<<"Delta Bumped: "<< ((bumpedSpot - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped: "<< ((bumpedVol - valoracionBase)/deltaPrice)<<endl;
        cout<<"Vega Bumped2: "<< ((bumpedVol2 - valoracionBase)/deltaPrice)<<endl;
        cout<<"Vega Bumped3: "<< ((bumpedVol3 - valoracionBase)/deltaPrice)<<endl;
        cout<<"Vega Bumped4: "<< ((bumpedVol4 - valoracionBase)/deltaPrice)<<endl;

    cout<<"Rho Bumped: "<< ((bumpedInt - valoracionBase)/deltaPrice)<<endl;
    cout<<"Valoracion Portfolio: "<<valoracionBase<<endl;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

}
BOOST_AUTO_TEST_CASE(Test_OptionGen_Opcion){
    BOOST_TEST_MESSAGE("Se ejecuta test de valoración para una única opción.");
    cout<<"///************TEST****************////"<<endl;

    double interes = 0.08;
    double spot = 305;
    double strike = 300;
    double sigma = 0.25;
    map<double,double> sigma_map;
    sigma_map[0.5] = 2.8;
    sigma_map[1.0] = 3.0;
    sigma_map[2.0] = 3.2;
    sigma_map[2.5] = 3.4;
    unsigned long paths = 10000;
    unsigned long samples = 365;

    ///Prueba calculo de sensibilidades con unica opcion

    Call<double> opcionCallVega(1,interes, strike, spot, sigma, 4.0 / 12.0);

    auto start = std::chrono::high_resolution_clock::now();
    double valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);

    double deltaPrice = 0.00001;

    double bumpedSpot = SimpleMonteCarlo2(opcionCallVega,spot + deltaPrice,sigma_map,interes,paths,samples);

    sigma_map[0.5]=sigma_map[0.5] + deltaPrice;
    double bumpedVol = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[0.5] = 2.8;
    sigma_map[1.0] = sigma_map[1.0] + deltaPrice;
    double bumpedVol2 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[1.0] = 3.0;
    sigma_map[2.0] = sigma_map[2.0] + deltaPrice;
    double bumpedVol3 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[2.0] = 3.2;
    sigma_map[2.5] = sigma_map[2.5] + deltaPrice;
    double bumpedVol4 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map [2.5] = 3.4;

    double bumpedInt = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes + deltaPrice,paths,samples);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds=end-start;
    cout<<"Delta Bumped: "<< ((bumpedSpot - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped: "<< ((bumpedVol - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped2: "<< ((bumpedVol2 - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped3: "<< ((bumpedVol3 - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped4: "<< ((bumpedVol4 - valoracionBase)/deltaPrice)<<endl;

    cout<<"Rho Bumped: "<< ((bumpedInt - valoracionBase)/deltaPrice)<<endl;
    cout<<"Valoracion Portfolio: "<<valoracionBase<<endl;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

}
BOOST_AUTO_TEST_CASE(Test_OptionGen_OpcionBS){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Black-Scholes, se comprueba valoración de BS y MC, así como sensibilidades análiticas y mediante DA");
    cout<<"///************TEST****************////"<<endl;
    double interes = 0.08;
    double spot = 305;
    double strike = 300;
    double sigma = 0.25;
    map<double,double> sigma_map;
    sigma_map[0.5] = sigma;
    sigma_map[1.0] = sigma;
    sigma_map[2.0] = sigma;
    sigma_map[2.5] = sigma;
    unsigned long paths = 10000;
    unsigned long samples = 365;

    ///Prueba calculo de sensibilidades con unica opcion

    Call<double> opcionCallVega(1,interes, strike, spot, sigma, 4.0 / 12.0);
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;
    cout<<"Rho: "<<opcionCallVega.griegas.rho()<<endl;
    OptionBS <double> opcionBS (call,interes,strike,spot,sigma, 4.0/12.0);
    cout<<"OptionBS Price:" <<opcionBS.price()<<endl;

    auto start = std::chrono::high_resolution_clock::now();
    double valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);

    double deltaPrice = 0.00001;

    double bumpedSpot = SimpleMonteCarlo2(opcionCallVega,spot + deltaPrice,sigma_map,interes,paths,samples);

    ///Alteracion del valor para cada una de las volatilidades del mapa.
    sigma_map[0.5] = sigma_map[0.5] + deltaPrice;
    double bumpedVol = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[0.5] = sigma;
    sigma_map[1.0] = sigma_map[1.0] + deltaPrice;
    double bumpedVol2 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[1.0] = sigma;
    sigma_map[2.0] = sigma_map[2.0] + deltaPrice;
    double bumpedVol3 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map[2.0] = sigma;
    sigma_map[2.5] = sigma_map[2.5] + deltaPrice;
    double bumpedVol4 = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    sigma_map [2.5] = sigma;

    double bumpedInt = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes + deltaPrice,paths,samples);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds=end-start;
    cout<<"Delta Bumped: "<< ((bumpedSpot - valoracionBase)/deltaPrice)<<endl;

    cout<<"Vega Bumped: "<< ((bumpedVol - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped2: "<< ((bumpedVol2 - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped3: "<< ((bumpedVol3 - valoracionBase)/deltaPrice)<<endl;
    cout<<"Vega Bumped4: "<< ((bumpedVol4 - valoracionBase)/deltaPrice)<<endl;

    cout<<"Rho Bumped: "<< ((bumpedInt - valoracionBase)/deltaPrice)<<endl;

    cout<<"Valoracion OptionBS por MC: "<<valoracionBase<<endl;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

}

BOOST_AUTO_TEST_CASE(Test_OptionGen_Asian){
    BOOST_TEST_MESSAGE("Se ejecuta test Asian");
    cout<<"///************TEST****************////"<<endl;
    double interes = 0.08;
    double spot = 305;
    double strike = 300;
    double sigma = 0.25;
    map<double,double> sigma_map;
    sigma_map[0.5] = sigma;
    sigma_map[1.0] = sigma;
    sigma_map[2.0] = sigma;
    sigma_map[2.5] = sigma;
    unsigned long paths = 10000;
    unsigned long samples = 365;

    ///Prueba calculo con opciones Asiaticas

  Asian<double> optionAsian1(avg_,new Call<double>(1,0.08,300.0,305.0,0.25,4));
  Asian<double> *optionAsian2 = new Asian<double>(max_,new Call<double>(1,0.08,300.0,305.0,0.25,4));
  Asian<double> *optionAsian3 = new Asian<double>(min_,new Call<double>(1,0.08,300.0,305.0,0.25,4));
    Composite<double> myOptions1;
    myOptions1.add(&optionAsian1);
    myOptions1.add(optionAsian2);
    myOptions1.add(optionAsian3);
    double valoracionAsian = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    cout<<"Valoracion Asiaticas: "<<valoracionAsian<<endl;

}