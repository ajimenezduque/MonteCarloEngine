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
    std::mt19937 gen(1234);

    std::normal_distribution<> normalDistribution(0,1);

    double maxExpiry = option.getExpiry();

    adouble dt = maxExpiry/NumberOfSamples;

    std::vector<adouble> underlying_values(NumberOfSamples + 1,Spot); //x[0]);

    map<adouble,adouble> runningSum{};

    adouble thisGaussian {};
    adouble diffusion{};
    adouble drift{};
    map<adouble, adouble> thisPayoff{};

    for (unsigned long i = 0; i < NumberOfPaths; i++) {
        adouble posIni = 0;
        for (unsigned long j = 1; j < underlying_values.size(); j++) {
             posIni = posIni + dt*j;
             if(posIni > (Vol.rbegin())->first){
                 posIni = Vol.rbegin()->first;
             }
             auto it = Vol.lower_bound(posIni);
           //  cout<<"IT: "<<posIni<<" Frist: "<<it->first <<" Second: "<<it->second<<endl;
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
        value+= it->second;
    }

    adouble NumberOfPaths1 = NumberOfPaths;
    return value/NumberOfPaths1;
}

BOOST_AUTO_TEST_CASE(Test_OptionGenAdept_Portfolio){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");
    cout<<"///************TEST****************////"<<endl;

    Stack stack;
    Composite <adouble> myOptions1;

    adouble interes = 0.08;
    adouble spot = 305.0;
    //adouble sigma = 0.25;
    map<adouble,adouble> sigma;
    sigma[0.5] = 2.8;
    sigma[1.0] = 3.0;
    sigma[2.0] = 3.2;
    sigma [2.5] = 3.4;
    unsigned long paths = 10000;
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

    stack.new_recording();

    auto start = std::chrono::high_resolution_clock::now();
    adouble valoracionBase = SimpleMonteCarlo2(myOptions1,spot,sigma,interes,paths,samples);
    valoracionBase.set_gradient(1.0);
    stack.compute_adjoint();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds = end-start;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

    cout<<"Delta: "<<spot.get_gradient()<<endl;
    cout<<"Vega: "<<sigma[0.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[1.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[2.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma[2.5].get_gradient()<<endl;
    cout<<"Rho: "<<interes.get_gradient()<<endl;

    cout<<"Valoracion Portfolio: "<<valoracionBase<<endl;

}
BOOST_AUTO_TEST_CASE(Test_OptionGenAdept_Opcion){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");
    cout<<"///************TEST****************////"<<endl;

    Stack stack;

    adouble interes = 0.08;
    adouble spot = 305.0;
    adouble sigma = 0.25;
    adouble strike = 300;
    map<adouble,adouble> sigma_map;
    sigma_map[0.5] = 2.8;
    sigma_map[1.0] = 3.0;
    sigma_map[2.0] = 3.2;
    sigma_map [2.5] = 3.4;
    unsigned long paths = 10000;
    unsigned long samples = 365;


    Call<adouble> opcionCallVega(1,interes, strike, spot,sigma, 4.0 / 12.0);
    stack.new_recording();

    auto start = std::chrono::high_resolution_clock::now();
    adouble valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    valoracionBase.set_gradient(1.0);
    stack.compute_adjoint();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds = end-start;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

    cout<<"Delta: "<<spot.get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[0.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[1.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.5].get_gradient()<<endl;
    cout<<"Rho: "<<interes.get_gradient()<<endl;

    cout<<"Valoracion Opcion: "<<valoracionBase<<endl;
}
BOOST_AUTO_TEST_CASE(Test_OptionGenAdept_OpcionBS){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Black-Scholes, se comprueba valoración de BS y MC, así como sensibilidades análiticas y mediante DA ");
    cout<<"///************TEST****************////"<<endl;
    Stack stack;

    adouble interes = 0.08;
    adouble spot = 305.0;
    adouble sigma = 0.25;
    adouble strike = 300;
    map<adouble,adouble> sigma_map;
    sigma_map[0.5] = sigma;
    sigma_map[1.0] = sigma;
    sigma_map[2.0] = sigma;
    sigma_map [2.5] = sigma;
    unsigned long paths = 10000;
    unsigned long samples = 365;


    Call<adouble> opcionCallVega(1,interes, strike, spot,sigma, 4.0 / 12.0);
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;
    cout<<"Rho: "<<opcionCallVega.griegas.rho()<<endl;
    OptionBS <adouble> opcionBS (call,interes,strike,spot,sigma, 4.0/12.0);
    cout<<"OptionBS Price:" <<opcionBS.price()<<endl;

    stack.new_recording();

    auto start = std::chrono::high_resolution_clock::now();
    adouble valoracionBase = SimpleMonteCarlo2(opcionCallVega,spot,sigma_map,interes,paths,samples);
    valoracionBase.set_gradient(1.0);
    stack.compute_adjoint();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_seconds = end-start;
    cout << "Tiempo de cálculo: "<<diff_seconds.count()<<" segundos"<< endl;

    cout<<"Delta: "<<spot.get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[0.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[1.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.5].get_gradient()<<endl;
    cout<<"Rho: "<<interes.get_gradient()<<endl;

    cout<<"Valoracion Opcion: "<<valoracionBase<<endl;
}
BOOST_AUTO_TEST_CASE(Test_OptionGenAdept_Asian){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Asiatica");
    cout<<"///************TEST****************////"<<endl;
    Stack stack;

    adouble interes = 0.08;
    adouble spot = 305.0;
    adouble sigma = 0.25;
    adouble strike = 300;
    map<adouble,adouble> sigma_map;
    sigma_map[0.5] = sigma;
    sigma_map[1.0] = sigma;
    sigma_map[2.0] = sigma;
    sigma_map [2.5] = sigma;
    unsigned long paths = 10000;
    unsigned long samples = 365;


    Asian<adouble> optionAsian1(avg_,new Call<adouble>(1,0.08,300.0,305.0,0.25,4));
    Asian<adouble> *optionAsian2 = new Asian<adouble>(max_,new Call<adouble>(1,0.08,300.0,305.0,0.25,4));
    Asian<adouble> *optionAsian3 = new Asian<adouble>(min_,new Call<adouble>(1,0.08,300.0,305.0,0.25,4));
    Composite<adouble> myOptions1;
    myOptions1.add(&optionAsian1);
    myOptions1.add(optionAsian2);
    myOptions1.add(optionAsian3);
    stack.new_recording();
    adouble valoracionAsian = SimpleMonteCarlo2(myOptions1,spot,sigma_map,interes,paths,samples);
    valoracionAsian.set_gradient(1.0);
    stack.compute_adjoint();
    cout<<"Valoracion Asiaticas: "<<valoracionAsian<<endl;

    cout<<"Delta: "<<spot.get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[0.5].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[1.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.0].get_gradient()<<endl;
    cout<<"Vega: "<<sigma_map[2.5].get_gradient()<<endl;
    cout<<"Rho: "<<interes.get_gradient()<<endl;

}
