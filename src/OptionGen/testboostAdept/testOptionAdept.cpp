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
//Spot, volatilidad y interes dependen del momento actual por lo que son inputs montcarlo
//double
adouble SimpleMonteCarlo2(
        OptionGen<adouble> &option,
        adouble Spot,//dato de mercado dejar como input
        adouble Vol,//input
        adouble r,//input
       // const adouble x[3],
        unsigned long NumberOfPaths,
        unsigned long NumberOfSamples)
{


    //use this construction for a random generator with a fixed random sequence
    std::mt19937 gen(1234);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> nomalDistribution(0,1);

    double maxExpiry = option.getExpiry();

    adouble dt = maxExpiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<adouble> underlying_values(NumberOfSamples + 1,Spot); //x[0]);

    adouble variance = Vol * Vol;//x[1] * x[1] ;
    adouble rootVariance = sqrt(variance * dt);

    adouble movedSpotFactor = exp((r - (0.5 * variance))*dt);

    map<adouble,adouble> runningSum{};

    for (unsigned long i = 0; i < NumberOfPaths; i++) {
        for (unsigned long j = 1; j < underlying_values.size(); j++) {
            adouble thisGaussian = nomalDistribution(gen);
            adouble diffusion = exp(rootVariance * thisGaussian);
            adouble drift = underlying_values[j - 1] * movedSpotFactor;
            underlying_values[j] = drift * diffusion;
        }

        map<adouble, adouble> thisPayoff{};

        thisPayoff = option.evaluate(underlying_values, NumberOfSamples/maxExpiry);

        for (auto it = thisPayoff.begin(); it != thisPayoff.end(); ++it) {
            auto f = runningSum.find(it->first);
            if (f != runningSum.end()) {
                //  cout<<"Entra autof : First: " <<it->first<<" second : "<<it->second<<" Fsecond : "<<f->second<<endl;
                adouble suma{};
                suma = it->second + f->second;
                runningSum.erase(f);
                runningSum.insert(make_pair(it->first, suma));
            } else {
                // cout<<"Entra no tiene valores: "<< it->first<<" Second : "<<it->second<<endl;
                runningSum.insert(make_pair(it->first, it->second));
            }

        }
    }

    adouble value{};
    for (auto it = runningSum.begin(); it != runningSum.end(); ++it) {
        // cout<< it->first << " => " << it->second << '\n';
        value += (exp(-r * it->first) * (it->second / NumberOfPaths));
    }

    return value;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");
   cout<<"Prueba"<<endl;
   Stack stack;
   adouble x[3];
    x[0] = 305.0;
    x[1] = 0.25;
    x[2] = 0.08;

    adouble interes = 0.08;
    adouble spot = 305.0;
    adouble sigma = 0.25;
    unsigned long paths = 100000;
    unsigned long samples = 12;

    Call<adouble> opcionCallDelta(0.01, 100, 100, 0.5, 4.0);
    Put<adouble> opcionPutTheta(0.08, 300.0, 305.0, 0.25, 4.0 / 12.0);
    Call<adouble> opcionCallVega(0.08, 300, 305, 0.25, 4.0 / 12.0);

    cout<<"Griegas: "<<endl;
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Theta: "<<opcionCallVega.griegas.theta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;
//
    Composite<adouble> Deltas;
    Deltas.add(&opcionCallDelta);

    Composite<adouble> Theta;
    Theta.add(&opcionPutTheta);

    Composite<adouble> Vega;
    Vega.add(&opcionCallVega);

    adouble valor;
    stack.new_recording();
    valor = SimpleMonteCarlo2(opcionCallVega, spot,sigma,interes, paths, samples);
    valor.set_gradient(1.0);
    stack.compute_adjoint();

    stack.reverse();
    std::cout << "Final list of gradients:\n";
    stack.print_gradients();
    cout<<endl;
    //for (int i=0;i<3;++i){
        cout<<"Gradiente "<<0<<": "<<spot.get_gradient()/100000<<endl;
        cout<<"Gradiente "<<1<<": "<<sigma.get_gradient()/100000<<endl;
        cout<<"Gradiente "<<2<<": "<<interes.get_gradient()/100000<<endl;
   // }
    cout<<"Vega "<<valor<<endl;




   /* valor = SimpleMonteCarlo2(Theta, 305.0, 0.25, 0.08, 100000, 12);
    cout<<"Theta "<<valor<<endl;
    valor = SimpleMonteCarlo2(Vega, 305.0, 0.25, 0.08, 100000, 12);
    cout<<"Vega "<<valor<<endl;*/


    ///Comprobacion Composite en release///
    Call<adouble> opcionCallDelta1 (0.08,100,305,0.25,4.0);
    Put<adouble> opcionPutTheta1 (0.08,300.0,305.0,0.25,4.0/12.0);
    Call<adouble> opcionCallVega1 (0.08,300,305,0.25,4.0/12.0);


    Composite<adouble> myOptions;
    myOptions.add(&opcionCallDelta1);
    myOptions.add(&opcionPutTheta1);
    myOptions.add(&opcionCallVega1);

    adouble result = SimpleMonteCarlo2(myOptions,spot,sigma,interes,100000,12);
    cout<<"Portfolio: "<<result<<endl;

}
