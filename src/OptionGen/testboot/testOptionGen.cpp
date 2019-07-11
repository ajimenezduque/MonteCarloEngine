//
// Created by alejandro on 10/07/19.
//
#define BOOST_TEST_MODULE pruebaOptionGen

#include <boost/test/included/unit_test.hpp>
#include <boost/test/included/unit_test_framework.hpp>
#include "../OptionGen.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <chrono>


using namespace std;
//Spot, volatilidad y interes dependen del momento actual por lo que son inputs montcarlo
vector<double> SimpleMonteCarlo2(Composite<double>  myOptions,                     //double Expiry,//puntero clase tipo Option, añadir metedo getExpiry
        //double Strike,//call y put
                                 double Expiry,
                                 double Spot,//dato de mercado dejar como input
                                 double Vol,//input
                                 double r,//input
                                 unsigned long NumberOfPaths,
                                 unsigned long NumberOfSamples)
{

    /*
    // use this construction for a real random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    */

    //use this construction for a random generator with a fixed random sequence
    std::mt19937 gen(1234);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> nomalDistribution(0,1);

    // 365 muestras
    double maxExpiry = 0.0;
    for (int i=0;i<myOptions.compSize();++i){
        maxExpiry = max(maxExpiry,myOptions.getOptionGen(i)->getExpiry());
    }

    //double dt = Expiry/NumberOfSamples;
    double dt = maxExpiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<double> underlying_values(NumberOfSamples + 1, Spot);

    auto variance = Vol * Vol ;
    auto rootVariance = std::sqrt(variance * dt);

    auto movedSpotFactor = std::exp((r - (0.5 * variance))*dt);

    //double runningSum = 0.0;
    vector<double> runningSum (myOptions.compSize(),0.0);
    for (unsigned long i = 0; i < NumberOfPaths; i++)
    {
        for (unsigned long j = 1; j < underlying_values.size(); j++)
        {
            auto thisGaussian = nomalDistribution(gen);
            auto diffusion = std::exp( rootVariance*thisGaussian);
            auto drift = underlying_values[j-1] * movedSpotFactor;
            underlying_values[j] = drift * diffusion;
        }
        //pasar valor maximo, media etc cualquier payoff que se ocurra
        //añadir funcion evalue
        //auto exerciseValue = underlying_values.back();
        //double thisPayoff{};
        //llamar al evaluate de la opcion->evaluate()

        vector <double> thisPayoff (myOptions.compSize(),0.0);
        //meter en composite
        for(int i=0;i<myOptions.compSize();++i){
            //pasar criterio de evaluacion a evaluate
            thisPayoff[i]=myOptions.getOptionGen(i)->evaluate(underlying_values);
        }
        //thisPayoff=myOptions.evaluate(underlying_values);
        //
        //thisPayoff = std::max(exerciseValue - Strike, 0.0);
        for(int i=0;i<myOptions.compSize();++i){
            runningSum[i] += thisPayoff[i];
        }

    }
    //mapa fecha, valor
    vector<double> mean(myOptions.compSize(),0.0);
    for (int i=0;i<myOptions.compSize();++i){
        mean[i] = runningSum[i] / NumberOfPaths;
        mean[i] *= exp(-r * myOptions.getOptionGen(i)->getExpiry());
    }
    //auto mean = runningSum / NumberOfPaths;
    // mean *= exp(-r*Expiry);

    return mean;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");

    Call<double> opcionCallDelta (0.01,100,100,0.5,4.0);
    Put<double> opcionPutTheta (0.08,300.0,305.0,0.25,4.0/12.0);
    Call<double> opcionCallVega (0.08,300,305,0.25,4.0/12.0);


    ///Comprobacion griegas con BS//
    BOOST_TEST(0.70542 == opcionCallDelta.griegas.delta(), boost::test_tools::tolerance(0.01));
    BOOST_TEST(-0.041109 == opcionPutTheta.griegas.theta()/365, boost::test_tools::tolerance(0.01));
    BOOST_TEST(65.56772 == opcionCallVega.griegas.vega(),boost::test_tools::tolerance(0.01));


    //OptionBS callBSDelta (call,0.01,100,100,0.5,4.0);
    //OptionBS putBSTheta (put,0.08,300.0,305.0,0.25,4.0/12.0);
//    OptionBS callBSVega (call,0.08,300,305,0.25,4.0/12.0);

    Composite<double> Deltas;
    Deltas.add(&opcionCallDelta);

    Composite<double> Theta;
    Theta.add(&opcionPutTheta);

    Composite<double> Vega;
    Vega.add(&opcionCallVega);

    ///Comprobacion Pricing BS vs Montecarlo 1 vs 1 ///
  /*  BOOST_TEST(SimpleMonteCarlo2(Deltas,4,100.0,0.5,0.01,100000,365).at(0) == callBSDelta.price(), boost::test_tools::tolerance(0.01));
    BOOST_TEST(SimpleMonteCarlo2(Theta,4.0/12.0,305.0,0.25,0.08,100000,365).at(0) == putBSTheta.price(), boost::test_tools::tolerance(0.01));
    BOOST_TEST(SimpleMonteCarlo2(Vega,4.0/12.0,305.0,0.25,0.08,100000,365).at(0) == callBSVega.price(), boost::test_tools::tolerance(0.01));*/


    ///Comprobacion ejecucion en release///
    Call<double> opcionCallDelta1 (0.08,100,305,0.25,4.0);
    Put<double> opcionPutTheta1 (0.08,300.0,305.0,0.25,4.0/12.0);
    Call<double> opcionCallVega1 (0.08,300,305,0.25,4.0/12.0);

    //OptionBS callBSDelta1 (call,0.08,100,305,0.25,4.0);
    //cout<<"Delta1 "<<callBSDelta1.price()<<endl;
    //OptionBS putBSTheta1 (put,0.08,300.0,305.0,0.25,4.0/12.0);
    //cout<<"Theta1 "<<putBSTheta1.price()<<endl;
    //OptionBS callBSVega1 (call,0.08,300,305,0.25,4.0/12.0);
    //cout<<"Vega1 "<<callBSVega1.price()<<endl;

    Composite<double> myOptions;
    myOptions.add(&opcionCallDelta1);
    myOptions.add(&opcionPutTheta1);
    myOptions.add(&opcionCallVega1);

    auto start = std::chrono::high_resolution_clock::now();
    //T,spot,sigma,interes,paths,samples
    vector<double> result = SimpleMonteCarlo2(myOptions,4.0,305,0.25,0.08,100000,365);
    auto end = std::chrono::high_resolution_clock::now();

    std::for_each(result.begin(), result.end(),   [](const double& i) { std::cout << "Result: " << i<<endl;});
    std::chrono::duration<double> diff=end-start;
    cout << "Tiempo de cálculo: "<<diff.count()<<" segundos"<< endl;


    ///////Griegas de forma Numerica////
    double deltaPrice = 0.001;

    //OptionBS callBSDeltaNum (call,0.01,100,100,0.5,4.0);
   // OptionBS callBSDeltaNum1 (call,0.01,100,100 + deltaPrice,0.5,4.0);

   //cout<<"BSDelta "<< opcionCallDelta.griegas.delta()<<endl;
    //cout<<"Numeric Delta "<<(callBSDeltaNum1.price() - callBSDeltaNum.price()) / deltaPrice<<endl ;


    //OptionBS callBSVegaNum (call,0.08,300,305,0.25,4.0/12.0);
    //OptionBS callBSVegaNum1 (call,0.08,300,305,0.25 + deltaPrice,4.0/12.0);

    //cout<<"BSVega "<< opcionCallVega.griegas.vega()<<endl;
    //cout<<"Numeric Vega "<<(callBSVegaNum1.price() - callBSVegaNum.price()) / deltaPrice<<endl ;

    //OptionBS putBSThetaNum (put,0.08,300.0,305.0,0.25,4.0/12.0);
    //OptionBS putBSThetaNum1  (put,0.08,300.0,305.0,0.25,(4.0/12.0)+deltaPrice);

    //cout<<"BSTheta "<< opcionPutTheta.griegas.theta()/365<<endl;
    //cout<<"Numeric Theta "<<-((putBSThetaNum1.price() - putBSThetaNum.price()) / deltaPrice)/365<<endl ;



    /*Put opcionPut (0.08,300.0,305.0,0.25,4);
    OptionBS putBS (put,0.08,300.0,305.0,0.25,4);

    Asian *optionAsian1 = new Asian(avg_,new Call(0.08,300.0,305.0,0.25,4));
    Asian *optionAsian2 = new Asian(max_,new Call(0.08,300.0,305.0,0.25,4));
    Asian *optionAsian3 = new Asian(min_,new Call(0.08,300.0,305.0,0.25,4));*/

}
