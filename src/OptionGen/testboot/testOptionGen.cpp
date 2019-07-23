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
            //runningSum[it->first] = it->second +
            auto f = runningSum.find(it->first);
            if (f != runningSum.end()) {
              //  cout<<"Entra autof : First: " <<it->first<<" second : "<<it->second<<" Fsecond : "<<f->second<<endl;
                double suma{};
                suma = it->second + f->second;
                runningSum.erase(f);
                runningSum.insert(make_pair(it->first, suma));
            } else {
               // cout<<"Entra no tiene valores: "<< it->first<<" Second : "<<it->second<<endl;
                runningSum.insert(make_pair(it->first, it->second));
            }

           // std::cout << it->first << " => " << it->second << '\n';
        }
    }

    double value{};
    for (auto it = runningSum.begin(); it != runningSum.end(); ++it) {
       // cout<< it->first << " => " << it->second << '\n';
        value += (exp(-r * it->first) * (it->second / NumberOfPaths));
    }

    return value;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen){
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");

    double interes = 0.08;
    double spot = 305.0;
    double sigma = 0.25;
    unsigned long paths = 100000;
    unsigned long samples = 12;
    double deltaPrice = 0.00001;
    Call<double> opcionCallDelta(1,0.01, 100, 100, 0.5, 4.0);

    Put<double> opcionPutTheta(-1,0.08, 300.0, 305.0, 0.25, 4.0 / 12.0);

    Call<double> opcionCallVega(1,0.08, 300, 305, 0.25, 4.0 / 12.0);

   // Call<double> opcionCallVegaNew(1,0.08, 300, 305, 0.25, 4.0 / 12.0);
    double spotM = SimpleMonteCarlo2(opcionCallVega, spot + deltaPrice,sigma,interes, paths, samples);
    double base = SimpleMonteCarlo2(opcionCallVega, spot,sigma,interes, paths, samples);
    double sigmaM = SimpleMonteCarlo2(opcionCallVega,spot,sigma + deltaPrice,interes,paths,samples);

    cout<<"Valoracion SpotM: "<<spotM<<endl;
    cout<<"Valoracion Spot1: "<<base<<endl;

    cout<<"Dif Spot: "<<((spotM-base)/deltaPrice)<<endl;
    cout<<"Dif Sigma: "<<((sigmaM-base)/deltaPrice)<<endl;

    cout<<"Griegas: "<<endl;
    cout<<"Delta: "<<opcionCallVega.griegas.delta()<<endl;
    cout<<"Theta: "<<opcionCallVega.griegas.theta()<<endl;
    cout<<"Vega: "<<opcionCallVega.griegas.vega()<<endl;
    cout<<"Fin griegas.--"<<endl;



    Composite<double> Deltas;
    Deltas.add(&opcionCallDelta);

    Composite<double> Theta;
    Theta.add(&opcionPutTheta);

    Composite<double> Vega;
    Vega.add(&opcionCallVega);

    double valor = SimpleMonteCarlo2(opcionCallDelta, 305.0, 0.25, 0.08, 100000, 12);
    cout<<"Delta "<<valor<<endl;
    valor = SimpleMonteCarlo2(Theta, 305.0, 0.25, 0.08, 100000, 12);
    cout<<"Theta "<<valor<<endl;
    valor = SimpleMonteCarlo2(Vega, 305.0, 0.25, 0.08, 100000, 12);
    cout<<"Vega "<<valor<<endl;

    ///Comprobacion ejecucion en release///
    Call<double> opcionCallDelta1 (1,0.08,100,305,0.25,4.0);
    Put<double> opcionPutTheta1 (-1,0.08,300.0,305.0,0.25,4.0/12.0);
    Call<double> opcionCallVega1 (-1,0.08,300,305,0.25,4.0/12.0);


    Composite<double> myOptions;
    myOptions.add(&opcionCallDelta1);
    myOptions.add(&opcionPutTheta1);
    myOptions.add(&opcionCallVega1);

    double result = SimpleMonteCarlo2(myOptions,305,0.25,0.08,100000,12);
    cout<<"Portfolio: "<<result<<endl;




















    /*for (int i = 0; i < result.size();i++){
        double sum{};
        double price{};
        // cout<<"Fecha: "<< get<0>(result[i])<<" Vector: ";
        for (int j =0; j < get<1>(result[i]).size();++j){
            sum += get<1>(result[i]).at(j);
            //cout<<get<1>(valor[i]).at(j)<<" ";
        }
        sum = sum / 100000;
        price = sum * exp(-0.08 * get<0>(result[i]));
        cout<<"Price: "<<i<<" "<<price<<endl;
    }*/













   /* map <double,vector<double>> result = SimpleMonteCarlo2(myOptions,305,0.25,0.08,100000,12);
    for (int i = 0; i< myOptions.compSize();i++) {
        double tau = myOptions.getOptionGen(i)->getExpiry();
        auto it = result.find(tau);
        double price{};
        if (it != result.end()) {
            double sum{};
            for (int i = 0; i < it->second.size(); ++i) {
                sum += it->second.at(i);
            }
            sum = sum / 100000;
            price = sum * exp(-0.08 * tau);
        }
        cout<<"First "<<tau<<" Price: "<<price<<endl;
    }*/

   /* for (auto &x : result){
        cout<<"First: "<< x.first<<" ";
        double sum{};
        double price{};
        double tau = myOptions.getOptionGen(it)->getExpiry();
        double exponencial = exp(-0.08 * tau);;
        if (tau == x.first){
            for(int i = 0;i <x.second.size();++i){
                sum += x.second.at(i);
            }
            sum = sum/100000;
            price = sum * exponencial;
            cout<<"Price: "<< price <<endl;
        }
        it++;
    }*/
  /*  auto start = std::chrono::high_resolution_clock::now();
    //T,spot,sigma,interes,paths,samples
   // vector<double> result = SimpleMonteCarlo2(myOptions,4.0,305,0.25,0.08,100000,365);

    map <double,vector<double>> result = SimpleMonteCarlo2(myOptions,4.0,305,0.25,0.08,100000,365);
    //vector <tuple<double,double>> result = SimpleMonteCarlo2(myOptions,305,0.25,0.08,100000,365);
    auto end = std::chrono::high_resolution_clock::now();

    //std::for_each(result.begin(), result.end(),   [](map<double, double>::iterator it) { std::cout << "index: " << it->first<< " Key: "<<it->second<<endl;});
    for (const auto& it : result) {
        std::cout << "key: "<<it.first << " value: " << it.second << std::endl;
    }*/
    /*for (const auto& it : result) {
        std::cout << "key: "<<get<0>(it) << " value: " << get<1>(it) << std::endl;
    }*/
   // std::chrono::duration<double> diff=end-start;
    //cout << "Tiempo de cálculo: "<<diff.count()<<" segundos"<< endl;


    ///////Griegas de forma Numerica////
   // double deltaPrice = 0.001;

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
