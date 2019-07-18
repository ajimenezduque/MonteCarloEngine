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
vector<tuple<double,vector<double>>> SimpleMonteCarlo2
(Composite<double>  myOptions,                     //double Expiry,//puntero clase tipo Option, añadir metedo getExpiry
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

    double maxExpiry = myOptions.getExpiry();

    double dt = maxExpiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<double> underlying_values(NumberOfSamples + 1, Spot);

    auto variance = Vol * Vol ;
    auto rootVariance = std::sqrt(variance * dt);

    auto movedSpotFactor = std::exp((r - (0.5 * variance))*dt);

    //vector<double> runningSum (myOptions.compSize(),0.0);
   // vector <tuple <double,vector<double> >> runningSum{};//(myOptions.compSize(),0.0);
    /*for (int i = 0;i<myOptions.compSize();i++){
        runningSum.emplace_back(make_pair(0.0,0.0));
    }*/
    vector <tuple <double,vector<double> >> runningSum{};
    vector<vector<double>> flujos{};
    for (int i = 0;i<myOptions.compSize();i++){
        vector<double> aux{};
        flujos.emplace_back(aux);
    }
    vector<double> fechas (myOptions.compSize(),0.0);
    for (unsigned long i = 0; i < NumberOfPaths; i++)
    {
        for (unsigned long j = 1; j < underlying_values.size(); j++)
        {
            auto thisGaussian = nomalDistribution(gen);
            auto diffusion = std::exp( rootVariance*thisGaussian);
            auto drift = underlying_values[j-1] * movedSpotFactor;
            underlying_values[j] = drift * diffusion;
        }

        //vector <double> thisPayoff (myOptions.compSize(),0.0);
        vector <tuple <double,double> > thisPayoff{};//(myOptions.compSize(),0.0);
        /*for (int i = 0;i<myOptions.compSize();i++){
            thisPayoff.emplace_back(make_pair(0.0,0.0));
        }*/


        thisPayoff = myOptions.evaluate(underlying_values, NumberOfSamples/maxExpiry);
        for (int i = 0; i<myOptions.compSize();i++){
            //guardar cada uno de las += thispayoff en vector y sumarlos despues...
                //get<1>(runningSum[i]) += get<1>(thisPayoff[i]);
                //get<0>(runningSum[i]) = get<0>(thisPayoff[i]);
                ///PRUEBA NUEVA VECTOR
                flujos[i].emplace_back(get<1>(thisPayoff[i]));
                fechas[i] = get<0>(thisPayoff[i]);
               // runningSum.insert(pair<double,vector<double>>(get<0>(thisPayoff[i]),flujos[i]));
        }
    }
    for (int i = 0;i<myOptions.compSize();++i){
        runningSum.emplace_back(make_tuple(fechas[i],flujos[i]));
    }

   /* for(int i =0;i<flujos.size();i++){
        cout<<"Vector: "<<i<<" Size: "<<flujos[i].size()<<endl;
        for (int j = 0; j<flujos[i].size();j++){
            cout<<" "<<flujos[i].at(j)<<endl;
        }
    }*/

    //map<double,vector<double>> mean{};
   // vector<double> mean(myOptions.compSize(),0.0);
    //vector<tuple<double,double>> mean2{};

    // double fecha{};
    //for (int i=0;i<myOptions.compSize();++i){
      //  fecha = get<0>(runningSum[i]);
        //flujos.emplace_back(get<1>(runningSum[i])/NumberOfPaths);
        //double media = get<1>(runningSum[i])/NumberOfPaths ;
       // mean.insert(pair<double,double>(get<0>(runningSum[i]),media * exp(-r * myOptions.getOptionGen(i)->getExpiry())));
      // mean2.emplace_back(make_pair(get<0>(runningSum[i]),media * exp(-r * myOptions.getOptionGen(i)->getExpiry())));
    //}
   // mean.insert(pair<double,vector<double>>(fecha,flujos));
  /* cout<<"Size of runningSum "<<runningSum.size()<<endl;
    for(int i =0;i<flujos.size();i++){
        cout<<"Vector: "<<i<<" Size: "<<flujos[i].size()<<endl;
        for (int j = 0; j<flujos[i].size();j++){
            cout<<" "<<flujos[i].at(j)<<endl;
        }
    }*/
    return runningSum;
}

BOOST_AUTO_TEST_CASE(Test_OptionGen) {
    BOOST_TEST_MESSAGE("Se ejecuta test opcion Put y call, comprobando pricing y griegas con BS y con MC");

    double interes = 0.08;
    double spot = 305.0;
    double sigma = 0.25;
    double paths = 100000;
    double samples = 12;
    Call<double> opcionCallDelta(0.01, 100, 100, 0.5, 4.0);
    Put<double> opcionPutTheta(0.08, 300.0, 305.0, 0.25, 4.0 / 12.0);
    Call<double> opcionCallVega(0.08, 300, 305, 0.25, 4.0 / 12.0);

    ///Comprobacion griegas con BS//
    BOOST_TEST(0.70542 == opcionCallDelta.griegas.delta(), boost::test_tools::tolerance(0.01));
    BOOST_TEST(-0.041109 == opcionPutTheta.griegas.theta() / 365, boost::test_tools::tolerance(0.01));
    BOOST_TEST(65.56772 == opcionCallVega.griegas.vega(), boost::test_tools::tolerance(0.01));


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

    /* cout<<"Deltas" <<SimpleMonteCarlo2(Deltas,4,305.0,0.25,0.08,100000,12).at(0)<<endl;
       cout<<"Theta" <<SimpleMonteCarlo2(Theta,4.0/12.0,305.0,0.25,0.08,100000,12).at(0)<<endl;
       cout<<"Vega"<< SimpleMonteCarlo2(Vega,4.0/12.0,305.0,0.25,0.08,100000,12).at(0)<<endl;*/

    //map<double,vector<double>> valor = SimpleMonteCarlo2(Deltas,305.0,0.25,0.08,100000,12);
    vector<tuple<double, vector<double>>> valor = SimpleMonteCarlo2(Deltas, 305.0, 0.25, 0.08, 100000, 12);
    /*   cout<<"Deltas" <<endl;
       for (auto &x : valor){
           cout<<"First: "<< x.first<<" ";
           double sum{};
           double price{};
           for(int i = 0;i <x.second.size();++i){
               sum += x.second.at(i);
           }
           sum = sum/100000;
           price = sum * exp(-0.08 * opcionCallDelta.getExpiry());
           cout<<"Price Delta "<< price <<endl;
       }
       valor = SimpleMonteCarlo2(Theta,305.0,0.25,0.08,100000,12);
       cout<<"Theta" <<endl;
       for (auto &x : valor){
           cout<<"First: "<< x.first<<" ";
           double sum{};
           double price{};
           for(int i = 0;i <x.second.size();++i){
               sum += x.second.at(i);
           }
           sum = sum/100000;
           price = sum * exp(-0.08 * opcionPutTheta.getExpiry());
           cout<<"Price Theta "<< price <<endl;
       }
       valor = SimpleMonteCarlo2(Vega,305.0,0.25,0.08,100000,12);
       cout<<"Vega"<<endl;
       for (auto &x : valor){
           cout<<"First: "<< x.first<<" ";
           double sum{};
           double price{};
           for(int i = 0;i <x.second.size();++i){
               sum += x.second.at(i);
           }
           sum = sum/100000;
           price = sum * exp(-0.08 * opcionCallVega.getExpiry());
           cout<<"Price Vega "<< price <<endl;
       }*/

   // vector<tuple<double, vector<double>>> valor = SimpleMonteCarlo2(Deltas, 305.0, 0.25, 0.08, 100000, 365);
    cout << "Deltas" << endl;
    for (int i = 0; i < valor.size();i++){
        double sum{};
        double price{};
       // cout<<"Fecha: "<< get<0>(valor[i])<<" Vector: ";
        for (int j =0; j < get<1>(valor[i]).size();++j){
            sum += get<1>(valor[i]).at(j);
            //cout<<get<1>(valor[i]).at(j)<<" ";
        }
        sum = sum / 100000;
        price = sum * exp(-0.08 * get<0>(valor[i]));
        cout<<"Price Delta: "<<price<<endl;
    }
    valor = SimpleMonteCarlo2(Theta, 305.0, 0.25, 0.08, 100000, 12);
    cout << "Theta" << endl;
    for (int i = 0; i < valor.size();i++){
        double sum{};
        double price{};
        // cout<<"Fecha: "<< get<0>(valor[i])<<" Vector: ";
        for (int j =0; j < get<1>(valor[i]).size();++j){
            sum += get<1>(valor[i]).at(j);
            //cout<<get<1>(valor[i]).at(j)<<" ";
        }
        sum = sum / 100000;
        price = sum * exp(-0.08 * get<0>(valor[i]));
        cout<<"Price Theta: "<<price<<endl;
    }
    valor = SimpleMonteCarlo2(Vega, 305.0, 0.25, 0.08, 100000, 12);
    cout << "Vega" << endl;
    for (int i = 0; i < valor.size();i++){
        double sum{};
        double price{};
        // cout<<"Fecha: "<< get<0>(valor[i])<<" Vector: ";
        for (int j =0; j < get<1>(valor[i]).size();++j){
            sum += get<1>(valor[i]).at(j);
            //cout<<get<1>(valor[i]).at(j)<<" ";
        }
        sum = sum / 100000;
        price = sum * exp(-0.08 * get<0>(valor[i]));
        cout<<"Price Vega: "<<price<<endl;
    }


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
    vector<tuple <double,vector<double>>> result = SimpleMonteCarlo2(myOptions,305,0.25,0.08,100000,12);

    for (int i = 0; i < result.size();i++){
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
    }
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
