#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include "Option.h"

using namespace std;

//pasar clase instrumento
//montecarlo input vector de opciones
//double Spot,//dato de mercado dejar como input
//                          double Vol,//input
//                          double r,//input
//                          unsigned long NumberOfPaths,
//                          unsigned long NumberOfSamples


double SimpleMonteCarlo2(Composite  myOptions,                     //double Expiry,//puntero clase tipo Option, añadir metedo getExpiry
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
    double dt = Expiry/NumberOfSamples;

    // vector size =  NumberOfSamples + Spot
    std::vector<double> underlying_values(NumberOfSamples + 1, Spot);

    auto variance = Vol * Vol ;
    auto rootVariance = std::sqrt(variance * dt);

    auto movedSpotFactor = std::exp((r - (0.5 * variance))*dt);

    double runningSum = 0.0;
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
        double thisPayoff{};
        //controlar para put
        //llamar al evaluate de la opcion->evaluate()

                //cout<<myOptions.evaluate(underlying_values)<<endl;
        thisPayoff=myOptions.evaluate(underlying_values);
             //thisPayoff = std::max(exerciseValue - Strike, 0.0);
        runningSum += thisPayoff;
    }

    auto mean = runningSum / NumberOfPaths;
    mean *= exp(-r*Expiry);

    return mean;
}
int main() {
   // optionType tipo, double interesAnual, double strike, double spot, double sigma, double tau

    Call option1 (0.08,300.0,305.0,0.25,(4.0/12.0));
    Put option2 (0.08,300.0,305.0,0.25,(4.0/12.0));
    Call  optionDelta(0.01,100.0,100.00001,0.5,4);
   // Asian asiatica(0.01,100.0,100.00001,0.5,4);
    vector<double> vec1 {1.0,2.0,380.0};
   // cout<<optionDelta.evaluate(vec)<<endl;

    Asian * callAsiatica = new Asian(avg_,new Call(0.08,300.0,305.0,0.25,(4.0/12.0)));
    Call *optionCall1 = new Call(0.08,300.0,305.0,0.25,(4.0/12.0));

    Call *optionCallDelta = new Call(0.01,100.0,100.00001,0.5,4);
    Put *optionPut2 = new Put(0.08,300.0,305.0,0.25,(4.0/12.0));


//std::auto_ptr<Option> myOption = *callAsiatica;
Composite myOptions;
//cout<<callAsiatica->evaluate(vec1)<<endl;
    myOptions.add(callAsiatica);
    myOptions.add(&option2);
    myOptions.add(&option1);
    myOptions.add(&optionDelta);

    double res = myOptions.evaluate(vec1);
    //cout<<myOptions.evaluate(vec1)<<endl;
    //cout<<res<<endl;
    //auto_ptr<Composite> my_options = &myOptions;
    //my_options->evaluate(vec1);

    //double resultMC = SimpleMonteCarlo2(myOptions,4,300,0.25,0.08,100000,10000);

    //cout<<resultMC<<endl;

















  /*  Asian a;
    std::shared_ptr<Option> p = std::make_shared<Call>(option1);
    std::shared_ptr<Option> p1 = std::make_shared<Put>(option2);
    std::shared_ptr<Option> p2 = std::make_shared<Asian>(a);
    Composite vec;
    vec.add(p);
    vec.add(p1);
    vec.add(p2);
    cout<<"Evaluate"<<endl;
    cout<<vec.evaluate()<<endl;

    double result = option1.price();
    cout<<"OptionDelta Delta:"<<optionDelta.delta()<<endl;
    cout<<"OptionDelta Vega: "<<optionDelta.vega()<<endl;
    cout<<"OptionDelta Theta:"<<optionDelta.theta()/365<<endl;
    std::cout << "Price: " << optionDelta.price() << std::endl;
    cout<<endl;
    cout<<"OptionVega Delta:"<<option1.delta()<<endl;
    cout<<"OptionVega Vega: "<<option1.vega()<<endl;
    cout<<"OptionVega Theta:"<<option1.theta()/365<<endl;
    std::cout << "Price: " << option1.price() << std::endl;
cout<<endl;
    cout<<"OptionTheta Delta: "<<option2.delta()<<endl;
    cout<<"OptionTheta Vega: "<<option2.vega()<<endl;
    cout<<"OptionTheta Theta: "<<option2.theta()/365<<endl;
    std::cout << "Price: " << option2.price() << std::endl;
    cout<<endl;

    double resultMC = SimpleMonteCarlo2(optionDelta.tau,optionDelta.strike,optionDelta.spot,optionDelta.sigma,optionDelta.interesAnual,100000,10000);

    std::cout << "Montecarlo:" << resultMC << std::endl;*/
    return 0;
}