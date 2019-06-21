#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include "Option.h"

using namespace std;
/*adouble algorithm_tfm(const adouble x[5]) {

    // adouble z = x[0]*x[1] + sin(x[0]);
    adouble z = log(x[0]) +  x[0]*x[1] - sin(x[1]);
    return z;
}*/
//pasar clase instrumento
double SimpleMonteCarlo2( double Expiry,//puntero clase tipo instrumento
                          double Strike,
                          double Spot,
                          double Vol,
                          double r,
                          unsigned long NumberOfPaths,
                          unsigned long NumberOfSamples)
{
    // create a mersenne_twister generator

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
        auto exerciseValue = underlying_values.back();
        auto thisPayoff = std::max(exerciseValue - Strike, 0.0);
        runningSum += thisPayoff;
    }

    auto mean = runningSum / NumberOfPaths;

    mean *= exp(-r*Expiry);

    return mean;
}
int main() {
   // optionType tipo, double interesAnual, double strike, double spot, double sigma, double tau


    Option option1 (call,0.08,300.0,305.0,0.25,(4.0/12.0));
    Option option2 (put,0.08,300.0,305.0,0.25,(4.0/12.0));
    Option  optionDelta(call,0.01,100.0,100.00001,0.5,4);
    /*Asian a;
    std::shared_ptr<Instrumento> p = std::make_shared<Option>(option1);
    std::shared_ptr<Instrumento> p1 = std::make_shared<Option>(option1);
    std::shared_ptr<Instrumento> p2 = std::make_shared<Asian>(a);
    Composite vec;
    vec.add(p);
    vec.add(p1);
    vec.add(p2);
    cout<<"Evaluate"<<endl;
    cout<<vec.evaluate()<<endl;*/

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
  //  double resultMC = SimpleMonteCarlo2(option1.tau,option1.strike,option1.spot,option1.sigma,option1.interesAnual,100000,10000);
    std::cout << "Price: " << result << std::endl;
  //  std::cout << "Montecarlo:" << resultMC << std::endl;
    return 0;
}