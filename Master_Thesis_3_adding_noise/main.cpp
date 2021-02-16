#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>

typedef boost::numeric::ublas::vector<double> vectord;
typedef boost::numeric::ublas::matrix<double> matrixd;
using namespace boost::numeric::odeint;
using namespace std;

#include "system.h"

int main(int argc, char **argv)
{
    /*
    string file = string(argv[1]);

    boost::property_tree::ptree desc;
    boost::property_tree::json_parser::read_json(file, desc);
    
    const int   n         = desc.get<int>("number of oscillators");
    const int   r_one     = desc.get<int>("range");
    float       om_loc    = desc.get<float>("omega");
    double      g         = desc.get<double>("g");
    const float a         = desc.get<float>("alpha");
    const int   layers    = desc.get<int>("layers");
    const double sim_time = desc.get<float>("simulation time");
    string      path      = desc.get<string>("path");
    string      model     = desc.get<string>("model");
    double      Dmax      = desc.get<double>("Dmax");
    double      Dstep     = desc.get<double>("Dstep");
    float       alpha_    = desc.get<float>("PSalpha");
    
    cout << "n: "         << n        << '\n' << "r: "        << r_one  << '\n' <<
            "omega loc: " << om_loc   << '\n' << "g: "        << g      << '\n' <<
            "Alpha: "     << a        << '\n' << "layers: "   << layers << '\n' <<
            "Maximum D: " << Dmax     << '\n' << "D step: "   << Dstep << '\n' <<
            "Sim time: "  << sim_time << '\n' << "path: "     << path << '\n';
     */
     
    clock_t tStart = clock();
    
    constexpr double sim_time   = 100;
    const size_t n = 256;
    int r_one = 35;
    float om_loc = 0.5;
    double g = 0.0;
    float a = 1.45;
    int layers = 2;
    string path = string("/Users/iwanphillips/C++/Master_Thesis_3_adding_noise/RK2_test_saturday") + "" + string(".txt");
    string model = "standard";
    float alpha_ = 0; // 1.95
    double Dmax = 0.05; double Dstep = 0.005;
    //double gamma = 0.4;
    
    float l_one = 0.085;
    float l_two = 0.01;
    
    vectord omega ( 2*n , 0.0 );
    
    matrixd x( n , 2 , 0. );
    
    boost::mt19937 rng;
    boost::cauchy_distribution<> cauchy( om_loc , g );
    boost::variate_generator< boost::mt19937&, boost::cauchy_distribution<> > gen( rng , cauchy );
    generate( omega.begin() , omega.end() , gen );
     
    static uniform_real_distribution<double> u_dist(-0.5,0.5);
    random_device rd;
    default_random_engine generator(rd());
    
    //
        int steps = (float) sim_time / dt;
        int n_pts = 16384;
    //
        while(steps > n_pts){n_pts = 2*n_pts;}
        float ha = 1;
        float Q = 10; // ha / (2*pow(6.283185, alpha_));
        float Q_d = 1; // Q / pow(dt, 1-alpha_);
        
        int c = 5;
        int max = c*n_pts;
        while(max < 2000000){max += steps; c += 1;}
        vector<float> X; //float X[max];
        for( int i=0 ; i<c ; ++i ){ //for( int i=0 ; i<50 ; ++i ){ //i<n
            float X_i[n_pts];
            memset(X_i, 0.0, sizeof X_i);
            long utime; utime=(long)time(NULL); long seed=utime;
            f_alpha(n_pts, X_i, Q_d, alpha_, &seed);
            //memcpy( X + i*n_pts, X_i, sizeof( X_i ) );
            vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
            X.insert(X.end(), Xi.begin(), Xi.end());
        }
    /*
       int c = 10; int max = c*steps;
       while(max < 2000000){max += steps; c += 1;}
       vector<float> X;
       for( int i=0 ; i<c ; ++i ){
           float X_i[steps];
           memset(X_i, 0.0, sizeof X_i);
           long utime; utime=(long)time(NULL); long seed=utime;
           cor_exp(steps, X_i, 0.01, &seed, dt, gamma);
           vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
           X.insert(X.end(), Xi.begin(), Xi.end());
       }
    //
        int c = 10; int max = c*steps;
        while(max < 800000){max += steps; c += 1;}
        float X[max];
        for( int i=0 ; i<c ; ++i ){
            float X_i[steps];
            memset(X_i, 0.0, sizeof X_i);
            exp_noise_second(0.4, steps, 1, X_i, dt);
            memcpy( X + i*steps, X_i, sizeof( X_i ) );
        }
    */
    
    string output_data = path;
    ofstream data_out(output_data);
    
    cout << "Starting integration..." << endl;
    
    for(size_t i = 0 ; i < n ; ++i )
    {
        double pos = i*2*M_PI/(n-1) - M_PI;
        double r1 = u_dist(generator); double r2 = u_dist(generator);
        x(i,0) = 6*r1*exp(-0.76*pos*pos);  if(layers == 2){x(i,1) = 6*r2*exp(-0.76*pos*pos);}
    }
    
    //double D = 0;
    //system ( sim_time , n , x , model , layers , omega , a , l_one , l_two , r_one , D , data_out , alpha_ );

    if(layers == 1){
        for( double D  = 0.0 ; D < Dmax ; D += Dstep ){
            for( double l_one = 0.0 ; l_one < 0.1 ; l_one += 0.01 ){
                system ( sim_time , n , x , model , layers , omega , a , l_one , l_two , r_one , D , data_out , alpha_ , alpha_ , X , n_pts , max );
            }
        }
    }
        
        //boost::mt19937 rng;
        
    if(layers == 2){
        for( double l_two = 0.0001 ; l_two < 0.041 ; l_two += 0.004 ){ // double gamma = 0.01 ; gamma < 3.1 ; gamma += 0.3
            for( double D  = 0.0 ; D < (0.1*l_two)+0.0001 ; D += (0.1*l_two)/10 ){
                
                /*
                cout << "c  = " << c << '\n';
                for( int i=0 ; i<c ; ++i ){
                    float X_i[steps];
                    memset(X_i, 0.0, sizeof X_i);
                    long utime; utime=(long)time(NULL); long seed=utime;
                    cor_exp(steps, X_i, 0.01, &seed, dt, gamma); // Q was 0.1
                    //for (int h = 0; h < steps; ++h){
                    //    cout << "h = " << h << '\t' << "X_i test = " << X_i[h] << '\n';
                    //}
                    vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
                    X.insert(X.end(), Xi.begin(), Xi.end());
                    Xi.clear(); // not necassary ?
                 } */
                
                //boost::cauchy_distribution<> cauchy( om_loc , g );
                //boost::variate_generator< boost::mt19937&, boost::cauchy_distribution<> > gen( rng , cauchy );
                //generate( omega.begin() , omega.end() , gen );
                
                system ( sim_time , n , x , model , layers , omega , a , l_one , l_two , r_one , D , data_out , alpha_ , alpha_ , X , n_pts , max );
                
                //X.clear();
            }
        }
    }
    
    data_out.close();
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}



