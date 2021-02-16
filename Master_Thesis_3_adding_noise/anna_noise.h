#ifndef anna_noise_h
#define anna_noise_h

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <vector>

using namespace std;

static normal_distribution<double> n_dist(0 , sqrt(dt));
static random_device rd;
static default_random_engine generator(rd());
 
static double dW ( double dt ) {
    return n_dist(generator);
}

void exp_noise_first(double gamma, int N, double c_sigma, float X[])
{
    double y_init = 1;
    double sum = 0;
    X[0] = y_init;
    
    for (int i = 1; i <= N; i++ ){
        double t = i * dt;
        sum += 15 * exp (gamma*t) * dW(dt);
        X[i] = (sqrt(dt) * 3) * exp (-gamma*t)*(y_init + sqrt(2*gamma)*sqrt(dt)*sum);
        }
    X[0] *= (sqrt(dt) * 3); // so that it is of similar magnitude to the coloured noise
}

void exp_noise_second(double gamma, int N, double c_sigma, float X[])
{
    double y_init = 1;
    double sum = 0;
    X[0] = y_init;
    
    double alpha = gamma/(1 - (0.5 * gamma));
    double sigma = sqrt(dt)/(1 - (0.5 * gamma));
    
    for (int i = 1; i <= N; i++ ){
        double t = i * dt;
        sum += 10 * exp (alpha*t) * dW(dt);
        X[i] = (sqrt(dt) * 3) * exp (-alpha*t)*(y_init + sigma*sum);
        }
    X[0] *= (sqrt(dt) * 3);
}

void exp_noise_EM(double gamma, int N, double c_sigma, float X[])
{
    double y = 1;
        
    for (int i = 0; i < N; i++ ){
        y += - gamma * y * dt + sqrt(2*gamma) * sqrt(dt) * dW(dt);
        X[i] = (sqrt(dt) * 3) * y;
        }
}

void exp_noise_anna(double gamma, int N, double c_sigma, float X[])
{
    X[0] = 1;
        
    for (int i = 0; i < N; i++ ){
        X[i+1] = (1 - gamma) * X[i] + 5 * sqrt(2*gamma) * dW(1);
        }
    for (int i = 0; i < N; i++ ){X[i] *= (sqrt(dt) * 3);}
}

#endif /* anna_noise_h */


/*
 
 int main()
 {
    int n_pts = 100000;
    int steps = n_pts * dt;
    double gamma = 0.4;
    double D = 15;
     
    cout << "steps = " << '\t' << steps << '\n';
     
    float X1[n_pts];
    memset(X1, 0.0, sizeof X1);
     
    float W[n_pts];
    for (int i = 0; i < n_pts; ++i){
        W[i] = dW(dt);
    }
     
     float W2[steps];
     for (int i = 0; i < steps; ++i){
         W2[i] = dW(1);
     }
     
     exp_noise_first(gamma, n_pts, D, X1, W);
     exp_noise_second(gamma, n_pts, D, X2, W);
     exp_noise_EM(gamma, n_pts, D, X3, W);
     exp_noise_anna(gamma, steps, D, X4, W2);
     
     ofstream myfile;
     myfile.open ("/Users/iwanphillips/C++/Master_Thesis_3_adding_noise/anna_noise.txt");
     
     for (int i = 0; i < n_pts; ++i){
         myfile << X1[i] << '\n';
     }
     
     myfile.close();
     
     return 0;
 }
*/
