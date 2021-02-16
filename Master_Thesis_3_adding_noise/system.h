#ifndef system_h
#define system_h

#include "noise.h"
#include "int_method.h"
#include "anna_noise.h"

// maybe I should use one function fv which produces a pair of (x first layer, x second layer) rather than creating two seperate functions
double fv1 ( double t, double x, matrixd xm, int i, int n, int r_one, double lamb1, double alpha, int layers, double lamb2, vectord omega )
{
    double loc_coup = 0;
    
    for( int k=0; k<n ; ++k ){
        float dist = abs(k-i);
        dist = abs(dist - round(dist/( (float) n ) ) * n );
        if(dist <= r_one && dist > 0){
            loc_coup += sin( x - xm(k,0) + alpha );
        }
    }
    
    double x2 =  omega[i] - (lamb1/(2*r_one + 1))*loc_coup;
    if(layers == 2){ x2 += (lamb2/2)*sin( xm(i,1) - x ); } // sin( x - xm(i,1) )
    
    return x2;
}

double fv2 ( double t, double x, matrixd xm, int i, int n, int r_one, double lamb1, double alpha, int layers, double lamb2, vectord omega )
{
    double loc_coup = 0.;
    
    for( int k=0; k<n ; ++k ){
        float dist = abs(k-i);
        dist = abs(dist - round(dist/( (float) n ) ) * n );
        if(dist <= r_one && dist > 0){
            loc_coup += sin( x - xm(k,1) + alpha );
        }
    }
    
    double x2 = omega[i+n] - (lamb1/(2*r_one + 1))*loc_coup + (lamb2/2)*sin( xm(i,0) - x ); // - not +
    
    return x2;
}


double gv ( double t, double x, double D, matrixd xm, int i, double xorig)
{
    double x2;
    if(xorig == xm(i,0)){x2 = D*sin( xm(i,1) - x );} // x = xm(i,0)
    else{x2 = D*sin( xm(i,0) - x );}             // x = xm(i,1)
    return x2;
}


double system ( int sim_time , int n , matrixd x , string model , int layers , vectord omega , double alpha ,
                double lamb1 , double lamb2 , int r_one , double c_sigma , ofstream &data_out, float alpha_, double gamma,  vector<float> X, int n_pts, int max ) // float X[]
{
    int steps = (float) sim_time / dt;
    
    double delta2 = 0.0;
    size_t count = 0;
    
    // matrixd loc_coup(n , 2, 0.) ;
    // matrixd int_coup(n , 2, 0.) ;
    
    if(model == "standard"){
        
        cout << "size of array: " << max << '\n';
        for (int h = 0; h < steps; ++h){
            cout << "h = " << h << '\t' << "X_h [i=0] = " << X[h] << '\n'; //  << '\t' << "Xe_h [i=0] = " << X_e[h]
            
            double t = h*dt;
            double delta = 0.0;
            
            for( int i=0 ; i<n ; ++i ){
            
                /*
                int_coup(i,0) = sin( x(i,0) - x(i,1) );
                loc_coup(i,0) = 0.; loc_coup(i,1) = 0.;
            
                for( int k=0; k<n ; ++k ){
                    float dist = abs(k-i);
                    dist = abs(dist - round(dist/( (float) n ) ) * n );
                    if(dist <= r_one && dist > 0){
                        loc_coup(i,0) += sin( x(i,0) - x(k,0) + alpha );
                        loc_coup(i,1) += sin( x(i,1) - x(k,1) + alpha );
                    }
                }
                
                x(i,0) += (omega(i) - (lamb1/(2*r_one + 1))*loc_coup(i,0))*dt; // + c_sigma * X[h]; // + c_sigma * dW(dt)
                if(layers == 2){
                    x(i,0) += (lamb2/2)*int_coup(i,0)*dt;
                    x(i,1) += (omega(i) - (lamb1/(2*r_one + 1))*loc_coup(i,1) - (lamb2/2)*int_coup(i,0))*dt;
                }
                */
                
                int utime; utime=(int)time(NULL); int seed=utime; //g1.size()
                double x1 = x(i,0); double x2 = x(i,1);
                int D1 = X[ (h + i*n_pts) % (X.size()-100) ]; int D2 = X[ (h + 40 + (i + 100)*n_pts) % (X.size()-70)  ];
                //int D1 = X[ (h + i*steps) % (max - 50) ]; int D2 = X[ (h + 10 + (i + 100)*steps) % (max - 40)   ]; //sizeof(X)/sizeof(*X)
                //int D1 = X_e[ (h + i*steps) % (sizeof(X_e)/sizeof(*X_e)) ]; int D2 = X_e[ (h + 10 + (i + 100)*steps) % (sizeof(X_e)/sizeof(*X_e)) ];
                x(i,0) = rk2_tv_step ( x1, t, dt, 1.0, x, i, n, D1, r_one, lamb1, alpha, c_sigma, layers, lamb2, omega, &seed, &fv1 );
                if(layers == 2){ x(i,1) = rk2_tv_step ( x2, t, dt, 1.0, x, i, n, D2, r_one, lamb1, alpha, c_sigma, layers, lamb2, omega, &seed, &fv2 ); }
                
                double diff = fmod(abs(x(i,1)-x(i,0)), 2*M_PI);
                if(diff > M_PI){diff = 2*M_PI - diff;}
                
                delta += diff * diff; // sin(x(i,1) - x(i,0)) * sin(x(i,1) - x(i,0));
                
                //if(t > sim_time-(3*dt) ){
                //    data_out << c_sigma  << '\t' << lamb1 << '\t' << lamb2 << '\t' << t << '\t' << i << '\t' << x(i,0) << '\t' << x(i,1) <<  '\t' << (x(i,1) - x(i,0)) * (x(i,1) - x(i,0)) << '\n';
                //}
            }
            if(t > sim_time-30){
                ++count;
                delta2 += sqrt( delta );
            }
        }
    }
            
    /*
    if(model == "strogatz"){
        for (int h = 0; h < steps; ++h){
            for( int i=0 ; i<n ; ++i ){
                
                int_coup(i,0) = sin( x(i,0) - x(i,1) );
                loc_coup(i,0) = 0.; loc_coup(i,1) = 0.;
                
                double a = i*2*M_PI/(n-1) - M_PI;
                for( int k=0; k<n ; ++k ){
                    double b = k*2*M_PI/(n-1) - M_PI;

                    loc_coup(i,0) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,0) - x(k,0) + alpha );
                    loc_coup(i,1) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,1) - x(k,1) + alpha );
                }
                x(i,0) = omega(i) - (2*M_PI*lamb1*loc_coup(i,0)/(double) n)*dt + (lamb2/2)*int_coup(i,0);
                x(i,1) = omega(i) - (2*M_PI*lamb1*loc_coup(i,1)/(double) n)*dt - (lamb2/2)*int_coup(i,0);
                
                cout << h << '\t' << "Strogatz x is " << '\t' << x(i,0) << '\t' << x(i,1) << '\n';
                data_out << h << '\t' << x(i,0) << '\t' << x(i,1) << '\n';
            }
        }
    }
    */

    if( count > 0 ){cout << "delta is " << '\t' << delta2 / double( count ) << '\n';}
    if( count > 0 ){data_out << c_sigma  << '\t' << lamb1 << '\t' << lamb2 << '\t' << 0 << '\t' << 0 << '\t' << 0  << '\t' << 0  << '\t' << delta2 / double( count ) << '\n';}
    return 0.0; // ( count != 0 ) ? delta2 / double( count ) : 0.0 ;
}


#endif /* system_h */
