#ifndef int_method_h
#define int_method_h

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

//double fv ( double t, double x, matrixd xm, int i, int n, int r_one, double lamb1, double alpha, int layers, double lamb2 );
double gv ( double t, double x, double D, matrixd xm, int i, double xorig);
double r8_normal_01 ( int *seed );
double r8_uniform_01 ( int *seed );
void timestamp ( );


double rk2_tv_step ( double x, double t, double h, double q, matrixd xm, int i, int n, double w, int r_one, double lamb1, double alpha, double D, int layers, double lamb2, vectord omega,
                    int *seed, double (*fv) ( double, double, matrixd, int, int, int, double, double, int, double, vectord ))
{
  double a21; double a31; double a32; double k1; double k2;
  double q1; double q2; double t1; double t2; // double w1;
  /*double w2;*/ double x1; double x2; double xstar;

  a21 =   1.0;
  a31 =   0.5;
  a32 =   0.5;

  q1 = 2.0;
  q2 = 2.0;

  t1 = t;
  x1 = x;
  //w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fv ( t1, x1, xm, i, n, r_one, lamb1, alpha, layers, lamb2, omega ) + h * gv ( t1, x1, D, xm, i, x ) * w;  // + h * gv ( t1, x1, D ) * w1

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  //w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fv ( t2, x2, xm, i, n, r_one, lamb1, alpha, layers, lamb2, omega ) + h * gv ( t2, x2, D, xm, i, x) * w;  // + h * gv ( t1, x1, D ) * w2

  xstar = x1 + a31 * k1 + a32 * k2;

  return xstar;
}

/*
double rk4_tv_step ( double x, double t, double h, double q, matrixd xm, int i, int n, double w, int r_one, double lamb1, double alpha, double D, int *seed )

//    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
//    Input, double H, the time step.
//    Input, double Q, the spectral density of the input white noise.
//    Input/output, int *SEED, a seed for the random number generator.
{
  double a21; double a31; double a32; double a41; double a42; double a43; double a51; double a52; double a53; double a54; double k1;
  double k2; double k3; double k4; double q1; double q2; double q3; double q4; double t1; double t2; double t3; double t4; double w1;
  double w2; double w3; double w4; double x1; double x2; double x3; double x4; double xstar;

  a21 =   0.66667754298442; a31 =   0.63493935027993; a32 =   0.00342761715422; a41 = - 2.32428921184321;
  a42 =   2.69723745129487; a43 =   0.29093673271592; a51 =   0.25001351164789; a52 =   0.67428574806272;
  a53 = - 0.00831795169360; a54 =   0.08401868181222;

  q1 = 3.99956364361748; q2 = 1.64524970733585; q3 = 1.59330355118722; q4 = 0.26330006501868;

  t1 = t;
  x1 = x;
  w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
  k1 = h * fv ( t1, x1, xm, i, n, r_one, lamb1, alpha ) + h * gv ( t1, x1, D ) * w1;

  t2 = t1 + a21 * h;
  x2 = x1 + a21 * k1;
  w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
  k2 = h * fv ( t2, x2, xm, i, n, r_one, lamb1, alpha ) + h * gv ( t2, x2, D ) * w2;

  t3 = t1 + a31 * h  + a32 * h;
  x3 = x1 + a31 * k1 + a32 * k2;
  w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
  k3 = h * fv ( t3, x3, xm, i, n, r_one, lamb1, alpha ) + h * gv ( t3, x3, D ) * w3;

  t4 = t1 + a41 * h  + a42 * h  + a43 * h;
  x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3;
  w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h );
  k4 = h * fv ( t4, x4, xm, i, n, r_one, lamb1, alpha ) + h * gv ( t4, x4, D ) * w4;

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;

  return xstar;
}
*/

























double r8_normal_01 ( int *seed )

//    R8_NORMAL_01 returns a unit pseudonormal R8.
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    Because this routine uses the Box Muller method, it requires pairs
//    of uniform random values to generate a pair of normal random values.
//    This means that on every other call, the code can use the second
//    value that it calculated.
//
//    However, if the user has changed the SEED value between calls,
//    the routine automatically resets itself and discards the saved data.
//
//  Parameters:
//    Input/output, int *SEED, a seed for the random number generator.
//    Output, double R8_NORMAL_01, a normally distributed random value.

{
# define R8_PI 3.141592653589793

  double r1;
  double r2;
  static int seed2 = 0;
  static int seed3 = 0;
  static int used = 0;
  double v1;
  static double v2 = 0.0;

//  If USED is odd, but the input SEED does not match
//  the output SEED on the previous call, then the user has changed
//  the seed.  Wipe out internal memory.

  if ( ( used % 2 ) == 1 )
  {
    if ( *seed != seed2 )
    {
      used = 0;
      seed2 = 0;
      seed3 = 0;
      v2 = 0.0;
    }
  }

//  If USED is even, generate two uniforms, create two normals,
//  return the first normal and its corresponding seed.

  if ( ( used % 2 ) == 0 )
  {
    r1 = r8_uniform_01 ( seed );

    if ( r1 == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_NORMAL_01 - Fatal error!\n";
      cerr << "  R8_UNIFORM_01 returned a value of 0.\n";
      exit ( 1 );
    }

    seed2 = *seed;
    r2 = r8_uniform_01 ( seed );
    seed3 = *seed;
    *seed = seed2;

    v1 = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * R8_PI * r2 );
    v2 = sqrt ( - 2.0 * log ( r1 ) ) * sin ( 2.0 * R8_PI * r2 );
  }

//  If USED is odd (and the input SEED matched the output value from
//  the previous call), return the second normal and its corresponding seed.

  else
  {
    v1 = v2;
    *seed = seed3;
  }

  used = used + 1;

  return v1;
# undef R8_PI
}


double r8_uniform_01 ( int *seed )
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}


void timestamp ( )
//  TIMESTAMP prints the current YMDHMS date as a time stamp.
//  Example: 31 May 2001 09:45:54 AM
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}




#endif /* int_method_h */
