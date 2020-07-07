#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <stdlib.h>
#include <chrono>

const double pi = atan(1)*4;
const double tolerance = 1e-9;
typedef std::complex<double> cdouble;

cdouble operator*(const cdouble& cd, int i) { return cd * std::complex<double>(i); }
cdouble operator*(int i, const cdouble& cd) { return cd * std::complex<double>(i); }
cdouble operator+(const cdouble& cd, int i) { return cd + std::complex<double>(i); }
cdouble operator+(int i, const cdouble& cd) { return cd + std::complex<double>(i); }
cdouble operator-(const cdouble& cd, int i) { return cd - std::complex<double>(i); }
cdouble operator-(int i, const cdouble& cd) { return std::complex<double>(i) - cd; }

double f_C7(const double &obs, const double &c, const double &b) {
// Calculate the C7 coefficient for a given combination of source and 
// observer radius and observer angle
    // Static allocate all the variables to avoid redundant allocations
    static cdouble s_obs, c_obs, c_obs2, s_obs2, c2, b2, C0, C1, C2, C3, C4, C5;

    // Do some precalculating to reduce the number of floating point operations required
    s_obs = sin(obs);
    c_obs = cos(obs);
    c_obs2 = std::pow(c_obs, 2);
    s_obs2 = std::pow(s_obs, 2);
    c2 = std::pow(c, 2);
    b2 = std::pow(b, 2);

    // Calculate the coefficients from the paper
    C0 = b2 + c2 + 2 * b * c * s_obs - 4;
    C1 = 24 * b * c_obs2 * (b - c * s_obs) - 48 * (-1 + c2) * c_obs2 + std::pow(C0, 2); 
    C2 = -432 * b2 * (-1 + c2) * std::pow(c_obs2, 2) + 432 * c_obs2 * std::pow(b - c * s_obs,2) + 72 * b * c_obs2 * (b - c * s_obs) * C0 + 
         288 * (-1 + c2) * c_obs2 * C0 + 2 * std::pow(C0, 3); 
    C3 = std::pow(C2 + std::sqrt(-4 * std::pow(C1, 3) + std::pow(C2, 2)), 1.0 / 3.0); 
    C4 = C1 / (6 * std::pow(2, 2.0 / 3.0) * C3) + C3 / (12 * std::cbrt(2)); 
    C5 = std::sqrt((b2 * c_obs2) / 4.0 - C0 / 6.0 + C4);

    // Return the C7 coefficient
    return (b2 * c_obs2 / 2.0 - C0 / 3.0 - C4 + (std::pow(b, 3) * std::pow(c_obs, 3) - 4 * c_obs * (b - c * s_obs) - b * c_obs * C0) / (4 * C5)).real();
}

double onefinite(const double &obs, const double &c) {
/* 
   Calculate the specular point for an infinite source and given
   observer radius and angle
        obs - observer angle in radians measured from positive x-axis in 
              the rotated reference frame (theta_obs)
 
        c   - c constant, radius of sphere divided by radius of source 
              (R_sph / R_src)
  
        b   - b constant, radius of sphere divided by radius of observer
              (R_sph / R_obs)
*/
    // Static allocate all the variables to avoid redundant allocations
    static cdouble s_obs, c_obs, c_obs2, s_obs2, c2, c2_4, c2_43, c4, b2, C0, C1, C2, C3, C4;
    static double p1, p2, l;

    // Do precalculations to reduce the number of floating point operations required
    s_obs = sin(2 * obs);
    c_obs = cos(2 * obs);
    c2 = std::pow(c, 2);
    c2_4 = c2 - 4;
    c2_43 = std::pow(c2_4, 3);
    c4 = std::pow(c, 4);

    // Calculated the one-finite coefficients from the paper
    C0 = 1152 * c2_4 * (-1 + c2 - c_obs + c2 * c_obs) + 864 * c2 * std::pow(s_obs, 2);
    C1 = std::sqrt(-4 * std::pow(160 - 128 * c2 + 4 * c4 + 96 * c_obs - 96 * c2 * c_obs, 3) + std::pow(-16 * c2_43 - C0, 2));
    C2 = std::pow(16 * c2_43 + C0 - C1, 1.0 / 3.0);
    C3 = (40 - 32 * c2 + c4 + 24 * c_obs - 24 * c2 * c_obs) / (3 * std::pow(2, 2.0 / 3.0) * C2) + C2 / (24 * std::cbrt(2));
    C4 = std::sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + C3);
    
    p1 = acos((0.5 * C4 - 0.5 * std::sqrt(-(c2_4 / 3.0) - C3 + (c * s_obs) / (2 * C4))).real());
    p2 = acos((0.5 * C4 + 0.5 * std::sqrt(-(c2_4 / 3.0) - C3 + (c * s_obs) / (2 * C4))).real());
    l = pi * (tan(1.0 / c) - obs) / (pi + 2 * tan(1.0 / c));

    // Chose the branch which is less close to the line from (-pi / 2, pi / 2) to ({l1}, 0)
    if (abs(p1 - l) >= abs(p2 - l)) return p1;
    else return p2;
}

double bothfinite(const double &obs, const double &c, const double &b) {
/* 
   Calculates the specular point for a given combination of source and
   observer radius and angle
        obs - observer angle in radians measured from positive x-axis in 
              the rotated reference frame (theta_obs)
  
        c   - c constant, radius of sphere divided by radius of source 
              (R_sph / R_src)
  
        b   - b constant, radius of sphere divided by radius of observer
              (R_sph / R_obs)
*/
    // Static allocate all the variables to avoid redundant allocations
    static cdouble s_obs, c_obs, c_obs2, s_obs2, c2, b2, C0, C1, C2, C3, C4, C5, C7;
    static double l1, l2, deriv;

    
    // Do precalculations to reduce the number of floating point operations required
    s_obs = sin(obs);
    c_obs = cos(obs);
    c_obs2 = std::pow(c_obs, 2);
    s_obs2 = std::pow(s_obs, 2);
    c2 = std::pow(c, 2);
    b2 = std::pow(b, 2);

    // Calculate the both-finite coefficients from the paper
    C0 = b2 + c2 + 2 * b * c * s_obs - 4;
    C1 = 24 * b * c_obs2 * (b - c * s_obs) - 48 * (-1 + c2) * c_obs2 + std::pow(C0, 2); 
    C2 = -432 * b2 * (-1 + c2) * std::pow(c_obs2, 2) + 432 * c_obs2 * std::pow(b - c * s_obs, 2) + 72 * b * c_obs2 * (b - c * s_obs) * C0 + 
                        288 * (-1 + c2) * c_obs2 * C0 + 2 * std::pow(C0, 3); 
    C3 = std::pow(C2 + std::sqrt(-4 * std::pow(C1, 3) + std::pow(C2, 2)), 1.0 / 3.0); 
    C4 = C1 / (6 * std::pow(2, 2.0 / 3.0) * C3) + C3 / (12 * std::cbrt(2)); 
    C5 = std::sqrt((b2 * c_obs2) / 4.0 - C0 / 6.0 + C4);
    C7 = b2 * c_obs2 / 2.0 - C0 / 3.0 - C4 + (std::pow(b, 3) * std::pow(c_obs, 3) - 4 * c_obs * (b - c * s_obs) - b * c_obs * C0) / (4 * C5);

    l1 = -std::atan2((c + b * std::sqrt(1 - b2 + c2)).real(), (b2 - c2).real());
    l2 = -std::atan2((c - b * std::sqrt(1 - b2 + c2)).real(), (b2 - c2).real());

    // If we are below the lower limit {l1} we know we are on the -+++ branch
    if (-pi / 2 <= obs && obs < l1) return -acos((b * c_obs) * 0.25 + C5 * 0.5 + std::sqrt(C7)*0.5).real();
	// Else, if we are above {l2} we know we have passed the intesection between the 
	// ++++ and ++-+ branches, and are on the ++-+ branch. 
    else if (obs > l2)  return acos((b * c_obs) * 0.25 + C5 * 0.5 - std::sqrt(C7)*0.5).real();
    else {
        // Otherwise, the branch deduction requires analysis of the derivative
        deriv = C7.real() - f_C7(obs - tolerance, c, b);
        // If the derivative less than 0 and we are below {l2} we are on the ++++ branch
        if (obs <= l2 && deriv <= 0) return acos((b * c_obs) * 0.25 + C5 * 0.5 + std::sqrt(C7)*0.5).real();
        // If the derivative is greater than zero or if we are above {l2} we must be on the ++-+ branch
        else if ((l2 < obs && obs < pi / 2) || deriv > 0) return acos((b * c_obs) * 0.25 + C5 * 0.5 - std::sqrt(C7)*0.5).real();
    }
}

// Calculate a pseudorandom floating point number between {a} and {b}
double randf(const double &a, const double &b) { return ((double)rand() / (double)RAND_MAX) * (b - a) + a; }

int main(int argn, char** argv) {
    if (argn < 3) {
    	// If too few arguments are specified, default to benchmarking
    	
        // One-finite
        auto start = std::chrono::high_resolution_clock::now();
        int percent = 0;
        long end = 1e7;
        std::cout << std::setw(3) << percent << std::flush;
        for (long ii = 0; ii < end; ii ++) {
            onefinite(randf(0, pi), randf(0, 1.0));
            if ((100 * ii) / end > percent) {
                percent = (100 * ii) / end;
                std::cout << '\r' << std::setw(3) << percent << std::flush;
            }
        }
        std::cout << "\r100" << std::flush;
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout << "\rOne finite  : ";
        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        std::cout << " us (" << std::setprecision(0) << std::scientific << (double)end << " iterations)" << std::endl;
        std::cout << std::defaultfloat;

        // Both-finite
        start = std::chrono::high_resolution_clock::now();
        percent = 0;
        end = 1e7;
        std::cout << std::setw(3) << percent << std::flush;
        for (long ii = 0; ii < end; ii ++) {
            bothfinite(randf(0, pi), randf(0, 1.0), randf(0, 1.0));
            if ((100 * ii) / end > percent) {
                percent = (100 * ii) / end;
                std::cout << '\r' << std::setw(3) << percent << std::flush;
            }
        }
        std::cout << "\r100" << std::flush;
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "\rBoth finite : ";
        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        std::cout << " us (" << std::setprecision(0) << std::scientific << (double)end << " iterations)" << std::endl;
        std::cout << std::defaultfloat;
    }
    else {
        std::cout << "Note: specify observer angle in degrees." << std::endl;
        double obs = atof(argv[1]) * pi / 180.0;
        double c = atof(argv[2]);
        if (argn > 3) {
            double b = atof(argv[3]);
            std::cout << bothfinite(obs, c, b) << std::endl;
        }
        else {
            std::cout << onefinite(obs, c) << std::endl;
        }
    }
}
