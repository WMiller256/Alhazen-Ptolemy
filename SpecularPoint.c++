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

double f_C7(const double &d, const double &c, const double &b) {
	static cdouble sd = sin(d);
	static cdouble cd = cos(d);
	static cdouble cd2 = std::pow(cd, 2);
	static cdouble sd2 = std::pow(sd, 2);
	static cdouble c2 = std::pow(c, 2);
	static cdouble b2 = std::pow(b, 2);

	static cdouble C0 = b2 + c2 + 2 * b * c * sd - 4;
	static cdouble C1 = 24 * b * cd2 * (b - c * sd) - 48 * (-1 + c2) * cd2 + std::pow(C0, 2); 
	static cdouble C2 = -432 * b2 * (-1 + c2) * std::pow(cd2, 2) + 432 * cd2 * std::pow(b - c * sd,2) + 72 * b * cd2 * (b - c * sd) * C0 + 
				 288 * (-1 + c2) * cd2 * C0 + 2 * std::pow(C0, 3); 
	static cdouble C3 = std::pow(C2 + std::sqrt(-4 * std::pow(C1, 3) + std::pow(C2, 2)), 1.0 / 3.0); 
	static cdouble C4 = C1 / (6 * std::pow(2, 2.0 / 3.0) * C3) + C3 / (12 * std::cbrt(2)); 
	static cdouble C5 = std::sqrt((b2 * cd2) / 4.0 - C0 / 6.0 + C4);
	return (b2 * cd2 / 2.0 - C0 / 3.0 - C4 + (std::pow(b, 3) * std::pow(cd, 3) - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5)).real();
}

double onefinite(const double &d, const double &c, bool print = false) {
	cdouble sd = sin(2 * d);
	cdouble cd = cos(2 * d);
	cdouble c2 = std::pow(c, 2);
	cdouble c2_4 = c2 - 4;
	cdouble c2_43 = std::pow(c2_4, 3);
	cdouble c4 = std::pow(c, 4);

	cdouble C0 = 1152 * c2_4 * (-1 + c2 - cd + c2 * cd) + 864 * c2 * std::pow(sd, 2);
	cdouble C1 = std::sqrt(-4 * std::pow(160 - 128 * c2 + 4 * c4 + 96 * cd - 96 * c2 * cd, 3) + std::pow(-16 * c2_43 - C0, 2));
	cdouble C2 = std::pow(16 * c2_43 + C0 - C1, 1.0 / 3.0);
	cdouble C3 = (40 - 32 * c2 + c4 + 24 * cd - 24 * c2 * cd) / (3 * std::pow(2, 2.0 / 3.0) * C2) + C2 / (24 * std::cbrt(2));
	cdouble C4 = std::sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + C3);
	
	double p1 = acos((0.5 * C4 - 0.5 * std::sqrt(-(c2_4 / 3.0) - C3 + (c * sd) / (2 * C4))).real());
	double p2 = acos((0.5 * C4 + 0.5 * std::sqrt(-(c2_4 / 3.0) - C3 + (c * sd) / (2 * C4))).real());
	double l1 = pi * (tan(1.0 / c) - d) / (pi + 2 * tan(1.0 / c));
	
	if (abs(p1 - l1) >= abs(p2 - l1)) return p1;
	else return p2;
}

double twofinite(const double &d, const double &c, const double &b) {
	cdouble sd = sin(d);
	cdouble cd = cos(d);
	cdouble cd2 = std::pow(cd, 2);
	cdouble sd2 = std::pow(sd, 2);
	cdouble c2 = std::pow(c, 2);
	cdouble b2 = std::pow(b, 2);

	cdouble C0 = b2 + c2 + 2 * b * c * sd - 4;
	cdouble C1 = 24 * b * cd2 * (b - c * sd) - 48 * (-1 + c2) * cd2 + std::pow(C0, 2); 
	cdouble C2 = -432 * b2 * (-1 + c2) * std::pow(cd2, 2) + 432 * cd2 * std::pow(b - c * sd,2) + 72 * b * cd2 * (b - c * sd) * C0 + 
				 288 * (-1 + c2) * cd2 * C0 + 2 * std::pow(C0, 3); 
	cdouble C3 = std::pow(C2 + std::sqrt(-4 * std::pow(C1, 3) + std::pow(C2, 2)), 1.0 / 3.0); 
	cdouble C4 = C1 / (6 * std::pow(2, 2.0 / 3.0) * C3) + C3 / (12 * std::cbrt(2)); 
	cdouble C5 = std::sqrt((b2 * cd2) / 4.0 - C0 / 6.0 + C4);
	cdouble C7 = b2 * cd2 / 2.0 - C0 / 3.0 - C4 + (std::pow(b, 3) * std::pow(cd, 3) - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5);

	double l1 = -std::atan2((c + b * sqrt(1 - b2 + c2)).real(), (b2 - c2).real());
	double l2 = -std::atan2((c - b * sqrt(1 - b2 + c2)).real(), (b2 - c2).real());

	if (-pi / 2 <= d && d < l1) {
		return -acos((b * cd) * 0.25 + C5 * 0.5 + std::sqrt(C7)*0.5).real();
	}
	else {
		double deriv = C7.real() - f_C7(d - tolerance, c, b);
		if (d <= l2 && deriv <= 0)                    return  acos((b * cd) * 0.25 + C5 * 0.5 + std::sqrt(C7)*0.5).real();
		else if ((l2 < d && d < pi / 2) || deriv > 0) return  acos((b * cd) * 0.25 + C5 * 0.5 - std::sqrt(C7)*0.5).real();
	}
}

double randf(double a, double b) { return (double(rand()) / double(RAND_MAX)) * (b - a) + a; }

int main(int argn, char** argv) {
	if (argn < 3) {
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
			twofinite(randf(0, pi), randf(0, 1.0), randf(0, 1.0));
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
		std::cout << "d is in degrees" << std::endl;
		double d = atof(argv[1]) * pi / 180.0;
		double c = atof(argv[2]);
		if (argn > 3) {
			double b = atof(argv[3]);
			std::cout << twofinite(d, c, b) << std::endl;
		}
		else {
			std::cout << onefinite(d, c, true) << std::endl;
		}
	}
}
