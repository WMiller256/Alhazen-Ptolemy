#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <stdlib.h>
#include <chrono>

const double pi = atan(1)*4;
typedef std::complex<double> cdouble;

cdouble operator*(const cdouble& cd, int i) { return cd * std::complex<double>(i); }
cdouble operator*(int i, const cdouble& cd) { return cd * std::complex<double>(i); }
cdouble operator+(const cdouble& cd, int i) { return cd + std::complex<double>(i); }
cdouble operator+(int i, const cdouble& cd) { return cd + std::complex<double>(i); }
cdouble operator-(const cdouble& cd, int i) { return cd - std::complex<double>(i); }
cdouble operator-(int i, const cdouble& cd) { return std::complex<double>(i) - cd; }

double onefinite(const double &d, const double &rs, bool print = false) {
	cdouble sd = sin(2 * d);
	cdouble cd = cos(2 * d);
	cdouble rs2 = std::pow(rs, 2);
	cdouble rs2_4 = rs2 - 4;
	cdouble rs2_43 = std::pow(rs2_4, 3);
	cdouble rs4 = std::pow(rs, 4);

	cdouble C0 = 1152 * rs2_4 * (-1 + rs2 - cd + rs2 * cd) + 864 * rs2 * std::pow(sd, 2);
	cdouble C1 = std::sqrt(-4 * std::pow(160 - 128 * rs2 + 4 * rs4 + 96 * cd - 96 * rs2 * cd, 3) + std::pow(-16 * rs2_43 - C0, 2));
	cdouble C2 = std::pow(16 * rs2_43 + C0 - C1, 1.0 / 3.0);
	cdouble C3 = (40 - 32 * rs2 + rs4 + 24 * cd - 24 * rs2 * cd) / (3 * std::pow(2, 2.0 / 3.0) * C2) + C2 / (24 * std::cbrt(2));
	cdouble C4 = std::sqrt(-0.25 * rs2_4 + (1.0 / 12.0) * rs2_4 + C3);
	double p = acos((0.5 * C4 - 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real());
 	if (print) {
		std::cout << "+--- " << acos((-0.5 * C4 - 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "+-+- " << acos((-0.5 * C4 + 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "++-+ " << p << std::endl;
		std::cout << "++++ " << acos((0.5 * C4 + 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "-+++ " << -acos((0.5 * C4 + 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "-+-+ " << -acos((0.5 * C4 - 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "--+- " << -acos((-0.5 * C4 + 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))).real()) << std::endl;
		std::cout << "---- " << -acos((-0.5 * C4 - 0.5 * std::sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))).real()) << std::endl;
	}
	return p;
}

double twofinite(const double &d, const double &rs, const double &rd, bool print = false) {
	cdouble sd = sin(d);
	cdouble cd = cos(d);
	cdouble cd2 = std::pow(cd, 2);
	cdouble sd2 = std::pow(sd, 2);
	cdouble rs2 = std::pow(rs, 2);
	cdouble rd2 = std::pow(rd, 2);

	cdouble C0 = rd2 + rs2 + 2 * rd * rs * sd - 4;
	cdouble C5 = 24 * rd * cd2 * (rd - rs * sd) - 48 * (-1 + rs2) * cd2 + std::pow(C0, 2); 
	cdouble C6 = -432 * rd2 * (-1 + rs2) * std::pow(cd2, 2) + 432 * cd2 * std::pow(rd - rs * sd,2) + 72 * rd * cd2 * (rd - rs * sd) * C0 + 288 * (-1 + rs2) * cd2 * C0 + 2 * std::pow(C0, 3); 
	cdouble C1 = std::pow(C6 + std::sqrt(-4 * std::pow(C5, 3) + std::pow(C6, 2)), 1.0 / 3.0); 
	cdouble C2 = (rd2 * cd2) / 4.0 - C0 / 6.0; 
	cdouble C3 = C5 / (6 * std::pow(2, 2.0 / 3.0) * C1) + C1 / (12 * std::cbrt(2)); 
	cdouble C4 = std::sqrt(C2 + C3);
	double p = acos((rd * cd) * 0.25 + C4 * 0.5 - std::sqrt(2 * C2 - C3 + (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real();
	if (print) {
		std::cout << "+--- " << acos((rd * cd) * 0.25 - C4 * 0.5 - std::sqrt(2 * C2 - C3 - (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
		std::cout << "+-+- " << acos((rd * cd) * 0.25 - C4 * 0.5 + std::sqrt(2 * C2 - C3 - (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
		std::cout << "++-+ " << p << std::endl;
		std::cout << "++++ " << acos((rd * cd) * 0.25 + C4 * 0.5 + std::sqrt(2 * C2 - C3 + (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
		std::cout << "-+++ " << -acos((rd * cd) * 0.25 + C4 * 0.5 + std::sqrt(2 * C2 - C3 + (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
		std::cout << "-+-+ " << -p << std::endl;
		std::cout << "--+- " << -acos((rd * cd) * 0.25 - C4 * 0.5 + std::sqrt(2 * C2 - C3 - (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
		std::cout << "---- " << -acos((rd * cd) * 0.25 - C4 * 0.5 - std::sqrt(2 * C2 - C3 - (std::pow(rd, 3) * std::pow(cd, 3) - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real() << std::endl;
	}

	return p;
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
		double rs = atof(argv[2]);
		if (argn > 3) {
			double rd = atof(argv[3]);
			std::cout << twofinite(d, rs, rd, true) << std::endl;
		}
		else {
			std::cout << onefinite(d, rs, true) << std::endl;
		}
	}
}
