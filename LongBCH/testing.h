//
// Created by ADMIN on 27-Dec-24.
//

#ifndef TESTING_H
#define TESTING_H
#include <vector>
#include <gmpxx.h>
#include <map>
#include "Polynomial.h"

using namespace std;

extern mpz_class n;
extern vector<mpz_class> pi;

Polynomial random_polynomial_factor_mono_degree(int degree, const mpz_class &q);
Polynomial random_polynomial_factor_equal_degree(int degree, const mpz_class &q);

void compare_three_algo(const std::string &filename, const mpz_class &q);
void analyze_factorization_algorithms(const std::string &output_filename, const mpz_class &q);
void compare_Ben_Or_Mignotte_Hybird(const std::string &output_filename, const mpz_class &q);
void get_statistics_attemps_in_node(const mpz_class &q);
void decoding_time_comparison(const std::string &output_filename, const mpz_class &q);
#endif //TESTING_H