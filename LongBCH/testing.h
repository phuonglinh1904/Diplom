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
struct factor_result {
    Polynomial poly;
    int recursion_count;
    int time_run;

    factor_result(Polynomial poly, int recursion_count, int time_run) : poly(poly), recursion_count(recursion_count),
                                                                        time_run(time_run) {}
};

void write_to_csv(const std::vector<factor_result> &results, const std::string &filename);
void write_result_2_to_csv(const std::vector<factor_result> &results, const std::string &filename);
void compare_three_algo(const std::string &filename, const mpz_class &q);
void analyze_factorization_algorithms(const std::string &output_filename, const mpz_class &q);
Polynomial random_polynomial_factor_1_degree(int degree);
Polynomial random_polynomial_factor_mono_degree(int degree, const mpz_class &q);
void compare_Ben_Or_Mignotte(const std::string &output_filename, const mpz_class &q);
#endif //TESTING_H