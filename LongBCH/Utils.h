//
// Created by ADMIN23 on 30.12.2024.
//

#ifndef LONGBCH_UTILS_H
#define LONGBCH_UTILS_H
#include "Polynomial.h"
#include <map>
#include <string>

vector<mpz_class> SieveOfEratosthenes();
vector<pair<mpz_class, mpz_class>> modified_trial_division(const vector<mpz_class>& P);
vector<mpz_class> calculate_pi(const vector<pair<mpz_class, mpz_class>>& factors);
#endif //LONGBCH_UTILS_H
