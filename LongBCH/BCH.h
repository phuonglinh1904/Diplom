//
// Created by ADMIN23 on 23.03.2025.
//

#ifndef LONGBCH_BCH_H
#define LONGBCH_BCH_H
#include "Polynomial.h"
#include "GFq.h"
#include  "Algorthms.h"
#endif //LONGBCH_BCH_H
map<mpz_class , vector<mpz_class>> find_cyclotomic_classes();
Polynomial compute_g();
std::vector<mpz_class> encode(const std::vector<mpz_class>& a, const Polynomial& g, const mpz_class& q);
vector<mpz_class> decode(const vector<mpz_class> &v, const mpz_class &q, bool is_hybrid, int &time);