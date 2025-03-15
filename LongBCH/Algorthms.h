//
// Created by ADMIN23 on 09.12.2024.
//

#ifndef LONGBCH_ALGORTHMS_H
#define LONGBCH_ALGORTHMS_H
#include "Polynomial.h"
#include <iostream>
#include <gmpxx.h>
using namespace std;
extern int cnt[2];
vector<pair<Polynomial, int>> factor(bool flag, const Polynomial &poly, const mpz_class &q, int &rec_count, int &time);
vector<pair<Polynomial, int>> factor_by_Shoup_Algorithm(const Polynomial &poly, const mpz_class &q, int &time);
vector<Polynomial> equal_degree_factorization_by_Ben_Or(const Polynomial &poly, int degree, const mpz_class &q, int &recursion_count);
std::vector<Polynomial> separate(const Polynomial &g, const Polynomial &s, int degree, const mpz_class &q, int &recursion_count);
vector<pair<Polynomial, mpz_class>> find_root_by_Mignotte(const Polynomial& f, const mpz_class& q, int &time);
vector<Polynomial> find_root_by_Ben_Or(const Polynomial &poly, const mpz_class &q, int &time);
vector<Polynomial> find_root_by_Cantor(const Polynomial &poly, const mpz_class &q, int &time);
vector<Polynomial> equal_degree_factorization_by_hybrid(const Polynomial &poly, int degree, const mpz_class &q, int &recursion_count);
vector<Polynomial> separate_hybrid(const Polynomial &g, const Polynomial &s, int degree, int &cnt_failed, const mpz_class &q, int &recursion_count);
vector<Polynomial> find_root_by_hybrid(const Polynomial &poly, const mpz_class &q, int &time);
#endif //LONGBCH_ALGORTHMS_H
