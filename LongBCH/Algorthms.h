//
// Created by ADMIN23 on 09.12.2024.
//

#ifndef LONGBCH_ALGORTHMS_H
#define LONGBCH_ALGORTHMS_H
#include "Polynomial.h"
#include <iostream>
#include <gmpxx.h>
#include <fstream>
using namespace std;
extern int cnt[2];
extern vector<mpz_class> pi;
vector<pair<Polynomial, int>> factor(bool flag, const Polynomial &g, const mpz_class &q, int &time);
vector<pair<Polynomial, int>> factor_by_Shoup_Algorithm(const Polynomial &g, const mpz_class &q, int &time);
vector<Polynomial> equal_degree_factorization_by_Ben_Or(const Polynomial &poly, int degree, const mpz_class &q);
std::vector<Polynomial> separate(const Polynomial &g, const Polynomial &s, int degree, const mpz_class &q);
vector<pair<Polynomial, mpz_class>> find_root_by_Mignotte(const Polynomial& g, const mpz_class& q, int &time);
vector<Polynomial> find_root_by_Ben_Or(const Polynomial &poly, const mpz_class &q, int &time);
vector<Polynomial> find_root_by_Cantor(const Polynomial &poly, const mpz_class &q, int &time);
std::vector<pair<Polynomial, mpz_class>> hybrid_factorization(const Polynomial &f, int degree, const mpz_class &q, int &time);

vector<Polynomial> equal_degree_factorization_by_hybird(const Polynomial &poly, int degree, const mpz_class &q);
vector<Polynomial> separate_hybird(const Polynomial &g, const Polynomial &s, int degree, int &cnt_failed, const mpz_class &q);
vector<Polynomial> find_root_by_hybrid(const Polynomial &poly, const mpz_class &q, int &time);
vector<mpz_class> find_root_by_chien(const Polynomial &poly, int &time);
vector<pair<Polynomial, mpz_class>> find_root_by_hybird_1(const Polynomial& g, const mpz_class& q, int &time);
#endif //LONGBCH_ALGORTHMS_H
