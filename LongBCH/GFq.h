//
// Created by ADMIN23 on 09.12.2024.
//

#ifndef LONGBCH_GFQ_H
#define LONGBCH_GFQ_H
#include <vector>
#include <string>
#include <map>
#include "Polynomial.h"
#include "Utils.h"

extern mpz_class m;
extern mpz_class primitive;
extern mpz_class n;
extern std::vector<mpz_class> deg_alpha;
extern std::vector<mpz_class> log_alpha;
extern std::vector<mpz_class> pi;
extern std::map<mpz_class, std::vector<mpz_class>> cyclotomic;
extern mpz_class delta;
std::vector<std::string> split_by(const std::string& s, const std::string& delimiter);
mpz_class mul_alpha(const mpz_class& a);
void initialize_tables_and_pi();
mpz_class gf_multiply(const mpz_class& a, const mpz_class& b);
mpz_class gf_inverse(const mpz_class& a);
mpz_class gf_add(const mpz_class &a, const mpz_class &b);
mpz_class gf_power(const mpz_class &base, mpz_class exponent);
#endif //LONGBCH_GFQ_H