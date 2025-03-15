//
// Created by ADMIN23 on 09.12.2024.
//

#ifndef LONGBCH_GFQ_H
#define LONGBCH_GFQ_H
#include <vector>
#include <string>
#include <map>
#include "Polynomial.h"

extern mpz_class m;
extern mpz_class primitive;
extern mpz_class n;
extern mpz_class epsilon;
extern std::vector<mpz_class> deg_alpha;
extern std::vector<mpz_class> log_alpha;
extern std::map<mpz_class, std::pair<std::vector<mpz_class>, mpz_class>> cyclotomic;
std::vector<std::string> split_by(const std::string& s, const std::string& delimiter);
mpz_class mul_alpha(const mpz_class& a);
void initialize_tables();
mpz_class gf_multiply(const mpz_class& a, const mpz_class& b);
mpz_class gf_inverse(const mpz_class& a);
mpz_class gf_add(const mpz_class &a, const mpz_class &b);
mpz_class gf_power(const mpz_class &base, mpz_class exponent);
std::map<mpz_class, std::pair<std::vector<mpz_class>, mpz_class>> find_cyclotomic_classes();
std::map<mpz_class, std::vector<std::vector<mpz_class>>> group_by_length_ordered(const std::map<mpz_class, std::pair<std::vector<mpz_class>, mpz_class>> &cyclotomic);
std::pair<mpz_class, mpz_class> random_two_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped);
std::tuple<mpz_class, mpz_class, mpz_class> random_three_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped);
std::tuple<mpz_class, mpz_class, mpz_class, mpz_class> random_four_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped);
#endif //LONGBCH_GFQ_H