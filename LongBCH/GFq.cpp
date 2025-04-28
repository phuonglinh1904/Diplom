//
// Created by ADMIN23 on 09.12.2024.
//

#include <vector>
#include <gmpxx.h>
#include <map>
#include <random>
#include "GFq.h"
#include <iostream>
#include <ctime>

using namespace std;


std::vector<std::string> split_by(const std::string &s, const std::string &delimiter) {
    std::vector<std::string> result;
    auto start = 0U;
    auto end = s.find(delimiter);
    while (end != std::string::npos) {
        result.push_back(s.substr(start, end - start));
        start = end + delimiter.length();
        end = s.find(delimiter, start);
    }
    result.push_back(s.substr(start));
    return result;
}

mpz_class mul_alpha(const mpz_class &a) {
    mpz_class c = a;
    c <<= 1;

    mpz_class temp = c >> m.get_ui();

    mpz_class bit_check;
    mpz_and(bit_check.get_mpz_t(), temp.get_mpz_t(), mpz_class(1).get_mpz_t());
    if (bit_check != 0) {
        c ^= primitive;
    }

    return c;
}

void initialize_tables_and_pi() {
    mpz_class id_tmp;
    log_alpha.resize(n.get_ui() + 1);
    deg_alpha.resize(2 * n.get_ui() + 1);
    deg_alpha[0] = mpz_class("1");
    for (mpz_class i = 1; i <= n; ++i) {
        deg_alpha[i.get_ui()] = mul_alpha(deg_alpha[i.get_ui() - 1]);
        id_tmp = n + i;
        deg_alpha[id_tmp.get_ui()] = deg_alpha[i.get_ui()];
        log_alpha[deg_alpha[i.get_ui()].get_ui()] = i;
    }

    vector<mpz_class> primes = SieveOfEratosthenes();
    vector<pair<mpz_class, mpz_class>> factors = modified_trial_division( primes);
    pi = calculate_pi(factors);
}

mpz_class gf_multiply(const mpz_class &a, const mpz_class &b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    mpz_class index = log_alpha[a.get_ui()] + log_alpha[b.get_ui()];
    return deg_alpha[index.get_ui()];
}

mpz_class gf_inverse(const mpz_class &a) {
    mpz_class index = n - log_alpha[a.get_ui()];
    return deg_alpha[index.get_ui()];
}

mpz_class gf_add(const mpz_class &a, const mpz_class &b) {
    return a ^ b;
}

mpz_class gf_power(const mpz_class &base, mpz_class exponent) {
    mpz_class res = 1;
    mpz_class a = base;

    while (exponent > 0) {
        if (exponent % 2 == 1) {
            res = gf_multiply(res, a);
        }
        a = gf_multiply(a, a);
        exponent >>= 1;
    }
    return res;
}

