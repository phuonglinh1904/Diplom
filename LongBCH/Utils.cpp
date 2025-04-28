//
// Created by ADMIN23 on 30.12.2024.
//
#include "Utils.h"
#include "GFq.h"
//Eratosthenes Sieve
vector<mpz_class> SieveOfEratosthenes() {
    unsigned long bound = n.get_ui();
    vector<bool> a(bound + 1, true);

    // Sàng Eratosthenes
    for (unsigned long i = 2; i * i <= bound; ++i) {
        if (a[i]) {
            for (unsigned long j = i * i; j <= bound; j += i) {
                a[j] = false;
            }
        }
    }
    size_t size = 0;
    for (unsigned long i = 2; i <= bound; ++i) {
        if (a[i]) ++size;
    }

    vector<mpz_class> primes(size);
    size_t index = 0;
    for (unsigned long i = 2; i <= bound; ++i) {
        if (a[i]) primes[index++] = i;
    }

    return primes;
}

// Hàm phân tích thừa số của q-1
vector<pair<mpz_class, mpz_class>> modified_trial_division(const vector<mpz_class>& P) {
    vector<pair<mpz_class, mpz_class>> factors;
    mpz_class temp = n;
    for (size_t i = 0; i < P.size(); i++) {
        mpz_class count = 0;
        while (temp % P[i] == 0) {
            temp = temp / P[i];
            count++;
        }
        if (count > 0) {
            factors.push_back(make_pair(P[i], count));
        }
    }
    if (temp > 1) {
        factors.push_back(make_pair(temp, 1));
    }

    return factors;
}
// Hàm tính các π_i từ phân tích thừa số
vector<mpz_class> calculate_pi(const vector<pair<mpz_class, mpz_class>>& factors) {
    vector<mpz_class> pi_list;

    for (const auto& factor_ : factors) {
        mpz_class prime = factor_.first;
        mpz_class power = factor_.second;
        for (mpz_class i = 0; i < power; ++i) {
            pi_list.push_back(prime);
        }
    }

    return pi_list;
}
