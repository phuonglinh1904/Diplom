//
// Created by ADMIN23 on 23.03.2025.
//
#include <random>
#include <chrono>
#include "BCH.h"

map<mpz_class, vector<mpz_class>> find_cyclotomic_classes() {
    for (mpz_class i = 0; i < n - 1; ++i) {
        bool is_in_class = false;
        for (const auto &[key, values]: cyclotomic) {
            if (find(values.begin(), values.end(), i) != values.end()) {
                is_in_class = true;
                break;
            }
        }
        if (is_in_class) {
            continue;
        }

        vector<mpz_class> new_class;
        mpz_class j = i;
        do {
            new_class.push_back(j);
            j = (j << 1) % n;
        } while (j != i);

        cyclotomic[i] = new_class;
    }
    return cyclotomic;
}

Polynomial compute_g() {
    vector<Polynomial> minimals;
    for (mpz_class i = 1; i < delta; ++i) {
        auto it = cyclotomic.find(i.get_ui());
        if (it != cyclotomic.end()) {
            Polynomial minimal = Polynomial(vector<mpz_class>{1});
            for (mpz_class idx: it->second) {
                vector<mpz_class> cur_multiplier = {deg_alpha[idx.get_ui()], 1};
                minimal = Polynomial::mul(minimal, cur_multiplier);
            }

            minimals.push_back(minimal);
        }
    }
    Polynomial g = Polynomial(vector<mpz_class>{1});
    for (const auto &minimal: minimals) {
        g = Polynomial::mul(g, minimal);
    }

    return g;
}

std::vector<mpz_class> encode(const std::vector<mpz_class>& a, const Polynomial& g, const mpz_class& q) {
    unsigned int n_int = n.get_ui();
    std::vector<mpz_class> c(g.get_degree() + 1);
    c[g.get_degree()] = 1;
    Polynomial c_poly(c);
    Polynomial a_poly(a);
    // c = x^(n-k) * a(x)
    auto res = Polynomial::mul(c_poly, a_poly);
    // c += (x^(n-k) * a(x)) mod g(x)
    res = Polynomial::add(res, Polynomial::mod(res, g, q));
    vector<mpz_class> result_coeffs = res.get_coeff();
    result_coeffs.resize(n_int);
    return result_coeffs;
}

//vector<mpz_class> decode(const vector<mpz_class> &v, const mpz_class &q, bool is_hybrid, int &time) {
//    Polynomial poly = Polynomial(v);
//    vector<mpz_class> syndromes(delta.get_ui() - 1);
//    for (int i = 0; i < syndromes.size(); ++i) {
//        syndromes[i] = Polynomial::calcPoly(poly, deg_alpha[i + 1]);
//    }
//    int mm = 0, l = 0;
//    Polynomial Lambda = Polynomial(vector<mpz_class>{1});
//    Polynomial B = Polynomial(vector<mpz_class>{1});
//    Polynomial T;
//    mpz_class dr;
//    int r = 1;
//    while (r < delta.get_ui()) {
//        dr = 0;
//        for (int j = 0; j <= l; ++j) {
//            dr ^= gf_multiply(Lambda.coeff[j], syndromes[r - 1 - j]);
//        }
//        if (dr != 0) {
//            T = B;
//            for (int i = 0; i < r - mm; ++i) {
//                T = Polynomial::mul(T, Polynomial("x"));
//            }
//            T  = Polynomial::mul_alpha(T, dr);
//            T = Polynomial::add(Lambda, T);
//            if (2 * l < r) {
//                B = Lambda;
//                B = Polynomial::mul_alpha(B, gf_inverse(dr));
//                l = r - l;
//                mm = r;
//            }
//            Lambda = T;
//
//        }
//        r += 2;
//    }
//    auto decoded = v;
//    if(is_hybrid){
//        auto time_hybrid = 0;
//        auto roots = find_root_by_hybrid(Lambda.normalize(), q, time_hybrid);
//        time = time_hybrid;
//        std::vector<unsigned long> log_vals(roots.size());
//        for (size_t i = 0; i < roots.size(); ++i) {
//            log_vals[i] = static_cast<unsigned long>(log_alpha[roots[i].coeff[0].get_ui()].get_ui());
//        }
//
//        for (size_t i = 0; i < roots.size(); ++i) {
//            auto log_val = log_vals[i];
//            mpz_class location = n - log_val;
//            unsigned long idx = location.get_ui();
//            decoded[idx] ^= 1;
//        }
//    } else{
//        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//        for (int i = 0; i < n; ++i) {
//            if (Polynomial::calcPoly(Lambda, gf_inverse(deg_alpha[i])) == 0) {
//                decoded[i] ^= 1;
//            }
//        }
//        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//        time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
//    }
//    return decoded;
//
//}
//
