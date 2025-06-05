//
// Created by ADMIN23 on 09.12.2024.
//
#include "Polynomial.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include "GFq.h"

using namespace std;

Polynomial::Polynomial(std::string str) {
    try {
        auto monoms_string = split_by(str, "+");
        vector<pair<mpz_class, int>> monoms(monoms_string.size());

        for (size_t i = 0; i < monoms.size(); i++) {
            auto ss = split_by(monoms_string[i], "x");

            mpz_class coeff;
            int pwr;

            if (ss.size() == 1) {
                pwr = 0;
            } else {
                if (ss[1] == "") {
                    pwr = 1;
                } else {
                    pwr = stoi(ss[1].substr(1, ss[1].size()));
                }
            }

            if (ss[0] == "") {
                coeff = 1;
            } else {
                coeff = mpz_class(ss[0]);
            }
            monoms[i] = {coeff, pwr};
        }

        int pwr = monoms[0].second;
        vector<mpz_class> result(pwr + 1);

        for (auto pr : monoms) {
            result[pr.second] = pr.first;
        }

        this->coeff = result;
    } catch (...) {
        cout << "Error during polynomial parse" << endl;
    }
}


void Polynomial::prune() {
    while (!this->coeff.empty() && this->coeff.back() == 0) {
        this->coeff.pop_back();
    }
}

int Polynomial::get_degree() const {
    return this->coeff.size() - 1;
}

std::string Polynomial::to_string(std::string default_variable_name) const {
    int degree = this->coeff.size() - 1;
    std::string result = "";

    for (int i = degree; i >= 0; i--) {
        if (this->coeff[i] == 0) continue;

        if (!result.empty()) {
            result += "+";
        }

        if (this->coeff[i] != 1 || i == 0) {
            result += this->coeff[i].get_str(); // Chuyển số lớn sang chuỗi
        }
        if (i != 0) {
            result += default_variable_name;
            if (i != 1) {
                result += "^" + std::to_string(i);; // Chuyển chỉ số lớn sang chuỗi
            }
        }
    }

    if (result.empty()) {
        return "0";
    }

    return result;
}

Polynomial Polynomial::diff() const {
    std::vector<mpz_class> v(this->coeff.size() - 1, 0);

    for (int i = 1; i < this->coeff.size(); ++i) {
        if (i % 2 == 1) {
            for (mpz_class j = 0; j < i; ++j) {
                v[i - 1] ^= this->coeff[i];
            }
        }
    }

    auto d = Polynomial(v);
    d.prune();

    return d;
}

Polynomial Polynomial::zip(const int degree, const mpz_class& q) const {
    vector<mpz_class> v(this->get_degree() / degree + 1);

    for (size_t i = 0; i <= this->get_degree(); i += degree) {
        v[i / degree] = gf_power(this->coeff[i], q.get_ui() / 2);
    }

    return Polynomial(v);
}

Polynomial Polynomial::add(const Polynomial &a, const Polynomial &b) {
    if (a.is_zero()) {
        return b;
    }
    if (b.is_zero()) {
        return a;
    }

    std::vector<mpz_class> v(std::max(a.coeff.size(), b.coeff.size()), 0);

    for (size_t i = 0; i < v.size(); ++i) {
        mpz_class ai = (i < a.coeff.size()) ? a.coeff[i] : 0;
        mpz_class bi = (i < b.coeff.size()) ? b.coeff[i] : 0;
        v[i] = ai ^ bi; // XOR hai số lớn
    }

    Polynomial p = Polynomial(v);
    p.prune();
    return p;
}

Polynomial Polynomial::mul(const Polynomial &a, const Polynomial &b) {
    if (a.is_zero() || b.is_zero()) {
        return Polynomial(vector<mpz_class>{});
    }

    vector<mpz_class> v(a.get_degree() + b.get_degree() + 1);

    for (size_t i = 0; i <= a.get_degree(); ++i) {
        for (size_t j = 0; j <= b.get_degree(); ++j) {
            v[i + j] ^= gf_multiply(a.coeff[i], b.coeff[j]);
        }
    }

    return Polynomial(v);
}

pair<Polynomial, Polynomial> Polynomial::div_internal(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    Polynomial r = a;
    r.prune();
    Polynomial normalized_b = b.normalize();
    normalized_b.prune();

    if (normalized_b.is_zero()) {
        throw invalid_argument("Division by zero polynomial in GF(2^n)");
    }

    mpz_class degree_of_result = r.get_degree() - normalized_b.get_degree() + 1;

    if (degree_of_result < 1) {
        return {Polynomial(std::vector<mpz_class>{}), r};
    }

    std::vector<mpz_class> quotient(degree_of_result.get_ui(), 0);

    while (r.get_degree() >= normalized_b.get_degree()) {
        size_t pos = (r.get_degree() - normalized_b.get_degree());
        mpz_class coeff_div = gf_inverse(normalized_b.coeff.back());
        coeff_div = gf_multiply(r.coeff.back(), coeff_div);
        quotient[pos] = coeff_div;

        for (size_t i = 0; i <= normalized_b.get_degree(); ++i) {
            r.coeff[pos + i] ^= gf_multiply(normalized_b.coeff[i], coeff_div);
        }

        r.prune();
    }

    return {Polynomial(quotient), r};
}


Polynomial Polynomial::div(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    return Polynomial::div_internal(a, b, q).first;
}


Polynomial Polynomial::mod(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    return Polynomial::div_internal(a, b, q).second;
}

Polynomial Polynomial::get_one() {
    return Polynomial(vector<mpz_class>{1});
}

bool operator==(const Polynomial &poly1, const Polynomial &poly2) {
    return poly1.coeff == poly2.coeff;
}


bool operator!=(const Polynomial &poly1, const Polynomial &poly2) {
    return !(poly1 == poly2);
}

bool operator<(const Polynomial &poly1, const Polynomial &poly2) {
    if (poly1.coeff.size() != poly2.coeff.size()) {
        return poly1.coeff.size() < poly2.coeff.size();
    }

    for (size_t i = poly1.coeff.size(); i-- > 0;) {
        if (poly1.coeff[i] != poly2.coeff[i]) {
            return poly1.coeff[i] < poly2.coeff[i];
        }
    }
    return false;
}

Polynomial Polynomial::normalize() const {
    mpz_class bl = this->coeff.back();
    mpz_class ib = gf_inverse(bl);

    vector<mpz_class> v = this->coeff;

    for (int i = 0; i < v.size(); i++) {
        v[i] = gf_multiply(v[i], ib);

    }

    return Polynomial(v);
}

Polynomial gcd(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    if (a.is_zero()) {
        return b;
    }
    if (b.is_zero()) {
        return a;
    }

    Polynomial temp_a = a.normalize();
    Polynomial temp_b = b.normalize();
    while (!temp_b.is_zero()) {
        Polynomial z = Polynomial::mod(temp_a, temp_b, q);;
        temp_a = temp_b;
        temp_b = z;
    }

    return temp_a.normalize();
}

Polynomial powmod(const Polynomial &a, mpz_class b, const Polynomial &mod, const mpz_class &q) {
    mpz_class power = b;
    Polynomial rez = Polynomial::get_one();
    Polynomial aa = a;
    while (power > 0) {
        if (power % 2 == 1) {
            rez = Polynomial::mul(rez, aa);
            rez = Polynomial::mod(rez, mod, q);
        }
        aa = Polynomial::mul(aa, aa);
        aa = Polynomial::mod(aa, mod, q);
        power /= 2;
    }
    return rez;
}

bool Polynomial::is_zero() const {
    return this->coeff.size() == 0;
}


bool Polynomial::is_one() const {
    return this->coeff.size() == 1 && this->coeff[0] == 1;
}

Polynomial Polynomial::get_random_polynomial(int max_degree, const mpz_class &q) {
    long long poly_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(poly_seed);

    uniform_int_distribution<int> dist(0, q.get_ui() - 1);
    uniform_int_distribution<int> dist_degree(0, max_degree - 1);

    auto degree = dist_degree(rng) + 1;

    vector<mpz_class> vr(degree + 1);
    vr[degree] = 1;

    for (size_t i = 0; i < degree; i++) {
        vr[i] = dist(rng);
    }

    Polynomial result(vr);

    result.prune();

    return result;
}

Polynomial Polynomial::get_random_polynomial_with_degree(int degree, const mpz_class &q) {
    long long poly_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 rng(poly_seed);

    std::uniform_int_distribution<int> dist(0, q.get_ui() - 1);

    std::vector<mpz_class> vr(degree + 1);
    vr[degree] = 1; // Đảm bảo hệ số bậc cao nhất khác 0

    for (size_t i = 0; i < degree; i++) {
        vr[i] = dist(rng);
    }

    return Polynomial(vr);
}

ostream &operator<<(ostream &strm, const Polynomial &poly) {
    return strm << poly.to_string();
}

Polynomial Polynomial::mul_alpha(const Polynomial &a, mpz_class lambda) {
    vector<mpz_class> b(a.get_degree() + 1);
    for (int i = 0; i <= a.get_degree(); i++) {
        if (a.coeff[i] != 0) {
            mpz_class index = log_alpha[a.coeff[i].get_ui()] + log_alpha[lambda.get_ui()];
            b[i] = deg_alpha[index.get_ui()];
        }
        else
            b[i] = 0;
    }
    return Polynomial(b);
}

std::vector<mpz_class> Polynomial::get_coeff() const {
    return coeff;
}

mpz_class Polynomial :: calcPoly(const Polynomial &a, mpz_class x) {
    auto res = a.coeff.back();
    for (int j = a.get_degree() - 1; j >= 0 && j <= a.get_degree(); --j) {
        res = a.coeff[j] ^ gf_multiply(x, res);
    }
    return res;
}