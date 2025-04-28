//
// Created by ADMIN23 on 09.12.2024.
//

#ifndef LONGBCH_POLYNOMIAL_H
#define LONGBCH_POLYNOMIAL_H

#include<iostream>
#include<vector>
#include <gmpxx.h>

using namespace std;

class Polynomial {
private:

    void prune();


    static std::pair<Polynomial, Polynomial> div_internal(const Polynomial &a, const Polynomial &b, const mpz_class& q);

public:

    vector<mpz_class> coeff;
    Polynomial() {};

    Polynomial(vector<mpz_class> coeff) : coeff(coeff) {};

    Polynomial(std::string s);

    Polynomial(const Polynomial &polynomial) : coeff(polynomial.coeff) {};

    Polynomial &operator=(const Polynomial &polynomial) = default;

    int get_degree() const;

    std::string to_string(std::string default_variable_name = "x") const;

    Polynomial diff() const;

    Polynomial zip(const int degree, const mpz_class& q) const;

    static Polynomial add(const Polynomial &a, const Polynomial &b);

    static Polynomial mul(const Polynomial &a, const Polynomial &b);

    static Polynomial div(const Polynomial &a, const Polynomial &b, const mpz_class& q);

    static Polynomial mod(const Polynomial &a, const Polynomial &b, const mpz_class& q);

    static Polynomial mul_alpha(const Polynomial &a, mpz_class lambda);

    // new function
    std::vector<mpz_class> get_coeff() const;

    static Polynomial get_one();

    static Polynomial get_random_polynomial(int max_degree, const mpz_class &q);

    static Polynomial get_random_polynomial_with_degree(int degree, const mpz_class &q);

    bool is_zero() const;

    bool is_one() const;

    friend bool operator==(const Polynomial &poly1, const Polynomial &poly2);

    friend bool operator<(const Polynomial &poly1, const Polynomial &poly2);

    friend bool operator!=(const Polynomial &poly1, const Polynomial &poly2);
    Polynomial normalize() const;
    static mpz_class calcPoly(const Polynomial &a, mpz_class x);
};


ostream &operator<<(ostream &Str, Polynomial const &v);

Polynomial gcd(const Polynomial &a, const Polynomial &b, const mpz_class& q);

Polynomial powmod(const Polynomial &a, mpz_class b, const Polynomial &mod, const mpz_class& q);



#endif //LONGBCH_POLYNOMIAL_H