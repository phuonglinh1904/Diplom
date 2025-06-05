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

vector<pair<mpz_class, mpz_class>> modified_trial_division(const vector<mpz_class> &P) {
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

vector<mpz_class> calculate_pi(const vector<pair<mpz_class, mpz_class>> &factors) {
    vector<mpz_class> pi_list;

    for (const auto &factor_: factors) {
        mpz_class prime = factor_.first;
        mpz_class power = factor_.second;
        for (mpz_class i = 0; i < power; ++i) {
            pi_list.push_back(prime);
        }
    }

    return pi_list;
}

mpz_class compute_T4(const mpz_class &D) {
    mpz_class trace = D;
    mpz_class max_i = (m - 2) / 2;
    mpz_class term = D;
    for (mpz_class i = 0; i < max_i; i++) {
        term = gf_power(term, 4); //  D^(2^(2*i))
        trace ^= term;
    }
    return trace;
}

mpz_class compute_a(const mpz_class &D, const mpz_class &q) {
    auto T4 = compute_T4(D);
    if (T4 == 0) {
        return 0;
    } else {

        return gf_power(2, (q - 1) / 3);
    }
}

mpz_class Tr(const mpz_class &y, const mpz_class &q) {
    mpz_class s = y;
    mpz_class current = y;

    for (mpz_class i = 0; i < q / 2; i++) {
        current = gf_power(current, 2); // y^(2^i)
        s ^= current;
    }
    return s;
}

mpz_class random_g(const mpz_class &q) {
    mpz_class g;
    do {
        g = 1 + rand() % (q.get_ui() - 1);
    } while (Tr(g, q) != 1);
    return g;
}

mpz_class compute_z1(mpz_class &D, const mpz_class &q) {
    if (m % 2 != 0) {
        mpz_class max_i = (m.get_ui() - 1) / 2;
        mpz_class result = D;
        mpz_class current = D;
        for (mpz_class i = 0; i < max_i; i ++) {
            current = gf_power(current, 4); // D^(2^i)
            result = gf_add(result, current);
        }
        return result;
    } else if (m % 4 == 2) {
        mpz_class D2 = gf_power(D, 2);         // D^2
        mpz_class sum_term = gf_power(D ^ D2, 4);
        mpz_class a = compute_a(D, q) ^ sum_term;
        mpz_class max_i = (m - 6) / 4;

        for (mpz_class i = 0; i < max_i; ++i) {
            sum_term = gf_power(sum_term, 16);  // (D + D^2)^{2^j}, j = 2 + 4i
            a ^= sum_term;
        }
        return a;
    } else if (m % 4 == 0) {
        auto T4 = compute_T4(D);
        mpz_class g, D1;
        if (T4 == 1) {
            g = 0;
            D1 = D;
        } else {
            g = random_g(q);
            mpz_class g2 = gf_multiply(g, g); // g^2
            D1 = D ^ g ^ g2;
        }
        // Bước 2: Tính a
        mpz_class a = 0;
        mpz_class max_j = (m / 4) - 1;

        for (mpz_class j = 1; j <= max_j; ++j) {
            mpz_class tmp1 = 2 * j - 1 + m / 2;
            mpz_class tmp2 = 2 * j - 2;

            mpz_class exp1, exp2;
            mpz_ui_pow_ui(exp1.get_mpz_t(), 2, tmp1.get_ui());
            mpz_ui_pow_ui(exp2.get_mpz_t(), 2, tmp2.get_ui());

            mpz_class c = exp1 + exp2;
            mpz_class term = gf_power(D1, c);
            a ^= term;
        }
        // Bước 3: Tính phần tổng theo I
        mpz_class sum_term = 1;
        for (mpz_class i = 0; i <= max_j; ++i) {
            mpz_class b = 2 * i + (m / 2);

            mpz_class exp;
            mpz_ui_pow_ui(exp.get_mpz_t(), 2, b.get_ui()); // exp = 2^b

            mpz_class term = gf_power(D1, exp);
            sum_term ^= term;
        }


        // Bước 4: D1^θ, với θ = 2^m - 1
        mpz_class theta = q.get_ui() - 1;
        mpz_class D1_theta = gf_power(D1, theta.get_ui());

        // Bước 5: Tính z1
        mpz_class z1 = g ^ a ^ gf_multiply(a, a); // g + a + a^2
        z1 ^= gf_multiply(D1_theta, sum_term);            // D1^θ * (1 + sum)

        return z1;

    } else {
        return -1;
    }

}
vector<Polynomial> solveDegree2(const Polynomial &poly, const mpz_class &q) {
    vector<Polynomial> roots;
    auto f_x = poly.normalize();
    auto D = gf_multiply(f_x.coeff[0], gf_inverse(gf_multiply(f_x.coeff[1], f_x.coeff[1])));
    auto z1 = compute_z1(D, q);
    auto z2 = z1 ^ 1;
    auto A = f_x.coeff[1];
    auto x1 = gf_multiply(A, z1);
    auto x2 = gf_multiply(A, z2);
    roots.push_back(Polynomial(vector<mpz_class>{x1, 1}));
    roots.push_back(Polynomial(vector<mpz_class>{x2, 1}));
    return roots;
}


mpz_class solve_cubic(const mpz_class &E, const mpz_class &q){
    auto D = gf_inverse(gf_multiply(E, E));
    auto w = compute_z1(D, q);
    auto v = gf_multiply(E, w);
    auto g0 = find_primitive_element();
    auto u = gf_cube_root(v, q, g0);
    if (u == -1) {
        w ^= 1;
        v = gf_multiply(E, w);
        u = gf_cube_root(v, q, g0);
    }
    mpz_class ui_inv = gf_inverse(u);
    mpz_class zi = gf_add(u, ui_inv);
    return zi;
}

vector<Polynomial> solveDegree3(const Polynomial &poly, const mpz_class &q){
    auto f_x = poly.normalize();
    auto A = f_x.coeff[2];
    auto B = f_x.coeff[1];
    auto C = f_x.coeff[0];
    auto E = gf_multiply(gf_add(C, gf_multiply(A, B)), gf_inverse(gf_square_root(gf_power(B ^ gf_multiply(A, A), 3), q)));
    auto z1 = solve_cubic(E, q);
    auto sqrt_a = gf_square_root(B ^ gf_multiply(A, A), q);
    auto y = gf_multiply(z1, sqrt_a);
    auto x = gf_add(A, y);
    Polynomial mono_root = Polynomial(vector<mpz_class>{x, 1});
    Polynomial poly_remain = Polynomial::div(poly, mono_root, q);
    auto roots_remain = solveDegree2(poly_remain, q);
    vector<Polynomial> roots;
    roots.push_back(mono_root);
    roots.push_back(roots_remain[0]);
    roots.push_back(roots_remain[1]);
    return roots;
}

vector<Polynomial> solveDegree4(const Polynomial &poly, const mpz_class &q){
    auto f_x = poly.normalize();
    auto A = f_x.coeff[3];
    auto B = f_x.coeff[2];
    auto C = f_x.coeff[1];
    auto D = f_x.coeff[0];
    auto a3 = D ^ gf_multiply(gf_multiply(B, C), gf_inverse(A)) ^ gf_power(gf_multiply(C, gf_inverse(A)),2)
                            ^ gf_multiply(C, gf_square_root(gf_multiply(C, gf_inverse(A)), q));
    auto a2 = B ^ gf_square_root(gf_multiply(A, C), q);
    auto a1 = A;
    auto E1 = gf_multiply(a1, gf_square_root(gf_multiply(a3, gf_inverse(gf_power(a2, 3))), q));
    auto E2 = gf_multiply(a3, gf_inverse(gf_power(a2, 2)));
    auto A1 = solve_cubic(E1, q);
    // auto A1 = A1_list[0];

    auto f_B1 = Polynomial(vector<mpz_class>{E2, 1 ^ gf_multiply(A1, A1), 1});
    auto B1_list = solveDegree2(f_B1, q);
    auto B1 = B1_list[0].get_coeff()[0];
    auto C1 = gf_multiply(E2, gf_inverse(B1));
    auto fx1 = Polynomial(vector<mpz_class>{B1, A1, 1});
    auto fx2 = Polynomial(vector<mpz_class>{C1, A1, 1});
    cout << "E1, E2: " << E1 << "; " << E2 << endl;
    cout << "Polinom is: " << Polynomial::mul(fx1, fx2) << endl;

    auto res1 = solveDegree2(fx1, q);
    auto res2 = solveDegree2(fx2, q);

    vector<Polynomial> roots;

    for (const auto& r1 : res1) {
        roots.push_back(r1);
    }

    for (const auto& r2 : res2) {
        roots.push_back(r2);
    }
    vector<Polynomial> roots_x;
    for (auto &root: roots) {
        mpz_class z = root.get_coeff()[0];
        mpz_class y = gf_multiply(z, gf_square_root(gf_multiply(a2, gf_inverse(a3)), q));
        mpz_class x = gf_inverse(y) ^ gf_square_root(gf_multiply(C, gf_inverse(A)), q);
        roots_x.push_back(Polynomial(vector<mpz_class>{x, 1}));
    }

    return roots_x;
}

