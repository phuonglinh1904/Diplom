//
// Created by ADMIN23 on 09.12.2024.
//
#include "Algorthms.h"
#include "GFq.h"
#include "Tree.h"
#include <vector>
#include <cmath>
#include <chrono>
#include "Utils.h"

using namespace std;
int cnt[2] = {0, 0};

vector<pair<Polynomial, int>> square_free_decomposition(const Polynomial &poly, const mpz_class &q) {
    vector<pair<Polynomial, int>> result;

    Polynomial a = poly;
    Polynomial c;
    int m = 1;

    do {
        auto ad = a.diff();
        c = gcd(a, ad, q);
        auto w = Polynomial::div(a, c, q);
        int i = 1;
        while (!w.is_one()) {
            auto y = gcd(w, c, q);
            auto qq = Polynomial::div(w, y, q);
            if (!qq.is_one()) {
                result.push_back({qq, i * m});
            }
            w = y;
            c = Polynomial::div(c, y, q);
            i++;
        }
        if (!c.is_one()) {
            a = c.zip(2, q);
            m = m * 2;
        }
    } while (!c.is_one());

    return result;
}

vector<pair<Polynomial, int>> distinct_degree_factorization(const Polynomial &poly, const mpz_class &q) {
    vector<pair<Polynomial, int>> result;

    Polynomial h("1x^1");
    Polynomial href = h;
    auto f = poly;
    int i = 1;

    while (!f.is_one()) {
        h = powmod(h, q, poly, q);
        auto sbs = Polynomial::add(h, href);
        auto g = gcd(sbs, f, q);
        if (!g.is_one()) {
            result.push_back({g, i});
            f = Polynomial::div(f, g, q);
        }
        i++;
    }

    return result;
}

Polynomial gcd_with_count(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    cnt[0]++;  // Đếm số lần gọi GCD
    return gcd(a, b, q);
}

// Hàm bọc quanh phép chia để đếm số lần gọi
Polynomial div_with_count(const Polynomial &a, const Polynomial &b, const mpz_class &q) {
    cnt[1]++;  // Đếm số lần gọi phép chia
    return Polynomial::div(a, b, q);
}

Polynomial trace(const Polynomial &y, const Polynomial &g, const int deg, const mpz_class q) {
    Polynomial s = y;
    Polynomial current = y;

    for (int i = 1; i < deg; i++) {
        current = powmod(current, q, g, q); // Tính y^(2^i) mod g y^(2^(m*deg - 1))
        s = Polynomial::mod(Polynomial::add(s, current), g, q); // Cộng vào truy vết
    }

    return s;
}

Polynomial T_m(const Polynomial &y, const Polynomial &g, const int deg, const mpz_class q) {
    Polynomial s = y;
    Polynomial current = y;

    for (int i = 1; i < deg; i++) {
        current = powmod(current, 2, g, q); // Tính y^(2^i) mod g y^(2^(m*deg - 1))
        s = Polynomial::mod(Polynomial::add(s, current), g, q); // Cộng vào truy vết
    }

    return s;
}

void Refine(vector<Polynomial> &U, const Polynomial &v, const Polynomial &f, const mpz_class &q) {
    vector<Polynomial> U_new;
    for (const auto &u: U) {
        // Gọi hàm GCD có đếm số lần gọi
        Polynomial gcd_uv = gcd_with_count(u, v, q);

        if (gcd_uv.is_one() || gcd_uv == u) { // Nếu gcd(u, v) = 1, giữ nguyên u
            U_new.push_back(u);
        } else {
            // Gọi hàm chia có đếm số lần gọi
            Polynomial u_div = div_with_count(u, gcd_uv, q);
            U_new.push_back(gcd_uv);
            U_new.push_back(u_div);
        }
    }
    U = U_new;
}


vector<vector<mpz_class>> create_polynomial_h(const Polynomial &g, int d, const mpz_class &q) {
    vector<vector<mpz_class>> h = {{0, 1},
                                   {1}};
    Polynomial lambda = Polynomial("x");
    for (int i = 1; i < d; i++) {
        // tính x^q
        lambda = powmod(lambda, q, g, q);
        vector<vector<mpz_class>> temp = h;
        h.insert(h.begin(), {0});
        for (int i = 0; i < temp.size(); i++) {
            temp[i] = Polynomial::mod(Polynomial::mul(lambda, Polynomial(temp[i])), g, q).get_coeff();
            h[i] = Polynomial::mod(Polynomial::add(h[i], temp[i]), g, q).get_coeff();
        }
    }
    return h;
}

vector<Polynomial> equal_degree_factorization_by_Shoup(const Polynomial &poly, int degree, const mpz_class &q) {
    vector<Polynomial> U;
    U.push_back(poly);  // Khởi tạo tập chứa f

    auto h = create_polynomial_h(poly, degree, q);
    h.pop_back();  // Loại bỏ phần tử cuối
    vector<Polynomial> separatingSet;
    for (auto &elem: h) {
        if (!Polynomial(elem).is_zero())
            separatingSet.push_back(Polynomial(elem));
    }

    const int max_attempts = 10000;
    int attempts_loop = 0;
    int z = 0;
    while (U.size() < poly.get_degree() / degree) {
        if (++attempts_loop > max_attempts) {
            std::cerr << "Shoup error! " << max_attempts << " attempts." << std::endl;
            return {};
        }
        for (size_t l = 0; l < separatingSet.size(); ++l) {
            Polynomial h_l = separatingSet[l];
            Polynomial x_z(vector<mpz_class>{z}); // Tạo đa thức bậc 0 từ z

            // Bước 1: Refine với h_l + x_z
            Polynomial refined_poly_1 = Polynomial::add(h_l, x_z);
            Refine(U, refined_poly_1, poly, q);

            // Bước 2: Refine với (h_l + x_z)^{(q-1)/2} - 1
            Polynomial exponentiated = powmod(refined_poly_1, (q - 1) / 2, poly, q);
            Polynomial refined_poly_2 = Polynomial::add(exponentiated, Polynomial(vector<mpz_class>{1}));
            Refine(U, refined_poly_2, poly, q);
        }
        z++; // Tăng giá trị của z sau mỗi vòng lặp
    }

    return U;
}




vector<Polynomial> factor_square_free_by_Shoup(const Polynomial &poly, const mpz_class &q) {
    vector<Polynomial> result;

    auto distDeg = distinct_degree_factorization(poly, q);

    for (auto t: distDeg) {

        auto r1 = equal_degree_factorization_by_Shoup(t.first, t.second, q);
        if (r1.empty()) {
            return {};
        }
        result.insert(result.end(), r1.begin(), r1.end());
    }

    return result;
}



vector<pair<Polynomial, int>> factor_by_Shoup_Algorithm(const Polynomial &poly, const mpz_class &q, int &time) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vector<pair<Polynomial, int>> result;

    vector<pair<Polynomial, int>> sqrfree = square_free_decomposition(poly, q);

    for (auto const &value: sqrfree) {
        auto r1 = factor_square_free_by_Shoup(value.first, q);
        if (r1.empty()) {
            return {};
        }
        for (size_t i = 0; i < r1.size(); i++) {
            result.push_back({r1[i], value.second});
        }
    }

    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // save time
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    cout << "Computation time = " << time << "[ms]" << std::endl;

    return result;
}


//=======================================================================================================
vector<pair<Polynomial, mpz_class>> find_root_by_Mignotte(const Polynomial& f, const mpz_class& q, int &time, size_t &total_memory) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mpz_class chi = q - 1;

    vector<Polynomial> hs;
    hs.push_back(Polynomial("x"));
    total_memory += hs[0].memory_usage();
    for (size_t i = 0; i < pi.size() - 1; ++i) {
        Polynomial hi = powmod(hs[i], pi[i], f, q);
        hs.push_back(hi);
        total_memory += hi.memory_usage();
    }
    vector<pair<Polynomial, mpz_class>> F = {make_pair(f, 0)};
    total_memory += f.memory_usage();

    for (int i = pi.size() - 1; i >= 0; i--) {
        vector<Polynomial> polynomials;
        for (const auto& pair : F) {
            polynomials.push_back(pair.first);
            total_memory += pair.first.memory_usage();
        }
        auto root = buildSubproductTree(polynomials, 0, polynomials.size() - 1);
        vector<Polynomial> remainders;
        computeRemainders(hs[i], root, remainders, q);
        total_memory += remainders.back().memory_usage() * remainders.size();
        vector<pair<Polynomial, mpz_class>> G;
        for (mpz_class j = 0; j < pi[i]; ++j) {
            int index = 0;
            for (const auto& pair : F) {
                const Polynomial& g = pair.first;

                const mpz_class& e = pair.second;
                mpz_class v = (e + j*chi) / pi[i];
                auto eps = deg_alpha[v.get_ui()];
                Polynomial term = Polynomial({eps});
                Polynomial subtracted = Polynomial::add(remainders[index], term);

                Polynomial gj = gcd(subtracted, g, q);
                if (gj.get_degree() >= 1) {
                    auto gj_ = gj.normalize();
                    G.push_back({gj_, v});
                    total_memory += gj_.memory_usage();
                }
                index++;
            }
        }
        size_t G_memory = 0;
        for (const auto& pair : G) {
            G_memory += pair.first.memory_usage();
        }
        total_memory += G_memory;
        F = G;
        if (F.size() == f.get_degree()) break;
    }

    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    cout << "Computation time = " << time << "[ms]" << std::endl;
    cout << "Total memory used in find_root_by_Mignotte: " << total_memory << " bytes" << std::endl;
    return F;
}
std::vector<Polynomial>
separate(const Polynomial &g, const Polynomial &s_trace, int degree, const mpz_class &q, size_t &total_memory) {
    if (g.get_degree() == degree) {
        return {g};
    }

    Polynomial s = Polynomial::mod(s_trace, g, q);
    total_memory += s.memory_usage();
    if (s.is_zero()) {
        auto result = equal_degree_factorization_by_Ben_Or(g, degree, q, total_memory);
        size_t result_memory = 0;
        for (const auto &res : result) {
            result_memory += res.memory_usage();
        }
        total_memory += result_memory;
        return result;
    }

    const int max_attempts = 20;
    int attempts_loop = 0;

    while (true) {
        if (++attempts_loop > max_attempts) {
            return equal_degree_factorization_by_Ben_Or(g, degree, q, total_memory);
        }
        mpz_class delta = 1 + rand() % (q.get_ui() - 1);
        Polynomial delta_s = Polynomial::mul_alpha(s, delta);
        total_memory += delta_s.memory_usage();
        // T(delta * s(x))
        Polynomial T = T_m(delta_s, g, m.get_ui(), q);
        total_memory += T.memory_usage();
        Polynomial g1 = gcd_with_count(g, T, q);
        total_memory += g1.memory_usage();
        if (!g1.is_one() && g1.get_degree() != g.get_degree()) {
            auto g2 = div_with_count(g, g1, q);
            total_memory += g2.memory_usage();
            auto result1 = separate(g1, s, degree, q, total_memory);
            if (result1.empty()) {
                return {};
            }
            auto result2 = separate(g2, s, degree, q, total_memory);
            if (result2.empty()) {
                return {};
            }
            result1.insert(result1.end(), result2.begin(), result2.end());
            size_t result_memory = 0;
            for (const auto &res : result1) {
                result_memory += res.memory_usage();
            }
            total_memory += result_memory;
            return result1;
        }
    }
}

vector<Polynomial>
equal_degree_factorization_by_Ben_Or(const Polynomial &poly, int degree, const mpz_class &q, size_t &total_memory) {
    if (poly.get_degree() == degree) {
        return vector<Polynomial>({poly});
    }

    Polynomial y = Polynomial::get_random_polynomial(poly.get_degree() - 1, q.get_ui());
    total_memory += y.memory_usage();
    Polynomial s = trace(y, poly, degree, q);
    total_memory += s.memory_usage();
    if (s.get_degree() == 0) {
        return equal_degree_factorization_by_Ben_Or(poly, degree, q, total_memory);
    }

    auto result = separate(poly, s, degree, q, total_memory);
    size_t result_memory = 0;
    for (const auto &res : result) {
        result_memory += res.memory_usage();
    }
    total_memory += result_memory;


    return result;
}

vector<Polynomial> find_root_by_Ben_Or(const Polynomial &poly, const mpz_class &q, int &time, size_t &total_memory) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    total_memory += poly.memory_usage();
    Polynomial h("x");
    total_memory += h.memory_usage();
    Polynomial href = h;
    total_memory += href.memory_usage();
    h = powmod(h, q, poly, q);
    total_memory += h.memory_usage();
    auto sbs = Polynomial::add(h, href);
    total_memory += sbs.memory_usage();
    auto poly_with_root = gcd(sbs, poly, q);
    total_memory += poly_with_root.memory_usage();
    vector<Polynomial> polys = {};
    if (poly_with_root.get_degree()) {
        try {
            polys = equal_degree_factorization_by_Ben_Or(poly_with_root, 1, q, total_memory);
            size_t polys_memory = 0;
            for (const auto &p : polys) {
                polys_memory += p.memory_usage();
            }
            total_memory += polys_memory;
        } catch (const std::exception &e) {
                cerr << "Exception occurred: " << e.what() << endl;
            polys.clear();
        }
    }

    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // save time
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    cout << "Computation time = " << time << "[ms]" << std::endl;
    cout << "Total memory used in find_root_by_Ben_Or: " << total_memory << " bytes" << endl;
    if (polys.empty()) return {};
    return polys;
}


// =========================================================================================================================== //
std::vector<Polynomial>
separate_hybird(const Polynomial &g, const Polynomial &s_trace, int degree, int &cnt_failed, const mpz_class &q, size_t &total_memory) {
    if (g.get_degree() == degree) {
        return {g};
    }

    Polynomial s = Polynomial::mod(s_trace, g, q);
    total_memory += s.memory_usage();
    if (s.is_zero()) {
        auto result = equal_degree_factorization_by_hybird(g, degree, q, total_memory);
        size_t result_memory = 0;
        for (const auto &res : result) {
            result_memory += res.memory_usage();
        }
        total_memory += result_memory;
        std::cout << "Total memory used in separate_hybird (equal_degree_factorization): " << total_memory << " bytes" << std::endl;
        return result;
    }

    const int max_failed = 5;
    const int max_attempts = 10;
    int attempts_loop = 0;

    if (cnt_failed < max_failed) {
        while (true) {
            if (++attempts_loop > max_attempts) {
                int time;
                auto results = find_root_by_Mignotte(g, q, time, total_memory);
                vector<Polynomial> polys = {};
                for (auto &res : results) {
                    polys.push_back(res.first);
                }
                size_t polys_memory = 0;
                for (const auto &p : polys) {
                    polys_memory += p.memory_usage();
                }
                total_memory += polys_memory;

                cnt_failed += 1;
                return polys;
            }
            mpz_class delta = 1 + rand() % (q.get_ui() - 1);
            Polynomial delta_s = Polynomial::mul_alpha(s, delta);
            total_memory += delta_s.memory_usage();
            // T(delta * s(x))
            Polynomial T = T_m(delta_s, g, m.get_ui(), q);
            total_memory += T.memory_usage();
            Polynomial g1 = gcd_with_count(g, T, q);
            total_memory += g1.memory_usage();
            if (!g1.is_one() && g1.get_degree() != g.get_degree()) {
                auto g2 = div_with_count(g, g1, q);
                total_memory += g2.memory_usage();
                auto result1 = separate_hybird(g1, s, degree, cnt_failed, q, total_memory);
                if (result1.empty()) {
                    return {};
                }
                auto result2 = separate_hybird(g2, s, degree, cnt_failed, q, total_memory);
                if (result2.empty()) {
                    return {};
                }
                result1.insert(result1.end(), result2.begin(), result2.end());
                size_t result_memory = 0;
                for (const auto &res : result1) {
                    result_memory += res.memory_usage();  // Đo bộ nhớ của kết quả
                }
                total_memory += result_memory;
                return result1;
            }
        }
    }
    else {
        int time;
        auto results = find_root_by_Mignotte(g, q, time, total_memory);
        vector<Polynomial> polys = {};
        for (auto &res : results) {
            polys.push_back(res.first);
        }
        size_t polys_memory = 0;
        for (const auto &p : polys) {
            polys_memory += p.memory_usage();
        }
        total_memory += polys_memory;
        return polys;
    }

}

vector<Polynomial>
equal_degree_factorization_by_hybird(const Polynomial &poly, int degree, const mpz_class &q, size_t &total_memory) {
    int cnt_failed = 0;
    if (poly.get_degree() == degree) {
        return vector<Polynomial>({poly});
    }

    Polynomial y = Polynomial::get_random_polynomial(poly.get_degree() - 1, q.get_ui());
    total_memory += y.memory_usage();

    Polynomial s = trace(y, poly, degree, q);
    total_memory += s.memory_usage();
    if (s.get_degree() == 0) {
        return equal_degree_factorization_by_hybird(poly, degree, q, total_memory);
    }

    vector<Polynomial> result = separate_hybird(poly, s, degree, cnt_failed, q, total_memory);
    size_t result_memory = 0;
    for (const auto& p : result) {
        result_memory += p.memory_usage();
    }
    total_memory += result_memory;
    return result;
}

vector<Polynomial> find_root_by_hybrid(const Polynomial &poly, const mpz_class &q, int &time, size_t &total_memory) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    vector<Polynomial> polys = {};
    total_memory += poly.memory_usage();
    if (poly.get_degree()) {
        try {
            polys = equal_degree_factorization_by_hybird(poly, 1, q, total_memory);
            // Đo bộ nhớ của mảng polys
            size_t polys_memory = 0;
            for (const auto& p : polys) {
                polys_memory += p.memory_usage();
            }
            total_memory += polys_memory;
        } catch (const std::exception &e) {
            cerr << "Exception occurred: " << e.what() << endl;
            polys.clear();
        }
    }

    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // save time
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    cout << "Computation time = " << time << "[ms]" << std::endl;
    cout << "Total memory used in find_root_by_hybrid: " << total_memory << " bytes" << endl;
    if (polys.empty()) return {};
    return polys;
}

//========================================================================
vector<mpz_class> find_root_by_chien(const Polynomial &poly, int &time) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    vector<mpz_class> roots;
    for (mpz_class i = 0; i < n; ++i) {
        if (Polynomial::calcPoly(poly, gf_inverse(deg_alpha[i.get_ui()])) == 0) {
           roots.push_back(i);
        }
    }
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    return roots;
}
//=======================================================================================================