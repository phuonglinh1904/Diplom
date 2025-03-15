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

void initialize_tables() {
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

std::map<mpz_class, std::pair<std::vector<mpz_class>, mpz_class>> find_cyclotomic_classes() {
    for (mpz_class i = 0; i <= n - 1; ++i) {
        bool is_in_class = false;

        // Kiểm tra xem i đã nằm trong một lớp chu kỳ nào đó chưa
        for (const auto &[key, values_length]: cyclotomic) {
            const auto &values = values_length.first;
            if (find(values.begin(), values.end(), i) != values.end()) {
                is_in_class = true;
                break;
            }
        }
        if (is_in_class) {
            continue;
        }

        // Tạo một lớp chu kỳ mới
        vector<mpz_class> new_class;
        mpz_class length = 0; // Độ dài của lớp chu kỳ
        mpz_class j = i;
        do {
            new_class.push_back(j);
            j = (j * 2) % n; // Lũy thừa theo modulo
            length++; // Tăng độ dài
        } while (j != i);

        // Lưu lớp chu kỳ và độ dài
        cyclotomic[i] = make_pair(new_class, length);
    }
    return cyclotomic;
}

// Hàm chia các lớp chu kỳ thành các nhóm có độ dài bằng nhau và lưu trữ theo thứ tự tăng dần
std::map<mpz_class, std::vector<std::vector<mpz_class>>> group_by_length_ordered(
    const std::map<mpz_class, std::pair<std::vector<mpz_class>, mpz_class>> &cyclotomic) {

    std::map<mpz_class, std::vector<std::vector<mpz_class>>> grouped;

    // Duyệt qua từng lớp chu kỳ
    for (const auto &[key, values_length] : cyclotomic) {
        const auto &values = values_length.first;  // Danh sách các giá trị trong lớp
        const auto &length = values_length.second; // Độ dài của lớp

        // Thêm lớp chu kỳ vào nhóm theo độ dài
        grouped[length].push_back(values);
    }

    return grouped;
}

// Hàm chọn hai phần tử ngẫu nhiên từ cùng một nhóm hoặc từ hai nhóm khác nhau
std::pair<mpz_class, mpz_class> random_two_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped) {
    std::srand(std::time(0)); // Khởi tạo seed ngẫu nhiên

    while (true) {
        // Ngẫu nhiên quyết định chọn cùng nhóm hoặc khác nhóm
        bool same_group = std::rand() % 2 == 0;

        if (same_group) {
            // Chọn một nhóm ngẫu nhiên
            auto it = grouped.begin();
            std::advance(it, std::rand() % grouped.size());
            const auto &group = it->second;

            // Chọn một lớp chu kỳ ngẫu nhiên trong nhóm
            const auto &cycle = group[std::rand() % group.size()];

            if (cycle.size() >= 2) {
                // Chọn hai phần tử ngẫu nhiên từ lớp chu kỳ
                size_t idx1 = std::rand() % cycle.size();
                size_t idx2;
                do {
                    idx2 = std::rand() % cycle.size();
                } while (idx1 == idx2);

                return {cycle[idx1], cycle[idx2]};
            }
        } else {
            // Chọn hai nhóm khác nhau ngẫu nhiên
            if (grouped.size() >= 2) {
                auto it1 = grouped.begin();
                auto it2 = grouped.begin();

                std::advance(it1, std::rand() % grouped.size());
                do {
                    it2 = grouped.begin();
                    std::advance(it2, std::rand() % grouped.size());
                } while (it1 == it2);

                const auto &group1 = it1->second;
                const auto &group2 = it2->second;

                const auto &cycle1 = group1[std::rand() % group1.size()];
                const auto &cycle2 = group2[std::rand() % group2.size()];

                mpz_class elem1 = cycle1[std::rand() % cycle1.size()];
                mpz_class elem2 = cycle2[std::rand() % cycle2.size()];

                return {elem1, elem2};
            }
        }
    }
}

// Hàm chọn ba phần tử ngẫu nhiên theo các trường hợp khác nhau
std::tuple<mpz_class, mpz_class, mpz_class> random_three_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped) {
    std::random_device rd;
    std::mt19937 gen(rd());

    while (true) {
        int case_type = gen() % 3; // Chọn ngẫu nhiên một trường hợp

        if (case_type == 0) { // Ba phần tử cùng một nhóm
            auto it = std::next(grouped.begin(), gen() % grouped.size());
            const auto &group = it->second;
            const auto &cycle = group[gen() % group.size()];

            if (cycle.size() >= 3) {
                std::vector<size_t> indices(cycle.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::shuffle(indices.begin(), indices.end(), gen);
                return {cycle[indices[0]], cycle[indices[1]], cycle[indices[2]]};
            }
        } else if (case_type == 1 && grouped.size() >= 3) { // Ba phần tử thuộc ba nhóm khác nhau
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            mpz_class elem1 = iters[0]->second[gen() % iters[0]->second.size()][0];
            mpz_class elem2 = iters[1]->second[gen() % iters[1]->second.size()][0];
            mpz_class elem3 = iters[2]->second[gen() % iters[2]->second.size()][0];

            return {elem1, elem2, elem3};
        } else if (case_type == 2 && grouped.size() >= 2) { // Hai phần tử cùng một nhóm, phần tử còn lại thuộc nhóm khác
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            const auto &group1 = iters[0]->second;
            const auto &group2 = iters[1]->second;
            const auto &cycle1 = group1[gen() % group1.size()];
            const auto &cycle2 = group2[gen() % group2.size()];

            if (cycle1.size() >= 2) {
                std::vector<size_t> indices(cycle1.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::shuffle(indices.begin(), indices.end(), gen);
                return {cycle1[indices[0]], cycle1[indices[1]], cycle2[gen() % cycle2.size()]};
            }
        }
    }
}

// Hàm chọn 4 phần tử ngẫu nhiên theo các trường hợp khác nhau
std::tuple<mpz_class, mpz_class, mpz_class, mpz_class> random_four_select(const std::map<mpz_class, std::vector<std::vector<mpz_class>>> &grouped) {
    std::random_device rd;
    std::mt19937 gen(rd());

    while (true) {
        int case_type = gen() % 5; // Chọn ngẫu nhiên một trường hợp

        if (case_type == 0 && grouped.size() >= 4) { // 4 phần tử thuộc 4 nhóm khác nhau
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            mpz_class elem1 = iters[0]->second[gen() % iters[0]->second.size()][0];
            mpz_class elem2 = iters[1]->second[gen() % iters[1]->second.size()][0];
            mpz_class elem3 = iters[2]->second[gen() % iters[2]->second.size()][0];
            mpz_class elem4 = iters[3]->second[gen() % iters[3]->second.size()][0];

            return {elem1, elem2, elem3, elem4};
        } else if (case_type == 1 && grouped.size() >= 3) { // 2 phần tử thuộc cùng 1 nhóm, 2 phần tử còn lại ở hai nhóm khác nhau
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            const auto &group1 = iters[0]->second;
            const auto &group2 = iters[1]->second;
            const auto &group3 = iters[2]->second;

            const auto &cycle1 = group1[gen() % group1.size()];
            const auto &cycle2 = group2[gen() % group2.size()];
            const auto &cycle3 = group3[gen() % group3.size()];

            if (cycle1.size() >= 2) {
                std::vector<size_t> indices1(cycle1.size());
                std::iota(indices1.begin(), indices1.end(), 0);
                std::shuffle(indices1.begin(), indices1.end(), gen);

                mpz_class elem1 = cycle1[indices1[0]];
                mpz_class elem2 = cycle1[indices1[1]];
                mpz_class elem3 = cycle2[gen() % cycle2.size()];
                mpz_class elem4 = cycle3[gen() % cycle3.size()];

                return {elem1, elem2, elem3, elem4};
            }
        } else if (case_type == 2 && grouped.size() >= 2) { // Hai phần tử thuộc cùng nhóm, 2 phần tử còn lại thuộc cùng nhóm khác
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            const auto &group1 = iters[0]->second;
            const auto &group2 = iters[1]->second;

            const auto &cycle1 = group1[gen() % group1.size()];
            const auto &cycle2 = group2[gen() % group2.size()];

            if (cycle1.size() >= 2 && cycle2.size() >= 2) {
                std::vector<size_t> indices1(cycle1.size());
                std::vector<size_t> indices2(cycle2.size());
                std::iota(indices1.begin(), indices1.end(), 0);
                std::iota(indices2.begin(), indices2.end(), 0);
                std::shuffle(indices1.begin(), indices1.end(), gen);
                std::shuffle(indices2.begin(), indices2.end(), gen);

                return {cycle1[indices1[0]], cycle1[indices1[1]], cycle2[indices2[0]], cycle2[indices2[1]]};
            }
        } else if (case_type == 3) { // 4 phần tử thuộc cùng 1 nhóm
            auto it = std::next(grouped.begin(), gen() % grouped.size());
            const auto &group = it->second;
            const auto &cycle = group[gen() % group.size()];

            if (cycle.size() >= 4) {
                std::vector<size_t> indices(cycle.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::shuffle(indices.begin(), indices.end(), gen);
                return {cycle[indices[0]], cycle[indices[1]], cycle[indices[2]], cycle[indices[3]]};
            }
        } else if (case_type == 4 && grouped.size() >= 2) { // 3 phần tử thuộc cùng nhóm, 1 phần tử thuộc nhóm khác
            std::vector<decltype(grouped.begin())> iters;
            for (auto it = grouped.begin(); it != grouped.end(); ++it) iters.push_back(it);
            std::shuffle(iters.begin(), iters.end(), gen);

            const auto &group1 = iters[0]->second;
            const auto &group2 = iters[1]->second;

            const auto &cycle1 = group1[gen() % group1.size()];
            const auto &cycle2 = group2[gen() % group2.size()];

            if (cycle1.size() >= 3) {
                std::vector<size_t> indices(cycle1.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::shuffle(indices.begin(), indices.end(), gen);
                return {cycle1[indices[0]], cycle1[indices[1]], cycle1[indices[2]], cycle2[gen() % cycle2.size()]};
            }
        }
    }
}



