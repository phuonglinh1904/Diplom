//
// Created by ADMIN on 27-Dec-24.
//
#include <fstream>
#include "testing.h"
#include "Algorthms.h"
#include "GFq.h"
#include "Utils.h"

Polynomial random_polynomial_factor_mono_degree(int degree, const mpz_class &q) {
    Polynomial poly_res = Polynomial::get_one();
    vector<Polynomial> polys;
    for (int i = 0; i < degree; i++) {
        Polynomial poly_degree_1 = Polynomial::get_random_polynomial_with_degree(1, q);
        // kiểm tra nếu chưa có trong polys thì push vào
        if (std::find(polys.begin(), polys.end(), poly_degree_1) == polys.end()) {
            polys.push_back(poly_degree_1);
            poly_res = Polynomial::mul(poly_res, poly_degree_1);
        }
        else {
            i--;
        }

    }
    return poly_res;
}

Polynomial random_polynomial_factor_equal_degree(int degree, const mpz_class &q) {
    Polynomial poly_res = Polynomial::get_one();
    for (int i = 0; i < degree / 5; i++) {
        Polynomial poly_degree_1 = Polynomial::get_random_polynomial_with_degree(1, q);
        poly_res = Polynomial::mul(poly_res, poly_degree_1);
    }
    int degree_remain = degree - poly_res.get_degree();
    Polynomial poly_remain = Polynomial::get_random_polynomial_with_degree(degree_remain, q);
    poly_res = Polynomial::mul(poly_res, poly_remain);
    return poly_res;
}


void compare_three_algo(const std::string &filename, const mpz_class &q) {
    std::ofstream file(filename);
    file << "n,Cantor,Ben-Or,Shoup\n";
    vector<int> degrees = {10,24,56,104,252};

    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_mono_degree(degrees[i], q);
        //poly = random_polynomial_factor_equal_degree(degrees[i], q);
        cout << "Polynomial with degree: " << degrees[i] << endl;
        int time_run_1 = 0, time_run_2 = 0, time_run_3 = 0;
        auto res_1 = factor(false, poly, q, time_run_1);
        if (res_1.empty()) {
            // i--;
            // continue;
            time_run_1 = 0;
            cout << "Cantor fail factor" << endl;
        }
        //cout << "OK1" << endl;
        auto res_2 = factor(true, poly, q, time_run_2);
        if (res_2.empty()) {
            // i--;
            // continue;
            time_run_2 = 0;
            cout << "Ben-Or fail factor" << endl;
        }
        //cout << "OK2" << endl;
        auto res_3 = factor_by_Shoup_Algorithm(poly, q, time_run_3);
        if (res_3.empty()) {
            // i--;
            // continue;
            time_run_3 = 0;
            cout << "Shoup fail factor" << endl;
        }
        //cout << "OK3" << endl;
        // write to csv
        file << degrees[i] << "," << time_run_1 << "," << time_run_2 << "," << time_run_3 << "\n";
    }
    cout << "Success" << endl;
}

void analyze_factorization_algorithms(const std::string &output_filename, const mpz_class &q) {
    std::ofstream output_file(output_filename);

    // Tiêu đề bảng: n, Số lần gọi GCD, Số lần chia, Thời gian chạy
    output_file << "n,Ben-Or_GCD,Shoup_GCD,Ben-Or_Div,Shoup_Div,Ben-Or_Time(ms),Shoup_Time(ms)\n";

    std::vector<int> degrees = {24, 56,  104, 252, 500};

    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_equal_degree(degrees[i], q);
        std::cout << "Testing polynomial with degree: " << degrees[i] << std::endl;
        int time_run_2 = 0, time_run_3 = 0;
        int benor_count[2] = {0}, shoup_count[2] = {0};

        // Chạy thuật toán Ben-Or
        std::fill(std::begin(cnt), std::end(cnt), 0);  // Reset bộ đếm
        auto res_2 = factor(true, poly, q, time_run_2);
        std::copy(std::begin(cnt), std::end(cnt), std::begin(benor_count));  // Lưu số lần gọi GCD & div
        if (res_2.empty() && degrees[i] < 500) {
            i--;
            continue;
            //time_run_2 = 0;
            std::cout << "Ben-Or factorization failed" << std::endl;
        }

        // Chạy thuật toán Shoup
        std::fill(std::begin(cnt), std::end(cnt), 0);  // Reset bộ đếm
        auto res_3 = factor_by_Shoup_Algorithm(poly, q, time_run_3);
        std::copy(std::begin(cnt), std::end(cnt), std::begin(shoup_count));  // Lưu số lần gọi GCD & div
        if (res_3.empty() && degrees[i] < 500) {
            i--;
            continue;
            //time_run_3 = 0;
            std::cout << "Shoup factorization failed" << std::endl;
        }

        // Ghi vào file
        output_file << degrees[i] << ","
                    << benor_count[0] << "," << shoup_count[0] << ","
                    << benor_count[1] << "," << shoup_count[1] << ","
                    << time_run_2 << "," << time_run_3 << "\n";
    }

    std::cout << "Analysis completed successfully!" << std::endl;
}

// viết hàm so sánh thời gian chạy hai thuật Ben-Or và Minhotte
void compare_Ben_Or_Mignotte_Hybird(const std::string &output_filename, const mpz_class &q) {
    std::ofstream output_file(output_filename);

    // Tiêu đề bảng: n, Số lần gọi GCD, Số lần chia, Thời gian chạy
//    output_file << "n,Ben-Or_Time(ms),Mignotte_Time(ms), Hybrid_algo_time(ms), Hybrid_1_algo_time(ms)\n";
    output_file << "n,Ben-Or_Time(ms),Mignotte_Time(ms), Hybrid_algo_time(ms)\n";

    std::vector<int> degrees = {10, 30, 50, 100, 150, 200, 250, 500};
    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_mono_degree(degrees[i], q);
        std::cout << "Testing polynomial with degree: " << degrees[i] << std::endl;
//        int time_run_1 = 0, time_run_2 = 0, time_run_3 = 0, time_run_4 = 0;
        int time_run_1 = 0, time_run_2 = 0, time_run_3 = 0;

        auto res_1 = find_root_by_Ben_Or(poly, q, time_run_1);
        if (res_1.empty()) {
            time_run_1 = 0;
            std::cout << "Ben-Or factorization failed" << std::endl;
        }

        auto res_2 = find_root_by_Mignotte(poly, q, time_run_2);

        auto res_3 = find_root_by_hybrid(poly, q, time_run_3);

//        auto res_4 = find_root_by_hybird_1(poly, q, time_run_4);

        output_file << degrees[i] << "," << time_run_1 << "," << time_run_2 << "," << time_run_3 << "\n";
    }

}
void decoding_time_comparison(const std::string &output_filename, const mpz_class &q) {
    std::ofstream output_file(output_filename);
    output_file << "n, Hybrid_algo_time(ms), Chien(ms)\n";

    std::vector<int> degrees = {10, 30, 50, 100, 150, 200, 250, 500};
    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_mono_degree(degrees[i], q);
        std::cout << "Testing polynomial with degree: " << degrees[i] << std::endl;
        int time_run_1 = 0, time_run_2 = 0;

        auto res_1 = find_root_by_hybrid(poly, q, time_run_1);
        auto res_2 = find_root_by_chien(poly, time_run_2);

        output_file << degrees[i] << "," << time_run_1 << "," << time_run_2  << "\n";
    }
}

void get_statistics_attemps_in_node(const mpz_class &q) {
    std::vector<int> degrees = {10, 30, 50, 100, 150, 200, 250, 500};
    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_mono_degree(degrees[i], q);
        int time = 0;
        auto res_1 = find_root_by_Ben_Or(poly, q, time);
        cout << "Ok" << endl;
    }
}

