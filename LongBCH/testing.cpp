//
// Created by ADMIN on 27-Dec-24.
//
#include <fstream>
#include "testing.h"
#include "Algorthms.h"
#include "GFq.h"
#include "Utils.h"


// testing with n random polynomials with degree d
pair<vector<factor_result>, vector<factor_result>> testing_random(int degree, int num) {
    vector<factor_result> factor_results_algo_optimal;
    vector<factor_result> factor_results_algo_no_optimal;
    int recursion_count, time_run;
    Polynomial poly;
    mpz_class q = n + 1;
    for (int i = 0; i < num; i++) {
        poly = Polynomial::get_random_polynomial_with_degree(degree, q);
        recursion_count = 0; // Reset
        time_run = 0;
        auto result = factor(true, poly, q, recursion_count, time_run);
        if (result.empty() || (time_run == 0 && recursion_count == 1)) {
            i--;
            //cout << "Failed to factor polynomial: " << poly << endl;
            continue;
        }
        factor_results_algo_optimal.push_back(factor_result(poly, recursion_count - 1, time_run));
        recursion_count = 0; // Reset
        time_run = 0;
        auto result2 = factor(false, poly, q, recursion_count, time_run);
        factor_results_algo_no_optimal.push_back(factor_result(poly, recursion_count - 1, time_run));
        cout << i << endl;
    }
    return {factor_results_algo_optimal, factor_results_algo_no_optimal};
}

void write_to_csv(const std::vector<factor_result> &results, const std::string &filename) {
    std::ofstream file(filename);
    file << "Numerical,Degree,Recursion_count,Time_run\n";
    for (int i = 0; i < results.size(); i++) {
        file << i + 1 << "," << results[i].poly.get_degree() << "," << results[i].recursion_count << "," << results[i].time_run << "\n";
    }
}


bool check_exist_poly(vector<factor_result> factor_results, const Polynomial  &poly) {
    for (int j = 0; j < factor_results.size(); j++) {
        if (poly == factor_results[j].poly) {
            return true;
        }
    }
    return false;
}


void write_result_2_to_csv(const std::vector<factor_result> &results, const std::string &filename) {
    std::ofstream file(filename);
    file << "Numerical,Polynomial,Recursion_count\n";
    for (int i = 0; i < results.size(); i++) {
        file << i + 1 << "," << results[i].poly << "," << results[i].recursion_count << "\n";
    }
}

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
        int time_run_1 = 0, time_run_2 = 0, time_run_3 = 0, rec_count;
        auto res_1 = factor(false, poly, q, rec_count, time_run_1);
        if (res_1.empty()) {
            // i--;
            // continue;
            time_run_1 = 0;
            cout << "Cantor fail factor" << endl;
        }
        //cout << "OK1" << endl;
        auto res_2 = factor(true, poly, q, rec_count, time_run_2);
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
        int time_run_2 = 0, time_run_3 = 0, rec_count;
        int benor_count[2] = {0}, shoup_count[2] = {0};

        // Chạy thuật toán Ben-Or
        std::fill(std::begin(cnt), std::end(cnt), 0);  // Reset bộ đếm
        auto res_2 = factor(true, poly, q, rec_count, time_run_2);
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
void compare_Ben_Or_Mignotte(const std::string &output_filename, const mpz_class &q) {
    std::ofstream output_file(output_filename);

    // Tiêu đề bảng: n, Số lần gọi GCD, Số lần chia, Thời gian chạy
    //output_file << "n,Ben-Or_Time(ms),Mignotte_Time(ms), Hybrid_algo_time(ms)\n";
    output_file << "n,Cantor_Time(ms),Ben-Or_Time(ms), Mignotte_Time(ms)\n";

    std::vector<int> degrees = {10, 30, 50, 100, 150, 200, 250};//, 500};
    //std::vector<int> degrees = {250};
    int max_ran = 0;
    for (int i = 0; i < degrees.size(); i++) {
        Polynomial poly = random_polynomial_factor_mono_degree(degrees[i], q);
        std::cout << "Testing polynomial with degree: " << degrees[i] << std::endl;
        int time_run_1 = 0, time_run_2 = 0, time_run_3 = 0;

        auto res_1 = find_root_by_Ben_Or(poly, q, time_run_1);
        if (res_1.empty()) {
            if (max_ran++ < 5) {
                i--;
                cout << "Error num: " << max_ran << endl;
                continue;
            }
            else {
                time_run_1 = 0;
                max_ran = 0;
                std::cout << "Ben-Or factorization failed" << std::endl;
            }
        }
        max_ran = 0;
        auto res_2 = find_root_by_Mignotte(poly, q, time_run_2);
        if (res_2.empty()) {
            if (max_ran++ < 5) {
                i--;
                cout << "Error num: " << max_ran << endl;
                continue;
            }
            else {
                time_run_2 = 0;
                max_ran = 0;
                std::cout << "Mignotte factorization failed" << std::endl;
            }
        }

        auto res_3 = find_root_by_hybrid(poly, q, time_run_3);

        output_file << degrees[i] << "," << time_run_1 << "," << time_run_2 << "," << time_run_3 << "\n";
    }

}

