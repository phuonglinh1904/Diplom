#include <gmpxx.h>
#include <map>
#include "GFq.h"
#include "testing.h"
#include "BCH.h"
#include <vector>
#include <random>

using namespace std;
//mpz_class m = 10; // GF(2^m)
//mpz_class primitive = 1135; // Primitive polynomial in GF(2^m)
//mpz_class m = 15; // GF(2^m)
//mpz_class primitive = 32771; // Primitive polynomial in GF(2^m)
mpz_class m = 16; // GF(2^m)
mpz_class primitive = 65581; // Primitive polynomial in GF(2^m)
//mpz_class m = 20;
//mpz_class primitive = 1050355; // Primitive polynomial in GF(2^m)
//mpz_class m = 25; // GF(2^m)
//mpz_class primitive = 33554441; // Primitive polynomial in GF(2^m)
//mpz_class m = 3; // GF(2^m)
//mpz_class primitive = 11; // Primitive polynomial in GF(2^m)
vector<mpz_class> deg_alpha;
vector<mpz_class> log_alpha;
mpz_class n = (mpz_class(1) << m.get_ui()) - 1;
vector<mpz_class> pi;
map<mpz_class, vector<mpz_class>> cyclotomic;
mpz_class delta = 1000;

int main() {
    mpz_class q = n + 1;
    initialize_tables_and_pi();
    std::mt19937 Generator;
    cout<<"Starting cyclotomic classes..."<<endl;
    auto cyclotomic = find_cyclotomic_classes();
    cout<<"Starting to calculate g..."<<endl;
    auto g = compute_g();
    cout << "Generating polynomial g: "<<g<<endl;
    cout<<"Starting BCH decoding simulation..."<<endl;
    std::vector<double> noises = {0.001, 0.002, 0.0025, 0.003, 0.0035};
    int trials_per_noise = 30;

    for (double noise : noises) {
        int success_count_hybrid = 0, success_count_chien = 0;
        int total_time_hybrid = 0, total_time_chien = 0;

        std::cout << "\n--- Noise level: " << noise << " ---\n";

        for (int i = 0; i < trials_per_noise; ++i) {
            std::uniform_real_distribution<> distribution(0, 1);
            std::uniform_int_distribution<> bits(0, 1);

            auto k = n.get_ui() - g.get_degree();
            std::vector<mpz_class> seq(k);
            for (auto &j: seq) j = bits(Generator);

            auto enc = encode(seq, g, q);
            auto noised = enc;

            for (auto &j: noised) {
                j ^= (distribution(Generator) < noise);
            }

            int time_hybrid = 0, time_chien = 0;
            bool ok_hybrid = true, ok_chien = true;
            auto dec_hybrid = decode(noised, q, true, time_hybrid);
            if (dec_hybrid != enc) {
                std::cout << "Decode by hybrid failed | Time: " << time_hybrid << " ms" << std::endl;
                ok_hybrid = false;
            }
            auto dec_chien = decode(noised, q, false, time_chien);
            if (dec_chien != enc) {
                std::cout << "Decode by chien failed | Time: " << time_chien << " ms" << std::endl;
                ok_chien = false;
            }
            if (ok_chien) {
                success_count_chien++;
                total_time_chien += time_chien;
            }
            if (ok_hybrid) {
                success_count_hybrid++;
                total_time_hybrid += time_hybrid;
            }
        }
        std::cout << "Success by chien: " << success_count_chien << "/" << trials_per_noise << std::endl;
        if (success_count_chien > 0) {
            std::cout << "Average decode time for chien: "
                      << (total_time_chien / success_count_chien) << " ms\n";
        } else {
            std::cout << "All decodings failed for chien at this noise level.\n";
        }
        std::cout << "Success by hybrid: " << success_count_hybrid << "/" << trials_per_noise << std::endl;
        if (success_count_hybrid > 0) {
            std::cout << "Average decode time for hybrid: "
                      << (total_time_hybrid / success_count_hybrid) << " ms\n";
        } else {
            std::cout << "All decodings failed for hybrid at this noise level.\n";
        }
    }
//    decoding_time_comparison("decoding_time_comparison_2^16.csv", q);
    return 0;
}
