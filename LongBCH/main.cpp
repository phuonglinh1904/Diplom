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
//mpz_class m = 16; // GF(2^m)
//mpz_class primitive = 65581; // Primitive polynomial in GF(2^m)
//mpz_class m = 20;
//mpz_class primitive = 1050355; // Primitive polynomial in GF(2^m)
mpz_class m = 25; // GF(2^m)
mpz_class primitive = 33554441; // Primitive polynomial in GF(2^m)
vector<mpz_class> deg_alpha;
vector<mpz_class> log_alpha;
mpz_class n = (mpz_class(1) << m.get_ui()) - 1;
vector<mpz_class> pi;
map<mpz_class, vector<mpz_class>> cyclotomic;
mpz_class delta = 1000;

int main() {
    mpz_class q = n + 1;
    initialize_tables_and_pi();
    std::string output_filename = "memory_comparison_results.csv";
    compare_memory_usage_Ben_Or_Mignotte_Hybrid(output_filename, q);
    return 0;
}
