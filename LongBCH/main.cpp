#include <gmpxx.h>
#include <map>
#include <iostream>
#include "GFq.h"
#include "testing.h"
#include "BCH.h"
#include "Polynomial.h"
#include "Utils.h"
#include <vector>

#include "Algorthms.h"

using namespace std;

// mpz_class m = 4; // GF(2^m)
// mpz_class primitive = 19; // Primitive polynomial in GF(2^m)
// mpz_class m = 8; // GF(2^m)
// mpz_class primitive = 285; // Primitive polynomial in GF(2^m)
mpz_class m = 18; // GF(2^m)
mpz_class primitive = 262183; // Primitive polynomial in GF(2^m)
// mpz_class m = 13; // GF(2^m)
// mpz_class primitive = 11135; // Primitive polynomial in GF(2^m)
// mpz_class m = 10; // GF(2^m)
// mpz_class primitive = 1135; // Primitive polynomial in GF(2^m)
// mpz_class m = 15; // GF(2^m)
// mpz_class primitive = 32771; // Primitive polynomial in GF(2^m)
// mpz_class m = 20; // GF(2^m)
// mpz_class primitive = 1050355; // Primitive polynomial in GF(2^m)
// mpz_class m = 25; // GF(2^m)
// mpz_class primitive = 33554441; // Primitive polynomial in GF(2^m)
//mpz_class m = 3; // GF(2^m)
//mpz_class primitive = 11; // Primitive polynomial in GF(2^m)
vector<mpz_class> deg_alpha;
vector<mpz_class> log_alpha;
mpz_class n = (mpz_class(1) << m.get_ui()) - 1;
vector<mpz_class> pi;
map<mpz_class, vector<mpz_class>> cyclotomic;
//mpz_class delta = 824;
mpz_class delta = 88;

int main() {
    mpz_class q = n + 1;
    initialize_tables_and_pi();
    Polynomial polys;
    int time;

    // Polynomial poly = random_polynomial_factor_mono_degree(100, q);
    // auto res = find_root_by_hybird_1(poly, q, time);
    // for (auto &r : res) {
    //     cout << r.first << endl;
    // }
    compare_Ben_Or_Mignotte_Hybird("output_2^18.csv", q);




    return 0;
}
