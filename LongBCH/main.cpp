#include <gmpxx.h>
#include <map>
#include "GFq.h"
#include "testing.h"
#include "Polynomial.h"
#include "Algorthms.h"
#include "Utils.h"
#include <vector>
#include <typeinfo>

using namespace std;
// mpz_class m = 10; // GF(2^m)
// mpz_class primitive = 1135; // Primitive polynomial in GF(2^m)
// mpz_class m = 15; // GF(2^m)
// mpz_class primitive = 32771; // Primitive polynomial in GF(2^m)
// mpz_class m = 20; // GF(2^m)
// mpz_class primitive = 1050355; // Primitive polynomial in GF(2^m)
mpz_class m = 25; // GF(2^m)
mpz_class primitive = 33554441; // Primitive polynomial in GF(2^m)
mpz_class epsilon = 2;
vector<mpz_class> deg_alpha;
vector<mpz_class> log_alpha;
mpz_class n = (mpz_class(1) << m.get_ui()) - 1;
map<mpz_class, pair<vector<mpz_class>, mpz_class>> cyclotomic;

int main() {
    mpz_class q = n + 1;
    initialize_tables();
    compare_Ben_Or_Mignotte("compare_2^25_algos.csv", q);


    return 0;
}
