#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <gmp.h>

#define BITS (4096)
#define SPACING 1
#define MAXTERMS (BITS/2-1)

void pow2(mpz_t out, int n) {
    mpz_set_ui(out, 1);
    mpz_mul_2exp(out, out, n);
}

uint64_t count = 0;
int total = 1;
int id = 0;

void recurse(mpz_t num, int min, int max, int bits) {
    if (bits == 0) {
        ++count;
        if ((count % total) != id) return;
        if ((count % 65536) == 0) printf(".");
        if (mpz_probab_prime_p(num, 1)) {
            mpz_t order;
            mpz_init(order);
            mpz_sub_ui(order, num, 1);
            mpz_fdiv_q_2exp(order, order, 1);
            if (mpz_probab_prime_p(order, 15) && mpz_probab_prime_p(num, 15)) {
                gmp_printf("\nFOUND: %Zx\n", num);
            }
            mpz_clear(order);
        }
        return;
    }
    mpz_t pow;
    mpz_init(pow);
    for (int pos = min + bits - 1; pos < max; ++pos) {
        pow2(pow, pos);
        mpz_sub(num, num, pow);
        recurse(num, min, pos, bits - 1);
        mpz_add(num, num, pow);
    }
    mpz_clear(pow);
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    total = argc > 1 ? strtoul(argv[1], NULL, 0) : 1;
    id = argc > 2 ? strtoul(argv[2], NULL, 0) : 0;
    printf("CPU %i/%i\n", id, total);

    mpz_t n;
    mpz_init(n);
    pow2(n, BITS);
    mpz_sub_ui(n, n, 1);
    for (int bits = 0; bits < MAXTERMS; ++bits) {
        recurse(n, 1, MAXTERMS, bits);
    }

    return 0;
}
