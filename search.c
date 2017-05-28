#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <gmp.h>

#define BITS 3072
#define SPACING 64
#define MAXTERMS 23

void pow2(mpz_t out, int n) {
    mpz_set_ui(out, 1);
    mpz_mul_2exp(out, out, n);
}

int main(int argc, char** argv) {

    int total = argc > 1 ? strtoul(argv[1], NULL, 0) : 1;
    int num = argc > 2 ? strtoul(argv[2], NULL, 0) : 0;
    printf("CPU %i/%i\n", num, total);

    mpz_t n, m, p;
    mpz_init(n);
    mpz_init(m);
    mpz_init(p);
    for (int bits = 0; bits <= MAXTERMS; ++bits) {
        int cnt = 0;
        for (int val = 0; val < (1 << MAXTERMS); ++val) {
            if (__builtin_popcount(val) != bits) continue;
            uint64_t hash = ((uint64_t)val) * 0x123456789ABCDEFULL;
            hash = (hash >> 32) | (hash << 32);
            hash *= 0x123456789ABCDEFULL;
            hash = (hash >> 32) | (hash << 32);
            hash *= 0x123456789ABCDEFULL;
            if ((hash % total) != num) continue;
            ++cnt;
            pow2(n, BITS);
            pow2(p, BITS - 1);
            mpz_sub_ui(n, n, 1);
            mpz_sub_ui(p, p, 1);
            for (int bit = 0; bit < MAXTERMS; ++bit) {
                if ((val >> bit) & 1) {
                    pow2(m, (bit + 1) * SPACING - 1);
                    mpz_sub(p, p, m);
                    mpz_mul_2exp(m, m, 1);
                    mpz_sub(n, n, m);
                }
            }
            if (mpz_probab_prime_p(p, 1) && mpz_probab_prime_p(n, 15) && mpz_probab_prime_p(p, 15)) {
                printf("2^%i", BITS);
                for (int bit = MAXTERMS - 1; bit >= 0; --bit) {
                    if ((val >> bit) & 1) {
                        printf(" - 2^%i", (bit + 1) * SPACING);
                    }
                }
                printf(" - 1\n");
            }
        }
    }

    return 0;
}
