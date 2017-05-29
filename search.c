#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <gmp.h>

#define BITS 256
#define SPACING 16
#define MAXTERMS 30

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
        for (int val = 1; val < (1ULL << MAXTERMS); val += 2) {
            if (__builtin_popcount(val) != bits) continue;
            for (int sign = 0; sign <= val ; ++sign) {
                if (sign & ~val) continue;
                if (__builtin_clz(sign) != __builtin_clz(val)) continue;
                uint64_t hash = ((uint64_t)val) * 0x123456789ABCDEFULL ^ sign;
                hash = (hash >> 32) | (hash << 32);
                hash *= 0x123456789ABCDEFULL;
                hash = (hash >> 32) | (hash << 32);
                hash *= 0x123456789ABCDEFULL;
                if ((hash % total) != num) continue;
                ++cnt;
                pow2(n, BITS);
                for (int bit = 0; bit < MAXTERMS; ++bit) {
                    if ((val >> bit) & 1) {
                        pow2(m, bit * SPACING);
                        if ((sign >> bit) & 1) {
                            mpz_sub(n, n, m);
                        } else {
                            mpz_add(n, n, m);
                        }
                    }
                }
                mpz_sub_ui(p, n, 1);
                mpz_fdiv_q_2exp(p, p, 1);
                if (mpz_probab_prime_p(p, 1) && mpz_probab_prime_p(n, 15) && mpz_probab_prime_p(p, 15)) {
                    printf("2^%i", BITS);
                    for (int bit = MAXTERMS; bit >= 0; --bit) {
                        if ((val >> bit) & 1) {
                            if (bit) {
                                printf(" %c 2^%i", (sign >> bit) & 1 ? '-' : '+', bit * SPACING);
                            } else {
                                printf(" %c 1", (sign >> bit) & 1 ? '-' : '+');
                            }
                        }
                    }
                    printf("\n");
                }
            }
        }
        printf("* %i bits: %i attempts\n", bits, cnt);
    }

    return 0;
}
