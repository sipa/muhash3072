#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <gmp.h>

typedef unsigned __int128 uint128_t;

/* [c0,c1,c2] += n * [d0,d1,d2]. c2 is 0 initially */
#define muladdn(c0,c1,c2,d0,d1,d2,n) { \
    uint128_t t = (uint128_t)d0 * n + c0; \
    c0 = t; \
    t >>= 64; \
    t += (uint128_t)d1 * n + c1; \
    c1 = t; \
    t >>= 64; \
    c2 = t + d2 * n; \
}

/* [c0,c1] *= n */
#define muln_fast(c0,c1,n) { \
    uint128_t t = (uint128_t)c0 * n; \
    c0 = t; \
    t >>= 64; \
    t += (uint128_t)c1 * n; \
    c1 = t; \
    t >>= 64; \
}

/** [c0,c1] = a * b */
#define mulset(c0,c1,a,b) { \
    uint128_t t = (uint128_t)a * b; \
    c2 = 0; \
    c1 = t >> 64; \
    c0 = t; \
}

/** Extract the lowest 64 bits of [c0,c1,c2] into n, and left shift the number 64 bits. */
#define extract(c0,c1,c2,n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = c2; \
    c2 = 0; \
}

/** Extract the lowest 64 bits of [c0,c1] into n, and left shift the number 64 bits. */
#define extract_fast(c0,c1,n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = 0; \
}

#ifdef X86_64_ASM

/** [c0,c1,c2] += a * b */
#define muladd(c0,c1,c2,a,b) { \
    uint64_t tl, th; \
    __asm__ ("mulq %3" : "=a"(tl), "=d"(th) : "a"(a), "g"(b) : "cc"); \
    __asm__ ("addq %3,%0; adcq %4,%1; adcq $0,%2" : "+r"(c0), "+r"(c1), "+r"(c2) : "a"(tl), "d"(th) : "cc"); \
}
/** [c0,c1,c2] += 2 * a * b */
#define muladd2(c0,c1,c2,a,b) { \
    uint64_t tl, th; \
    __asm__ ("mulq %3" : "=a"(tl), "=d"(th) : "a"(a), "g"(b) : "cc"); \
    __asm__ ("addq %3,%0; adcq %4,%1; adcq $0,%2" : "+r"(c0), "+r"(c1), "+r"(c2) : "a"(tl), "d"(th) : "cc"); \
    __asm__ ("addq %3,%0; adcq %4,%1; adcq $0,%2" : "+r"(c0), "+r"(c1), "+r"(c2) : "a"(tl), "d"(th) : "cc"); \
}
/** [c0,c1] += a */
#define sumadd_fast(c0,c1,a) { \
    __asm__ ("add %2,%0; adc $0,%1" : "+r"(c0), "+r"(c1) : "r"(a) : "cc"); \
}

#else

/** [c0,c1,c2] += a * b */
#define muladd(c0,c1,c2,a,b) { \
    uint64_t tl, th; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = t >> 64;         /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
}

/** [c0,c1,c2] += 2 * a * b */
#define muladd2(c0,c1,c2,a,b) { \
    uint64_t tl, th; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = t >> 64;         /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
}

/** [c0,c1] += a */
#define sumadd_fast(c0,c1,a) { \
    c0 += (a);                 /* overflow is handled on the next line */ \
    c1 += (c0 < (a)) ? 1 : 0;  /* never overflows by contract (verified the next line) */ \
}

#endif


struct num3072 {
    uint64_t d[48];
};

static void mul3072(struct num3072* out, const struct num3072* a, const struct num3072* b) {
    uint64_t c0 = 0, c1 = 0;
    struct num3072 tmp;

    /* Compute limbs 0..46 of a*b into tmp, including one reduction. */
    for (int j = 0; j < 47; ++j) {
        uint64_t d0 = 0, d1 = 0, d2 = 0, c2 = 0;
        mulset(d0, d1, a->d[1 + j], b->d[48 + j - (1 + j)]);
        for (int i = 2 + j; i < 48; ++i) muladd(d0, d1, d2, a->d[i], b->d[48 + j - i]);
        muladdn(c0, c1, c2, d0, d1, d2, 1103717);
        for (int i = 0; i <= j; ++i) muladd(c0, c1, c2, a->d[i], b->d[j - i]);
        extract(c0, c1, c2, tmp.d[j]);
    }
    /* Compute limb 47 of a*b into tmp (unaffected by the first reduction). */
    {
        uint64_t c2 = 0;
        for (int i = 0; i <= 47; ++i) muladd(c0, c1, c2, a->d[i], b->d[47 - i]);
        extract(c0, c1, c2, tmp.d[47]);
    }
    /* Perform a second reduction.*/
    muln_fast(c0, c1, 1103717);
    for (int j = 0; j < 48; ++j) {
        sumadd_fast(c0, c1, tmp.d[j]);
        extract_fast(c0, c1, out->d[j]);
    }
}

static void sqr3072(struct num3072* out, const struct num3072* a) {
    uint64_t c0 = 0, c1 = 0;
    struct num3072 tmp;

    /* Compute limbs 0..47 of a*a into tmp, including one reduction. */
    for (int j = 0; j < 47; ++j) {
        uint64_t d0 = 0, d1 = 0, d2 = 0, c2 = 0;
        for (int i = 0; i < (47 - j) / 2; ++i) muladd2(d0, d1, d2, a->d[i + j + 1], a->d[47 - i]);
        if ((47 - j) & 1) muladd(d0, d1, d2, a->d[(47 - j) / 2 + j + 1], a->d[47 - (47 - j)/2]);
        muladdn(c0, c1, c2, d0, d1, d2, 1103717);
        for (int i = 0; i < (j + 1) / 2; ++i) muladd2(c0, c1, c2, a->d[i], a->d[j - i]);
        if ((j + 1) & 1) muladd(c0, c1, c2, a->d[(j + 1) / 2], a->d[j - (j + 1) / 2]);
        extract(c0, c1, c2, tmp.d[j]);
    }
    /* Compute limb 47 of a*b into tmp (unaffected by the first reduction). */
    {
        uint64_t c2 = 0;
        for (int i = 0; i < 24; ++i) muladd2(c0, c1, c2, a->d[i], a->d[47 - i]);
        extract(c0, c1, c2, tmp.d[47]);
    }
    /* Perform a second reduction. */
    muln_fast(c0, c1, 1103717);
    for (int j = 0; j < 48; ++j) {
        sumadd_fast(c0, c1, tmp.d[j]);
        extract_fast(c0, c1, out->d[j]);
    }
}


void inv3072(struct num3072* out, const struct num3072* a) {
    struct num3072 p[12]; // p[i] = a^(2^(2^i)-1)
    struct num3072 x;

    p[0] = *a;

    for (int i = 0; i < 11; ++i) {
        p[i + 1] = p[i];
        for (int j = 0; j < (1 << i); ++j) sqr3072(&p[i + 1], &p[i + 1]);
        mul3072(&p[i + 1], &p[i + 1], &p[i]);
    }

    x = p[11];

    for (int j = 0; j < 512; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[9]);
    for (int j = 0; j < 256; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[8]);
    for (int j = 0; j < 128; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[7]);
    for (int j = 0; j < 64; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[6]);
    for (int j = 0; j < 32; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[5]);
    for (int j = 0; j < 8; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[3]);
    for (int j = 0; j < 2; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[1]);
    for (int j = 0; j < 1; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[0]);
    for (int j = 0; j < 5; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[2]);
    for (int j = 0; j < 3; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[0]);
    for (int j = 0; j < 2; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[0]);
    for (int j = 0; j < 4; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[0]);
    for (int j = 0; j < 4; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[1]);
    for (int j = 0; j < 3; ++j) sqr3072(&x, &x);
    mul3072(&x, &x, &p[0]);

    *out = x;
}

int main(void) {
    unsigned char ab[384], bb[384];
    FILE* f = fopen("/dev/urandom", "r");
    fread(ab, 384, 1, f);
    fread(bb, 384, 1, f);
    fclose(f);

/*
    struct num3072 an, bn;
    memcpy(&an, ab, 384);
    memcpy(&bn, bb, 384);
    printf("begin: ");
    for (int i = 0; i < 48; ++i) {
        printf("%016lx", (unsigned long)an.d[47 - i]);
    }
    printf("\n");
    for (int i = 0; i < 1000; ++i) {
        mul3072(&an, &an, &bn);
        inv3072(&an, &an);
    }

    printf("end: ");
    for (int i = 0; i < 48; ++i) {
        printf("%016lx", (unsigned long)an.d[47 - i]);
    }
    printf("\n");
*/

    mpz_t ag, bg, t, m;
    mpz_init(ag);
    mpz_init(bg);
    mpz_init(t);
    mpz_init(m);
    mpz_set_ui(m, 2);
    mpz_pow_ui(m, m, 3072);
    mpz_sub_ui(m, m, 1103717);
    mpz_import(ag, 384, -1, 1, -1, 0, ab);
    mpz_import(bg, 384, -1, 1, -1, 0, bb);
    gmp_printf("begin: %Zx\n", ag);
    for (int i = 0; i < 1000000; ++i) {
/*        mpz_tdiv_q_2exp(ag, t, 3072);
        mpz_tdiv_r_2exp(t, t, 3072);
        mpz_addmul_ui(t, ag, 1103717);
        mpz_tdiv_r_2exp(ag, t, 3072);
        mpz_tdiv_q_2exp(t, t, 3072);
        mpz_addmul_ui(ag, t, 1103717);*/
        mpz_invert(ag, ag, m);
    }
    gmp_printf("end: %Zx\n", ag);

    return 0;
}
