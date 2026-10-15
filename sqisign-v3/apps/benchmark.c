// SPDX-License-Identifier: Apache-2.0

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <api.h>
#include <sig.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>
#if defined(TARGET_BIG_ENDIAN)
#include <tutil.h>
#endif

typedef int (*sign_fn_t)(unsigned char *,
                         unsigned long long *,
                         const unsigned char *,
                         unsigned long long,
                         const unsigned char *);
typedef int (*open_fn_t)(unsigned char *,
                         unsigned long long *,
                         const unsigned char *,
                         unsigned long long,
                         const unsigned char *);

typedef struct
{
    uint64_t keypair;
    uint64_t sign;
    uint64_t verify;
} bench_result_t;

bench_result_t
bench(size_t runs, bool use_new)
{
    const size_t m_len = 32;
    const size_t sm_len = CRYPTO_BYTES + m_len;

    unsigned char *pkbuf = calloc(runs, CRYPTO_PUBLICKEYBYTES);
    unsigned char *skbuf = calloc(runs, CRYPTO_SECRETKEYBYTES);
    unsigned char *smbuf = calloc(runs, sm_len);
    unsigned char *mbuf = calloc(runs, m_len);

    unsigned char *pk[runs], *sk[runs], *sm[runs], *m[runs];
    bench_result_t result;
    for (size_t i = 0; i < runs; ++i) {
        pk[i] = pkbuf + i * CRYPTO_PUBLICKEYBYTES;
        sk[i] = skbuf + i * CRYPTO_SECRETKEYBYTES;
        sm[i] = smbuf + i * sm_len;
        m[i] = mbuf + i * m_len;
        if (randombytes(m[i], m_len))
            abort();
    }

    unsigned long long len;
    sign_fn_t sign = use_new ? sqisign_sign_new : crypto_sign;
    open_fn_t open = use_new ? sqisign_open_new : crypto_sign_open;

    printf("\n\n%s (%zu iterations, %s implementation)\n", CRYPTO_ALGNAME, runs, use_new ? "new" : "original");

    BENCH_CODE_1(runs);
    crypto_sign_keypair(pk[i], sk[i]);
    BENCH_CODE_2("keypair", result.keypair);

    BENCH_CODE_1(runs);
    len = sm_len;
    sign(sm[i], &len, m[i], m_len, sk[i]);
    if (len != sm_len)
        abort();
    BENCH_CODE_2("sign", result.sign);

    int ret;
    BENCH_CODE_1(runs);
    len = m_len;
    ret = open(m[i], &len, sm[i], sm_len, pk[i]);
    if (ret)
        abort();
    BENCH_CODE_2("verify", result.verify);

    free(pkbuf);
    free(skbuf);
    free(smbuf);
    free(mbuf);

    return result;
}

static void
print_comparison(const char *name, uint64_t original, uint64_t new)
{
    printf("  %-10s | %9" BENCH_FMT " | %9" BENCH_FMT " | %7.2f%% | %6.2fx\n",
           name,
           BENCH_TRUNC(original),
           BENCH_TRUNC(new),
           100. * ((double)new - (double)original) / (double)original,
           (double)original / (double)new);
}

int
main(int argc, char *argv[])
{
    enum
    {
        METHOD_ORI,
        METHOD_NEW,
        METHOD_BOTH
    };
    uint32_t seed[12] = { 0 };
    int iterations = SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int method = METHOD_BOTH;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }

        if (strncmp(argv[i], "--method=", 9) == 0) {
            const char *value = argv[i] + 9;
            if (strcmp(value, "ori") == 0) {
                method = METHOD_ORI;
            } else if (strcmp(value, "new") == 0) {
                method = METHOD_NEW;
            } else if (strcmp(value, "both") == 0) {
                method = METHOD_BOTH;
            } else {
                fprintf(stderr, "Invalid method: %s\n", value);
                return 1;
            }
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>] [--method=<ori|new|both>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        printf("Where <method> selects the protocol implementation; default is both\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    cpucycles_init();

    bench_result_t original = { 0 };
    bench_result_t new = { 0 };

    if (method == METHOD_ORI || method == METHOD_BOTH) {
        original = bench(iterations, false);
    }
    if (method == METHOD_NEW || method == METHOD_BOTH) {
        new = bench(iterations, true);
    }

    if (method == METHOD_BOTH) {
        printf("\n\nTime comparison (original / new, %s):\n", BENCH_UNIT);
        printf("  %-10s | %9s | %9s | %8s | %6s\n", "", "original", "new", "change", "speedup");
        print_comparison("keypair", original.keypair, new.keypair);
        print_comparison("sign", original.sign, new.sign);
        print_comparison("verify", original.verify, new.verify);
        // print_comparison(
        //     "total", original.keypair + original.sign + original.verify, new.keypair + new.sign + new.verify);
    }

    return 0;
}
