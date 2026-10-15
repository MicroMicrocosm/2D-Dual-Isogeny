/**
 * @file bench_signature.c
 *
 * @note Modified by Sai-Yu Wang on 2026/09/18.
 * Add enum `signature_method_t` to represent the signature method to be used.
 * Add parameter `--method=<ori|new|both>` to specify the signature method to be used.
 * Add parameter `method` to the function `bench_sqisign` to specify the signature method to be used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include <verification.h>
#include <signature.h>

#include <tools.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

typedef int (*sign_fn_t)(signature_t *, const public_key_t *, secret_key_t *, const unsigned char *, size_t);
typedef int (*verify_fn_t)(signature_t *, const public_key_t *, const unsigned char *, size_t);

typedef struct
{
    uint64_t keygen_cycles;
    uint64_t sign_cycles;
    uint64_t verify_cycles;
    double keygen_ms;
    double sign_ms;
    double verify_ms;
} bench_result_t;

static __inline__ uint64_t
rdtsc(void)
{
    return (uint64_t)cpucycles();
}

bench_result_t
bench_sqisign(uint64_t bench, const bool use_new_sign)
{
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1;
    clock_t t;
    float ms;
    bench_result_t result;

    public_key_t pks[bench];
    secret_key_t sks[bench];
    signature_t sigs[bench];
    unsigned char msg[32] = { 0 };

    sign_fn_t sign = use_new_sign ? protocols_sign_new : protocols_sign;
    verify_fn_t verify = use_new_sign ? protocols_verify_new : protocols_verify;

    for (uint64_t i = 0; i < bench; i++) {
        public_key_init(&pks[i]);
        secret_key_init(&sks[i]);
    }

    printf("\n\nBenchmarking signatures for " STRINGIFY(SQISIGN_VARIANT) " (%llu iterations, %s implementation):\n",
           (unsigned long long)bench,
           use_new_sign ? "new" : "original");

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        if (!protocols_keygen(&pks[i], &sks[i])) {
            printf("keygen failed\n");
            abort();
        }
    }
    t1 = rdtsc();
    result.keygen_cycles = (t1 - t0) / bench;
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    result.keygen_ms = (double)ms / bench;
    printf("  Average keygen time [%.2f ms]\n", (float)(ms / bench));
    printf("  \x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", result.keygen_cycles);

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        sign(&(sigs[i]), &(pks[i]), &(sks[i]), msg, sizeof(msg) / sizeof(*msg));
    }
    t1 = rdtsc();
    result.sign_cycles = (t1 - t0) / bench;
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    result.sign_ms = (double)ms / bench;
    printf("  Average signature with %s implementation time [%.2f ms]\n",
           use_new_sign ? "new" : "original",
           (float)(ms / bench));
    printf("  \x1b[34mAvg signature with %s implementation: %'" PRIu64 " cycles\x1b[0m\n",
           use_new_sign ? "new" : "original",
           result.sign_cycles);

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        int check = verify(&(sigs[i]), &(pks[i]), msg, sizeof(msg) / sizeof(*msg));
        if (!check) {
            printf("verify with %s implementation failed ! \n", use_new_sign ? "new" : "original");
        }
    }
    t1 = rdtsc();
    result.verify_cycles = (t1 - t0) / bench;
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    result.verify_ms = (double)ms / bench;
    printf("  Average verification with %s implementation time [%.2f ms]\n",
           use_new_sign ? "new" : "original",
           (float)(ms / bench));
    printf("  \x1b[34mAvg verification with %s implementation: %'" PRIu64 " cycles\x1b[0m\n",
           use_new_sign ? "new" : "original",
           result.verify_cycles);
    return result;
}

static void
print_comparison_cycles(const char *name, uint64_t original, uint64_t new)
{
    printf("  %-10s | %9" BENCH_FMT " | %9" BENCH_FMT " | %7.2f%% | %6.2fx\n",
           name,
           BENCH_TRUNC(original),
           BENCH_TRUNC(new),
           100. * ((double)original - new) / original,
           (double)original / new);
}

static void
print_comparison_ms(const char *name, double original, double new)
{
    printf("  %-10s | %9.2f | %9.2f | %7.2f%% | %6.2fx\n",
           name,
           original,
           new,
           100. * (new - original) / original,
           original / new);
}

// run all tests in module
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

    if ((method == METHOD_ORI) || (method == METHOD_BOTH)) {
        original = bench_sqisign(iterations, false);
    }
    if ((method == METHOD_NEW) || (method == METHOD_BOTH)) {
        new = bench_sqisign(iterations, true);
    }

    if (method == METHOD_BOTH) {
        printf("\n\nTime comparison (original / new, %s):\n", BENCH_UNIT);
        printf("  %-10s | %9s | %9s | %8s | %6s\n", "", "original", "new", "change", "speedup");
        print_comparison_cycles("keygen", original.keygen_cycles, new.keygen_cycles);
        print_comparison_cycles("sign", original.sign_cycles, new.sign_cycles);
        print_comparison_cycles("verify", original.verify_cycles, new.verify_cycles);
        // print_comparison_cycles("total",
        //                         original.keygen_cycles + original.sign_cycles + original.verify_cycles,
        //                         new.keygen_cycles + new.sign_cycles + new.verify_cycles);

        printf("\nTime comparison (original / new, ms):\n");
        printf("  %-10s | %9s | %9s | %8s | %6s\n", "", "original", "new", "change", "speedup");
        print_comparison_ms("keygen", original.keygen_ms, new.keygen_ms);
        print_comparison_ms("sign", original.sign_ms, new.sign_ms);
        print_comparison_ms("verify", original.verify_ms, new.verify_ms);
        // print_comparison_ms("total",
        //                     original.keygen_ms + original.sign_ms + original.verify_ms,
        //                     new.keygen_ms + new.sign_ms + new.verify_ms);
    }

    return 0;
}
