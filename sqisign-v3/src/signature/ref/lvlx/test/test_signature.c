/**
 * @file test_signature.c
 *
 * @note Modified by Sai-Yu Wang on 2026/09/18.
 * Add enum `signature_method_t` to represent the signature method to be used.
 * Add parameter `--method=<ori|new|both>` to specify the signature method to be used.
 * Add parameter `method` to the function `test_sqisign` to specify the signature method to be used.
 */
#include <stdio.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include <verification.h>
#include <signature.h>

#include <tools.h>
#include <rng.h>
#include <bench_test_arguments.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

typedef int (*sign_fn_t)(signature_t *, const public_key_t *, secret_key_t *, const unsigned char *, size_t);
typedef int (*verify_fn_t)(signature_t *, const public_key_t *, const unsigned char *, size_t);

int
test_sqisign(int repeat, const bool use_new_sign)
{
    int res = 1;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    const unsigned char msg[32] = { 0 };

    sign_fn_t sign = use_new_sign ? protocols_sign_new : protocols_sign;
    verify_fn_t verify = use_new_sign ? protocols_verify_new : protocols_verify;

    public_key_init(&pk);
    secret_key_init(&sk);

    printf("\n\nTesting signatures for " STRINGIFY(SQISIGN_VARIANT) " with %s implementation:\n\n",
           use_new_sign ? "new" : "original");

    for (int i = 0; i < repeat; ++i) {
        printf("Key #%d \n", i);

        if (!protocols_keygen(&pk, &sk)) {
            printf("keygen failed ! \n");
            res = 0;
            continue;
        }
        for (int j = 0; j < repeat; j++) {
            printf("  Sig #%d \n", j);
            int check;

            check = sign(&sig, &pk, &sk, msg, 32);
            if (!check) {
                printf("sign failed ! \n");
                res = 0;
            }

            check = verify(&sig, &pk, msg, 32);
            if (!check) {
                printf("verify failed ! \n");
                res = 0;
            }
        }
    }

    return res;
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
    int res = 0;
    int method = METHOD_BOTH;

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
        printf("Where <iterations> is the number of iterations used for testing; if not "
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

    if ((method == METHOD_ORI) || (method == METHOD_BOTH)) {
        res = test_sqisign(iterations, false);

        if (!res) {
            printf("\nSome tests with original implementation failed!\n");
        } else {
            printf("All tests with original implementation passed!\n");
        }
    }

    if ((method == METHOD_NEW) || (method == METHOD_BOTH)) {
        res = test_sqisign(iterations, true);

        if (!res) {
            printf("\nSome tests with new implementation failed!\n");
        } else {
            printf("All tests with new implementation passed!\n");
        }
    }
    return (!res);
}
