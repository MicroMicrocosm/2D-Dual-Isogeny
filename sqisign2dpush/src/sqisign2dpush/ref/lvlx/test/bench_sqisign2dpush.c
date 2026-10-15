#include <sqisign2dpush.h>
#include <rng.h>

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(__i386__) || defined(__x86_64__)
#include <x86intrin.h>
#define HAVE_CPU_CYCLES 1
#endif

typedef struct
{
    double *values;
    int count;
} samples_t;

typedef struct
{
    double isogeny_ms;
    double elliptic_curve_ms;
    double ideal_ms;
    double other_ms;
} breakdown_ms_t;

typedef struct
{
    double keygen_cycles;
    double sign_cycles;
    double verify_cycles;
    double keygen_ms;
    double sign_ms;
    double verify_ms;
    breakdown_ms_t keygen_breakdown;
    breakdown_ms_t sign_breakdown;
    breakdown_ms_t verify_breakdown;
} benchmark_result_t;

static unsigned long long
read_cycles(void)
{
#if HAVE_CPU_CYCLES
    _mm_lfence();
    unsigned long long cycles = __rdtsc();
    _mm_lfence();
    return cycles;
#else
    return (unsigned long long)clock();
#endif
}

static double
elapsed_ms(clock_t start, clock_t end)
{
    return 1000.0 * (double)(end - start) / (double)CLOCKS_PER_SEC;
}

static double
monotonic_ms(void)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return 1000.0 * (double)now.tv_sec + (double)now.tv_nsec / 1.0e6;
}

static breakdown_ms_t
keygen_breakdown(const keygen_timings_t *timings)
{
    breakdown_ms_t result = { 0 };
    result.isogeny_ms = timings->ms_keygen_ec_isog;
    result.elliptic_curve_ms = timings->ms_keygen_ec_non_isog;
    result.ideal_ms = timings->ms_keygen_quat;
    return result;
}

static breakdown_ms_t
sign_breakdown(const sign_timings_t *timings)
{
    breakdown_ms_t result = { 0 };
    result.isogeny_ms = timings->ms_sign_isog;
    result.elliptic_curve_ms = timings->ms_sign_ec_non_isog;
    result.ideal_ms = timings->ms_sign_quat;
    return result;
}

static breakdown_ms_t
verify_breakdown(const verify_timings_t *timings)
{
    breakdown_ms_t result = { 0 };

    result.isogeny_ms = timings->ms_ver_ec_isog;
    result.elliptic_curve_ms = timings->ms_ver_ec_non_isog;
    result.ideal_ms = timings->ms_ver_quat;
    return result;
}

/* Categories are disjoint. Keep the signed residual visible and reject
 * overlap beyond clock granularity and float accumulator rounding. */
static int
finish_breakdown(breakdown_ms_t *value, double total_ms, const char *phase, int iteration)
{
    value->other_ms = total_ms - value->isogeny_ms - value->elliptic_curve_ms - value->ideal_ms;
    double tolerance_ms = 4.0 * 1000.0 / CLOCKS_PER_SEC + total_ms * 1e-5;
    if (!isfinite(value->other_ms) || value->other_ms < -tolerance_ms) {
        fprintf(stderr, "%s timing overlap on iteration %d: residual %.6f ms\n", phase, iteration, value->other_ms);
        return 0;
    }
    return 1;
}

static double
average(const samples_t *samples)
{
    double total = 0.0;
    for (int i = 0; i < samples->count; ++i)
        total += samples->values[i];
    return total / samples->count;
}

static int
compare_double(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return (a > b) - (a < b);
}

static double
median(const samples_t *samples)
{
    double *sorted = malloc((size_t)samples->count * sizeof(*sorted));
    if (!sorted)
        return 0.0;
    memcpy(sorted, samples->values, (size_t)samples->count * sizeof(*sorted));
    qsort(sorted, (size_t)samples->count, sizeof(*sorted), compare_double);
    double result = sorted[samples->count / 2];
    if (samples->count % 2 == 0)
        result = (result + sorted[samples->count / 2 - 1]) / 2.0;
    free(sorted);
    return result;
}

static double
stddev(const samples_t *samples, double mean)
{
    double variance = 0.0;
    for (int i = 0; i < samples->count; ++i) {
        double delta = samples->values[i] - mean;
        variance += delta * delta;
    }
    return sqrt(variance / samples->count);
}

static void
accumulate_breakdown(breakdown_ms_t *total, const breakdown_ms_t *current)
{
    total->isogeny_ms += current->isogeny_ms;
    total->elliptic_curve_ms += current->elliptic_curve_ms;
    total->ideal_ms += current->ideal_ms;
    total->other_ms += current->other_ms;
}

static void
print_stats(const char *name, const samples_t *samples)
{
    double mean = average(samples);
    double min = samples->values[0];
    double max = samples->values[0];
    for (int i = 1; i < samples->count; ++i) {
        if (samples->values[i] < min)
            min = samples->values[i];
        if (samples->values[i] > max)
            max = samples->values[i];
    }
    printf("  %-7s average %.3f, stddev %.3f, median %.3f, min %.3f, max %.3f MCycles\n",
           name,
           mean,
           stddev(samples, mean),
           median(samples),
           min,
           max);
}

static void
print_time_stats(const char *name, const samples_t *samples)
{
    double min = samples->values[0];
    double max = samples->values[0];
    for (int i = 1; i < samples->count; ++i) {
        if (samples->values[i] < min)
            min = samples->values[i];
        if (samples->values[i] > max)
            max = samples->values[i];
    }
    double mean = average(samples);
    printf("  %-7s average %.3f, stddev %.3f, median %.3f, min %.3f, max %.3f ms\n",
           name,
           mean,
           stddev(samples, mean),
           median(samples),
           min,
           max);
}

static void
print_breakdown(const char *name, const breakdown_ms_t *breakdown, double total_ms, int iterations)
{
    double isogeny_ms = breakdown->isogeny_ms / iterations;
    double elliptic_curve_ms = breakdown->elliptic_curve_ms / iterations;
    double ideal_ms = breakdown->ideal_ms / iterations;
    double other_ms = breakdown->other_ms / iterations;
    double percent_scale = total_ms > 0.0 ? 100.0 / total_ms : 0.0;

    printf("  %-7s isogeny %.3f ms (%.1f%%), "
           "elliptic-curve %.3f ms (%.1f%%), ideal/quaternion %.3f ms (%.1f%%), "
           "other %.3f ms (%.1f%%)\n",
           name,
           isogeny_ms,
           isogeny_ms * percent_scale,
           elliptic_curve_ms,
           elliptic_curve_ms * percent_scale,
           ideal_ms,
           ideal_ms * percent_scale,
           other_ms,
           other_ms * percent_scale);
}

static void
print_comparison(const benchmark_result_t *original, const benchmark_result_t *new_result)
{
    const char *phases[] = { "keygen", "sign", "verify" };
    const double original_cycles[] = { original->keygen_cycles, original->sign_cycles, original->verify_cycles };
    const double new_cycles[] = { new_result->keygen_cycles, new_result->sign_cycles, new_result->verify_cycles };
    const double original_ms[] = { original->keygen_ms, original->sign_ms, original->verify_ms };
    const double new_ms[] = { new_result->keygen_ms, new_result->sign_ms, new_result->verify_ms };

    printf("Comparison (negative change means faster):\n");
    printf("  MCycles\n");
    printf("  %-7s %12s %12s %10s %10s\n", "phase", "original", "new", "change", "speedup");
    for (size_t i = 0; i < sizeof(phases) / sizeof(phases[0]); ++i) {
        double cycles_change = 100.0 * (new_cycles[i] - original_cycles[i]) / original_cycles[i];
        printf("  %-7s %12.3f %12.3f %9.2f%% %9.3fx\n",
               phases[i],
               original_cycles[i],
               new_cycles[i],
               cycles_change,
               original_cycles[i] / new_cycles[i]);
    }
    printf("\n");
    printf("  CPU ms\n");
    printf("  %-7s %12s %12s %10s %10s\n", "phase", "original", "new", "change", "speedup");
    for (size_t i = 0; i < sizeof(phases) / sizeof(phases[0]); ++i) {
        double ms_change = 100.0 * (new_ms[i] - original_ms[i]) / original_ms[i];
        printf("  %-7s %12.3f %12.3f %9.2f%% %9.3fx\n",
               phases[i],
               original_ms[i],
               new_ms[i],
               ms_change,
               original_ms[i] / new_ms[i]);
    }

    const char *breakdown_names[] = { "isogeny", "elliptic-curve", "ideal/quaternion", "other" };
    const breakdown_ms_t *original_breakdowns[] = { &original->keygen_breakdown,
                                                    &original->sign_breakdown,
                                                    &original->verify_breakdown };
    const breakdown_ms_t *new_breakdowns[] = { &new_result->keygen_breakdown,
                                               &new_result->sign_breakdown,
                                               &new_result->verify_breakdown };
    printf("\n  Breakdown averages (CPU ms)\n");
    printf("  %-7s %-16s %12s %12s %10s\n", "phase", "category", "original", "new", "change");
    printf("  --------\n");
    for (size_t phase = 0; phase < sizeof(phases) / sizeof(phases[0]); ++phase) {
        for (size_t category = 0; category < sizeof(breakdown_names) / sizeof(breakdown_names[0]); ++category) {
            const breakdown_ms_t *original_breakdown = original_breakdowns[phase];
            const breakdown_ms_t *new_breakdown = new_breakdowns[phase];
            double original_value;
            double new_value;
            switch (category) {
                case 0:
                    original_value = original_breakdown->isogeny_ms;
                    new_value = new_breakdown->isogeny_ms;
                    break;
                case 1:
                    original_value = original_breakdown->elliptic_curve_ms;
                    new_value = new_breakdown->elliptic_curve_ms;
                    break;
                case 2:
                    original_value = original_breakdown->ideal_ms;
                    new_value = new_breakdown->ideal_ms;
                    break;
                default:
                    original_value = original_breakdown->other_ms;
                    new_value = new_breakdown->other_ms;
                    break;
            }
            double change = original_value == 0.0 ? 0.0 : 100.0 * (new_value - original_value) / original_value;
            printf("  %-7s %-16s %12.3f %12.3f %9.2f%%\n",
                   phases[phase],
                   breakdown_names[category],
                   original_value,
                   new_value,
                   change);
        }
        printf("  --------\n");
    }
}

int
bench(const int iterations, const bool using_dual, benchmark_result_t *benchmark_result)
{
    int result = 0;
    unsigned char message[32] = { 0 };
    public_key_t *pks = NULL;
    secret_key_t *sks = NULL;
    signature_t *sigs = NULL;
    double *keygen_cycles = NULL;
    double *sign_cycles = NULL;
    double *verify_cycles = NULL;
    double *keygen_ms_samples = NULL;
    double *sign_ms_samples = NULL;
    double *verify_ms_samples = NULL;
    double keygen_seconds_total = 0.0;
    double sign_seconds_total = 0.0;
    double verify_seconds_total = 0.0;
    breakdown_ms_t keygen_breakdown_total = { 0 };
    breakdown_ms_t sign_breakdown_total = { 0 };
    breakdown_ms_t verify_breakdown_total = { 0 };

    keygen_cycles = calloc((size_t)iterations, sizeof(*keygen_cycles));
    sign_cycles = calloc((size_t)iterations, sizeof(*sign_cycles));
    verify_cycles = calloc((size_t)iterations, sizeof(*verify_cycles));
    keygen_ms_samples = calloc((size_t)iterations, sizeof(*keygen_ms_samples));
    sign_ms_samples = calloc((size_t)iterations, sizeof(*sign_ms_samples));
    verify_ms_samples = calloc((size_t)iterations, sizeof(*verify_ms_samples));
    pks = calloc((size_t)iterations, sizeof(*pks));
    sks = calloc((size_t)iterations, sizeof(*sks));
    sigs = calloc((size_t)iterations, sizeof(*sigs));
    if (!keygen_cycles || !sign_cycles || !verify_cycles || !keygen_ms_samples || !sign_ms_samples ||
        !verify_ms_samples || !pks || !sks || !sigs) {
        fprintf(stderr, "failed to allocate benchmark samples\n");
        free(keygen_cycles);
        free(sign_cycles);
        free(verify_cycles);
        free(keygen_ms_samples);
        free(sign_ms_samples);
        free(verify_ms_samples);
        free(pks);
        free(sks);
        free(sigs);
        return 1;
    }

    randombytes_init((unsigned char *)"some", (unsigned char *)"string", 128);
    for (int i = 0; i < iterations; ++i) {
        secret_key_init(&sks[i]);
        secret_sig_init(&sigs[i]);
    }

    /* Generate all key pairs before starting the signing phase. */
    for (int i = 0; i < iterations; ++i) {
        keygen_timings_t keygen_timings;
        unsigned long long cycle_start = read_cycles();
        double wall_start = monotonic_ms();
        clock_t cpu_start = clock();
        protocols_keygen(&pks[i], &sks[i], &keygen_timings);
        unsigned long long cycle_end = read_cycles();
        keygen_ms_samples[i] = monotonic_ms() - wall_start;
        double keygen_ms = elapsed_ms(cpu_start, clock());
        keygen_seconds_total += keygen_ms / 1000.0;
        keygen_cycles[i] = (double)(cycle_end - cycle_start) / 1.0e6;
        breakdown_ms_t current = keygen_breakdown(&keygen_timings);
        if (!finish_breakdown(&current, keygen_ms, "keygen", i)) {
            result = 1;
            goto cleanup;
        }
        accumulate_breakdown(&keygen_breakdown_total, &current);
    }

    printf("keygen completed\n");

    /* Sign all messages after key generation has completed. */
    for (int i = 0; i < iterations; ++i) {
        int ret;
        sign_timings_t sign_timings;
        unsigned long long cycle_start = read_cycles();
        unsigned long long cycle_end;
        double wall_start = monotonic_ms();
        clock_t cpu_start = clock();
        if (using_dual) {
            ret = protocols_sign_new(&sigs[i], &pks[i], &sks[i], message, sizeof(message), 0, &sign_timings);
        } else {
            ret = protocols_sign(&sigs[i], &pks[i], &sks[i], message, sizeof(message), 0, &sign_timings);
        }
        if (ret) {
            fprintf(stderr, "signing failed on iteration %d\n", i);
            result = 1;
            goto cleanup;
        }
        cycle_end = read_cycles();
        sign_ms_samples[i] = monotonic_ms() - wall_start;
        sign_cycles[i] = (double)(cycle_end - cycle_start) / 1.0e6;
        double sign_ms = elapsed_ms(cpu_start, clock());
        sign_seconds_total += sign_ms / 1000.0;
        breakdown_ms_t sign_current = sign_breakdown(&sign_timings);
        if (!finish_breakdown(&sign_current, sign_ms, "sign", i)) {
            result = 1;
            goto cleanup;
        }
        accumulate_breakdown(&sign_breakdown_total, &sign_current);
    }

    printf("sign completed\n");

    /* Verify all signatures after the signing phase has completed. */
    for (int i = 0; i < iterations; ++i) {
        verify_timings_t verify_timings;
        unsigned long long cycle_start = read_cycles();
        unsigned long long cycle_end;
        double wall_start = monotonic_ms();
        clock_t cpu_start = clock();
        if (!protocols_verify(&sigs[i], &pks[i], message, sizeof(message), &verify_timings)) {
            fprintf(stderr, "verification failed on iteration %d\n", i);
            result = 1;
            goto cleanup;
        }
        cycle_end = read_cycles();
        verify_ms_samples[i] = monotonic_ms() - wall_start;
        verify_cycles[i] = (double)(cycle_end - cycle_start) / 1.0e6;
        double verify_ms = elapsed_ms(cpu_start, clock());
        verify_seconds_total += verify_ms / 1000.0;
        breakdown_ms_t verify_current = verify_breakdown(&verify_timings);
        if (!finish_breakdown(&verify_current, verify_ms, "verify", i)) {
            result = 1;
            goto cleanup;
        }
        accumulate_breakdown(&verify_breakdown_total, &verify_current);
    }

    printf("verify completed\n\n");

    printf("SQIsign2DPush benchmark (%d iterations)\n", iterations);
    printf("CPU cycles: %s\n", HAVE_CPU_CYCLES ? "rdtsc" : "clock fallback");
    print_stats("keygen", &(samples_t){ keygen_cycles, iterations });
    print_stats("sign", &(samples_t){ sign_cycles, iterations });
    print_stats("verify", &(samples_t){ verify_cycles, iterations });
    printf("Wall-clock time:\n");
    print_time_stats("keygen", &(samples_t){ keygen_ms_samples, iterations });
    print_time_stats("sign", &(samples_t){ sign_ms_samples, iterations });
    print_time_stats("verify", &(samples_t){ verify_ms_samples, iterations });
    double average_keygen_ms = 1000.0 * keygen_seconds_total / iterations;
    double average_sign_ms = 1000.0 * sign_seconds_total / iterations;
    double average_verify_ms = 1000.0 * verify_seconds_total / iterations;
    printf("  average time: keygen %.3f ms, sign %.3f ms, "
           "verify %.3f ms\n",
           average_keygen_ms,
           average_sign_ms,
           average_verify_ms);
    printf("\n  breakdown averages [CPU ms (total percentage)]\n");
    printf("  other = unclassified time, hashing, management and measurement overhead\n");
    print_breakdown("keygen", &keygen_breakdown_total, average_keygen_ms, iterations);
    print_breakdown("sign", &sign_breakdown_total, average_sign_ms, iterations);
    print_breakdown("verify", &verify_breakdown_total, average_verify_ms, iterations);
    benchmark_result->keygen_cycles = average(&(samples_t){ keygen_cycles, iterations });
    benchmark_result->sign_cycles = average(&(samples_t){ sign_cycles, iterations });
    benchmark_result->verify_cycles = average(&(samples_t){ verify_cycles, iterations });
    benchmark_result->keygen_ms = average_keygen_ms;
    benchmark_result->sign_ms = average_sign_ms;
    benchmark_result->verify_ms = average_verify_ms;
    benchmark_result->keygen_breakdown = (breakdown_ms_t){ keygen_breakdown_total.isogeny_ms / iterations,
                                                           keygen_breakdown_total.elliptic_curve_ms / iterations,
                                                           keygen_breakdown_total.ideal_ms / iterations,
                                                           keygen_breakdown_total.other_ms / iterations };
    benchmark_result->sign_breakdown = (breakdown_ms_t){ sign_breakdown_total.isogeny_ms / iterations,
                                                         sign_breakdown_total.elliptic_curve_ms / iterations,
                                                         sign_breakdown_total.ideal_ms / iterations,
                                                         sign_breakdown_total.other_ms / iterations };
    benchmark_result->verify_breakdown = (breakdown_ms_t){ verify_breakdown_total.isogeny_ms / iterations,
                                                           verify_breakdown_total.elliptic_curve_ms / iterations,
                                                           verify_breakdown_total.ideal_ms / iterations,
                                                           verify_breakdown_total.other_ms / iterations };
cleanup:
    for (int i = 0; i < iterations; ++i) {
        secret_key_finalize(&sks[i]);
        secret_sig_finalize(&sigs[i]);
    }
    free(pks);
    free(sks);
    free(sigs);
    free(keygen_cycles);
    free(sign_cycles);
    free(verify_cycles);
    free(keygen_ms_samples);
    free(sign_ms_samples);
    free(verify_ms_samples);
    return result;
}

int
main(int argc, char **argv)
{
    int iterations = 100;
    const char *method = "both";
    bool method_set = false;
    bool iterations_set = false;
    bool run_original = true;
    bool run_new = true;
    benchmark_result_t original_result = { 0 };
    benchmark_result_t new_result = { 0 };

    for (int arg = 1; arg < argc; ++arg) {
        if (!strncmp(argv[arg], "--iterations=", strlen("--iterations="))) {
            char *end = NULL;
            long parsed_iterations;

            if (iterations_set) {
                fprintf(stderr, "usage: %s [--iterations=num] [--method=ori|new|both]\n", argv[0]);
                return 1;
            }
            errno = 0;
            parsed_iterations = strtol(argv[arg] + strlen("--iterations="), &end, 10);
            if (errno == 0 && end != argv[arg] + strlen("--iterations=") && *end == '\0' && parsed_iterations > 0 &&
                parsed_iterations <= INT_MAX) {
                iterations = (int)parsed_iterations;
                iterations_set = true;
                continue;
            }
        }

        if (!strncmp(argv[arg], "--method=", strlen("--method="))) {
            if (method_set) {
                fprintf(stderr, "usage: %s [--iterations=num] [--method=ori|new|both]\n", argv[0]);
                return 1;
            }
            method = argv[arg] + strlen("--method=");
            if (!strcmp(method, "ori") || !strcmp(method, "new") || !strcmp(method, "both")) {
                method_set = true;
                continue;
            }
        }

        fprintf(stderr, "usage: %s [--iterations=num] [--method=ori|new|both]\n", argv[0]);
        return 1;
    }

    run_original = !strcmp(method, "ori") || !strcmp(method, "both");
    run_new = !strcmp(method, "new") || !strcmp(method, "both");

    if (run_original) {
        printf("------------------------------------------------\n");
        printf("Running SQIsign2DPush benchmark for original algorithm with %d iterations\n", iterations);
        if (bench(iterations, false, &original_result)) {
            fprintf(stderr, "Benchmark for original algorithm failed\n");
            return 1;
        }
        printf("------------------------------------------------\n\n");
    }

    if (run_new) {
        printf("------------------------------------------------\n");
        printf("Running SQIsign2DPush benchmark for new algorithm with %d iterations\n", iterations);
        if (bench(iterations, true, &new_result)) {
            fprintf(stderr, "Benchmark for new algorithm failed\n");
            return 1;
        }
        printf("------------------------------------------------\n\n");
    }

    if (run_original && run_new) {
        printf("------------------------------------------------\n");
        print_comparison(&original_result, &new_result);
        printf("------------------------------------------------\n\n");
    }

    return 0;
}