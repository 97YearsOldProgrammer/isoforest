#ifndef HMM_LOG_MATH_H
#define HMM_LOG_MATH_H

#include <math.h>
#include "constants.h"

/* --------------- Log-Space Helper Macros --------------- */
// Check if a log-space value represents zero probability
#define IS_LOG_ZERO(x) ((x) <= LOG_ZERO + 1.0)

// Safe log: returns LOG_ZERO for non-positive values
#define SAFE_LOG(x) ((x) > 0.0 ? log(x) : LOG_ZERO)

// Safe exp: prevents overflow/underflow
#define SAFE_EXP(x) ((x) > LOG_EPSILON ? exp(x) : 0.0)

// Log-space multiplication: simply addition in log space
#define LOG_MUL(a, b) ((IS_LOG_ZERO(a) || IS_LOG_ZERO(b)) ? LOG_ZERO : ((a) + (b)))

/* --------------- Log-Space Inline Functions --------------- */
// Log-space addition: log(exp(a) + exp(b)) with numerical stability
static inline double log_add(double a, double b) {
    if (IS_LOG_ZERO(a)) return b;
    if (IS_LOG_ZERO(b)) return a;
    if (a > b) return a + log1p(exp(b - a));
    else return b + log1p(exp(a - b));
}

/* --------------- Utility Functions --------------- */
int    power(int base, int exp);
int    base4_to_int(int *array, int beg, int length);
double total_prob(double *array, int length);
double log_sum_exp(double *logs, int n);
void   tolerance_checker(double *array, int len, const double epsilon);
void   log_space_converter(double *array, int len);

/* --------------- K-mer Cache Functions --------------- */
#include "types.h"
#include "model_structs.h"

KmerCache *allocate_kmer_cache(int T, int exon_k, int intron_k, int donor_k, int acceptor_k);
void       init_kmer_cache(KmerCache *cache, int *numerical_sequence);
void       free_kmer_cache(KmerCache *cache);

/* Inline cache lookup - O(1) instead of O(k) */
static inline int get_exon_kmer(KmerCache *cache, int pos) {
    int start = pos - cache->exon_kmer_len + 1;
    if (start < 0 || pos >= cache->T) return -1;
    return cache->exon[pos];
}

static inline int get_intron_kmer(KmerCache *cache, int pos) {
    int start = pos - cache->intron_kmer_len + 1;
    if (start < 0 || pos >= cache->T) return -1;
    return cache->intron[pos];
}

static inline int get_donor_kmer(KmerCache *cache, int pos) {
    if (pos < 0 || pos + cache->donor_kmer_len > cache->T) return -1;
    return cache->donor[pos];
}

static inline int get_acceptor_kmer(KmerCache *cache, int pos) {
    int start = pos - cache->acceptor_kmer_len + 1;
    if (start < 0 || pos >= cache->T) return -1;
    return cache->acceptor[pos];
}

#endif
