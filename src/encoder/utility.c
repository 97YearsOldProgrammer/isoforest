#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* --------------- Sequence Transcription --------------- */

void numerical_transcription(Observed_events *info, const char *seq) {
    if(DEBUG)  printf("Start transforming original sequence into base4:\n");
    // turns original sequence into int
    size_t len  = strlen(seq);
    info->T     = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for (size_t i = 0; i < len; i++) {
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }

    if(DEBUG)  printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    if(DEBUG)  printf("\tFinished\n");
    if(DEBUG)  printf("\n");
}

/* --------------- Math Utilities --------------- */

int power(int base, int exp) {
    int result = 1;
    for (int i = 0 ; i < exp ; i++) {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int beg, int length) {
    int values = 0;
    int value;
    for (int i = 0; i < length; i++) {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    return values;
}

double log_sum_exp(double *array, int n) {
    if(!array){
        printf("Something wrong with log_sum_exp trick. Invalid Input\n");
        return LOG_ZERO;
    }
    if(n <= 0) return LOG_ZERO;

    // Find max, skipping LOG_ZERO values
    double max = LOG_ZERO;
    int valid_count = 0;
    for (int i = 0; i < n; i++) {
        if (!IS_LOG_ZERO(array[i])) {
            if (valid_count == 0 || array[i] > max) {
                max = array[i];
            }
            valid_count++;
        }
    }

    if (valid_count == 0) return LOG_ZERO;

    // Compute sum in numerically stable way
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (!IS_LOG_ZERO(array[i])) {
            sum += exp(array[i] - max);
        }
    }
    return max + log(sum);
}

/* --------------- Validators --------------- */

void tolerance_checker(double *array, int len, const double epsilon) {
    if(!array){
        printf("This is not a valid input for float20.0\n");
        return;
    }

    for(int i = 0 ; i < len ; i++) {
        if( array[i] < epsilon )    array[i] = 0;
    }
}

void log_space_converter(double *array, int len) {
    if(!array){
        printf("This is not a valid input for log_space converter\n");
        return;
    }

    for (int i = 0; i < len ; i++) {
        if (array[i] <= PROB_EPSILON)   array[i] = LOG_ZERO;
        else                            array[i] = log(array[i]);
    }
}

/* --------------- K-mer Cache Implementation --------------- */
/*
 * Sliding window algorithm for O(T) precomputation instead of O(T*k)
 *
 * For a k-mer starting at position i:
 *   index = seq[i]*4^(k-1) + seq[i+1]*4^(k-2) + ... + seq[i+k-1]*4^0
 *
 * Sliding: index[i+1] = ((index[i] & mask) << 2) | seq[i+k]
 *   where mask = 4^(k-1) - 1
 */

KmerCache *allocate_kmer_cache(int T, int exon_k, int intron_k, int donor_k, int acceptor_k) {
    KmerCache *cache = malloc(sizeof(KmerCache));
    if (!cache) return NULL;

    cache->T = T;
    cache->exon_kmer_len = exon_k;
    cache->intron_kmer_len = intron_k;
    cache->donor_kmer_len = donor_k;
    cache->acceptor_kmer_len = acceptor_k;

    cache->exon = malloc(T * sizeof(int));
    cache->intron = malloc(T * sizeof(int));
    cache->donor = malloc(T * sizeof(int));
    cache->acceptor = malloc(T * sizeof(int));

    if (!cache->exon || !cache->intron || !cache->donor || !cache->acceptor) {
        free_kmer_cache(cache);
        return NULL;
    }

    // Initialize to -1 (invalid)
    for (int i = 0; i < T; i++) {
        cache->exon[i] = -1;
        cache->intron[i] = -1;
        cache->donor[i] = -1;
        cache->acceptor[i] = -1;
    }

    return cache;
}

/* Compute k-mer indices using sliding window - O(T) total */
static void compute_kmer_indices_sliding(int *cache, int *seq, int T, int k) {
    if (k <= 0 || T < k) return;

    int mask = (1 << (2 * (k - 1))) - 1;  // 4^(k-1) - 1

    // Compute first k-mer using base4_to_int
    int idx = base4_to_int(seq, 0, k);
    cache[k - 1] = idx;  // Store at position of last base in k-mer

    // Slide through sequence - O(1) per position
    for (int i = k; i < T; i++) {
        idx = ((idx & mask) << 2) | seq[i];
        cache[i] = idx;
    }
}

/* Compute k-mer indices for "start at position" style (donor) */
static void compute_kmer_indices_start(int *cache, int *seq, int T, int k) {
    if (k <= 0 || T < k) return;

    int mask = (1 << (2 * (k - 1))) - 1;

    // Compute first k-mer
    int idx = base4_to_int(seq, 0, k);
    cache[0] = idx;  // Store at start position

    // Slide through sequence
    for (int i = 1; i <= T - k; i++) {
        idx = ((idx & mask) << 2) | seq[i + k - 1];
        cache[i] = idx;
    }
}

void init_kmer_cache(KmerCache *cache, int *numerical_sequence) {
    if (!cache || !numerical_sequence) return;

    int T = cache->T;

    // Exon/intron: indexed by END position (bps-okmer+1 to bps)
    compute_kmer_indices_sliding(cache->exon, numerical_sequence, T, cache->exon_kmer_len);
    compute_kmer_indices_sliding(cache->intron, numerical_sequence, T, cache->intron_kmer_len);

    // Donor: indexed by START position (bps to bps+dkmer-1)
    compute_kmer_indices_start(cache->donor, numerical_sequence, T, cache->donor_kmer_len);

    // Acceptor: indexed by END position (bps-akmer+1 to bps) - same as sliding
    compute_kmer_indices_sliding(cache->acceptor, numerical_sequence, T, cache->acceptor_kmer_len);

    if (DEBUG) {
        printf("K-mer cache initialized: T=%d, exon_k=%d, intron_k=%d, donor_k=%d, acc_k=%d\n",
               T, cache->exon_kmer_len, cache->intron_kmer_len,
               cache->donor_kmer_len, cache->acceptor_kmer_len);
    }
}

void free_kmer_cache(KmerCache *cache) {
    if (!cache) return;
    if (cache->exon) free(cache->exon);
    if (cache->intron) free(cache->intron);
    if (cache->donor) free(cache->donor);
    if (cache->acceptor) free(cache->acceptor);
    free(cache);
}
