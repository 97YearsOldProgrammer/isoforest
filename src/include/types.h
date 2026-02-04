#ifndef HMM_TYPES_H
#define HMM_TYPES_H

/* --------------- K-mer Index Cache --------------- */
/*
 * Precomputed k-mer indices using sliding window algorithm.
 * Instead of calling base4_to_int() O(T * k) times in hot loops,
 * we precompute all indices once in O(T) total time.
 *
 * Sliding window formula:
 *   new_index = ((old_index & mask) << 2) | new_base
 *   where mask = 4^(k-1) - 1
 */
typedef struct {
    int *exon;          // indices for exon k-mer at each position
    int *intron;        // indices for intron k-mer at each position
    int *donor;         // indices for donor k-mer at each position
    int *acceptor;      // indices for acceptor k-mer at each position
    int exon_kmer_len;
    int intron_kmer_len;
    int donor_kmer_len;
    int acceptor_kmer_len;
    int T;              // sequence length
} KmerCache;

/* --------------- Core Input Structure --------------- */
typedef struct
{
    char *original_sequence;            // org seq
    int T;                              // seq len
    int *numerical_sequence;            // org seq to num seq
    int flank;                          // flank size if provided
    KmerCache *kmer_cache;              // precomputed k-mer indices
} Observed_events;

#endif
