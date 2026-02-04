#ifndef BATCH_H
#define BATCH_H

#include "model.h"
#include "decoder/mcts.h"

/* --------------- Batch Processing Configuration --------------- */
typedef struct {
    char *don_emission;
    char *acc_emission;
    char *exon_emission;
    char *intron_emission;
    char *ped_exon;
    char *ped_intron;
    char *model_file;
    int flank_size;
    int use_mcts;
    int n_isoforms;
    double mcts_explore_c;
    int output_json;
} BatchConfig;

/* --------------- Batch Processing Results --------------- */
typedef struct {
    char *sequence_file;
    int n_donors;
    int n_acceptors;
    int n_isoforms;
    Locus *locus;           // NULL if not using random forest
    int success;            // 1 = success, 0 = failure
    char *error_message;    // NULL if success
} BatchResult;

/* --------------- Batch Processing API --------------- */

// Process a single sequence file (thread-safe, can be called from OpenMP)
BatchResult* process_single_sequence(const char *seq_file, BatchConfig *config);

// Process multiple sequences in parallel
BatchResult** process_batch_parallel(char **seq_files, int n_files, BatchConfig *config, int n_threads);

// Read sequence file list from batch file (one path per line)
char** read_batch_file(const char *batch_file, int *n_files);

// Free batch results
void free_batch_result(BatchResult *result);
void free_batch_results(BatchResult **results, int n_results);

// Print batch results summary
void print_batch_summary(BatchResult **results, int n_results);

#endif
