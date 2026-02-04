#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "batch.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* --------------- Single Sequence Processing --------------- */

BatchResult* process_single_sequence(const char *seq_file, BatchConfig *config) {
    BatchResult *result = calloc(1, sizeof(BatchResult));
    result->sequence_file = strdup(seq_file);
    result->success = 0;

    // Initialize data structures
    Observed_events info;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;
    Pos_prob pos;

    memset(&info, 0, sizeof(Observed_events));
    info.flank = config->flank_size;
    memset(&ed, 0, sizeof(Explicit_duration));
    ed.min_len_exon = -1;
    ed.min_len_intron = -1;
    memset(&l, 0, sizeof(Lambda));

    // Parse sequence
    read_sequence_file(seq_file, &info);
    if (!info.original_sequence) {
        result->error_message = strdup("Failed to read sequence file");
        return result;
    }
    numerical_transcription(&info, info.original_sequence);

    // Flank size check
    if (info.flank * 2 >= info.T) {
        result->error_message = strdup("Flank size too large for sequence");
        free(info.original_sequence);
        free(info.numerical_sequence);
        return result;
    }

    // Parse model
    if (config->model_file) {
        if (parse_json_model(config->model_file, &l, &ed) != 0) {
            result->error_message = strdup("Failed to parse model file");
            free(info.original_sequence);
            free(info.numerical_sequence);
            return result;
        }
    } else {
        donor_parser(&l, config->don_emission);
        acceptor_parser(&l, config->acc_emission);
        exon_intron_parser(&l, config->exon_emission, 0);
        exon_intron_parser(&l, config->intron_emission, 1);
        explicit_duration_probability(&ed, config->ped_exon, 0);
        explicit_duration_probability(&ed, config->ped_intron, 1);
    }

    // Compute transition matrices
    if (config->model_file) {
        compute_transition_matrices_log(&l);
    } else {
        compute_transition_matrices(&l);
    }

    l.log_values_len = (ed.max_len_exon > ed.max_len_intron) ? ed.max_len_exon : ed.max_len_intron;
    l.log_values = calloc(l.log_values_len, sizeof(double));

    // Log space conversion (if not already in log space)
    if (!config->model_file) {
        int don_size = power(4, l.B.don_kmer_len);
        int acc_size = power(4, l.B.acc_kmer_len);
        int exon_size = power(4, l.B.exon_kmer_len);
        int intron_size = power(4, l.B.intron_kmer_len);

        tolerance_checker(ed.exon, ed.exon_len, 1e-15);
        tolerance_checker(ed.intron, ed.intron_len, 1e-15);
        tolerance_checker(l.A.dons, don_size, 1e-15);
        tolerance_checker(l.A.accs, acc_size, 1e-15);

        log_space_converter(ed.exon, ed.exon_len);
        log_space_converter(ed.intron, ed.intron_len);
        log_space_converter(l.A.dons, don_size);
        log_space_converter(l.A.accs, acc_size);
        log_space_converter(l.B.exon, exon_size);
        log_space_converter(l.B.intron, intron_size);
    }

    // Initialize k-mer cache
    info.kmer_cache = allocate_kmer_cache(
        info.T,
        l.B.exon_kmer_len,
        l.B.intron_kmer_len,
        l.B.don_kmer_len,
        l.B.acc_kmer_len
    );
    if (!info.kmer_cache) {
        result->error_message = strdup("Failed to allocate k-mer cache");
        goto cleanup;
    }
    init_kmer_cache(info.kmer_cache, info.numerical_sequence);

    // Allocate algorithm structures
    allocate_fw(&info, &fw, &ed);
    allocate_bw(&bw, &ed, &info);
    allocate_pos(&pos, &info);

    // Forward-backward algorithm
    basis_fw_algo(&l, &ed, &fw, &info);
    fw_algo(&l, &fw, &info, &ed);
    basis_bw_algo(&l, &bw, &info, &ed);
    bw_algo(&l, &bw, &info, &ed);

    // Posterior probabilities
    pos_prob(&bw, &fw, &info, &pos);
    parse_splice_sites(&pos, &info);

    result->n_donors = pos.dons;
    result->n_acceptors = pos.accs;

    // MCTS isoform generation (optional)
    if (config->use_mcts && pos.dons > 0 && pos.accs > 0) {
        Locus *loc = create_locus(config->n_isoforms);

        MCTSTree *tree = create_mcts_tree(&pos, &ed, &info, config->mcts_explore_c);
        int max_iterations = config->n_isoforms * 10;
        generate_isoforms_mcts(tree, loc, max_iterations);

        result->n_isoforms = loc->n_isoforms;
        result->locus = loc;  // Transfer ownership

        free_mcts_tree(tree);
    }

    result->success = 1;

cleanup:
    // Cleanup
    free_splice_sites(&pos);
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free_pos(&pos, &info);
    free(info.original_sequence);
    free(info.numerical_sequence);
    free_kmer_cache(info.kmer_cache);
    free_lambda(&l);
    free_explicit_duration(&ed);

    return result;
}

/* --------------- Parallel Batch Processing --------------- */

BatchResult** process_batch_parallel(char **seq_files, int n_files, BatchConfig *config, int n_threads) {
    BatchResult **results = malloc(n_files * sizeof(BatchResult*));

#ifdef _OPENMP
    if (n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
    printf("Processing %d sequences with %d threads...\n", n_files, omp_get_max_threads());

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_files; i++) {
        int tid = omp_get_thread_num();
        printf("  [Thread %d] Processing: %s\n", tid, seq_files[i]);
        results[i] = process_single_sequence(seq_files[i], config);
    }
#else
    printf("Processing %d sequences (OpenMP not available, single-threaded)...\n", n_files);
    for (int i = 0; i < n_files; i++) {
        printf("  Processing: %s\n", seq_files[i]);
        results[i] = process_single_sequence(seq_files[i], config);
    }
#endif

    return results;
}

/* --------------- Batch File Reading --------------- */

char** read_batch_file(const char *batch_file, int *n_files) {
    FILE *fp = fopen(batch_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open batch file: %s\n", batch_file);
        *n_files = 0;
        return NULL;
    }

    // Count lines first
    int count = 0;
    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        // Skip empty lines and comments
        if (line[0] != '\n' && line[0] != '#' && line[0] != '\0') {
            count++;
        }
    }

    if (count == 0) {
        fclose(fp);
        *n_files = 0;
        return NULL;
    }

    // Allocate and read
    char **files = malloc(count * sizeof(char*));
    rewind(fp);

    int idx = 0;
    while (fgets(line, sizeof(line), fp) && idx < count) {
        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0') continue;

        // Remove trailing newline
        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            line[len-1] = '\0';
        }

        files[idx++] = strdup(line);
    }

    fclose(fp);
    *n_files = idx;
    return files;
}

/* --------------- Cleanup --------------- */

void free_batch_result(BatchResult *result) {
    if (!result) return;
    if (result->sequence_file) free(result->sequence_file);
    if (result->error_message) free(result->error_message);
    if (result->locus) free_locus(result->locus);
    free(result);
}

void free_batch_results(BatchResult **results, int n_results) {
    if (!results) return;
    for (int i = 0; i < n_results; i++) {
        free_batch_result(results[i]);
    }
    free(results);
}

/* --------------- Summary Output --------------- */

void print_batch_summary(BatchResult **results, int n_results) {
    printf("\n=== Batch Processing Summary ===\n");

    int success_count = 0;
    int total_donors = 0;
    int total_acceptors = 0;
    int total_isoforms = 0;

    for (int i = 0; i < n_results; i++) {
        BatchResult *r = results[i];
        if (r->success) {
            success_count++;
            total_donors += r->n_donors;
            total_acceptors += r->n_acceptors;
            total_isoforms += r->n_isoforms;
            printf("  [OK] %s: %d donors, %d acceptors",
                   r->sequence_file, r->n_donors, r->n_acceptors);
            if (r->n_isoforms > 0) {
                printf(", %d isoforms", r->n_isoforms);
            }
            printf("\n");
        } else {
            printf("  [FAIL] %s: %s\n", r->sequence_file,
                   r->error_message ? r->error_message : "Unknown error");
        }
    }

    printf("\n--- Totals ---\n");
    printf("  Sequences: %d/%d succeeded\n", success_count, n_results);
    printf("  Donors: %d\n", total_donors);
    printf("  Acceptors: %d\n", total_acceptors);
    if (total_isoforms > 0) {
        printf("  Isoforms: %d\n", total_isoforms);
    }
    printf("==============================\n");
}
