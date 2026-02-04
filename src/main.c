#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "model.h"
#include "decoder/mcts.h"
#include "batch.h"

int     DEBUG                   = 0;
int     use_mcts                = 0;
int     n_isoforms              = 10000;
int     model_in_log_space      = 0;

void print_usage(const char *program_name) {
    printf("RFHMM - Random Forest Hidden Markov Model for gene prediction\n");
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    printf("Required options (use ONE):\n");
    printf("  -s, --sequence FILE           Input sequence file (single mode)\n");
    printf("  -B, --batch FILE              Batch file with sequence paths (parallel mode)\n");
    printf("\nModel options (use ONE of the following):\n");
    printf("  -M, --model FILE              JSON model file (.splicemodel)\n");
    printf("  OR individual model files:\n");
    printf("  -d, --don_emission FILE       Donor emission file (default: ../models/don.pwm)\n");
    printf("  -a, --acc_emission FILE       Acceptor emission file (default: ../models/acc.pwm)\n");
    printf("  -e, --exon_emission FILE      Exon emission file (default: ../models/exon.mm)\n");
    printf("  -i, --intron_emission FILE    Intron emission file (default: ../models/intron.mm)\n");
    printf("  -x, --ped_exon FILE           Exon length distribution file (default: ../models/exon.len)\n");
    printf("  -n, --ped_intron FILE         Intron length distribution file (default: ../models/intron.len)\n");
    printf("\nAnalysis parameters:\n");
    printf("  -f, --flank NUM               Flank size for analysis (default: 99)\n");
    printf("\nMCTS isoform generation options:\n");
    printf("  -T, --mcts                    Enable Monte Carlo Tree Search (MCTS) algorithm\n");
    printf("  -N, --n_isoforms NUM          Maximum isoform capacity (default: 10000)\n");
    printf("  -C, --explore FLOAT           MCTS exploration constant (default: 1.4)\n");
    printf("\nParallel processing (with --batch):\n");
    printf("  -t, --threads NUM             Number of threads (default: all available cores)\n");
    printf("\nOutput control:\n");
    printf("  -j, --json                    Emit isoform locus information as JSON\n");
    printf("  -v, --verbose                 Show debug and progress information\n");
    printf("  -h, --help                    Show this help message\n");
    printf("\nExamples:\n");
    printf("  %s -s input.fasta --mcts\n", program_name);
    printf("  %s -s input.fasta --mcts -N 100\n", program_name);
    printf("  %s -s input.fasta --mcts -N 500 -C 2.0\n", program_name);
    printf("  %s --batch sequences.txt --threads 8 --mcts\n", program_name);
    printf("\nBatch file format (one sequence path per line):\n");
    printf("  # This is a comment\n");
    printf("  /path/to/sequence1.fasta\n");
    printf("  /path/to/sequence2.fasta\n");
    printf("\nMCTS Algorithm:\n");
    printf("  Monte Carlo Tree Search explores isoform combinations using UCB1.\n");
    printf("  It learns which splice site combinations produce valid isoforms\n");
    printf("  and focuses exploration on promising paths.\n");
}

int main(int argc, char *argv[])
{
    // Hard Code Path for Default(shall be deleted)
    char *default_don_emission      = "../models/don.pwm";
    char *default_acc_emission      = "../models/acc.pwm";
    char *default_exon_emission     = "../models/exon.mm";
    char *default_intron_emission   = "../models/intron.mm";
    char *default_Ped_exon          = "../models/exon.len";
    char *default_Ped_intron        = "../models/intron.len";

    // Variables for command-line inputs
    char *don_emission          = default_don_emission;
    char *acc_emission          = default_acc_emission;
    char *exon_emission         = default_exon_emission;
    char *intron_emission       = default_intron_emission;
    char *Ped_exon              = default_Ped_exon;
    char *Ped_intron            = default_Ped_intron;
    char *seq_input             = NULL;
    char *batch_file            = NULL;
    char *model_file            = NULL;
    int output_json             = 0;
    int flank_size              = DEFAULT_FLANK;
    int n_threads               = 0;            // 0 = use all available
    double mcts_explore_c       = 1.4;          // UCB1 exploration constant

    static struct option long_options[] = {
        {"sequence",        required_argument, 0, 's'},
        {"batch",           required_argument, 0, 'B'},
        {"model",           required_argument, 0, 'M'},
        {"don_emission",    required_argument, 0, 'd'},
        {"acc_emission",    required_argument, 0, 'a'},
        {"exon_emission",   required_argument, 0, 'e'},
        {"intron_emission", required_argument, 0, 'i'},
        {"ped_exon",        required_argument, 0, 'x'},
        {"ped_intron",      required_argument, 0, 'n'},
        {"flank",           required_argument, 0, 'f'},
        {"n_isoforms",      required_argument, 0, 'N'},
        {"threads",         required_argument, 0, 't'},
        {"explore",         required_argument, 0, 'C'},
        {"json",            no_argument,       0, 'j'},
        {"mcts",            no_argument,       0, 'T'},
        {"verbose",         no_argument,       0, 'v'},
        {"help",            no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "s:B:M:d:a:e:i:x:n:f:N:t:C:jTvh", long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                seq_input = optarg;
                break;
            case 'B':
                batch_file = optarg;
                break;
            case 't':
                n_threads = atoi(optarg);
                if (n_threads < 0) {
                    fprintf(stderr, "Error: threads must be non-negative\n");
                    return 1;
                }
                break;
            case 'M':
                model_file = optarg;
                break;
            case 'd':
                don_emission = optarg;
                break;
            case 'a':
                acc_emission = optarg;
                break;
            case 'e':
                exon_emission = optarg;
                break;
            case 'i':
                intron_emission = optarg;
                break;
            case 'x':
                Ped_exon = optarg;
                break;
            case 'n':
                Ped_intron = optarg;
                break;
            case 'f':
                flank_size = atoi(optarg);
                if (flank_size < 0) {
                    fprintf(stderr, "Error: flank size must be non-negative\n");
                    return 1;
                }
                break;
            case 'v':
                DEBUG = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'N':
                n_isoforms = atoi(optarg);
                if (n_isoforms < 1) {
                    fprintf(stderr, "Error: n_isoforms must be at least 1\n");
                    return 1;
                }
                break;
            case 'j':
                output_json = 1;
                break;
            case 'T':
                use_mcts = 1;
                break;
            case 'C':
                mcts_explore_c = atof(optarg);
                if (mcts_explore_c <= 0.0) {
                    fprintf(stderr, "Error: exploration constant must be positive\n");
                    return 1;
                }
                break;
            case '?':
                // getopt_long already printed an error message
                print_usage(argv[0]);
                return 1;
            default:
                abort();
        }
    }

    if (seq_input == NULL && batch_file == NULL) {
        fprintf(stderr, "Error: --sequence or --batch is required\n");
        print_usage(argv[0]);
        return 1;
    }

    if (seq_input != NULL && batch_file != NULL) {
        fprintf(stderr, "Error: Cannot use both --sequence and --batch\n");
        return 1;
    }

    /* --------------- Batch Processing Mode --------------- */
    if (batch_file != NULL) {
        BatchConfig config = {
            .don_emission = don_emission,
            .acc_emission = acc_emission,
            .exon_emission = exon_emission,
            .intron_emission = intron_emission,
            .ped_exon = Ped_exon,
            .ped_intron = Ped_intron,
            .model_file = model_file,
            .flank_size = flank_size,
            .use_mcts = use_mcts,
            .n_isoforms = n_isoforms,
            .mcts_explore_c = mcts_explore_c,
            .output_json = output_json
        };

        int n_files;
        char **seq_files = read_batch_file(batch_file, &n_files);
        if (!seq_files || n_files == 0) {
            fprintf(stderr, "Error: No sequences found in batch file\n");
            return 1;
        }

        BatchResult **results = process_batch_parallel(seq_files, n_files, &config, n_threads);
        print_batch_summary(results, n_files);

        // Cleanup
        free_batch_results(results, n_files);
        for (int i = 0; i < n_files; i++) {
            free(seq_files[i]);
        }
        free(seq_files);

        return 0;
    }

    /* --------------- Single Sequence Mode --------------- */
    /* --------------- Initialize Data Structure --------------- */
    Observed_events info;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;
    Pos_prob pos;

    memset(&info, 0, sizeof(Observed_events));
    info.flank = flank_size;
    memset(&ed, 0, sizeof(Explicit_duration));
    ed.min_len_exon     = -1;
    ed.min_len_intron   = -1;
    ed.max_len_exon     = 0;
    ed.max_len_intron   = 0;
    memset(&l, 0, sizeof(Lambda));

    if (DEBUG) printf("\n=== Starting RFHMM Analysis ===\n");

    /* --------------- Parse Sequence --------------- */
    if (DEBUG) printf("\n--- Phase 1: Loading Sequence ---\n");
    read_sequence_file(seq_input, &info);
    numerical_transcription(&info, info.original_sequence);
    
    // FLANK Size Check
    if (info.flank * 2 >= info.T) {
        fprintf(stderr, "Error: Flank size (%d) is too large for sequence length (%d)\n", 
                info.flank, info.T);
        fprintf(stderr, "Flank size should be less than %d\n", info.T / 2);
        return 1;
    }

    /* --------------- Parse HMM Model Input --------------- */
    if (DEBUG) printf("\n--- Phase 2: Parsing Model Files ---\n");

    if (model_file) {
        // Use unified JSON model file (values already in log-space)
        if (parse_json_model(model_file, &l, &ed) != 0) {
            fprintf(stderr, "Error: Failed to parse model file %s\n", model_file);
            return 1;
        }
        model_in_log_space = 1;
    } else {
        // Parse individual model files
        donor_parser(&l, don_emission);
        acceptor_parser(&l, acc_emission);
        exon_intron_parser(&l, exon_emission, 0);
        exon_intron_parser(&l, intron_emission, 1);

        explicit_duration_probability(&ed, Ped_exon, 0);
        explicit_duration_probability(&ed, Ped_intron, 1);
    }

    if (DEBUG) printf("\n--- Phase 3: Computing Transition Matrices ---\n");
    if (model_in_log_space) {
        compute_transition_matrices_log(&l);
    } else {
        compute_transition_matrices(&l);
    }
    
    l.log_values_len = (ed.max_len_exon > ed.max_len_intron) ? ed.max_len_exon : ed.max_len_intron;
    l.log_values     = calloc(l.log_values_len, sizeof(double));

    /* --------------- Validation Check --------------- */
    if (DEBUG) {
        printf("\n--- Phase 4: Validation ---\n");
        printf("Parameters loaded:\n");
        printf("  Sequence length: %d bp\n", info.T);
        printf("  Flank size: %d bp\n", info.flank);
        printf("  Exon length range: %d-%d bp\n", ed.min_len_exon, ed.max_len_exon);
        printf("  Intron length range: %d-%d bp\n", ed.min_len_intron, ed.max_len_intron);
        printf("  Analysis range: %d to %d\n", 
               info.flank+ed.min_len_exon, info.T-info.flank-ed.min_len_exon);
        
        print_transition_matrices_summary(&l);
        print_duration_summary(&ed);
    }

    /* --------------- Log Space Conversion --------------- */
    int don_size    = power(4, l.B.don_kmer_len);
    int acc_size    = power(4, l.B.acc_kmer_len);
    int exon_size   = power(4, l.B.exon_kmer_len);
    int intron_size = power(4, l.B.intron_kmer_len);

    if (model_in_log_space) {
        // JSON model already in log-space, skip conversion
        if (DEBUG) printf("\n--- Phase 5: Model already in log space ---\n");
    } else {
        if (DEBUG) printf("\n--- Phase 5: Converting to Log Space ---\n");

        // Check for numerical issues
        tolerance_checker(ed.exon, ed.exon_len, 1e-15);
        tolerance_checker(ed.intron, ed.intron_len, 1e-15);
        tolerance_checker(l.A.dons, don_size, 1e-15);
        tolerance_checker(l.A.accs, acc_size, 1e-15);

        // Convert to log space for numerical stability
        log_space_converter(ed.exon, ed.exon_len);
        log_space_converter(ed.intron, ed.intron_len);
        log_space_converter(l.A.dons, don_size);
        log_space_converter(l.A.accs, acc_size);
        log_space_converter(l.B.exon, exon_size);
        log_space_converter(l.B.intron, intron_size);
    }

    /* --------------- Initialize K-mer Cache --------------- */
    if (DEBUG) printf("\n--- Phase 5.5: Initializing K-mer Cache ---\n");

    info.kmer_cache = allocate_kmer_cache(
        info.T,
        l.B.exon_kmer_len,
        l.B.intron_kmer_len,
        l.B.don_kmer_len,
        l.B.acc_kmer_len
    );
    if (!info.kmer_cache) {
        fprintf(stderr, "Error: Failed to allocate k-mer cache\n");
        return 1;
    }
    init_kmer_cache(info.kmer_cache, info.numerical_sequence);

    /* --------------- Exe Forward Backward Algorithm --------------- */
    if (DEBUG) printf("\n--- Phase 6: Forward-Backward Algorithm ---\n");

    // Allocate memory for algorithms
    allocate_fw(&info, &fw, &ed);
    allocate_bw(&bw, &ed, &info);
    allocate_pos(&pos, &info);

    // Run forward algorithm
    basis_fw_algo(&l, &ed, &fw, &info);
    if (DEBUG) printf("  Forward basis complete\n");
    fw_algo(&l, &fw, &info, &ed);
    if (DEBUG) printf("  Forward algorithm complete\n");
    
    // Run backward algorithm
    basis_bw_algo(&l, &bw, &info, &ed);
    if (DEBUG) printf("  Backward basis complete\n");
    bw_algo(&l, &bw, &info, &ed);
    if (DEBUG) printf("  Backward algorithm complete\n");
    
    // Calculate posterior probabilities
    pos_prob(&bw, &fw, &info, &pos);
    if (DEBUG) printf("  Posterior probabilities calculated\n");
    
    // Parse splice sites from posterior probabilities
    parse_splice_sites(&pos, &info);
    
    if (DEBUG) {
        printf("\n=== Results Summary ===\n");
        printf("Found %d donor sites and %d acceptor sites\n", pos.dons, pos.accs);
    }

    /* --------------- For MCTS Isoform Generation --------------- */
    if (use_mcts) {
        if (DEBUG) printf("\n--- Phase 7: MCTS Isoform Generation ---\n");

        if (pos.dons == 0 || pos.accs == 0) {
            printf("Warning: No splice sites found. Cannot run MCTS.\n");
        } else {
            Locus *loc = create_locus(n_isoforms);

            if (DEBUG) {
                printf("Generating isoforms using MCTS:\n");
                printf("  Locus capacity: %d\n", n_isoforms);
                printf("  Exploration constant (C): %.2f\n", mcts_explore_c);
                printf("  Donor sites: %d\n", pos.dons);
                printf("  Acceptor sites: %d\n", pos.accs);
            }

            // Create MCTS tree
            MCTSTree *tree = create_mcts_tree(&pos, &ed, &info, mcts_explore_c);

            // Run MCTS - iterations scale with problem size
            int max_iterations = n_isoforms * 10;  // More iterations for better coverage
            generate_isoforms_mcts(tree, loc, max_iterations);

            if (DEBUG) printf("Unique isoforms found: %d\n", loc->n_isoforms);
            if (output_json) {
                print_locus_json(loc, &info, stdout);
            } else {
                print_locus(loc, &info);
            }

            // Clean up
            free_mcts_tree(tree);
            free_locus(loc);
        }
    }


    /* --------------- Memory Cleanup --------------- */
    if (DEBUG) printf("\n--- Phase 8: Cleanup ---\n");

    free_splice_sites(&pos);
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free_pos(&pos, &info);
    free(info.original_sequence);
    free(info.numerical_sequence);
    free_kmer_cache(info.kmer_cache);
    free_lambda(&l);
    free_explicit_duration(&ed);
    
    if (DEBUG) printf("\n=== RFHMM Analysis Complete ===\n");
    return 0;
}