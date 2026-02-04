#ifndef MCTS_H
#define MCTS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "model.h"
#include "randomf.h"

/* --------------- MCTS Node Structure --------------- */

typedef struct MCTSNode {
    // Splice site decision at this node
    int site_idx;               // Index into sorted splice sites (-1 for root)
    int is_donor;               // 1 if this node chose a donor, 0 if acceptor

    // MCTS statistics
    int visits;                 // n: number of times this node was visited
    int wins;                   // w: number of valid isoforms through this path

    // Tree structure
    struct MCTSNode *parent;
    struct MCTSNode **children;
    int n_children;
    int children_capacity;

    // Partial isoform state at this node
    int *donors;                // Donors selected so far
    int *acceptors;             // Acceptors selected so far
    int n_introns;              // Number of complete introns
    int last_acceptor_pos;      // Position of last acceptor (for constraint checking)
    int expecting_acceptor;     // 1 if next choice must be acceptor, 0 if donor
} MCTSNode;

typedef struct {
    MCTSNode *root;

    // Sorted splice sites for efficient lookup
    int *donor_positions;       // Sorted donor positions
    double *donor_values;       // Corresponding posterior values
    int n_donors;

    int *acceptor_positions;    // Sorted acceptor positions
    double *acceptor_values;    // Corresponding posterior values
    int n_acceptors;

    // Constraints
    int min_exon_len;
    int max_exon_len;
    int min_intron_len;
    int max_intron_len;
    int gene_start;
    int gene_end;

    // UCB1 exploration constant
    double exploration_c;

    // Statistics
    int total_iterations;
    int valid_isoforms_found;
} MCTSTree;

/* --------------- Function Declarations --------------- */

/* Tree Creation/Destruction */
MCTSTree* create_mcts_tree(Pos_prob *pos, Explicit_duration *ed,
                           Observed_events *info, double exploration_c);
void free_mcts_tree(MCTSTree *tree);

/* Node Operations */
MCTSNode* create_mcts_node(MCTSNode *parent, int site_idx, int is_donor);
void free_mcts_node(MCTSNode *node);

/* Core MCTS Algorithm */
void mcts_iteration(MCTSTree *tree);
MCTSNode* mcts_select(MCTSTree *tree, MCTSNode *node);
MCTSNode* mcts_expand(MCTSTree *tree, MCTSNode *node);
int mcts_simulate(MCTSTree *tree, MCTSNode *node);
void mcts_backpropagate(MCTSNode *node, int reward);

/* UCB1 Calculation */
double compute_ucb1(MCTSNode *node, int parent_visits, double c);

/* Isoform Generation */
void generate_isoforms_mcts(MCTSTree *tree, Locus *loc, int max_iterations);
Isoform* extract_isoform_from_node(MCTSNode *node, MCTSTree *tree);

/* Constraint Checking */
int is_valid_donor_choice(MCTSTree *tree, MCTSNode *node, int donor_idx);
int is_valid_acceptor_choice(MCTSTree *tree, MCTSNode *node, int acceptor_idx, int pending_donor_idx);

/* Debug/Stats */
void print_mcts_stats(MCTSTree *tree);

#endif
