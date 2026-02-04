#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mcts.h"

/* --------------- Node Creation/Destruction --------------- */

MCTSNode* create_mcts_node(MCTSNode *parent, int site_idx, int is_donor) {
    MCTSNode *node = calloc(1, sizeof(MCTSNode));

    node->site_idx = site_idx;
    node->is_donor = is_donor;
    node->visits = 0;
    node->wins = 0;
    node->parent = parent;
    node->children = NULL;
    node->n_children = 0;
    node->children_capacity = 0;

    // Initialize partial isoform state
    if (parent == NULL) {
        // Root node
        node->donors = NULL;
        node->acceptors = NULL;
        node->n_introns = 0;
        node->last_acceptor_pos = -1;
        node->expecting_acceptor = 0;  // Start by choosing a donor
    } else {
        // Copy parent state and extend
        int parent_introns = parent->n_introns;
        int new_introns = parent_introns;

        if (is_donor) {
            // Adding a donor - intron not complete yet
            node->donors = malloc((parent_introns + 1) * sizeof(int));
            node->acceptors = malloc((parent_introns + 1) * sizeof(int));

            if (parent->donors) {
                memcpy(node->donors, parent->donors, parent_introns * sizeof(int));
            }
            if (parent->acceptors) {
                memcpy(node->acceptors, parent->acceptors, parent_introns * sizeof(int));
            }
            node->donors[parent_introns] = site_idx;  // Will be converted to position later
            node->n_introns = parent_introns;  // Not complete until acceptor added
            node->last_acceptor_pos = parent->last_acceptor_pos;
            node->expecting_acceptor = 1;  // Must choose acceptor next
        } else {
            // Adding an acceptor - completes an intron
            new_introns = parent_introns + 1;
            node->donors = malloc(new_introns * sizeof(int));
            node->acceptors = malloc(new_introns * sizeof(int));

            if (parent->donors) {
                memcpy(node->donors, parent->donors, parent_introns * sizeof(int));
                // The pending donor from parent
                node->donors[parent_introns] = parent->donors[parent_introns];
            }
            if (parent->acceptors) {
                memcpy(node->acceptors, parent->acceptors, parent_introns * sizeof(int));
            }
            node->acceptors[parent_introns] = site_idx;  // Will be converted to position later
            node->n_introns = new_introns;
            node->last_acceptor_pos = site_idx;  // Will be converted to position later
            node->expecting_acceptor = 0;  // Can choose donor next (or terminate)
        }
    }

    return node;
}

void free_mcts_node(MCTSNode *node) {
    if (!node) return;

    // Recursively free children
    for (int i = 0; i < node->n_children; i++) {
        free_mcts_node(node->children[i]);
    }
    free(node->children);
    free(node->donors);
    free(node->acceptors);
    free(node);
}

/* --------------- Tree Creation/Destruction --------------- */

static int compare_int(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

MCTSTree* create_mcts_tree(Pos_prob *pos, Explicit_duration *ed,
                           Observed_events *info, double exploration_c) {
    MCTSTree *tree = calloc(1, sizeof(MCTSTree));

    int FLANK = (info->flank != 0) ? info->flank : 99;

    // Copy and sort donor positions
    tree->n_donors = pos->dons;
    tree->donor_positions = malloc(tree->n_donors * sizeof(int));
    tree->donor_values = malloc(tree->n_donors * sizeof(double));
    memcpy(tree->donor_positions, pos->dons_bps, tree->n_donors * sizeof(int));
    memcpy(tree->donor_values, pos->dons_val, tree->n_donors * sizeof(double));

    // Sort donors by position (keep values aligned)
    // Simple bubble sort for small arrays - could optimize
    for (int i = 0; i < tree->n_donors - 1; i++) {
        for (int j = 0; j < tree->n_donors - i - 1; j++) {
            if (tree->donor_positions[j] > tree->donor_positions[j+1]) {
                int tmp_pos = tree->donor_positions[j];
                tree->donor_positions[j] = tree->donor_positions[j+1];
                tree->donor_positions[j+1] = tmp_pos;

                double tmp_val = tree->donor_values[j];
                tree->donor_values[j] = tree->donor_values[j+1];
                tree->donor_values[j+1] = tmp_val;
            }
        }
    }

    // Copy and sort acceptor positions
    tree->n_acceptors = pos->accs;
    tree->acceptor_positions = malloc(tree->n_acceptors * sizeof(int));
    tree->acceptor_values = malloc(tree->n_acceptors * sizeof(double));
    memcpy(tree->acceptor_positions, pos->accs_bps, tree->n_acceptors * sizeof(int));
    memcpy(tree->acceptor_values, pos->accs_val, tree->n_acceptors * sizeof(double));

    // Sort acceptors
    for (int i = 0; i < tree->n_acceptors - 1; i++) {
        for (int j = 0; j < tree->n_acceptors - i - 1; j++) {
            if (tree->acceptor_positions[j] > tree->acceptor_positions[j+1]) {
                int tmp_pos = tree->acceptor_positions[j];
                tree->acceptor_positions[j] = tree->acceptor_positions[j+1];
                tree->acceptor_positions[j+1] = tmp_pos;

                double tmp_val = tree->acceptor_values[j];
                tree->acceptor_values[j] = tree->acceptor_values[j+1];
                tree->acceptor_values[j+1] = tmp_val;
            }
        }
    }

    // Set constraints
    tree->min_exon_len = ed->min_len_exon;
    tree->max_exon_len = ed->max_len_exon;
    tree->min_intron_len = ed->min_len_intron;
    tree->max_intron_len = ed->max_len_intron;
    tree->gene_start = FLANK;
    tree->gene_end = info->T - FLANK - 1;

    // UCB1 constant
    tree->exploration_c = exploration_c;

    // Create root node
    tree->root = create_mcts_node(NULL, -1, 0);

    // Statistics
    tree->total_iterations = 0;
    tree->valid_isoforms_found = 0;

    srand(time(NULL));

    return tree;
}

void free_mcts_tree(MCTSTree *tree) {
    if (!tree) return;

    free_mcts_node(tree->root);
    free(tree->donor_positions);
    free(tree->donor_values);
    free(tree->acceptor_positions);
    free(tree->acceptor_values);
    free(tree);
}

/* --------------- Constraint Checking --------------- */

int is_valid_donor_choice(MCTSTree *tree, MCTSNode *node, int donor_idx) {
    int donor_pos = tree->donor_positions[donor_idx];

    // Get the position to measure exon from
    int exon_start;
    if (node->n_introns == 0 && !node->expecting_acceptor) {
        // First donor - measure from gene start
        exon_start = tree->gene_start;
    } else if (node->last_acceptor_pos >= 0) {
        // Measure from last acceptor
        exon_start = tree->acceptor_positions[node->last_acceptor_pos];
    } else {
        exon_start = tree->gene_start;
    }

    int exon_len = donor_pos - exon_start;

    // Check exon length constraints
    if (exon_len < tree->min_exon_len) return 0;
    if (tree->max_exon_len > 0 && exon_len > tree->max_exon_len) return 0;

    // Check that donor is before gene end (leaving room for intron + final exon)
    if (donor_pos + tree->min_intron_len + tree->min_exon_len > tree->gene_end) return 0;

    return 1;
}

int is_valid_acceptor_choice(MCTSTree *tree, MCTSNode *node, int acceptor_idx, int pending_donor_idx) {
    int acceptor_pos = tree->acceptor_positions[acceptor_idx];
    int donor_pos = tree->donor_positions[pending_donor_idx];

    // Acceptor must come after donor
    if (acceptor_pos <= donor_pos) return 0;

    int intron_len = acceptor_pos - donor_pos;

    // Check intron length constraints
    if (intron_len < tree->min_intron_len) return 0;
    if (tree->max_intron_len > 0 && intron_len > tree->max_intron_len) return 0;

    // Check that there's room for final exon
    int final_exon_len = tree->gene_end - acceptor_pos;
    if (final_exon_len < tree->min_exon_len) return 0;

    return 1;
}

/* --------------- UCB1 Calculation --------------- */

double compute_ucb1(MCTSNode *node, int parent_visits, double c) {
    if (node->visits == 0) {
        return INFINITY;  // Unvisited nodes have highest priority
    }

    double exploitation = (double)node->wins / node->visits;
    double exploration = c * sqrt(log(parent_visits) / node->visits);

    return exploitation + exploration;
}

/* --------------- MCTS Selection --------------- */

MCTSNode* mcts_select(MCTSTree *tree, MCTSNode *node) {
    // Walk down the tree using UCB1 until we find an expandable node
    while (node->n_children > 0) {
        // Check if node is fully expanded
        int can_expand = 0;

        if (node->expecting_acceptor) {
            // Need to check if there are unexpanded acceptor choices
            int pending_donor = node->donors[node->n_introns];
            for (int i = 0; i < tree->n_acceptors; i++) {
                if (is_valid_acceptor_choice(tree, node, i, pending_donor)) {
                    // Check if this acceptor is already a child
                    int found = 0;
                    for (int j = 0; j < node->n_children; j++) {
                        if (!node->children[j]->is_donor && node->children[j]->site_idx == i) {
                            found = 1;
                            break;
                        }
                    }
                    if (!found) {
                        can_expand = 1;
                        break;
                    }
                }
            }
        } else {
            // Check for unexpanded donor choices
            for (int i = 0; i < tree->n_donors; i++) {
                if (is_valid_donor_choice(tree, node, i)) {
                    int found = 0;
                    for (int j = 0; j < node->n_children; j++) {
                        if (node->children[j]->is_donor && node->children[j]->site_idx == i) {
                            found = 1;
                            break;
                        }
                    }
                    if (!found) {
                        can_expand = 1;
                        break;
                    }
                }
            }
        }

        if (can_expand) {
            // This node can still be expanded
            return node;
        }

        // Node is fully expanded, select best child using UCB1
        MCTSNode *best_child = NULL;
        double best_ucb1 = -INFINITY;

        for (int i = 0; i < node->n_children; i++) {
            double ucb1 = compute_ucb1(node->children[i], node->visits, tree->exploration_c);
            if (ucb1 > best_ucb1) {
                best_ucb1 = ucb1;
                best_child = node->children[i];
            }
        }

        if (!best_child) {
            // No valid children - this is a terminal node
            return node;
        }

        node = best_child;
    }

    return node;
}

/* --------------- MCTS Expansion --------------- */

MCTSNode* mcts_expand(MCTSTree *tree, MCTSNode *node) {
    // Find a valid unexpanded action and create a new child

    if (node->expecting_acceptor) {
        // Must choose an acceptor
        int pending_donor = node->donors[node->n_introns];

        // Collect valid unexpanded acceptors
        int *valid_acceptors = malloc(tree->n_acceptors * sizeof(int));
        int n_valid = 0;

        for (int i = 0; i < tree->n_acceptors; i++) {
            if (is_valid_acceptor_choice(tree, node, i, pending_donor)) {
                // Check if already expanded
                int found = 0;
                for (int j = 0; j < node->n_children; j++) {
                    if (!node->children[j]->is_donor && node->children[j]->site_idx == i) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    valid_acceptors[n_valid++] = i;
                }
            }
        }

        if (n_valid == 0) {
            free(valid_acceptors);
            return node;  // No expansion possible
        }

        // Randomly select one (could weight by posterior value)
        int chosen_idx = valid_acceptors[rand() % n_valid];
        free(valid_acceptors);

        // Create new child node
        MCTSNode *child = create_mcts_node(node, chosen_idx, 0);

        // Store actual position in the node's acceptors array
        child->acceptors[child->n_introns - 1] = chosen_idx;
        child->last_acceptor_pos = chosen_idx;

        // Add to parent's children
        if (node->n_children >= node->children_capacity) {
            node->children_capacity = node->children_capacity == 0 ? 8 : node->children_capacity * 2;
            node->children = realloc(node->children, node->children_capacity * sizeof(MCTSNode*));
        }
        node->children[node->n_children++] = child;

        return child;

    } else {
        // Can choose a donor (or terminate if we have at least one intron)
        int *valid_donors = malloc(tree->n_donors * sizeof(int));
        int n_valid = 0;

        for (int i = 0; i < tree->n_donors; i++) {
            if (is_valid_donor_choice(tree, node, i)) {
                int found = 0;
                for (int j = 0; j < node->n_children; j++) {
                    if (node->children[j]->is_donor && node->children[j]->site_idx == i) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    valid_donors[n_valid++] = i;
                }
            }
        }

        if (n_valid == 0) {
            free(valid_donors);
            return node;  // No expansion possible - terminal node
        }

        // Randomly select one
        int chosen_idx = valid_donors[rand() % n_valid];
        free(valid_donors);

        // Create new child node
        MCTSNode *child = create_mcts_node(node, chosen_idx, 1);

        // Store actual position
        child->donors[child->n_introns] = chosen_idx;

        // Add to parent's children
        if (node->n_children >= node->children_capacity) {
            node->children_capacity = node->children_capacity == 0 ? 8 : node->children_capacity * 2;
            node->children = realloc(node->children, node->children_capacity * sizeof(MCTSNode*));
        }
        node->children[node->n_children++] = child;

        return child;
    }
}

/* --------------- MCTS Simulation (Random Rollout) --------------- */

int mcts_simulate(MCTSTree *tree, MCTSNode *node) {
    // Random rollout from current node to terminal state
    // Returns 1 if valid isoform, 0 otherwise

    // Copy current state
    int max_introns = tree->n_donors < tree->n_acceptors ? tree->n_donors : tree->n_acceptors;
    int *sim_donors = malloc(max_introns * sizeof(int));
    int *sim_acceptors = malloc(max_introns * sizeof(int));
    int sim_n_introns = node->n_introns;
    int sim_expecting_acceptor = node->expecting_acceptor;
    int sim_last_acceptor = node->last_acceptor_pos;
    int pending_donor = -1;

    // Copy existing introns
    for (int i = 0; i < sim_n_introns; i++) {
        sim_donors[i] = node->donors[i];
        sim_acceptors[i] = node->acceptors[i];
    }

    // If we're expecting an acceptor, there's a pending donor
    if (sim_expecting_acceptor && node->n_introns >= 0) {
        pending_donor = node->donors[node->n_introns];
    }

    // Random rollout
    int max_steps = 2 * max_introns;
    for (int step = 0; step < max_steps; step++) {
        if (sim_expecting_acceptor) {
            // Must choose an acceptor
            int *valid = malloc(tree->n_acceptors * sizeof(int));
            int n_valid = 0;

            for (int i = 0; i < tree->n_acceptors; i++) {
                int acc_pos = tree->acceptor_positions[i];
                int don_pos = tree->donor_positions[pending_donor];

                if (acc_pos <= don_pos) continue;

                int intron_len = acc_pos - don_pos;
                if (intron_len < tree->min_intron_len) continue;
                if (tree->max_intron_len > 0 && intron_len > tree->max_intron_len) continue;

                int final_exon = tree->gene_end - acc_pos;
                if (final_exon < tree->min_exon_len) continue;

                valid[n_valid++] = i;
            }

            if (n_valid == 0) {
                // Can't complete - invalid
                free(valid);
                free(sim_donors);
                free(sim_acceptors);
                return 0;
            }

            int chosen = valid[rand() % n_valid];
            free(valid);

            sim_acceptors[sim_n_introns] = chosen;
            sim_last_acceptor = chosen;
            sim_n_introns++;
            sim_expecting_acceptor = 0;

        } else {
            // Can choose a donor or terminate

            // Check if we can terminate (valid final exon)
            int can_terminate = 0;
            if (sim_n_introns > 0) {
                int last_acc_pos = tree->acceptor_positions[sim_last_acceptor];
                int final_exon = tree->gene_end - last_acc_pos;
                if (final_exon >= tree->min_exon_len &&
                    (tree->max_exon_len <= 0 || final_exon <= tree->max_exon_len)) {
                    can_terminate = 1;
                }
            }

            // Find valid donors
            int *valid = malloc(tree->n_donors * sizeof(int));
            int n_valid = 0;

            int exon_start = (sim_n_introns == 0) ? tree->gene_start
                                                   : tree->acceptor_positions[sim_last_acceptor];

            for (int i = 0; i < tree->n_donors; i++) {
                int don_pos = tree->donor_positions[i];

                if (don_pos <= exon_start) continue;

                int exon_len = don_pos - exon_start;
                if (exon_len < tree->min_exon_len) continue;
                if (tree->max_exon_len > 0 && exon_len > tree->max_exon_len) continue;

                if (don_pos + tree->min_intron_len + tree->min_exon_len > tree->gene_end) continue;

                valid[n_valid++] = i;
            }

            if (n_valid == 0 && !can_terminate) {
                // Dead end
                free(valid);
                free(sim_donors);
                free(sim_acceptors);
                return (sim_n_introns > 0) ? 1 : 0;  // Valid if we have at least one intron
            }

            if (n_valid == 0 || (can_terminate && rand() % 3 == 0)) {
                // Terminate
                free(valid);
                break;
            }

            // Choose a donor
            int chosen = valid[rand() % n_valid];
            free(valid);

            sim_donors[sim_n_introns] = chosen;
            pending_donor = chosen;
            sim_expecting_acceptor = 1;
        }
    }

    // Check if we ended in a valid state
    int valid = (sim_n_introns > 0 && !sim_expecting_acceptor);

    free(sim_donors);
    free(sim_acceptors);

    return valid ? 1 : 0;
}

/* --------------- MCTS Backpropagation --------------- */

void mcts_backpropagate(MCTSNode *node, int reward) {
    while (node != NULL) {
        node->visits++;
        node->wins += reward;
        node = node->parent;
    }
}

/* --------------- Main MCTS Iteration --------------- */

void mcts_iteration(MCTSTree *tree) {
    // 1. Selection
    MCTSNode *node = mcts_select(tree, tree->root);

    // 2. Expansion
    if (node->visits > 0 || node == tree->root) {
        MCTSNode *expanded = mcts_expand(tree, node);
        if (expanded != node) {
            node = expanded;
        }
    }

    // 3. Simulation
    int reward = mcts_simulate(tree, node);

    // 4. Backpropagation
    mcts_backpropagate(node, reward);

    tree->total_iterations++;
    if (reward) {
        tree->valid_isoforms_found++;
    }
}

/* --------------- Extract Isoform from Node --------------- */

Isoform* extract_isoform_from_node(MCTSNode *node, MCTSTree *tree) {
    if (node->n_introns == 0 || node->expecting_acceptor) {
        return NULL;  // Incomplete isoform
    }

    Isoform *iso = malloc(sizeof(Isoform));
    iso->beg = tree->gene_start;
    iso->end = tree->gene_end;
    iso->n_introns = node->n_introns;
    iso->val = 0.0;

    if (node->n_introns > 0) {
        iso->dons = malloc(node->n_introns * sizeof(int));
        iso->accs = malloc(node->n_introns * sizeof(int));

        for (int i = 0; i < node->n_introns; i++) {
            iso->dons[i] = tree->donor_positions[node->donors[i]];
            iso->accs[i] = tree->acceptor_positions[node->acceptors[i]];
            // Accumulate posterior values
            iso->val += tree->donor_values[node->donors[i]];
            iso->val += tree->acceptor_values[node->acceptors[i]];
        }
    } else {
        iso->dons = NULL;
        iso->accs = NULL;
    }

    return iso;
}

/* --------------- Collect All Terminal Isoforms --------------- */

static void collect_isoforms_recursive(MCTSNode *node, MCTSTree *tree, Locus *loc, int min_visits) {
    // If this is a terminal node (can represent a complete isoform)
    if (!node->expecting_acceptor && node->n_introns > 0 && node->visits >= min_visits) {
        if (loc->n_isoforms < loc->capacity) {
            Isoform *iso = extract_isoform_from_node(node, tree);
            if (iso) {
                loc->isoforms[loc->n_isoforms++] = iso;
            }
        }
    }

    // Recurse to children
    for (int i = 0; i < node->n_children; i++) {
        if (loc->n_isoforms >= loc->capacity) break;
        collect_isoforms_recursive(node->children[i], tree, loc, min_visits);
    }
}

/* --------------- Main Generation Function --------------- */

void generate_isoforms_mcts(MCTSTree *tree, Locus *loc, int max_iterations) {
    // Run MCTS iterations
    for (int i = 0; i < max_iterations; i++) {
        mcts_iteration(tree);

        if (DEBUG && i % 100 == 0) {
            printf("MCTS iteration %d: %d valid paths found\n", i, tree->valid_isoforms_found);
        }
    }

    // Collect isoforms from well-explored paths (visited at least 2 times)
    collect_isoforms_recursive(tree->root, tree, loc, 2);

    // Also add single-exon isoform (no introns)
    if (loc->n_isoforms < loc->capacity) {
        Isoform *single_exon = malloc(sizeof(Isoform));
        single_exon->beg = tree->gene_start;
        single_exon->end = tree->gene_end;
        single_exon->n_introns = 0;
        single_exon->dons = NULL;
        single_exon->accs = NULL;
        single_exon->val = 0.0;
        loc->isoforms[loc->n_isoforms++] = single_exon;
    }

    if (DEBUG) {
        print_mcts_stats(tree);
    }
}

/* --------------- Debug Statistics --------------- */

static void count_nodes_recursive(MCTSNode *node, int *total, int *terminal) {
    (*total)++;
    if (!node->expecting_acceptor && node->n_introns > 0) {
        (*terminal)++;
    }
    for (int i = 0; i < node->n_children; i++) {
        count_nodes_recursive(node->children[i], total, terminal);
    }
}

void print_mcts_stats(MCTSTree *tree) {
    int total_nodes = 0;
    int terminal_nodes = 0;
    count_nodes_recursive(tree->root, &total_nodes, &terminal_nodes);

    printf("\nMCTS Statistics:\n");
    printf("  Total iterations: %d\n", tree->total_iterations);
    printf("  Valid isoforms found: %d\n", tree->valid_isoforms_found);
    printf("  Tree nodes: %d\n", total_nodes);
    printf("  Terminal nodes (potential isoforms): %d\n", terminal_nodes);
    printf("  Root visits: %d\n", tree->root->visits);
    printf("  Root win rate: %.2f%%\n",
           tree->root->visits > 0 ? 100.0 * tree->root->wins / tree->root->visits : 0.0);
}
