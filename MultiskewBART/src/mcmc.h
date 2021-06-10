/*

  This file contains functions which are used to carry out the MCMC. This
  includes functions for computing the growth probability, probability of
  selecting nodes to split, and code for carrying out the MCMC.

 */

#ifndef MCMC_H
#define MCMC_H

#include "tree_utils.h"
#include "opts.h"

// Metropolis hastings
bool do_mh(double loglik_new, double loglik_old,
           double new_to_old, double old_to_new);

// Various tree probabilities
double growth_prior(int leaf_depth, const Hypers& hypers);



// Sampling and selecting nodes
Node* birth_node(Node* tree, double* leaf_node_probability);
double probability_node_birth(Node* tree);
Node* death_node(Node* tree, double* p_not_grand);


// MCMC on the trees
void birth_death(Node* tree, const arma::mat& X, const arma::mat& Y,
                 const Hypers& hypers);
void node_birth(Node* tree, const arma::mat& X, const arma::mat& Y,
                const Hypers& hypers);
void node_death(Node* tree, const arma::mat& X, const arma::mat& Y,
                const Hypers& hypers);
Node* draw_prior(Node* tree, const arma::mat& X, const arma::mat& Y, const Hypers& hypers);


// Functions for the perturb algorithm
double calc_cutpoint_likelihood(Node* node);
std::vector<double> get_perturb_limits(Node* branch);
void perturb_decision_rule(Node* tree,
                           const arma::mat& X,
                           const arma::mat& Y,
                           const Hypers& hypers);

// The Bayesian backfitting algorithm 
void TreeBackfit(Nodevec& forest, arma::mat& Y_hat,
                 const Hypers& hypers, const arma::mat& X, const arma::mat& Y,
                 const Opts& opts);

#endif
