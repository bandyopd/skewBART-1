/*

  This file contains functions for manipulating basic tree structures,
  extracting leaves/branches, counting variables, and so forth. Also, we have
  code here for extracting quantities such as the means.

 */

#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include "node.h"


std::vector<Node*> not_grand_branches(Node* tree);
void not_grand_branches(std::vector<Node*>& ngb, Node* node);
void branches(Node* n, Nodevec& branch_vec);
std::vector<Node*> branches(Node* root);
Node* rand(Nodevec& ngb);


// Getting the counts
arma::uvec get_var_counts(std::vector<Node*>& forest, const Hypers& hypers);
void get_var_counts(arma::uvec& counts, Node* node, const Hypers& hypers);


// Functions for collecting quantities from the forest
/* arma::mat get_means(Nodevec& forest); */
/* void get_means(Node* node, std::vector<double>& means); */


// Predictions
arma::mat predict(const Nodevec& forest, const arma::mat& X,
                  const Hypers& hypers);
arma::mat predict(Node* n, const arma::mat& X, const Hypers& hypers);
arma::vec predict(Node* n, const arma::vec& x, const Hypers& hypers);



#endif
