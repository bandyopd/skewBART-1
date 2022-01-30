#include "tree_utils.h"

using namespace arma;
using namespace Rcpp;

std::vector<Node*> not_grand_branches(Node* tree) {
  std::vector<Node*> ngb(0);
  not_grand_branches(ngb, tree);
  return ngb;
}

void not_grand_branches(Nodevec& ngb, Node* node) {
  if(!node->is_leaf) {
    bool left_is_leaf = node->left->is_leaf;
    bool right_is_leaf = node->right->is_leaf;
    if(left_is_leaf && right_is_leaf) {
      ngb.push_back(node);
    }
    else {
      not_grand_branches(ngb, node->left);
      not_grand_branches(ngb, node->right);
    }
  }
}
void branches(Node* n, Nodevec& branch_vec) {
  if(!(n->is_leaf)) {
    branch_vec.push_back(n);
    branches(n->left, branch_vec);
    branches(n->right, branch_vec);
  }
}

std::vector<Node*> branches(Node* root) {
  std::vector<Node*> branch_vec;
  branch_vec.resize(0);
  branches(root, branch_vec);
  return branch_vec;
}

Node* rand(Nodevec& ngb) {

  int N = ngb.size();
  arma::vec p = ones<vec>(N) / ((double)(N));
  int i = sample_class(p);
  return ngb[i];
}


arma::uvec get_var_counts(Nodevec& forest, const Hypers& hypers) {
  arma::uvec counts = zeros<uvec>(hypers.group.n_cols);
  int num_tree = forest.size();
  for(int t = 0; t < num_tree; t++) {
    get_var_counts(counts, forest[t], hypers);
  }
  return counts;
}

void get_var_counts(arma::uvec& counts, Node* node, const Hypers& hypers) {
  if(!node->is_leaf) {
    counts(node->var) = counts(node->var) + 1;
    get_var_counts(counts, node->left, hypers);
    get_var_counts(counts, node->right, hypers);
  }
}

// arma::uvec get_var_counts(Nodevec& forest, const Hypers& hypers) {
//   arma::uvec counts = zeros<uvec>(hypers.s.size());
//   int num_tree = forest.size();
//   for(int t = 0; t < num_tree; t++) {
//     get_var_counts(counts, forest[t], hypers);
//   }
//   return counts;
// }

// void get_var_counts(arma::uvec& counts, Node* node, const Hypers& hypers) {
//   if(!node->is_leaf) {
//     int group_idx = hypers.group(node->var);
//     counts(group_idx) = counts(group_idx) + 1;
//     get_var_counts(counts, node->left, hypers);
//     get_var_counts(counts, node->right, hypers);
//   }
// }

// arma::vec get_means(Nodevec& forest) {
//   std::vector<double> means(0);
//   int num_tree = forest.size();
//   for(int t = 0; t < num_tree; t++) {
//     get_means(forest[t], means);
//   }

//   // Convert std::vector to armadillo vector, deep copy
//   vec out(&(means[0]), means.size());
//   return out;
// }
//
// void get_means(Node* node, std::vector<double>& means) {
//
//   if(node->is_leaf) {
//     means.push_back(node->mu);
//   }
//   else {
//     get_means(node->left, means);
//     get_means(node->right, means);
//   }
// }


arma::mat predict(const Nodevec& forest, const arma::mat& X,
                  const Hypers& hypers) {

  mat out = zeros<mat>(X.n_rows, hypers.Sigma.n_rows);
  int num_tree = forest.size();

  for(int t = 0; t < num_tree; t++) {
    out = out + predict(forest[t], X, hypers);
  }
  return out;
}

// LEAVE OFF HERE
arma::mat predict(Node* n, const arma::mat& X, const Hypers& hypers) {
  mat out = zeros<mat>(X.n_rows, hypers.Sigma.n_rows);
  int N = X.n_rows;
  for(int i = 0; i < N; i++) {
    vec x = trans(X.row(i));
    out.row(i) = trans(predict(n, x, hypers));
  }
  return out;
}

arma::vec predict(Node* n, const arma::vec& x, const Hypers& hypers) {
  if(n->is_leaf) return n->mu;

  bool go_left = x(n->var) < n->val;
  if(go_left) {
    return predict(n->left, x, hypers);
  }
  else {
    return predict(n->right, x, hypers);
  }
}

