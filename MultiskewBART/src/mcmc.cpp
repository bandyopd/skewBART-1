#include "mcmc.h"

using namespace arma;
using namespace Rcpp;

bool do_mh(double loglik_new, double loglik_old,
           double new_to_old, double old_to_new) {

  double cutoff = loglik_new + new_to_old - loglik_old - old_to_new;

  return log(unif_rand()) < cutoff ? true : false;

}

double growth_prior(int leaf_depth, const Hypers& hypers) {
  return hypers.gamma * pow(1.0 + leaf_depth, -hypers.beta);
}


Node* birth_node(Node* tree, double* leaf_node_probability) {
  std::vector<Node*> leafs = leaves(tree);
  Node* leaf = rand(leafs);
  *leaf_node_probability = 1.0 / ((double)leafs.size());

  return leaf;
}

double probability_node_birth(Node* tree) {
  return tree->is_leaf ? 1.0 : 0.5;
}


Node* death_node(Node* tree, double* p_not_grand) {
  std::vector<Node*> ngb = not_grand_branches(tree);
  Node* branch = rand(ngb);
  *p_not_grand = 1.0 / ((double)ngb.size());

  return branch;
}

void birth_death(Node* tree, const arma::mat& X, const arma::mat& Y,
                 const Hypers& hypers) {


  double p_birth = probability_node_birth(tree);

  if(unif_rand() < p_birth) {
    node_birth(tree, X, Y, hypers);
  }
  else {
    node_death(tree, X, Y, hypers);
 }
}


void node_birth(Node* tree, const arma::mat& X, const arma::mat& Y,
                const Hypers& hypers) {

  // Rcout << "Sample leaf";
  double leaf_probability = 0.0;
  Node* leaf = birth_node(tree, &leaf_probability);

  // Rcout << "Compute prior";
  int leaf_depth = leaf->depth;
  double leaf_prior = growth_prior(leaf_depth, hypers);

  // Get likelihood of current state
  // Rcout << "Current likelihood";
  double ll_before = LogLT(tree, Y, X, hypers);
  ll_before += log(1.0 - leaf_prior);

  // Get transition probability
  // Rcout << "Transistion";
  double p_forward = log(probability_node_birth(tree) * leaf_probability);

  // Birth new leaves
  // Rcout << "Birth";
  leaf->BirthLeaves(hypers);

  // Get likelihood after
  // Rcout << "New Likelihood";
  double ll_after = LogLT(tree, Y, X, hypers);
  ll_after += log(leaf_prior) +
    log(1.0 - growth_prior(leaf_depth + 1, hypers)) +
    log(1.0 - growth_prior(leaf_depth + 1, hypers));

  // Get Probability of reverse transition
  // Rcout << "Reverse";
  std::vector<Node*> ngb = not_grand_branches(tree);
  double p_not_grand = 1.0 / ((double)(ngb.size()));
  double p_backward = log((1.0 - probability_node_birth(tree)) * p_not_grand);

  // Do MH
  double log_trans_prob = ll_after + p_backward - ll_before - p_forward;
  if(log(unif_rand()) > log_trans_prob) {
    leaf->DeleteLeaves();
  }
  else {
    // Rcout << "Accept!";
  }
}

void node_death(Node* tree, const arma::mat& X, const arma::mat& Y,
                const Hypers& hypers) {

  // Select branch to kill Children
  double p_not_grand = 0.0;
  Node* branch = death_node(tree, &p_not_grand);

  // Compute before likelihood
  int leaf_depth = branch->left->depth;
  double leaf_prob = growth_prior(leaf_depth - 1, hypers);
  double left_prior = growth_prior(leaf_depth, hypers);
  double right_prior = growth_prior(leaf_depth, hypers);
  double ll_before = LogLT(tree, Y, X, hypers) +
    log(1.0 - left_prior) + log(1.0 - right_prior) + log(leaf_prob);

  // Compute forward transition prob
  double p_forward = log(p_not_grand * (1.0 - probability_node_birth(tree)));

  // Save old leafs, do not delete (they are dangling, need to be handled by the end)
  Node* left = branch->left;
  Node* right = branch->right;
  branch->left = branch;
  branch->right = branch;
  branch->is_leaf = true;

  // Compute likelihood after
  double ll_after = LogLT(tree, Y, X, hypers) + log(1.0 - leaf_prob);

  // Compute backwards transition
  std::vector<Node*> leafs = leaves(tree);
  double p_backwards = log(1.0 / ((double)(leafs.size())) * probability_node_birth(tree));

  // Do MH and fix dangles
  branch->left = left;
  branch->right = right;
  double log_trans_prob = ll_after + p_backwards - ll_before - p_forward;
  if(log(unif_rand()) > log_trans_prob) {
    branch->is_leaf = false;
  }
  else {
    branch->DeleteLeaves();
  }
}

Node* draw_prior(Node* tree, const arma::mat& X, const arma::mat& Y, const Hypers& hypers) {

  // Compute loglik before
  Node* tree_0 = tree;
  double loglik_before = LogLT(tree_0, Y, X, hypers);

  // Make new tree and compute loglik after
  Node* tree_1 = new Node(hypers);
  tree_1->GenBelow(hypers);
  double loglik_after = LogLT(tree_1, Y, X, hypers);

  // Do MH
  if(log(unif_rand()) < loglik_after - loglik_before) {
    delete tree_0;
    tree = tree_1;
  }
  else {
    delete tree_1;
  }
  return tree;
}

double calc_cutpoint_likelihood(Node* node) {
  if(node->is_leaf) return 1;

  double out = 1.0 / (node->uppers(node->var) - node->lowers(node->var));
  out = out * calc_cutpoint_likelihood(node->left);
  out = out * calc_cutpoint_likelihood(node->right);

  return out;
}

/*
  get_perturb_limits:

  Input: a branch node pointer
  Output: a 2-d vector out, such that (out[0], out[1]) represents the interval in
          which branch->val can vary without contradicting any of the branch
          nodes further down the tree. Note: THIS DOES NOT MODIFY THE TREE,
          IT ONLY COMPUTES THE LIMITS

  Algorithm: First, manually traverse backwards up the tree, checking to see
             if we encounter any nodes splitting on the same var. If we do, we
             update min and max. Next, we collect the branches below the
             current node. If any left ancestor splits on the current node,
             then we modify the min; otherwise, we modify the max.

*/

std::vector<double> get_perturb_limits(Node* branch) {
  double min = 0.0;
  double max = 1.0;

  Node* n = branch;
  while(!(n->is_root)) {
    if(n->is_left) {
      n = n->parent;
      if(n->var == branch->var) {
        if(n->val > min) {
          min = n->val;
        }
      }
    }
    else {
      n = n->parent;
      if(n->var == branch->var) {
        if(n->val < max) {
          max = n->val;
        }
      }
    }
  }
  std::vector<Node*> left_branches = branches(n->left);
  std::vector<Node*> right_branches = branches(n->right);
  for(int i = 0; i < left_branches.size(); i++) {
    if(left_branches[i]->var == branch->var) {
      if(left_branches[i]->val > min)
        min = left_branches[i]->val;
    }
  }
  for(int i = 0; i < right_branches.size(); i++) {
    if(right_branches[i]->var == branch->var) {
      if(right_branches[i]->val < max) {
        max = right_branches[i]->val;
      }
    }
  }

  std::vector<double> out; out.push_back(min); out.push_back(max);
  return out;
}


void perturb_decision_rule(Node* tree,
                           const arma::mat& X,
                           const arma::mat& Y,
                           const Hypers& hypers) {

  // Randomly choose a branch; if no branches, we automatically reject
  std::vector<Node*> bbranches = branches(tree);
  if(bbranches.size() == 0)
    return;

  // Select the branch
  Node* branch = rand(bbranches);

  // Calculuate tree likelihood before proposal
  double ll_before = LogLT(tree, Y, X, hypers);

  // Calculate product of all 1/(B - A) here
  double cutpoint_likelihood = calc_cutpoint_likelihood(tree);

  // Calculate backward transition density
  std::vector<double> lims = get_perturb_limits(branch);
  double backward_trans = 1.0/(lims[1] - lims[0]);

  // save old split
  int old_feature = branch->var;
  double old_value = branch->val;

  // Modify the branch
  branch->var = hypers.SampleVar();
  lims = get_perturb_limits(branch);
  branch->val = lims[0] + (lims[1] - lims[0]) * unif_rand();
  tree->get_limits_below();

  // Calculate likelihood after proposal
  double ll_after = LogLT(tree, Y, X, hypers);

  // Calculate product of all 1/(B-A)
  double cutpoint_likelihood_after = calc_cutpoint_likelihood(tree);

  // Calculate forward transition density
  double forward_trans = 1.0/(lims[1] - lims[0]);

  // Do MH
  double log_trans_prob =
    ll_after + log(cutpoint_likelihood_after) + log(backward_trans)
    - ll_before - log(cutpoint_likelihood) - log(forward_trans);

  if(log(unif_rand()) > log_trans_prob) {
    branch->var = old_feature;
    branch->val = old_value;
    tree->get_limits_below();
  }
}


void TreeBackfit(Nodevec& forest, arma::mat & Y_hat,
                 const Hypers& hypers, const arma::mat& X, const arma::mat& Y,
                 const Opts& opts) {

  double MH_BD = 0.7;
  double MH_PRIOR = 0.4;

  int num_tree = hypers.num_tree;
  for(int t = 0; t < num_tree; t++) {
    // Rcout << "Getting backfit quantities";
    arma::mat Y_star = Y_hat - predict(forest[t], X, hypers);
    arma::mat res = Y - Y_star;

    if(unif_rand() < MH_PRIOR) {
      // Rcout << "Draw Prior";
      forest[t] = draw_prior(forest[t], X, res, hypers);
      // Rcout << "Done";
    }
    if(forest[t]->is_leaf || unif_rand() < MH_BD) {
      // Rcout << "BD step";
      birth_death(forest[t], X, res, hypers);
      // Rcout << "Done";
    }
    else {
      // Rcout << "Change step";
      perturb_decision_rule(forest[t], X, res, hypers);
      // Rcout << "Done";
    }
    // Rcout << "Update mu";
    forest[t]->UpdateMu(res, X, hypers);
    // Rcout << "Refit";
    Y_hat = Y_star + predict(forest[t], X, hypers);
    // Rcout << "Done";
  }
}




