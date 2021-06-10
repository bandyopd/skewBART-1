#include "node.h"

using namespace arma;
using namespace Rcpp;

Node::Node(const Hypers& hypers) {

  is_leaf = true;
  is_root = true;
  is_left = true;
  left = this;
  right = this;
  parent = this;

  depth = 0;
  var = 0;
  val = 0.0;
  lowers = zeros<sp_mat>(hypers.group.size(), 1);
  uppers = zeros<sp_mat>(hypers.group.size(), 1);
  parent_split_indicators = zeros<sp_umat>(hypers.group.size(), 1);

  mu = zeros<vec>(hypers.Sigma.n_rows);
  ss = new SuffStats(hypers.Sigma.n_rows);
}

Node::Node(Node* parent, const Hypers& hypers, bool is_left) {
  parent->is_leaf = false;

  is_leaf = true;
  is_root = false;
  left = this;
  right = this;
  this->parent = parent;
  this->is_left = is_left;

  depth = parent->depth + 1;
  var = 0;
  val = 0.0;
  lowers = parent->lowers;
  uppers = parent->uppers;
  parent_split_indicators = parent->parent_split_indicators;
  mu = zeros<vec>(hypers.Sigma.n_rows);
  ss = new SuffStats(hypers.Sigma.n_rows);

  // Update the bounds
  int pvar = parent->var;
  if(is_left) {
    uppers(pvar)= parent->val;
  }
  else {
    lowers(pvar) = parent->val;
  }


}

Node::~Node() {
  if(!is_leaf) {
    delete left;
    delete right;
  }
  delete ss;
}

void Node::BirthLeaves(const Hypers& hypers) {
  if(is_leaf) {
    // Rcout << "Get Vars";
    var = hypers.SampleVar();
    // Rcout << "OK Vars";
    if(parent_split_indicators(var) == 0) {
      parent_split_indicators(var) = 1;
      uppers(var) = 1.0;
      lowers(var) = 0.0;
    }
    // Rcout << "Sampling val";
    val = (uppers(var) - lowers(var)) * unif_rand() + lowers(var);
    // Rcout << "Make leftright";
    left = new Node(this, hypers, true);
    right = new Node(this, hypers, false);
  }
}

void Node::GenBelow(const Hypers& hypers) {
  double grow_prob = SplitProb(this, hypers);
  double u = unif_rand();
  if(u < grow_prob) {
    // Rcout << "BL";
    BirthLeaves(hypers);
    // Rcout << "Grow left";
    left->GenBelow(hypers);
    right->GenBelow(hypers);
  }
}

void Node::get_limits_below() {

  if(is_root) {
    lowers = zeros<sp_mat>(lowers.size(), 1);
    uppers = zeros<sp_mat>(uppers.size(), 1);
    parent_split_indicators = zeros<sp_umat>(lowers.size(), 1);
  }
  lowers = parent->lowers;
  uppers = parent->uppers;
  parent_split_indicators = parent->parent_split_indicators;
  if(!is_root) {
    if(is_left) {
      uppers(parent->var) = parent->val;
    }
    else {
      lowers(parent->var) = parent->val;
    }
  }
  if(!is_leaf) {
    if(parent_split_indicators(var) == 0) {
      parent_split_indicators(var) = 1;
      uppers(var) = 1.0;
      lowers(var) = 0.0;
    }
    left->get_limits_below();
    right->get_limits_below();
  }
}


void Node::DeleteLeaves() {
  delete left;
  delete right;
  left = this;
  right = this;
  is_leaf = true;
  if(is_root || parent->parent_split_indicators(var) == 0) {
    parent_split_indicators(var) = 0;
    uppers(var) = 0.0;
    lowers(var) = 0.0;
  }
  var = 0;
  val = 0.0;
}


void Node::UpdateMu(const arma::mat& Y, const arma::mat& X,
                     const Hypers& hypers) {

  GetSuffStats(this, Y, X, hypers);
  Nodevec leafs = leaves(this);
  int num_leaves = leafs.size();
  for(int l = 0; l < num_leaves; l++) {
    mat n_A_plus_B = leafs[l]->ss->n * hypers.A + hypers.B;
    vec Y_tilde = solve(n_A_plus_B, hypers.A * leafs[l]->ss->sum_Y);
    leafs[l]->mu = rmvnorm(Y_tilde, n_A_plus_B);
  }

}

void GetSuffStats(Node* n, const arma::mat& Y,
                  const arma::mat& X,
                  const Hypers& hypers) {

  n->ResetSuffStat();
  int N = Y.n_rows;
  for(int i = 0; i < N; i++) {
    vec y = trans(Y.row(i));
    mat yyt = y * y.t();
    vec x = trans(X.row(i));
    GetSuffStats(n, y, yyt, x, hypers);
  }
}

void GetSuffStats(Node* n, const arma::vec& y, const arma::mat& yyt,
                  const arma::vec& x, const Hypers& hypers) {

  n->ss->n += 1.0;
  n->ss->sum_Y = n->ss->sum_Y + y;
  n->ss->sum_YYt = n->ss->sum_YYt + yyt;
  if(!(n->is_leaf)) {
    if(x(n->var) < n->val) {
      GetSuffStats(n->left, y, yyt, x, hypers);
    }
    else {
      GetSuffStats(n->right, y, yyt, x, hypers);
    }
  }
}

// RECALL: A = inv(Sigma) and B = inv(Sigma_mu)
double LogLT(Node* n, const arma::mat& Y,
             const arma::mat& X, const Hypers& hypers) {

  GetSuffStats(n, Y, X, hypers);

  // Rcout << "Leaves ";
  std::vector<Node*> leafs = leaves(n);
  int num_leaves = leafs.size();

  double out = 0.0;

  double log_det_B, sign;
  log_det(log_det_B, sign, hypers.B);

  for(int l = 0; l < num_leaves; l++) {
    mat n_A_plus_B = leafs[l]->ss->n * hypers.A + hypers.B;
    vec Y_tilde = solve(n_A_plus_B, hypers.A * leafs[l]->ss->sum_Y);
    out -= 0.5 * trace(hypers.A * leafs[l]->ss->sum_YYt);
    out += 0.5 * as_scalar(Y_tilde.t() * n_A_plus_B * Y_tilde);
    out += 0.5 * log_det_B;

    double log_det_naplusb, sign2; log_det(log_det_naplusb, sign2, n_A_plus_B);
    out -= 0.5 * log_det_naplusb;
  }

  return out;
}

double SplitProb(Node* node, const Hypers& hypers) {
  return hypers.gamma * pow(1.0 + node->depth, -hypers.beta);
}

void leaves(Node* x, Nodevec& leafs) {
  if(x->is_leaf) {
    leafs.push_back(x);
  }
  else {
    leaves(x->left, leafs);
    leaves(x->right, leafs);
  }
}

std::vector<Node*> leaves(Node* x) {
  std::vector<Node*> leafs(0);
  leaves(x, leafs);
  return leafs;
}
