#ifndef NODE_H
#define NODE_H

#include <RcppArmadillo.h>
#include "hypers.h"

struct Node;
typedef std::vector<Node*> Nodevec;

class SuffStats {
 public:

  arma::vec sum_Y;
  arma::mat sum_YYt;
  double n;
  int J;

 SuffStats(int JJ) : J(JJ) {
    n = 0.0;
    sum_Y = arma::zeros<arma::vec>(J);
    sum_YYt = arma::zeros<arma::mat>(J);
  }

};

class Node {

  public:

  bool is_leaf;
  bool is_root;
  bool is_left;
  Node* left;
  Node* right;
  Node* parent;

  // Branch params
  int depth;
  int var;
  double val;
  arma::sp_vec lowers;
  arma::sp_vec uppers;
  arma::sp_uvec parent_split_indicators;

  // Sufficient Statistics
  SuffStats* ss;
  void ResetSuffStat() {
    ss->n = 0.0;
    ss->sum_Y = arma::zeros<arma::vec>(mu.size());
    ss->sum_YYt = arma::zeros<arma::mat>(mu.size(), mu.size());
    if(!is_leaf) {
      left->ResetSuffStat();
      right->ResetSuffStat();
    }
  }


  // Leaf parameters
  arma::vec mu;

  // Constructor / Destructor
  Node(const Hypers& hypers);
  Node(Node* parent, const Hypers& hypers, bool is_left);
  ~Node();

  // Updates and such
  void BirthLeaves(const Hypers& hypers);
  void GenBelow(const Hypers& hypers);
  void get_limits_below();
  void DeleteLeaves();
  void UpdateMu(const arma::mat& Y, const arma::mat& X, const Hypers& hypers);

};


// Functions for computing with the forest
void GetSuffStats(Node* n, const arma::mat& y,
                  const arma::mat& X, const Hypers& hypers);
void GetSuffStats(Node* n, const arma::vec& y, const arma::mat& yyt,
                    const arma::vec& x, const Hypers& hypers);

double LogLT(Node* n, const arma::mat& Y,
             const arma::mat& X, const Hypers& hypers);


// Split probability, used in growing nodes
double SplitProb(Node* node, const Hypers& hypers);

// Other node-specific functions
void leaves(Node* x, Nodevec& leafs);
std::vector<Node*> leaves(Node* x);

#endif
