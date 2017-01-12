#include <algorithm>
#include <exception>
#include <vector>

#include <Rcpp.h>

namespace {
  inline void daxpy(const int num,
		    const double a,
		    const double* x,
		    double* y) {
    if (a==0) return;
    for (int i=0; i < num; ++i) {
      y[i] += a * x[i];
    }
  }
}

//' Simulate simple multiple testing strategy for single binomial event
//'
//' @param p event probability
//' @param sizes monotonically increasing vector of sample sizes
//' @param crits vector below which sampling stops
//' @return a vector giving probability of stopping at each sample
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gsProbs_2(const double p,
			      const Rcpp::IntegerVector& sizes,
			      const Rcpp::IntegerVector& crits) {
  if (sizes.size() != crits.size()) {
    throw std::runtime_error("Size mismatch");
  }

  if (0 == sizes.size()) {
    throw std::runtime_error("Huh?");
  }
  
  if (p < 0 || p > 1) {
    throw std::runtime_error("Your probability should be between 0 and 1!");
  }

  if ((sizes[0] < 0) || (crits[0] < 0)) {
    throw std::runtime_error("Sizes and critical values must be non-negative");
  }
  
  // fill the probabilities
  std::vector<double> probs(1 + sizes[0]);
  for (int i=0; i<=sizes[0]; ++i) {
    probs[i] = R::dbinom(i, sizes[0], p, 0);
  }
  
  // here's where the result is allocated
  Rcpp::NumericVector result(sizes.size());

  // and after the first look, which is easy
  result[0] = std::accumulate(probs.begin(),
			      probs.begin() + 1 + crits[0],
			      double(0));

  /* 
     Right, now lets go look at the other looks.
  */
  
  for (int i=1; i < sizes.size(); ++i) {
    if (sizes[i] < sizes[i - 1]) {
      throw std::runtime_error("Your sizes should be non-decreasing!");
    }
    if (crits[i] < 0) {
      throw std::runtime_error("Critical values must be non-negative");
    }

    /* 
       The indexing here is:
         (a) correct
         (b) hairy
    */

    // the probability of being in each state at this interim
    std::vector<double> newProb(1 + sizes[i]);

    /*
      How many bernoulli moves might have been made between this and
      the last interim?

      The 1+ includes the possibility of no move.
    */
    
    const int num_steps = 1 + sizes[i] - sizes[i - 1];
    for (int j=0; j < num_steps; ++j) {
      const double step_prob = R::dbinom(j,
					 sizes[i] - sizes[i - 1],
					 p,
					 0);
      
      const int toDo = std::min(sizes[i-1] - crits[i - 1],
				sizes[i] - crits[i-1] - j);
      daxpy(toDo, step_prob,
	    &probs[1 + crits[i - 1]],
	    &newProb[1 + crits[i - 1] + j]);
    }

    result[i] = std::accumulate(newProb.begin(),
				newProb.begin() + 1 + crits[i],
				double(0));
    std::swap(newProb, probs);
  }
  
  return result;
}
