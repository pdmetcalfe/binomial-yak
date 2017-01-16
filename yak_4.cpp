#include <algorithm>
#include <queue>
#include <complex>
#include <exception>
#include <vector>

#include <Rcpp.h>

#include <fftw3.h>

namespace {
  inline int nextGoodSize(int size) {
    std::priority_queue<int, std::vector<int>, std::greater<int> > hamming;
    hamming.push(1);
    while (!hamming.empty()) {
      int current = hamming.top();
      if (current > size) return current;
      // take off the head
      hamming.pop();
      
      // take off duplicates
      while ((!hamming.empty()) && (current == hamming.top())) {
	hamming.pop();
      }
      
      // and now push on the next multiples
      hamming.push(2 * current);
      hamming.push(3 * current);
      hamming.push(5 * current);
      hamming.push(7 * current);
    }

    return size;
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
Rcpp::NumericVector gsProbs_4(const double p,
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
  std::vector<double> probs(nextGoodSize(1 + sizes[0]));
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

    std::vector<double> newProb(nextGoodSize(1 + sizes[i]));
    std::vector< std::complex<double> > probOut(newProb.size() / 2 + 1);
    std::vector<double> stepVector(newProb.size());
    std::vector< std::complex<double> > stepOut(probOut.size());

    fftw_plan probForward =
      fftw_plan_dft_r2c_1d(newProb.size(),
			   &newProb[0],
			   reinterpret_cast<fftw_complex*>(&probOut[0]),
			   FFTW_ESTIMATE);
    fftw_plan stepForward =
      fftw_plan_dft_r2c_1d(stepVector.size(),
			   &stepVector[0],
			   reinterpret_cast<fftw_complex*>(&stepOut[0]),
			   FFTW_ESTIMATE
			   );
    fftw_plan goBack =
      fftw_plan_dft_c2r_1d(newProb.size(),
			   reinterpret_cast<fftw_complex*>(&probOut[0]),
			   &newProb[0],
			   FFTW_ESTIMATE);
    
    std::copy(probs.begin() + 1 + crits[i - 1],
	      probs.begin() + 1 + sizes[i - 1],
	      newProb.begin() + 1 + crits[i - 1]);

    for (int k=0; k<=(sizes[i] - sizes[i-1]);++k) {
      stepVector[k] = R::dbinom(k, sizes[i] - sizes[i-1], p, 0);
    }

    fftw_execute(probForward);
    fftw_execute(stepForward);
    
    for (int k=0;k<stepOut.size();++k) {
      probOut[k] *= stepOut[k] / double(stepVector.size());
    }
    fftw_execute(goBack);

    fftw_destroy_plan(goBack);
    fftw_destroy_plan(probForward);
    fftw_destroy_plan(stepForward);
    
    result[i] = std::accumulate(newProb.begin(),
				newProb.begin() + 1 + crits[i],
				double(0));
    std::swap(newProb, probs);
  }
  
  return result;
}
