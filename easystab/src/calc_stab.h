#ifndef _CALC_STAB_H_
#define _CALC_STAB_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

#include "simple_rng.h"

using namespace std;

struct IdxPair {
  double value;
  size_t index;

  bool operator<(const IdxPair& p) const {
    return value < p.value;
  }
};

// This is bad coding, but it's the quick way to do this so that
// things are maximally optimized.  
template <typename Array1, typename Array2, int K>
void _calc_stability_matrix_fixed(Array1& dest, const Array2& X, 
				  size_t n, double beta) {

  if(K == 1) {
    for(size_t i = 0; i < n; ++i)
      dest[i] = 1;
    return;
  }

  const double beta_0d = max(0.0, exp(beta) - 1);

  IdxPair d[K];
  double B[K], C[K], D[K];

  for(size_t i = 0; i < n; ++i) {

    // First map them over
    double min_value = d[0].value;

    for(size_t k = 0; k < K; ++k) {
      double dv = X[i*K + k];
      d[k].value = max(dv, double(0));
      d[k].index = k;

      min_value = min(min_value, d[k].value);
    }

    if(min_value == 0) {
      size_t zero_count = 0;
      for(size_t k = 0; k < K; ++k)
	zero_count += ((d[k].value == 0) ? 1 : 0);

      for(size_t k = 0; k < K; ++k)
	dest[i*K + k] = (d[k].value == 0) ? (1.0 / zero_count) : 0;

      continue;
    }
    
    for(size_t k = 0; k < K; ++k)
      d[k].value /= (min_value + 1e-32);
 
    sort(d, d + K); // sorts in descending order by value

    B[0] = 1./d[0].value;
    C[0] = 1.0;

    for(size_t k = 1; k < K; ++k) {
      B[k] = B[k-1] + 1.0 / d[k].value;
      C[k] = exp(beta_0d * (k - d[k].value * B[k-1]));
    }

    D[K-1] = 0;
    
    for(size_t k = K-1; k >= 1; --k)
      D[k - 1] = D[k] + C[k] / (B[k-1]*(B[k-1]*d[k].value + 1) );

    for(size_t k = 0; k < K; ++k)
      dest[i*K + d[k].index] = (C[k] / B[k] - D[k]) / d[k].value;
  }
}

template <typename Array1, typename Array2>
void _calc_stability_matrix_var(Array1& dest, const Array2& X, 
				size_t n, size_t K, double beta) {

  const double beta_0d = max(0.0, exp(beta) - 1);

  IdxPair *d = new IdxPair[K];
  double *B = new double[3*K];
  double *C = B + K;
  double *D = B + 2*K;

  for(size_t i = 0; i < n; ++i) {

    // First map them over
    double min_value = d[0].value;

    for(size_t k = 0; k < K; ++k) {
      double dv = X[i*K + k];
      d[k].value = max(dv, double(0));
      d[k].index = k;

      min_value = min(min_value, d[k].value);
    }

    if(min_value == 0) {
      size_t zero_count = 0;
      for(size_t k = 0; k < K; ++k)
	zero_count += ((d[k].value == 0) ? 1 : 0);

      for(size_t k = 0; k < K; ++k)
	dest[i*K + k] = (d[k].value == 0) ? (1.0 / zero_count) : 0;

      continue;
    }
    
    for(size_t k = 0; k < K; ++k)
      d[k].value /= (min_value + 1e-32);
 
    sort(d, d + K); // sorts in descending order by value

    B[0] = 1./d[0].value;
    C[0] = 1.0;

    for(size_t k = 1; k < K; ++k) {
      B[k] = B[k-1] + 1.0 / d[k].value;
      C[k] = exp(beta_0d * (k - d[k].value * B[k-1]));
    }

    D[K-1] = 0;
    
    for(size_t k = K-1; k >= 1; --k)
      D[k - 1] = D[k] + C[k] / (B[k-1]*(B[k-1]*d[k].value + 1) );

    for(size_t k = 0; k < K; ++k)
      dest[i*K + d[k].index] = (C[k] / B[k] - D[k]) / d[k].value;
  }
  
  delete[] d;
  delete[] B;
}

template <typename Array1, typename Array2>
void _calc_stability_matrix(Array1& dest, const Array2& X, 
			    size_t n, size_t K, double beta) {

  switch(K) {
  case 1: _calc_stability_matrix_fixed<Array1, Array2, 1>(dest, X, n, beta); return;
  case 2: _calc_stability_matrix_fixed<Array1, Array2, 2>(dest, X, n, beta); return;
  case 3: _calc_stability_matrix_fixed<Array1, Array2, 3>(dest, X, n, beta); return;
  case 4: _calc_stability_matrix_fixed<Array1, Array2, 4>(dest, X, n, beta); return;
  case 5: _calc_stability_matrix_fixed<Array1, Array2, 5>(dest, X, n, beta); return;
  case 6: _calc_stability_matrix_fixed<Array1, Array2, 6>(dest, X, n, beta); return;
  case 7: _calc_stability_matrix_fixed<Array1, Array2, 7>(dest, X, n, beta); return;
  case 8: _calc_stability_matrix_fixed<Array1, Array2, 8>(dest, X, n, beta); return;
  case 9: _calc_stability_matrix_fixed<Array1, Array2, 9>(dest, X, n, beta); return;
  case 10: _calc_stability_matrix_fixed<Array1, Array2, 10>(dest, X, n, beta); return;
  case 11: _calc_stability_matrix_fixed<Array1, Array2, 11>(dest, X, n, beta); return;
  case 12: _calc_stability_matrix_fixed<Array1, Array2, 12>(dest, X, n, beta); return;
  case 13: _calc_stability_matrix_fixed<Array1, Array2, 13>(dest, X, n, beta); return;
  case 14: _calc_stability_matrix_fixed<Array1, Array2, 14>(dest, X, n, beta); return;
  case 15: _calc_stability_matrix_fixed<Array1, Array2, 15>(dest, X, n, beta); return;
  case 16: _calc_stability_matrix_fixed<Array1, Array2, 16>(dest, X, n, beta); return;
  default: _calc_stability_matrix_var<Array1, Array2>(dest, X, n, K, beta); return;
  }
}

struct SortElement {
  size_t K;
  double value;
  size_t index;

  bool operator<(const SortElement& e2) const {
    return ((K == e2.K) 
	    ? (value > e2.value) 
	    : (K < e2.K) );
  }
};

class KmapSorter_1 {
public: 
  KmapSorter_1(const vector<size_t>& _counts)
    : counts(_counts)
  {}

  bool operator()(size_t i, size_t j) const {
    return counts[i] > counts[j];
  }

private:
  const vector<size_t>& counts;

};

class KmapSorter_2 {
public: 
  KmapSorter_2(const vector<size_t>& _counts, const vector<size_t>& _avg_index)
    : counts(_counts), avg_index(_avg_index)
  {}

  bool operator()(size_t i, size_t j) const{
    return avg_index[i]*counts[j] < avg_index[j]*counts[i];
  }

private:
  const vector<size_t>& counts, avg_index;
};

void stability_matrix(double *dest, double *X, size_t n, size_t K, double beta) {

  _calc_stability_matrix(dest, X, n, K, beta);
}

void sort_stability_matrix(double *dest, int *indexes, int *K_map, 
			   double *src_stab_matrix, int *labels,
			   size_t n, size_t K, int Kmap_mode)
{
  vector<size_t> Km(K);
  {
    vector<size_t> avg_index(K, 0);
    vector<size_t> counts(K, 0);
    vector<size_t> K_mapping(K);

    for(size_t i = 0; i < n; ++i) {
      ++counts[labels[i]];
      avg_index[labels[i]] += i;
    }

    for(size_t k = 0; k < K; ++k) 
      K_mapping[k] = k;

    switch(Kmap_mode) {
    case 1:
      sort(K_mapping.begin(), K_mapping.end(), KmapSorter_1(counts));
      break;
    case 2:
      sort(K_mapping.begin(), K_mapping.end(), KmapSorter_2(counts,avg_index));
    default: break;
    }

    for(size_t k = 0; k < K; ++k)
      Km[K_mapping[k]] = k;

    for(size_t k = 0; k < K; ++k)
      K_map[k] = int(Km[k]);
  }

  // Get the best version for the mapping
  vector<SortElement> sorting(n);

  for(size_t i = 0; i < n; ++i) {
    sorting[i].K = Km[labels[i]];
    sorting[i].value = src_stab_matrix[i*K + labels[i]];
    sorting[i].index = i;
  }

  sort(sorting.begin(), sorting.end());
  
  for(size_t i = 0; i < n; ++i) {
    indexes[i] = sorting[i].index;
    for(size_t k = 0; k < K; ++k) 
      dest[i*K + Km[k]] = src_stab_matrix[sorting[i].index * K + k];
  }
}

double log_score(const vector<double>& src, size_t n, size_t K) {

  double total = 0;

  for(size_t i = 0; i < n; ++i)
    total += *max_element(src.begin() + i*K, src.begin() + (i+1)*K);

  return total / n;
}

template <typename Array1, typename Array2>
double log_score(const Array1& src, const Array2& labels, size_t n, size_t K) {

  double total = 0;

  for(size_t i = 0; i < n; ++i)
    total += src[i*K + labels[i]];

  return total / n;
}

template <typename Array>
static inline void fill_baseline_matrix(Array& dest, const double *src,
					size_t n, size_t K, size_t seed, 
					bool use_permutations, 
					bool by_dimension) {

  switch(2*int(use_permutations) + int(by_dimension)) {
    
  case 2*0 + 0:
    {
      CheapRNG uniform_int(seed);
  
      for(size_t i = 0; i < n*K; ++i) {
	size_t idx = uniform_int(n*K - 1);
	dest[i] = src[idx];
      }
    
      return;
    }
  case 2*0 + 1:
    {
      CheapRNG uniform_int(seed);
  
      for(size_t i = 0; i < n; ++i) {
	for(size_t k = 0; k < K; ++k) {

	  size_t idx = uniform_int(n - 1);

	  dest[i*K + k] = src[idx*K + k];
	}
      }
    
    return;
    }
  case 2*1 + 0:
    {
      Shuffler sh(seed, n*K);

      for(size_t i = 0; i < n*K; ++i)
	dest[i] = src[sh[i] ];

      return;
    }
  case 2*1 + 1:
    {
      for(size_t k = 0; k < K; ++k) {
	Shuffler sh(seed, n);

	for(size_t i = 0; i < n; ++i) {
	  size_t idx = sh[i];
	  dest[i*K + k] = src[idx*K + k];
	}
      }

      return;
    }
  }
}

template <typename Array>
void calculateScores(double * scores, double * confusion_matrix, double * _stability_matrix, 
		     const Array& src, int* labels, size_t n, size_t K, size_t seed, 
		     size_t n_baselines, double beta, bool use_permutations, bool by_dimension) {

  vector<double> data_buffer(n*K);

  double *stability_matrix = (_stability_matrix == NULL) ? new double[n*K] : _stability_matrix;

  // Get the first one
  _calc_stability_matrix(stability_matrix, src, n, K, beta);
  const double dist_score = log_score(stability_matrix, labels, n, K);

  if(confusion_matrix != NULL) {
    // Now set up the confusion matrix
    vector<size_t> cl_counts(K, 0);
    fill(confusion_matrix, confusion_matrix + K*K, 0);

    for(size_t i = 0; i < n; ++i) {
      size_t k = size_t(labels[i]); 
      ++cl_counts[k];

      for(size_t ki = 0; ki < K; ++ki)
	confusion_matrix[k*K + ki] += stability_matrix[i*K + ki];
    }

    for(size_t k = 0; k < K; ++k) {
      for(size_t ki = 0; ki < K; ++ki) {
	confusion_matrix[k*K + ki] /= (1e-32 + cl_counts[k]);
      }
    }
  }

  // Set up the rest of the baselines
  CheapRNG uniform_int(seed);

  vector<size_t> seeds(n_baselines);
  vector<double> stab_buffer(n*K);

  generate(seeds.begin(), seeds.end(), uniform_int);

  for(size_t i = 0; i < n_baselines; ++i) {
    fill_baseline_matrix(data_buffer, src, n, K, seeds[i], use_permutations, by_dimension);
    _calc_stability_matrix(stab_buffer, data_buffer, n, K, beta);
    scores[i] = dist_score - log_score(stab_buffer, n, K);
  }

  if(_stability_matrix == NULL) 
    delete[] stability_matrix;
}

template <typename Array1>
double score(const Array1& dist, int* labels, size_t n, size_t K, size_t seed,
	     size_t n_baselines, double beta, bool use_permutations, bool by_dimension)
{
  double *scores = new double[n_baselines];

  calculateScores(scores, (double*)NULL, (double*)NULL, dist,labels, 
		  n, K, seed, n_baselines, beta, use_permutations, by_dimension);

  double result =  accumulate(scores, scores + n_baselines, double(0)) / n_baselines;

  delete[] scores;
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Now for the silhouette distances

static inline double calculateSilhouette(double *silhouettes,  double *silhouette_distances, 
					 int *labels, size_t n, size_t K) {

  double silhouette_total = 0;

  for(size_t i = 0; i < n; ++i) {
    
    size_t second_min_index = (labels[i] == 0) ? 1 : 0;
    double second_best_value = silhouette_distances[i*K + second_min_index];
    
    for(size_t j = second_min_index + 1; j < K; ++j) {
      if(j == size_t(labels[i]))
	continue;

      if(silhouette_distances[i*K + j] < second_best_value) {
	second_best_value = silhouette_distances[i*K + j];
	second_min_index = j;
      }
    }

    double best_value = silhouette_distances[i*K + labels[i] ];

    double s = (second_best_value - best_value) / max(best_value, second_best_value);

    if(silhouettes != NULL) 
      silhouettes[i] = s;

    silhouette_total += s;
  }

  return silhouette_total / n;
}

static inline double sqr(double x) {return x*x;}

void make_stability_image(double *stab_image, size_t image_nx, size_t image_ny, 
			  double image_x_lower, double image_x_upper, 
			  double image_y_lower, double image_y_upper, 
			  double *centroids,  size_t K, // centroids is a K x 2 matrix
			  double beta, double *xvec, double *yvec, double edge_buffer)
{
  if(K == 0) {
    fill(stab_image, stab_image + image_nx*image_ny, 0);
    return;
  } else if(K == 1) {
    fill(stab_image, stab_image + image_nx*image_ny, 1);
    return;
  }

  if(image_x_upper < image_x_lower)
    swap(image_x_upper, image_x_lower);

  if(image_y_upper < image_y_lower)
    swap(image_y_upper, image_y_lower);

  if (image_x_upper == 0 && image_x_lower == 0
      && image_y_upper == 0 && image_y_lower == 0) {
    
    image_x_lower = centroids[0 + 0];
    image_x_upper = centroids[0 + 0];
    image_y_lower = centroids[0 + 1];
    image_y_upper = centroids[0 + 1];

    for(size_t k = 1; k < K; ++k) {
      image_x_lower = min(image_x_lower, centroids[2*k + 0]);
      image_x_upper = max(image_x_upper, centroids[2*k + 0]);
      image_y_lower = min(image_y_lower, centroids[2*k + 1]);
      image_y_upper = max(image_y_upper, centroids[2*k + 1]);
    }

    // Add in a buffer region
    double buffer_x = edge_buffer * (image_x_upper - image_x_lower);
    double buffer_y = edge_buffer * (image_y_upper - image_y_lower);

    image_x_lower -= buffer_x;
    image_x_upper += buffer_x;
    image_y_lower -= buffer_y;
    image_y_upper += buffer_y;
  }


  // do the calculation at the center of every pixel, not a corner.
  double vx      = (image_x_upper - image_x_lower) / (image_nx);
  double start_x = image_x_lower + vx / 2;
    
  double vy      = (image_y_upper - image_y_lower) / (image_ny);
  double start_y = image_y_lower + vy / 2;

  double *X = new double[K];
  double *s = new double[K];

  for(size_t yi = 0; yi < image_ny; ++yi)
    yvec[yi] = start_y + vy * yi;

  for(size_t xi = 0; xi < image_nx; ++xi) 
    xvec[xi] = start_x + vx * xi;

  for(size_t yi = 0; yi < image_ny; ++yi) {
    for(size_t xi = 0; xi < image_nx; ++xi) {
      
      for(size_t k = 0; k < K; ++k)
	X[k] = sqrt( sqr(xvec[xi] - centroids[2*k + 0]) + sqr(yvec[yi] - centroids[2*k + 1]));

      _calc_stability_matrix(s, X, 1, K, beta);

      stab_image[yi*image_nx + xi] = *max_element(s, s + K);
    }
  }

  // Now, go through and set the pixels of the centroids to 1.

  for(size_t k = 0; k < K; ++k) {
    long xi = round( (centroids[2*k + 0] - start_x) / vx);
    long yi = round( (centroids[2*k + 1] - start_y) / vy);

    if(xi >= 0 && xi < long(image_nx) 
       && yi >= 0 && yi < long(image_ny) ) {

      stab_image[yi*image_nx + xi] = 1;
    }
  }

  delete[] s;
  delete[] X;
}

#endif /* _CALC_STAB_H_ */
