#ifndef _CALC_STAB_H_
#define _CALC_STAB_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <cmath>
#include <iostream>

#include "simple_rng.hpp"

using namespace std;

struct IdxPair {
  double value;
  size_t index;

  bool operator<(const IdxPair& p) const {
    return value > p.value;
  }
};

struct Buffers {
  Buffers(size_t K) : B(K), C(K), D(K), E(K), G(K), d(K) {}
  vector<double> B, C, D, E, G;
  vector<IdxPair> d;
};

template <typename Array1, typename Array2>
void _calc_stability_matrix(Array1& dest, const Array2& X, 
			    size_t n, size_t K, double beta, 
			    Buffers& buffer) {

  const double eps = 1e-8;

  for(size_t i = 0; i < n; ++i) {
    vector<IdxPair>& d = buffer.d;
    vector<double>& B = buffer.B;
    vector<double>& C = buffer.C;
    vector<double>& D = buffer.D;
    vector<double>& E = buffer.E;
    vector<double>& G = buffer.G;

    // First map them over
    double norm = 0;

    for(size_t k = 0; k < K; ++k) {
      double dv = X[i*K + k];
      norm += dv*dv;
      d[k].value = max(dv, double(0));
      d[k].index = k;
    }

    if(norm == 0) {
      for(size_t k = 0; k < K; ++k) {
	dest[i*K + k] = 1.0 / K;
      }
      continue;
    }

    for(vector<IdxPair>::iterator it = d.begin(); it != d.end(); ++it) {
      it->value /= sqrt(norm);
    }
 
    sort(d.begin(), d.end()); // sorts in descending order by value
    
    if(d[0].value == 0 && d[1].value != 0) {
      for(size_t k = 0; k < K; ++k)
	dest[i*K + k] = 0;
    
      dest[i*K + d[0].index] = 1;
      continue;
    }    

    for(vector<IdxPair>::iterator it = d.begin(); it != d.end(); ++it) {
      it->value += eps;
    }


    B[0] = 1./(d[0].value + eps);
    D[0] = 1;
    G[0] = d[0].value;

    for(size_t k = 1; k < K; ++k)  {
      B[k] = B[k-1] + 1. / (d[k].value + eps);
      C[k] = k - d[k].value * B[k-1];
      D[k] = exp(-exp(beta)*C[k]);
      G[k] = D[k] / B[k];
    }

    E[K-1] = 0;

    for(size_t k = K-1; k > 0; --k)
      E[k-1] = E[k] + D[k] / (B[k-1]*(d[k].value * B[k-1] + 1.));

    for(size_t k = 0; k < K; ++k)
      dest[i*K + d[k].index] = min(double(1.0), 
				   max(double(0), 
				       double((G[k] - E[k]) / (d[k].value + eps))));
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

  Buffers b(K);
  _calc_stability_matrix(dest, X, n, K, beta, b);
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

void sorted_stability_matrix(double *dest, int *indexes, int *K_map, 
			     double *X, int *labels,
			     size_t n, size_t K, double beta, int Kmap_mode) {

  double *stab_matrix = new double[n*K];

  stability_matrix(stab_matrix, X, n, K, beta);

  sort_stability_matrix(dest, indexes, K_map, stab_matrix, 
			labels, n, K, Kmap_mode);

  delete[] stab_matrix;
}

double log_score(const vector<double>& src, size_t n, size_t K) {

  double total = 0;

  for(size_t i = 0; i < n; ++i)
    total += *max_element(src.begin() + i*K, src.begin() + (i+1)*K);

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
void calculateScores(double * scores, const Array& src, size_t n, size_t K, size_t seed, 
		     size_t n_baselines, double beta, bool use_permutations, bool by_dimension) {

  
  const size_t n_threads = omp_get_max_threads();

  vector<vector<double> > data_buffers(n_threads);
  vector<vector<double> > stab_buffers(n_threads);
  vector<Buffers> buffers(n_threads, Buffers(K));

  typedef vector<vector<double> >::iterator buf_iter;

  for(buf_iter it = data_buffers.begin(); it != data_buffers.end(); ++it)
    it->resize(n*K);

  for(buf_iter it = stab_buffers.begin(); it != stab_buffers.end(); ++it)
    it->resize(n*K);

  // Get the first one
  _calc_stability_matrix(stab_buffers[0], src, n, K, beta, buffers[0]);
  const double dist_score = log_score(stab_buffers[0], n, K);

  // Set up the rest of the baselines
  CheapRNG uniform_int(seed);

  vector<size_t> seeds(n_baselines);

  generate(seeds.begin(), seeds.end(), uniform_int);

#pragma omp parallel for shared(data_buffers, stab_buffers, seeds, buffers, scores, src, beta, n, K) 
  for(size_t i = 0; i < n_baselines; ++i) {
    
    size_t nt = omp_get_thread_num();

    fill_baseline_matrix(data_buffers[nt], src, n, K, seeds[i], use_permutations, by_dimension);
    _calc_stability_matrix(stab_buffers[nt], data_buffers[nt], n, K, beta, buffers[nt]);
    scores[i] = dist_score - log_score(stab_buffers[nt], n, K);
  }
}

template <typename Array1>
double score(const Array1& dist, size_t n, size_t K, size_t seed,
	     size_t n_baselines, double beta, bool use_permutations, bool by_dimension)
{
  double *scores = new double[n_baselines];

  calculateScores(scores, dist, n, K, seed, n_baselines, beta, use_permutations, by_dimension);

  double result =  accumulate(scores, scores + n_baselines, double(0)) / n_baselines;

  delete[] scores;
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// Now for the silhouette distances

static inline double calculateSilhouette(double *silhouettes,  double *silhouette_distances, 
				       int *labels, size_t n, size_t K) {

  double silhouette_total = 0;

#pragma omp parallel for reduction(+:silhouette_total)
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
			  double beta, double *xvec, double *yvec)
{
  /*printf("centroids\n");
  for(int i = 0; i<K; i++){
    printf("%f, %f\n", centroids[i*2+0], centroids[i*2+1]);
  }
  */
  if(K == 0) {
    fill(stab_image, stab_image + image_nx*image_ny, 0);
    return;
  }

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
    double buffer_x = 0.1 * (image_x_upper - image_x_lower);
    double buffer_y = 0.1 * (image_y_upper - image_y_lower);

    image_x_lower -= buffer_x;
    image_x_upper += buffer_x;
    image_y_lower -= buffer_y;
    image_y_upper += buffer_y;
  }

  double *X = new double[image_ny * image_nx * K];

  double vx      = (image_x_upper - image_x_lower) / (image_nx - 1);
  double start_x = image_x_lower;
    
  double vy      = (image_y_upper - image_y_lower) / (image_ny - 1);
  double start_y = image_y_lower;

  for(size_t yi = 0; yi < image_ny; ++yi) {
    for(size_t xi = 0; xi < image_nx; ++xi) {
      double x = start_x + vx * xi;
      double y = start_y + vy * yi;

      if(yi == 0)
	xvec[xi] = x;
      if(xi == 0)
	yvec[yi] = y;

      for(size_t k = 0; k < K; ++k)
	X[yi * (image_nx*K) + xi*K + k] 
	  = sqrt( sqr(x - centroids[2*k + 0]) + sqr(y - centroids[2*k + 1]));
    }
  }

  Buffers buffer(K);

  double *stab_matrix = new double[image_ny * image_nx * K];  

  _calc_stability_matrix(stab_matrix, X, image_nx*image_ny, K, beta, buffer);

  for(size_t i = 0; i < image_nx * image_ny; ++i)
    stab_image[i] = *max_element(stab_matrix + i*K, stab_matrix + (i+1)*K);

  delete[] stab_matrix;
  delete[] X;
}

#endif /* _CALC_STAB_H_ */
