#ifndef _CLUST_UTILS_H_
#define _CLUST_UTILS_H_

#include <algorithm>
#include <vector>

using namespace std;

static inline double D(const double *dist, size_t n, size_t _i, size_t _j) {

  size_t i = min(_i, _j);
  size_t j = max(_i, _j);

  //return (i == j) ? 0 : dist[n*(i-1) - i*(i-1)/2 + j-i];
  return (i == j) ? 0 : dist[n*i - (i+1)*i/2 + j - i - 1];
}

void calculateAverageLinkageDistances(double *dists, int *labels,
				      size_t n, size_t K, double *src_dists) {
  
  // First, calculate the nearest clusters
  vector<size_t> cluster_counts(K, 0);
  for(size_t i = 0; i < n; ++i)
    ++cluster_counts[labels[i]];

  fill(dists, dists + n*K, 0);

  for(size_t i = 0; i < n; ++i) {

    for(size_t j = 0; j < n; ++j) {
      if(i == j) continue;
      dists[i*K + labels[j]] += D(src_dists, n, i, j);
    }
    
    for(size_t k = 0; k < K; ++k) 
      dists[i*K + k] 
	/= max(size_t(1), cluster_counts[k] - (k == size_t(labels[i]) ? 1 : 0) );
  }
}

void calculateRepresentativeDistances(double *rep_distances, int *labels, 
				      size_t n, size_t K, double *src_dists) {

  // First, calculate the nearest clusters
  vector<vector<size_t> > clusters(K);

  for(size_t i = 0; i < n; ++i)
    clusters[labels[i]].push_back(i);

  vector<size_t> best_so_far_index(K, 0);
  vector<double> best_so_far_value(K, 1e16);

  for(size_t i = 0; i < n; ++i) {

    const size_t k = labels[i];
    const vector<size_t>& idx_vect = clusters[k];

    double d_total = 0;
    for(size_t j = 0; j < idx_vect.size(); ++j) {
      d_total += D(src_dists, n, i, idx_vect[j]);
    }

    if(d_total < best_so_far_value[k]) {
      
      {
	if(d_total < best_so_far_value[k]) {
	  
	  best_so_far_value[k] = d_total;
	  best_so_far_index[k] = i;
	}
      }
    }
  }

  for(size_t i = 0; i < n; ++i) {
     for(size_t k = 0; k < K; ++k)
       rep_distances[i*K + k] = 
	 D(src_dists, n, i, best_so_far_index[k]);
  }
}
  


#endif /* _CLUST_UTILS_H_ */
