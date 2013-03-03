#include "calc_stab.h"
#include "dist_utils.h"

#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <stdio.h>

extern "C" SEXP _make_stability_image(SEXP stab_image_, SEXP image_nx_, SEXP image_ny_,
				      SEXP image_x_lower_, SEXP image_x_upper_,
				      SEXP image_y_lower_, SEXP image_y_upper_,
				      SEXP centroids_, SEXP K_, SEXP beta_, 
				      SEXP xvec_, SEXP yvec_, SEXP buffer_){

  double * stab_image = REAL(stab_image_);
  size_t image_nx = INTEGER(image_nx_)[0];
  size_t image_ny = INTEGER(image_ny_)[0];
  double image_x_lower = REAL(image_x_lower_)[0];
  double image_x_upper = REAL(image_x_upper_)[0];
  double image_y_lower = REAL(image_y_lower_)[0];
  double image_y_upper = REAL(image_y_upper_)[0];
  double * centroids = REAL(centroids_);
  size_t K = INTEGER(K_)[0];
  double beta = REAL(beta_)[0];
  double * xvec = REAL(xvec_);
  double * yvec = REAL(yvec_);
  double buffer = REAL(buffer_)[0];

  make_stability_image(stab_image, image_nx, image_ny, 
		       image_x_lower, image_x_upper, image_y_lower, image_y_upper,
		       centroids, K, beta, xvec, yvec, buffer);

  SEXP retval;
  PROTECT(retval = NEW_INTEGER(1));
  INTEGER(retval)[0] = 0;
  UNPROTECT(1);

  return retval;
}

extern "C" SEXP _sort_stability_matrix(SEXP dest_, SEXP index_map_, SEXP K_map_, SEXP src_, 
				       SEXP labels_, SEXP n_, SEXP K_, SEXP Kmap_mode_){
  double * dest = REAL(dest_);
  int * index_map = INTEGER(index_map_);
  int * K_map = INTEGER(K_map_);
  double * src = REAL(src_);
  int * labels = INTEGER(labels_);
  int n = INTEGER(n_)[0];
  int K = INTEGER(K_)[0];
  int Kmap_mode = INTEGER(Kmap_mode_)[0];

  // R array start with index 1, c++ starts with 0
  for(int i=0; i<n; i++)
    labels[i] --;

  sort_stability_matrix(dest, index_map, K_map, src, labels, n, K, Kmap_mode);

  for(int i=0; i<n; i++){
    labels[i] ++;
    index_map[i] ++;
  }

  for(int k=0; k<K; k++){
    K_map[k]++;
  }

  SEXP retval;
  PROTECT(retval = NEW_INTEGER(1));
  INTEGER(retval)[0] = 0;
  UNPROTECT(1);

  return retval;
}

extern "C" SEXP _score(SEXP dist_, SEXP labels_, SEXP n_, SEXP K_, SEXP seed_,
		       SEXP n_baselines_, SEXP beta_, SEXP use_permutations_, SEXP by_dimension_){
  double * dist = REAL(dist_);
  int * labels = INTEGER(labels_);
  int n = INTEGER(n_)[0];
  int K = INTEGER(K_)[0];
  int seed = INTEGER(seed_)[0];
  int n_baselines = INTEGER(n_baselines_)[0];
  double beta = REAL(beta_)[0];
  bool use_permutations = LOGICAL(use_permutations_)[0] != 0;
  bool by_dimension = LOGICAL(by_dimension_)[0] != 0;
  
  SEXP retval;

  for(int i=0; i<n; i++)
    labels[i] --;

  double res = score(dist, labels, n, K, seed, n_baselines, beta, 
		     use_permutations, by_dimension);

  for(int i=0; i<n; i++)
    labels[i] ++;

  PROTECT(retval = NEW_NUMERIC(1));
  REAL(retval)[0] = res;
  UNPROTECT(1);

  return retval;
}

extern "C" SEXP _calculateScores(SEXP scores_, SEXP confusion_, SEXP stability_matrix_,
				 SEXP src_, SEXP labels_, 
				 SEXP n_, SEXP K_, SEXP seed_, SEXP n_baselines_,
				 SEXP beta_, SEXP use_permutations_, SEXP by_dimension_){
  double* scores = REAL(scores_);
  double* confusion = REAL(confusion_);
  double* stability_matrix = REAL(stability_matrix_);
  double* src = REAL(src_);
  int * labels = INTEGER(labels_);
  size_t n = INTEGER(n_)[0];
  size_t K = INTEGER(K_)[0];
  size_t seed = INTEGER(seed_)[0];
  size_t n_baselines = INTEGER(n_baselines_)[0];
  double beta = REAL(beta_)[0];
  bool use_permutations = LOGICAL(use_permutations_)[0] != 0;
  bool by_dimension = LOGICAL(by_dimension_)[0] != 0;

  for(int i=0; i<n; i++)
    labels[i] --;

  calculateScores(scores, confusion, stability_matrix, src, labels, 
		  n, K, seed, n_baselines, 
		  beta, use_permutations, by_dimension);

  for(int i=0; i<n; i++){
    labels[i] ++;
  }

  SEXP retval;
  PROTECT(retval = NEW_INTEGER(1));
  INTEGER(retval)[0] = 0;
  UNPROTECT(1);

  return retval;

}

extern "C" SEXP _calculateAverageLinkageDistances(SEXP _dists, SEXP _labels, SEXP _n, SEXP _K, SEXP _src_dists){
  double *dists = REAL(_dists);
  size_t n = INTEGER(_n)[0];

  int * labels = INTEGER(_labels);

  size_t K = INTEGER(_K)[0];

  double * src_dists = REAL(_src_dists);

  for(int i=0; i<n; i++)
    labels[i] --;

  calculateAverageLinkageDistances(dists, labels, n, K, src_dists);

  for(int i=0; i<n; i++){
    labels[i] ++;
  }

  SEXP retval;
  PROTECT(retval = NEW_INTEGER(1));
  INTEGER(retval)[0] = 0;
  UNPROTECT(1);

  return retval;

}

extern "C" SEXP _calculateRepresentativeDistances(SEXP _rep_distances, SEXP _labels, 
						  SEXP _n, SEXP _K, SEXP _src_dists) {
  double * rep_distances = REAL(_rep_distances);
  size_t n = INTEGER(_n)[0];
  size_t K = INTEGER(_K)[0];
  double *src_dists = REAL(_src_dists);

  int * labels = new int[n];

  for (int i = 0; i<n; i++)
    labels[i] = INTEGER(_labels)[i]-1;

  calculateRepresentativeDistances(rep_distances, labels, n, K, src_dists);

  delete [] labels;

  SEXP retval;
  PROTECT(retval = NEW_INTEGER(1));
  INTEGER(retval)[0] = 0;
  UNPROTECT(1);

  return retval;
}

