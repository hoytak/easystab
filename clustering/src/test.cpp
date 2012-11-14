#include <R.h>
/*#include <Rinternals.h>*/
#include <stdio.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

int* test_function(int* a, int *b, int size, bool s);

extern "C" SEXP print_matrix(SEXP s_matrix) {
  int size = INTEGER(GET_DIM(s_matrix))[0] * INTEGER(GET_DIM(s_matrix))[1];
  printf("number of points: %d\n", size);
  for(int i = 0; i<size; i++)
    printf("%d ", INTEGER(s_matrix)[i]);
  printf("\nfinished output\n");

  SEXP retval;
  PROTECT(retval = NEW_NUMERIC(1));
  DOUBLE_DATA(retval)[0] = (size) ;
  UNPROTECT(1);

  return retval;
}

extern "C" SEXP pretty_matrix(SEXP s_matrix) {
  int *matrix, height, width;
  SEXP retval;
  if (isMatrix(s_matrix) && isInteger(s_matrix)) {
    matrix = INTEGER(s_matrix);
    width  = INTEGER(GET_DIM(s_matrix))[1];
    height = INTEGER(GET_DIM(s_matrix))[0];
  }
  else {
    printf("invalid matrix.\n");
    return R_NilValue;
  }
  PROTECT(retval = NEW_NUMERIC(1));
  DOUBLE_DATA(retval)[0] = (4 * height + 1) ;
  UNPROTECT(1);

  return retval;
}

extern "C" SEXP test_function_wrapper(SEXP a_, SEXP b_, SEXP size_, SEXP s_){
  int * a = INTEGER(a_);
  int * b = INTEGER(b_);
  int * size = INTEGER(size_);
  int * s = LOGICAL(s_);
  int* res = test_function(a, b, *size, *s==0);
  SEXP retval;
  PROTECT(retval = NEW_INTEGER(*size));
  for(int i=0; i<*size; i++)
    INTEGER_DATA(retval)[i] = res[i];
  UNPROTECT(1);
  return retval;
}

int* test_function(int* a, int *b, int size, bool s){
  int *sum = new int[size];
  for (int i = 0; i<size; i++){
    sum[i] = a[i]+b[i];
  }
  if (s)
    for (int i=0; i<size; i++)
      sum[i] = 99;
  return sum;
}

//
//R_CallMethodDef callMethods[] =
//  {
//    {"pretty_matrix", (DL_FUNC)&pretty_matrix, 1},
//    {NULL,NULL, 0}
//  };
//
//void R_init_pretty_matrix(DllInfo *dll)
//{
//  R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
//}

