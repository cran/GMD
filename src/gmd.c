#define USING_R 1

/*

  Compute GM-Distance between discrete distributions, using method from:

  Systematic Clustering of Transcription Start Site Landscapes
  Author(s): Zhao et al, 2011
  Source: PLoS ONE 6(8): e23409. doi:10.1371/journal.pone.0023409 
  
  URL: http://dx.plos.org/10.1371/journal.pone.0023409
       http://www.plosone.org/article/info:doi/10.1371/journal.pone.0023409
  

  (c) Xiaobei Zhao 2011 <xiaobei@binf.ku.dk>

*/


#ifdef USING_R
# include <R.h>
# include <Rdefines.h>
# include <Rinternals.h>
# include <Rmath.h>
#endif
#include "gmd.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// ------------------------------------------------------------------------
// SEXP conversion, for R and C integration
// ------------------------------------------------------------------------
SEXP mdpa(SEXP v1, SEXP v2){
  int nProtected = 0;

  PROTECT(v1 = AS_NUMERIC(v1)); //a vector of real numbers
  nProtected++;

  PROTECT(v2 = AS_NUMERIC(v2));
  nProtected++;

  int s1 = LENGTH(v1);
  int s2 = LENGTH(v2);
	
  double* p_v1;
  double* p_v2;
  p_v1 = NUMERIC_POINTER(v1);
  p_v2 = NUMERIC_POINTER(v2);

  ARRAY a1={p_v1,s1};
  ARRAY a2={p_v2,s2};

  //
  SEXP ret;
  PROTECT(ret = NEW_NUMERIC(1));
  nProtected++;

  NUMERIC_POINTER(ret)[0] = __mdpa(a1, a2);
  
  UNPROTECT(nProtected);
  return(ret);
}



SEXP gmd0(SEXP v1, SEXP v2, SEXP pseudocount){
  int nProtected = 0;

  PROTECT(v1 = AS_NUMERIC(v1)); //a vector of real numbers
  nProtected++;

  PROTECT(v2 = AS_NUMERIC(v2));
  nProtected++;

  int i_pseudocount;
  i_pseudocount = INTEGER_VALUE(pseudocount);


  int s1 = LENGTH(v1);
  int s2 = LENGTH(v2);
	
  double* p_v1;
  double* p_v2;
  p_v1 = NUMERIC_POINTER(v1);
  p_v2 = NUMERIC_POINTER(v2);

  ARRAY a1={p_v1,s1};
  ARRAY a2={p_v2,s2};

  //
  SEXP ret;
  PROTECT(ret = NEW_NUMERIC(1));
  nProtected++;

  NUMERIC_POINTER(ret)[0] = __gmd0(a1, a2, i_pseudocount);
  
  UNPROTECT(nProtected);
  return(ret);
}




SEXP gmd(SEXP v1, SEXP v2, SEXP pseudocount){
  int j;
  int nProtected = 0;

  PROTECT(v1 = AS_NUMERIC(v1)); //a vector of real numbers
  nProtected++;

  PROTECT(v2 = AS_NUMERIC(v2));
  nProtected++;

  int i_pseudocount;
  i_pseudocount = INTEGER_VALUE(pseudocount);

  //
  int s1 = LENGTH(v1);
  int s2 = LENGTH(v2);
  double* p_v1;
  double* p_v2;
  p_v1 = NUMERIC_POINTER(v1);
  p_v2 = NUMERIC_POINTER(v2);


  ARRAY a1={p_v1,s1};
  ARRAY a2={p_v2,s2};

  // position
  SEXP v_position;
  int *p_v_position;   
  int s3 = s1+1+s2;

  PROTECT(v_position=NEW_INTEGER(s3));
  nProtected++;
  p_v_position=INTEGER_POINTER(v_position);

  ARRAYINT position;
  position.p=p_v_position;
  position.size=s3;

  // res
  SEXP res;
  PROTECT(res = NEW_NUMERIC(1));
  nProtected++;

  NUMERIC_POINTER(res)[0] = __gmd(a1, a2, i_pseudocount, position);

  for (j=0;j<s3;j++){
    INTEGER_POINTER(v_position)[j]=position.p[j];
  }

  //ret
  SEXP ret;
  PROTECT(ret = allocVector(VECSXP, 2)); // Creating a list with 2 vector elements
  nProtected++;
  SET_VECTOR_ELT(ret, 0, res);         
  SET_VECTOR_ELT(ret, 1, v_position);      

  // names
  char *names[2] = {"distance", "position"};
  SEXP ret_names;
  PROTECT(ret_names = allocVector(STRSXP, 2));
  nProtected++;
  int i;
  for(i = 0; i < 2; i++){
    SET_STRING_ELT(ret_names, i,  mkChar(names[i]));
  }
  setAttrib(ret, R_NamesSymbol, ret_names); //and attaching the vector names
  

  // free
  UNPROTECT(nProtected);
  return ret;
}

