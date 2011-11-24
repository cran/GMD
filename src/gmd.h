
#include "base.h"

// ------------------------------------------------------------------------
// Function Declarations (Function Prototypes)
// ------------------------------------------------------------------------
void get_sliderList_with_padding(ARRAY a, double pad, int wstep, ARRAY_2D res);
double __mdpa(ARRAY a1, ARRAY a2);
double __gmd0 (ARRAY P, ARRAY Q, int pseudocount);
double __gmd (ARRAY HX, ARRAY HY, int pseudocount, ARRAYINT position);



// ------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------


void get_sliderList_with_padding(ARRAY a, double pad, int wstep, ARRAY_2D res){
  /* 
     get a list of sub-list by sliding and padding given wstep.
     @wstep: width of a step
  
     res.row=finalLen-sliderLen+1; sliderLen=a.size
     res.col=finalLen
  */

  int i,j,k;
  double a0[a.size];// change pointer and its szie into an array container
  double* a_p=a.p;
  for (k=0;k<a.size;k++){
    a0[k]=*a_p++;
  }

  //calc the first index of each row, and save it in an array
  int r[res.row]; //first index of each row
  for (i=0,j=0;i<res.row;i+=wstep,j++){
    r[j]=i;
  }

  //assign element to res
  for (i=0;i<res.row;i++){
    for (j=0;j<res.col;j++){
      if (j<r[i] || j>(r[i]+a.size-1)){
        res.p[i][j]=pad;
      }
      else{
        res.p[i][j]=a0[j-r[i]];
      }
    }
  }
  
}




double __mdpa(ARRAY a1, ARRAY a2){
  /*
    Calc distance by MDPA (minimum difference of pair assignments)
    reference: (Cha and Srihari 2002) On measuring the distance between histograms.
    
    Get distance between two histograms with 
    same lengths (number of bins) and the same overall mass.
  */

  int size=a1.size;
  // int size2=a2.size;

  // // check if equal length
  // if (size2 != size){ 
  //   fprintf(stderr, "The `mdpa` function should take two equal-length vectors as input!\n"); 
  //   exit(1);
  // }
  // //

  double* p_t1;
  p_t1 = (double*)malloc(sizeof(double)*size); //... ... malloc
  double* p_t2;
  p_t2 = (double*)malloc(sizeof(double)*size); //... ... malloc
  double* p_t3;
  p_t3 = (double*)malloc(sizeof(double)*size); //... ... malloc

  ARRAY t1={p_t1,size};
  ARRAY t2={p_t2,size};
  ARRAY t3={p_t3,size};

  double res;
  //res of subtraction stored in t1
  subArrayWithArray(a1,a2,t1);
  //printf("t1=a1-a2: \n");
  //ARRAYPrint(t1);

  //res of prefixsumArray stored in t2
  prefixsumArray(t1,t2);
  //printf("t2=prefix-sum(t1): \n");
  //ARRAYPrint(t2);

  //abs()
  absArray(t2,t3);
  //printf("t3=abs(t2): \n");
  //ARRAYPrint(t3);

  //sum()
  res=sumArray(t3);

  free(p_t1);                                  //... ... free
  free(p_t2);                                  //... ... free
  free(p_t3);                                  //... ... free

  return res;
}


double __gmd0 (ARRAY P, ARRAY Q, int pseudocount){
  /*
    base version of `__gmd`, 
    a version without sliding.
    P and Q should be equal length 
    and will be normalized
  */
  
  // add pseudocount
  addArrayWithInt(P,pseudocount);
  addArrayWithInt(Q,pseudocount);
  
  // normalize
  normalizeArray_inplace(P); 
  normalizeArray_inplace(Q);

  double res;
  res=__mdpa(P, Q);

  return res;

}





double __gmd (ARRAY HX, ARRAY HY, int pseudocount, ARRAYINT position){
  /*
    Calc distance by GM-Distance (Generalized Minimum distance of distributions)
    reference: Zhao et al. (2011) Systematic Clustering of Transcription Start Site Landscapes. 
    Source: PLoS ONE 6(8): e23409. doi:10.1371/journal.pone.0023409
    
    Get distance between two histograms with 
    same lengths (number of bins) and the same overall mass.

    @HX,@HY: `ARRAY` of histograms
    @position
    initialize:
    ARRAYINT position;
    position.size=HX.size+1+HY.size;
    position.p=(int*)malloc(sizeof(int)*position.size);             //... ... malloc


    Return minimum distance and the shiftted positions of 
    the shorter (P) against the longer (Q) that produce it.
  */

  int i;
  int m=HX.size;
  int n=HY.size;

  //
  ARRAY P; 
  ARRAY Q;
  ARRAY tmp; 
  ARRAY padding; 

  // print
  //printARRAY(HX);
  //printARRAY(HY);

  // add pseudocount
  addArrayWithInt(HX,pseudocount);
  addArrayWithInt(HY,pseudocount);

  // save the shorter one as P, 
  // and save the longer one as `tmp`
  if (m <= n) {
    P=HX;
    tmp=HY;
    padding.size=m;
  }
  else{
    P=HY;
    tmp=HX;
    padding.size=n;
  }

  // malloc and initialize `padding`
  padding.p=(double*)malloc(sizeof(double)*padding.size); //... ... malloc
  initializeArray(padding, 0.0);
  
  // save padding+P+padding as Q
  Q.size=padding.size+tmp.size+padding.size;
  Q.p=(double*)malloc(sizeof(double)*Q.size);             //... ... malloc
  concatenateArray3(padding, tmp, padding, &Q);
  //printARRAY(Q);


  // P_res
  // padding P while sliding along Q, 
  // get a set of padded sliders, save as P_res
  int wstep=1;
  double pad=0.0;
  int finalLen=Q.size;
  int sliderLen=P.size;
  ARRAY_2D P_res;
  P_res.row=finalLen-sliderLen+1;
  P_res.col=finalLen;


  // memory allocation
  P_res.p=(double**)malloc(P_res.row*sizeof(double*));      //... ... malloc
  double* tmp_block;
  tmp_block = (double*)malloc(P_res.row * P_res.col * sizeof(double));//... ... malloc
  for (i=0;i<P_res.row;i++){
    P_res.p[i] = &tmp_block[i * P_res.col];
  }

  get_sliderList_with_padding(P, pad, wstep, P_res);
  // //printARRAY_2D(P_res);


  // normalize Q, in place
  normalizeArray_inplace(Q); 


  // normalize P, in place, via P_tmp
  ARRAY P_tmp;
  P_tmp.size=P_res.col;


  // //printf("P_tmp.size=%d\n",P_tmp.size); //??
  for (i=0;i<P_res.row;i++){
    P_tmp.p=P_res.p[i];
    normalizeArray_inplace(P_tmp);
  }


  // get a set of distances
  ARRAY a_dist;
  a_dist.size=P_res.row;
  a_dist.p=(double*)malloc(sizeof(double)*a_dist.size);             //... ... malloc


  for (i=0;i<P_res.row;i++){
    P_tmp.p=P_res.p[i];
    a_dist.p[i]=__mdpa(P_tmp, Q);
  }
  //printARRAY(a_dist);

  double res;
  res=minAarray(a_dist);

  // positions (in place)  
  // mark hit as `1`
  for (i=0;i<position.size;i++){ //position.size == P_res.row
    // //printf("a_dist.p[%d]=%lg\n",i,a_dist.p[i]);
    if (a_dist.p[i]==res){
      // //printf("min: a_dist.p[%d]=%lg\n",i,a_dist.p[i]);
      position.p[i]=1;
    }
    else {
      position.p[i]=0;
    }
  }
  

  // ------------------------------------------------------------------------
  // free  
  // free after sliding/padding

  //
  free(a_dist.p);                                               //... ... free

  // 
  free(P_res.p);                                                //... ... free
  free(tmp_block);                                              //... ... free

  //
  free(Q.p);                                                    //... ... free
  free(padding.p);                                              //... ... free


  //return
  //printf("gmd|res=%.7f\n",res);
  return res;

}

