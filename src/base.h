#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ------------------------------------------------------------------------
// Structure declarations
// ------------------------------------------------------------------------


typedef struct arrayint{
  /*  
      `ARRAYINT`, 
      with pionter and size   
  
  */
  int* p; //point to first element
  int size;  //length of array
} ARRAYINT;



typedef struct array{
  /*  
      `ARRAY`, 
      with pionter and size   
  
  */
  double* p; //point to first element
  int size;  //length of array
} ARRAY;




typedef struct array_2d{
  /*  
      `ARRAY_2D`, 
      with pionter and size (nrow, ncol)
  
  */
  double** p;
  int row;
  int col;
} ARRAY_2D;




// ------------------------------------------------------------------------
// Function declarations
// ------------------------------------------------------------------------
void absArray(ARRAY a, ARRAY res);
void addArrayWithInt(ARRAY a, int myint);
void concatenateArray3(ARRAY a1, ARRAY a2, ARRAY a3, ARRAY* pres);
void initializeArray(ARRAY a, double d);
double minAarray(ARRAY a);
void normalizeArray_inplace(ARRAY a);
void prefixsumArray(ARRAY a, ARRAY res);
// void printARRAYINT(ARRAYINT a);
// void printARRAY(ARRAY a);
// void printARRAY_2D(ARRAY_2D a);
void subArrayWithArray(ARRAY a1, ARRAY a2, ARRAY res);
double sumArray(ARRAY a);



// ------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------


void absArray(ARRAY a, ARRAY res){
  /*
    get abs() of an (double)array element-wise
  */
  double* a_p = a.p;
  double* res_p = res.p;
  int i;
  for(i=0;i<a.size;i++){
    *res_p++ = fabs(*a_p++);
  }
}


void addArrayWithInt(ARRAY a, int myint){
  /*
    (in place) increment each element of an arrary by the value of myint
  */
  double* a_p = a.p;
  int i;
  for(i=0;i<a.size;i++){
    *a_p++ += myint;
  }
}


void concatenateArray3(ARRAY a1, ARRAY a2, ARRAY a3, ARRAY* pres){
  /*
    Concatenate THREE (double)arrays and store to res
    Do "malloc" before calling, and "free" after.
  */
  double* a1_p=a1.p;
  double* a2_p=a2.p;
  double* a3_p=a3.p;
  //double* res_p=(*pres).p;
  double* res_p=pres->p;
  //
  int i;
  for (i=0;i<a1.size;i++){
    *res_p++ = *a1_p++;
  }
  for (i=0;i<a2.size;i++){
    *res_p++ = *a2_p++;
  }
  for (i=0;i<a3.size;i++){
    *res_p++ = *a3_p++;
  }
}




void initializeArray(ARRAY a, double d){
  /*
  (in place) Do "malloc" before calling, and "free" after.
  */
  double* a_p=a.p;
  int i;
  for(i=0;i<a.size;i++){
    *a_p++=d;
  }
}


double minAarray(ARRAY a){
  /*
    min()
  */
  int j,k;
  double res;

  double a0[a.size];// change pointer and its szie into an array container
  double* a_p=a.p;
  for (k=0;k<a.size;k++){
    a0[k]=*a_p++;
  }
  res=a0[0];//initialize
  for (j=1;j<a.size;j++){
    res = (res<=a0[j]) ? res : a0[j];
  }
  return res;
}


void normalizeArray_inplace(ARRAY a){ 
  /*
    (in place) arrayNormalize an (double)array
  */
  double* a_p = a.p;
  int i;
  double tmp_sum;
  tmp_sum=sumArray(a);

  for(i=0;i<a.size;i++){
    *a_p /= tmp_sum;
    a_p++;
  }
}



void prefixsumArray(ARRAY a, ARRAY res){
  /*
    get prefixsum() of an (double)array
    e.g.
    (1,2,3)->(1,3,6)
  */
  double* a_p = a.p;
  double* res_p = res.p;
  int i;
  double tmp=0.0;
  for(i=0;i<a.size;i++){
    tmp += *a_p++;
    *res_p++ = tmp;
  }
}


// void printARRAYINT(ARRAYINT a){
//   /* 
//      print an `ARRAYINT` of "int" elements, 
//      seperate by sep
//   */
//   char* sep=", "; //default argument
//   int i;

//   // start with "["
//   printf("%c",'[');

//   // the vector
//   for(i=0;i<a.size;i++){
//     printf("%d%s",*a.p++,sep);
//   }

//   // end with "]"
//   printf("%c%c",']','\n');
//   printf(", with legnth: %d",a.size);
//   printf("\n");
// }


// void printARRAY(ARRAY a){
//   /* 
//      print an `ARRAY` of "double" elements, 
//      seperate by sep
//   */
//   char* sep=", "; //default argument
//   double* a_p=a.p;
//   int i;
//   printf("%c",'[');
//   for(i=0;i<a.size;i++){
//     printf("%lg",*a_p++);
//     printf("%s",sep);
//   }
//   printf("%c%c",']','\n');
//   printf(", with legnth: %d",a.size);
//   printf("\n");
// }


// void printARRAY_2D(ARRAY_2D a){
//   // print an 2D array of "double" elements, seperate by sep
//   //char* sep=" "; //default argument
//   char* sep=", "; //default argument
//   int i,j;
//   printf("[\n");
//   for (i=0;i<a.row;i++){
//     printf("[");
//     for (j=0;j<a.col;j++){
//       //intPrint(i);
//       //intPrint(j);
//       printf("%lg",a.p[i][j]);
//     printf("%s",sep);
//     }
//     printf("], \n");
//   }
//   printf("]\n");
//   printf(", with row by col: %d x %d",a.row,a.col);
//   printf("\n");
// }


void subArrayWithArray(ARRAY a1, ARRAY a2, ARRAY res) {
  /* 
     subtract two arrays (of the same length) element-wise, 
     res=a1-a2.
  */
  int i;
  double* a1_p = a1.p;
  double* a2_p = a2.p;
  double* res_p = res.p;
  for(i=0;i<a1.size;i++){
    *res_p++ = (*a1_p++) - (*a2_p++);
  }
}



double sumArray(ARRAY a){
  /* 
     get sum() of an (double)array
  */
  double* a_p = a.p;
  int i;
  double res=0.0; //initialization
  for(i=0;i<a.size;i++){
    res += *a_p++;
  }
  return res;
}




// ------------------------------------------------------------------------
// tmp
// ------------------------------------------------------------------------





