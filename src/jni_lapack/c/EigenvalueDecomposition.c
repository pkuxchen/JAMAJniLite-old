#include <jni.h>
#include <assert.h>
#include <stdlib.h>

#define jniRowMajor 101
#define jniColMajor 102

#define jniNoTrans 111
#define jniTrans   112
#define jniConjTrans    113

#define jniUpper   121
#define jniLower   122

#define jniNonUnit 131
#define jniUnit    132

#define jniLeft    141
#define jniRight   142

/* two functions that deal with matrix layout */
void CRswitch (double *in, double *out, int m, int n, int ldin, int ldout);

void RCswitch (double *in, double *out, int m, int n, int ldin, int ldout);

double *create_vectord (int dim);

/* Calling fortran lapack from liblapack */
extern void dtrmv_(char *uplo, char *trans, char *diag, int *n, double *A,
                   int *lda, double *x, int *incx);

extern int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
                  double *wr, double *wi, double *vl, int *ldvl, double *vr,
                  int *ldvr, double *work, int *lwork, int *info);

extern int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                  double *w, double *work, int *lwork, int *info);

extern void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
                   int *n, double *alpha, double *A, int *lda, double *B,
                   int *ldb);

//from blas
JNIEXPORT void Java_JAMAJniLite_EigenvalueDecomposition_dtrmv
(JNIEnv *env, jclass klass, jint Layout, jint Uplo, jint Trans, jint Diag,
 jint n, jdoubleArray A, jdoubleArray x, jint incx){
    
    /*  DTRMV  performs one of the matrix-vector operations
     x := A*x,   or   x := A**T*x,
     where x is an n element vector and  A is an n by n unit, or non-unit,
     upper or lower triangular matrix. */
    
    // LDA specifies the first dimension of A; lda = n
    
    char Ts, uplo, diag;
    
    if (Diag == jniNonUnit) {diag = 'N';}
    else if (Diag == jniUnit) {diag = 'U';}
    else {fprintf(stderr, "** Illegal Diag setting \n"); return;}
    
    if (Layout == jniRowMajor){
        
        if (Trans == jniNoTrans) {Ts = 'T';}
        else if (Trans == jniTrans) {Ts = 'N';}
        else if (Trans == jniConjTrans) {Ts = 'N';}
        else {fprintf(stderr, "** Illegal Trans setting \n"); return;}
        
        if (Uplo == jniUpper) {uplo = 'L';}                // A is an upper triangular matrix.
        else if (Uplo == jniLower) {uplo = 'U';}           // A is a lower triangular matrix
        else {fprintf(stderr, "** Illegal Uplo setting \n"); return;}
        
    }
    else if(Layout == jniColMajor){
        
        if (Trans == jniTrans) {Ts = 'T';}
        else if (Trans == jniNoTrans) {Ts = 'N';}
        else if (Trans == jniConjTrans) {Ts = 'C';}
        else {fprintf(stderr, "** Illegal Trans setting \n"); return;}
        
        if (Uplo == jniUpper) {uplo = 'U';}                // A is an upper triangular matrix.
        else if (Uplo == jniLower) {uplo = 'L';}           // A is a lower triangular matrix
        else {fprintf(stderr, "** Illegal Uplo setting \n"); return;}
    
    }
    else{fprintf(stderr, "** Illegal Matrix_Layout setting \n"); return;}
    
    double *AElems, *xElems;
    
    AElems = (*env)-> GetDoubleArrayElements (env, A, NULL);
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    
    assert(AElems && xElems);
    
    dtrmv_(&uplo, &Ts, &diag, &n, AElems, &n, xElems, &incx);
    
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, AElems, JNI_ABORT);
}

JNIEXPORT jint JNICALL Java_JAMAJniLite_EigenvalueDecomposition_dgeev
  (JNIEnv *env, jclass obj, jint layout, jchar jjobvl, jchar jjobvr, jint n,
   jdoubleArray ja, jint lda, jdoubleArray jwr, jdoubleArray jwi,
   jdoubleArray jvl, jint ldvl, jdoubleArray jvr, jint ldvr,
   jdoubleArray jwork, jint lwork, jintArray jinfo)
{
      
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *wr = (*env)-> GetDoubleArrayElements (env, jwr, NULL);
    double *wi = (*env)-> GetDoubleArrayElements (env, jwi, NULL);
    double *vl = (*env)-> GetDoubleArrayElements (env, jvl, NULL);
    double *vr = (*env)-> GetDoubleArrayElements (env, jvr, NULL);
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char jobvr = (char) jjobvr;
    char jobvl = (char) jjobvl;
    int result;
  
    if (layout == jniColMajor){
        result = dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
                        work, &lwork, info);
    } else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        } if (ldvl < n) {
            fprintf(stderr, "** Illegal value of ldvl for row-major layout\n");
            return -1;
        } if (ldvr < n) {
            fprintf(stderr, "** Illegal value of ldvr for row-major layout\n");
            return -1;
        } if ((jobvl !='N') && (jobvl != 'n') && (jobvl != 'V') && (jobvl != 'v')){
            fprintf(stderr, "** Illegal value of jobvl\n");
            return -1;
        } if ((jobvr !='N') && (jobvr != 'n') && (jobvr != 'V') && (jobvr != 'v')){
            fprintf(stderr, "** Illegal value of jobvr\n");
            return -1;
        }
          
        double *a_t = create_vectord(n*n);
        double *vl_t, *vr_t;
        int ldvl_t = n, ldvr_t = n, lda_t = n;
        RCswitch(a, a_t, n, n, lda, lda_t);
        if (jobvl == 'V' || jobvl == 'v') vl_t = create_vectord(n*n);
        if (jobvr == 'V' || jobvr == 'v') vr_t = create_vectord(n*n);
        result = dgeev_(&jobvl, &jobvr, &n, a_t, &lda_t, wr, wi, vl_t, &ldvl_t,
                          vr_t, &ldvr_t, work, &lwork, info);
        CRswitch(a_t, a, n, n, lda_t, lda);
        CRswitch(vl_t, vl, n, n, ldvl_t, ldvl);
        CRswitch(vr_t, vr, n, n, ldvr_t, ldvr);
        free(a_t);
        if (jobvl == 'V' || jobvl == 'v') free(vl_t);
        if (jobvr == 'V' || jobvr == 'v') free(vr_t);
    } else {fprintf(stderr, "** Illegal layout setting\n");  return -1;}

    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwr, wr, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwi, wi, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jvl, vl, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jvr, vr, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
      
    return result;
}

JNIEXPORT void JNICALL Java_JAMAJniLite_EigenvalueDecomposition_dsyev
  (JNIEnv *env, jclass obj, jint layout, jchar jjobz, jchar juplo, jint n,
   jdoubleArray ja, jint lda, jdoubleArray jw, jdoubleArray jwork,
   jint lwork, jintArray jinfo)
{
      
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *w = (*env)-> GetDoubleArrayElements (env, jw, NULL);
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char jobz = (char) jjobz;
    char uplo = (char) juplo;
  
    if (layout == jniColMajor){
        dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, info);
    } else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return;
        }
          
        int lda_t = n;
        double* a_t = create_vectord(lda_t*n);
        RCswitch(a, a_t, n, n, lda, lda_t);
        dsyev_(&jobz, &uplo, &n, a_t, &lda_t, w, work, &lwork, info);
        CRswitch(a_t, a, n, n, lda_t, lda);
        free(a_t);
    } else {
        fprintf(stderr, "** Illegal layout setting\n");
        return;
    }
      
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jw, w, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
      
    return;
}

//from blas
JNIEXPORT void Java_JAMAJniLite_EigenvalueDecomposition_dtrsm
(JNIEnv *env, jclass klass, jint Layout, jint Side, jint Uplo, jint TransA,
 jint Diag, jint m, jint n, jdouble alpha, jdoubleArray  A, jdoubleArray B){
    
    /* dtrmm: performs one of the matrix-matrix operations
     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
     where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     op( A ) = A   or   op( A ) = A**T.*/
    
    double *aElems, *bElems;
    
    int lda;
    char side, uplo, TA, diag;
    
    if (TransA == jniNoTrans) {TA = 'N';}
    else if (TransA == jniTrans) {TA = 'T';}
    else if (TransA == jniConjTrans) {TA = 'C';}
    else {fprintf(stderr, "** Illegal TransA setting \n"); return;}
    
    if (Diag == jniNonUnit) {diag = 'N';}
    else if (Diag == jniUnit) {diag = 'U';}
    else {fprintf(stderr, "** Illegal Diag setting \n"); return;}
    
    if (Layout == jniRowMajor){
        
        if (Side == jniLeft){ side = 'R'; lda = m;}        //
        else if(Side == jniRight){side = 'L'; lda = n;}    //
        else{fprintf(stderr, "** Illegal Side setting \n"); return;}
        
        if (Uplo == jniUpper){ uplo = 'L';}                //
        else if(Uplo == jniLower){ uplo = 'U';}
        else{fprintf(stderr, "** Illegal Uplo setting \n"); return;}
        
        aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
        bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);
        
        assert(aElems && bElems);
        
        dtrsm_(&side, &uplo, &TA, &diag, &n, &m, &alpha, aElems, &lda,
               bElems, &n);
        
    }
    else if(Layout == jniColMajor){
        
        if (Side == jniLeft){ side = 'L'; lda = m;}        //
        else if(Side == jniRight){side = 'R'; lda = n;}    //
        else{fprintf(stderr, "** Illegal Side setting \n"); return;}
        
        if (Uplo == jniUpper){ uplo = 'U';}                //
        else if(Uplo == jniLower){ uplo = 'L';}
        else{fprintf(stderr, "** Illegal Uplo setting \n"); return;}
        
        aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
        bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);
        
        assert(aElems && bElems);
        
        dtrsm_(&side, &uplo, &TA, &diag, &m, &n, &alpha, aElems, &lda,
               bElems, &m);
        
    }
    else{fprintf(stderr, "** Illegal Matrix_Layout setting \n"); return;}
    
    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
    
}



/* switch row-major and col-major*/
void RCswitch (double *in, double *out, int m, int n, int ldin, int ldout)
{
    int i,j;
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            out[j * ldout + i] = in[i * ldin + j];
        }
    }
    
}


void CRswitch (double *in, double *out, int m, int n, int ldin, int ldout)
{
    int i, j;

    for(i = 0; i < m; i++)
    {
        for( j = 0; j < n; j++)
        {
            out[i * ldout + j] = in[j * ldin + i];
        }
    }

}


double *create_vectord (int dim)
{
  double *vector;
  
  vector = (double *) malloc (dim * sizeof (double));
  assert(vector);

  return vector;
}
