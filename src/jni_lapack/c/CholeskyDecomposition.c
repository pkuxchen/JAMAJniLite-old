#include <jni.h>
#include <assert.h>
#include <stdlib.h>

#define jniRowMajor 101
#define jniColMajor 102

/* two functions that deal with matrix layout */
void CRswitch (double *in, double *out, int m, int n, int ldin, int ldout);

void RCswitch (double *in, double *out, int m, int n, int ldin, int ldout);

double *create_vectord (int dim);

/* Calling fortran lapack from liblapack */
extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);

extern int dpotri_(char *uplo, int *n, double *a, int *lda, int *info); 

extern int dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                   double *b, int *ldb, int *info);


JNIEXPORT void JNICALL Java_JAMAJniLite_CholeskyDecomposition_dpotrf
  (JNIEnv *env, jclass obj, jint layout, jchar juplo, jint n, jdoubleArray ja,
   jint lda, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char uplo = (char) juplo;
  
    if (layout == jniColMajor) dpotrf_(&uplo, &n, a, &lda, info);
    else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return;
        }
        
        double *a_t = create_vectord(n*n);
        int lda_t = n;
        RCswitch(a, a_t, n, n, lda, lda_t);
        dpotrf_(&uplo, &n, a_t, &lda_t, info);
        CRswitch(a_t, a, n, n, lda_t, lda);
        free(a_t);
    } else { fprintf(stderr, "** Illegal layout setting \n");  return;}

    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
      
    return;
}


JNIEXPORT jint JNICALL Java_JAMAJniLite_CholeskyDecomposition_dpotri
  (JNIEnv *env, jclass obj, jint layout, jchar juplo, jint n, jdoubleArray ja,
   jint lda, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char uplo = (char) juplo;
    int result;
  
    if (layout == jniColMajor) result = dpotri_(&uplo, &n, a, &lda, info);
    else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        }
        double *a_t = create_vectord(n*n);
        int lda_t = n;
        RCswitch(a, a_t, n, n, lda, lda_t);
        result = dpotri_(&uplo, &n, a_t, &lda_t, info);
        CRswitch(a_t, a, n, n, lda_t, lda);
        free(a_t);
    } else {
        fprintf(stderr, "** Illegal layout setting\n");
        return -1;
    }
      
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
      
    return result;
}


JNIEXPORT jint JNICALL Java_JAMAJniLite_CholeskyDecomposition_dpotrs
  (JNIEnv *env, jclass obj, jint layout, jchar juplo, jint n, jint nrhs,
   jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *b = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char uplo = (char) juplo;
    int result;
      
    if (layout == jniColMajor){
        result = dpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
    }
    else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        }if (ldb < nrhs) {
            fprintf(stderr, "** Illegal value of ldb for row-major layout\n");
            return -1;
        }
          
        int lda_t = n, ldb_t = n;
        double *a_t = create_vectord(lda_t*n);
        double *b_t = create_vectord(ldb_t*nrhs);
        RCswitch(a, a_t, n, n, lda, lda_t);
        RCswitch(b, b_t, n, nrhs, ldb, ldb_t);
        result = dpotrs_(&uplo, &n, &nrhs, a_t, &lda_t, b_t, &ldb_t, info);
        CRswitch(b_t, b, n, nrhs, ldb_t, ldb);
        free(a_t);
        free(b_t);
    }
    else { fprintf(stderr, "** Illegal layout setting\n");  return -1;}

    (*env)-> ReleaseDoubleArrayElements (env, ja, a, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jb, b, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);

    return result;
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
