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
extern int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

extern void dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);

/* LU */
JNIEXPORT jint JNICALL Java_JAMAJniLite_LUDecomposition_dgetrf
  (JNIEnv *env, jclass obj, jint layout, jint m, jint n, jdoubleArray ja,
   jint lda, jintArray jipiv, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    int *ipiv = (*env)-> GetIntArrayElements(env,jipiv,NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    int result;

    if (layout == jniColMajor) {result = dgetrf_(&m, &n, a, &lda, ipiv, info);}
    else if (layout == jniRowMajor) {
        if (lda < n){
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        }
        double *a_t = create_vectord(m*n);
        int lda_t = m;
        RCswitch(a, a_t, m, n, lda, lda_t);
        result = dgetrf_(&m, &n, a_t, &lda_t, ipiv, info);
        CRswitch(a_t, a, m, n, lda_t, lda);
        free(a_t);
    } else {fprintf(stderr, "** Illegal Layout setting \n");  return -1;}

    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseIntArrayElements (env, jipiv, ipiv, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);

    return result;
}


JNIEXPORT void JNICALL Java_JAMAJniLite_LUDecomposition_dgetrs
  (JNIEnv *env, jclass obj, jint layout, jchar jtrans, jint n, jint nrhs,
   jdoubleArray ja, jint lda, jintArray jipiv, jdoubleArray jb, jint ldb,
   jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    int *ipiv = (*env)-> GetIntArrayElements(env,jipiv,NULL);
    double *b = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char trans = (char) jtrans;

    if (layout == jniColMajor) dgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
    else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return;
        } if (ldb < nrhs) {
            fprintf(stderr, "** Illegal value of ldb for row-major layout\n");
            return;
        }
        double *a_t = create_vectord(n*n);
        double *b_t = create_vectord(n*nrhs);
        int lda_t = n, ldb_t = n;
        RCswitch(a, a_t, n, n, lda, lda_t);
        RCswitch(b, b_t, n, nrhs, ldb, ldb_t);
        dgetrs_(&trans, &n, &nrhs, a_t, &lda_t, ipiv, b_t, &ldb_t, info);
        CRswitch(b_t, b, n, nrhs, ldb_t, ldb);
        free(a_t);
        free(b_t);
    } else { fprintf(stderr, "** Illegal layout setting \n");  return;}
  
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, JNI_ABORT);
    (*env)-> ReleaseIntArrayElements (env, jipiv, ipiv, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jb, b, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
      
    return;
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
