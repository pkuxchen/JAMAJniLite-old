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
extern int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
                   double *work, int *lwork, int *info);

extern int dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau,
                   double *work, int *lwork, int *info);

extern int dormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);


JNIEXPORT jint JNICALL Java_JAMAJniLite_QRDecomposition_dgeqrf
  (JNIEnv *env, jclass obj, jint layout, jint m, jint n, jdoubleArray ja,
   jint lda, jdoubleArray jtau, jdoubleArray jwork, jint lwork,
   jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *tau = (*env)-> GetDoubleArrayElements (env, jtau, NULL);
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    int result;
  
    if (layout == jniColMajor){
        result = dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
    } else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        }
          
        double *a_t = create_vectord(m*n);
        int lda_t = m;
        RCswitch(a, a_t, m, n, lda, lda_t);
        result = dgeqrf_(&m, &n, a_t, &lda_t, tau, work, &lwork, info);
        CRswitch(a_t, a, m, n, lda_t, lda);
        free(a_t);
    } else {
        fprintf(stderr, "** Illegal layout setting\n");
        return -1;
    }
      
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jtau, tau, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);

    return result;
}


JNIEXPORT jint JNICALL Java_JAMAJniLite_QRDecomposition_dorgqr
  (JNIEnv *env, jclass obj, jint layout, jint m, jint n, jint k,
   jdoubleArray ja, jint lda, jdoubleArray jtau, jdoubleArray jwork,
   jint lwork, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *tau = (*env)-> GetDoubleArrayElements (env, jtau, NULL);
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    int result;
  
    if (layout == jniColMajor) {
        result = dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
    } else if (layout == jniRowMajor) {
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        }
          
        double *a_t = create_vectord(m*n);
        int lda_t = m;
        RCswitch(a, a_t, m, n, lda, lda_t);
        result = dorgqr_(&m, &n, &k, a_t, &lda_t, tau, work, &lwork, info);
        CRswitch(a_t, a, m, n, lda_t, lda);
        free(a_t);
    } else {
        fprintf(stderr, "** Illegal layout setting\n");
        return -1;
    }

    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jtau, tau, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
    (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);

    return result;
}


JNIEXPORT jint Java_JAMAJniLite_QRDecomposition_dormqr
(JNIEnv *env, jclass klass, jint layout, jchar jside, jchar jtrans,
 jint m, jint n, jint k, jdoubleArray a, jint lda, jdoubleArray tau,
 jdoubleArray c, jint ldc, jdoubleArray jwork, jint lwork, jintArray jinfo){
    
    double *aElems, *tauElems, *cElems;
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char side = (char) jside;
    char trans = (char) jtrans;
    int result;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    cElems = (*env)-> GetDoubleArrayElements (env, c, NULL);
    
    assert(aElems && tauElems && cElems);
    
    result = dormqr_ (&side, &trans, &m, &n, &k, aElems, &lda, tauElems, cElems, &ldc, work, &lwork, info);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, c, cElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
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
