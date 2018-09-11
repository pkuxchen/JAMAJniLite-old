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
extern int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
                  double *wr, double *wi, double *vl, int *ldvl, double *vr,
                  int *ldvr, double *work, int *lwork, int *info);

extern int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
                   int *lda, double *s, double *u, int *ldu, double *vt,
                   int *ldvt, double *work, int *lwork, int *info);

extern int dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s,
                   double *u, int *ldu, double *vt, int *ldvt, double *work,
                   int *lwork, int *iwork, int *info);



JNIEXPORT jint JNICALL Java_JAMAJniLite_SingularValueDecomposition_dgeev
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


JNIEXPORT jint JNICALL Java_JAMAJniLite_SingularValueDecomposition_dgesvd
  (JNIEnv *env, jclass obj, jint layout, jchar jjobu, jchar jjobvt, jint m,
   jint n, jdoubleArray ja, jint lda, jdoubleArray js, jdoubleArray ju,
   jint ldu, jdoubleArray jvt, jint ldvt, jdoubleArray jwork, jint lwork,
   jintArray jinfo){
      
      double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
      double *s = (*env)-> GetDoubleArrayElements (env, js, NULL);
      double *u = (*env)-> GetDoubleArrayElements (env, ju, NULL);
      double *vt = (*env)-> GetDoubleArrayElements (env, jvt, NULL);
      double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
      int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
      char jobu = (char) jjobu;
      char jobvt = (char) jjobvt;
      int result;
      int min = m;

      if (min > n) min = n;
  
      if (layout == jniColMajor) {
          result = dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt,
                           &ldvt, work, &lwork, info);
      } else if (layout == jniRowMajor) {
          if (lda < n) {
              fprintf(stderr, "** Illegal value of lda for row-major layout\n");
              return -1;
          } if ((jobu != 'A') && (jobu != 'a') && (jobu != 'S') &&
                (jobu != 's') && (jobu != 'O') && (jobu != 'o') &&
                (jobu != 'N') && (jobu != 'n')) {
              fprintf(stderr, "** Illegal value of jobu\n");
              return -1;
          } if ((jobvt != 'A') && (jobvt != 'a') && (jobvt != 'S') &&
                (jobvt != 's') && (jobvt != 'O') && (jobvt != 'o') &&
                (jobvt != 'N') && (jobvt != 'n')) {
              fprintf(stderr, "** Illegal value of jobvt\n");
              return -1;
          }
          
          int nrows_u = ( jobu == 'A' || jobu == 'S' ||
                         jobu == 'a' || jobu == 's') ? m : 1;
          int ncols_u = (jobu == 'A' || jobu == 'a') ? m :
                        ((jobu == 'S' || jobu == 's') ? min : 1);
          int nrows_vt = (jobvt == 'A' || jobvt == 'a') ? n :
                        ((jobvt == 'S' || jobvt == 's') ? min : 1);
          int lda_t = m;
          int ldu_t = nrows_u;
          int ldvt_t = nrows_vt;
          
          if (ldu < ncols_u) {
              fprintf(stderr, "** Illegal value of ldu for row-major layout\n");
              return -1;
          } if (ldvt < n) {
              fprintf(stderr, "** Illegal value of ldvt for row-major layout\n");
              return -1;
          }
          
          double *a_t = create_vectord(lda_t*n);
          double *u_t = NULL, *vt_t = NULL;
          RCswitch(a, a_t, m, n, lda, lda_t);
          if (jobu == 'A' || jobu == 'S' || jobu == 'a' || jobu == 's') {
              u_t = create_vectord(ldu_t*ncols_u);
          } if (jobvt == 'A' || jobvt == 'S' || jobvt == 'a' || jobvt == 's') {
              vt_t = create_vectord(ldvt_t*n);
          }
          result = dgesvd_(&jobu, &jobvt, &m, &n, a_t, &lda_t, s, u_t, &ldu_t,
                           vt_t, &ldvt_t, work, &lwork, info);
          CRswitch(a_t, a, m, n, lda_t, lda);
          if (jobu == 'A' || jobu == 'S' || jobu == 'a' || jobu == 's') {
              CRswitch(u_t, u, nrows_u, ncols_u, ldu_t, ldu);
              free(u_t);
          } if (jobvt == 'A' || jobvt == 'S' || jobvt == 'a' || jobvt == 's'){
              CRswitch(vt_t, vt, nrows_vt, n, ldvt_t, ldvt);
              free(vt_t);
          }
          free(a_t);
      }
      else {fprintf(stderr, "** Illegal layout setting\n");  return -1;}

      (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
      (*env)-> ReleaseDoubleArrayElements (env, js, s, 0);
      (*env)-> ReleaseDoubleArrayElements (env, ju, u, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jvt, vt, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
      (*env)-> ReleaseIntArrayElements (env, jinfo, info, 0);
  
      return result;
  }


JNIEXPORT jint JNICALL Java_JAMAJniLite_SingularValueDecomposition_dgesdd
  (JNIEnv *env, jclass obj, jint layout, jchar jjobz, jint m, jint n,
   jdoubleArray ja, jint lda, jdoubleArray js, jdoubleArray ju, jint ldu,
   jdoubleArray jvt, jint ldvt, jdoubleArray jwork, jint lwork,
   jintArray jiwork, jintArray jinfo)
{
    double *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *s = (*env)-> GetDoubleArrayElements (env, js, NULL);
    double *u = (*env)-> GetDoubleArrayElements (env, ju, NULL);
    double *vt = (*env)-> GetDoubleArrayElements (env, jvt, NULL);
    double *work = (*env)-> GetDoubleArrayElements (env, jwork, NULL);
    int *iwork = (*env)-> GetIntArrayElements (env, jiwork, NULL);
    int *info = (*env)-> GetIntArrayElements(env,jinfo,NULL);
    char jobz = (char) jjobz;
    int result;
    int min = m;

    if (min > n) min = n;
    if (layout == jniColMajor) {
        result = dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
                        work, &lwork, iwork, info);
    } else if (layout == jniRowMajor) {
        if ((jobz != 'A') && (jobz != 'a') && (jobz != 'S') &&
            (jobz != 's') && (jobz != 'O') && (jobz != 'o') &&
            (jobz != 'N') && (jobz != 'n')) {
            fprintf(stderr, "** Illegal value of jobz\n");  return -1;
        }
        int nrows_u = ((jobz == 'a' || jobz == 'A') ||
                        (jobz == 's' || jobz == 'S') ||
                        ((jobz == 'o' || jobz == 'O') && m<n) ) ? m : 1;
        int ncols_u = ((jobz == 'a' || jobz == 'A') ||
                    ((jobz == 'o' || jobz == 'O') && m<n) ) ? m :
                    ((jobz == 's' || jobz == 'S') ? min : 1);
        int nrows_vt = ((jobz == 'a' || jobz == 'A') ||
                    ((jobz == 'o' || jobz == 'O') && m>=n) ) ? n :
                    ((jobz == 's' || jobz == 'S') ? min : 1);
        int lda_t = m;
        int ldu_t = nrows_u;
        int ldvt_t = nrows_vt;
        double* a_t = create_vectord(lda_t*n);
        double* u_t = NULL;
        double* vt_t = NULL;
          
        if (lda < n) {
            fprintf(stderr, "** Illegal value of lda for row-major layout\n");
            return -1;
        } if (ldu < ncols_u) {
            fprintf(stderr, "** Illegal value of ldu for row-major layout\n");
            return -1;
        } if (ldvt < n) {
            fprintf(stderr, "** Illegal value of ldvt for row-major layout\n");
            return -1;
        }
          
        RCswitch(a, a_t, m, n, lda, lda_t);
        if((jobz == 'a' || jobz == 'A') || (jobz == 's' || jobz == 'S') ||
            ((jobz == 'o' || jobz == 'O') && (m<n) )) {
            u_t = create_vectord(ldu_t*ncols_u);
        } if((jobz == 'a' || jobz == 'A') || (jobz == 's' || jobz == 'S') ||
            ((jobz == 'o' || jobz == 'O') && (m>=n) )) {
            vt_t = create_vectord(ldvt_t*n);
        }
        result = dgesdd_(&jobz, &m, &n, a_t, &lda_t, s, u_t, &ldu_t,
                        vt_t, &ldvt_t, work, &lwork, iwork, info);
        CRswitch(a_t, a, m, n, lda_t, lda);
        if((jobz == 'a' || jobz == 'A') || (jobz == 's' || jobz == 'S') ||
            ((jobz == 'o' || jobz == 'O') && (m<n) )) {
            CRswitch (u_t, u, nrows_u, ncols_u, ldu_t, ldu);
            free(u_t);
        } if((jobz == 'a' || jobz == 'A') || (jobz == 's' || jobz == 'S') ||
            ((jobz == 'o' || jobz == 'O') && (m>=n) )) {
            CRswitch (vt_t, vt, nrows_vt, n, ldvt_t, ldvt);
            free(vt_t);
        }
        free(a_t);
    } else {
        fprintf(stderr, "** Illegal layout setting\n");
        return -1;}
      
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, js, s, 0);
    (*env)-> ReleaseDoubleArrayElements (env, ju, u, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jvt, vt, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jwork, work, 0);
    (*env)-> ReleaseIntArrayElements (env, jiwork, iwork, 0);
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
