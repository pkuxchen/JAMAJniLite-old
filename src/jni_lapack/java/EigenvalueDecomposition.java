package JAMAJniLite;

import java.lang.*;

/** Eigenvalues and eigenvectors of a real matrix.
 <P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and
    V.times(V.transpose()) equals the identity matrix.
 <P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
 **/

public class EigenvalueDecomposition  implements java.io.Serializable {
 static {
    /* load library (which will contain blas functions.)*/
    System.loadLibrary("lapack_EigenvalueDecomposition");
     
 }

   /* ------------------------
    * Class variables
    * ------------------------ */
    /** Row and column dimension (square matrix).*/
    private int n;

    /** Symmetry flag.*/
    private boolean issymmetric;

    /** Arrays for internal storage of eigenvalues.
     @serial internal storage of eigenvalues.
     */
    private double[] wr, wi;
    
    /**Array for elements storage of A */
    private double[] a;
    
    /** Array for internal storage of eigenvectors.
     @serial internal storage of eigenvectors.
     */
    private double[] VR;
    
    private int[] info = new int[]{0};
    private double[] work = new double[1000];
    private int lwork = 3;


    /* ------------------------
     * Constructor
     * ------------------------ */
    
    /** Check for symmetry, then construct the eigenvalue decomposition
     Structure to access D and V.
     @param Arg    Square matrix
     */
   
    public EigenvalueDecomposition (Matrix Arg) {
        double[][] A = Arg.getArray();
        a = Arg.getColumnPackedCopy();
        n = Arg.getColumnDimension();
        VR = new double[n * n];
        wr = new double[n];
        wi = new double[n];
        
        issymmetric = true;
        for (int j = 0; (j < n) & issymmetric; j++) {
            for (int i = 0; (i < n) & issymmetric; i++) {
                issymmetric = (A[i][j] == A[j][i]);
            }
        }
        
        int matrix_layout =  EigenvalueDecomposition.LAYOUT.ColMajor;
        int lda = n;
        if (issymmetric) {
            char jobz = EigenvalueDecomposition.JOBV.Compute;
            char uplo = EigenvalueDecomposition.UPLO.Upper;
            lwork = 50;
            dsyev(matrix_layout, jobz, uplo, n, a, lda, wr, work, lwork, info);
        } else{
            char jobvl = EigenvalueDecomposition.JOBV.NoCompute;
            char jobvr = EigenvalueDecomposition.JOBV.Compute;
            double[] vl = new double[1];
            int ldvl = 1;
            int ldvr = n;
            lwork = 30;
            dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, VR, ldvr, work, lwork, info);
        }
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */
    /** Return the eigenvector matrix */
    public Matrix getV () {
        if (issymmetric){
            return new Matrix(a, n);
        } else{
            return new Matrix(VR, n);
        }
        
    }
    
    /** Return the real parts of the eigenvalues
     @return     real(diag(D))
     */
    
    public double[] getRealEigenvalues () {
        return wr;
    }
    
    /** Return the imaginary parts of the eigenvalues
     @return     imag(diag(D))
     */
    
    public double[] getImagEigenvalues () {
        return wi;
    }
    
    /** Return the block diagonal eigenvalue matrix
     @return     D
     */
    public Matrix getD () {
        Matrix X = new Matrix(n,n);
        double[][] D = X.getArray();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D[i][j] = 0.0;
            }
            D[i][i] = wr[i];
            if (wi[i] > 0) {
                D[i][i+1] = wi[i];
            } else if (wi[i] < 0) {
                D[i][i-1] = wi[i];
            }
        }
        return X;
    }


    public final static class LAYOUT {
        private LAYOUT() {}
        public final static int RowMajor = 101;
        public final static int ColMajor = 102;
    }
    
    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static char NoTrans = 'N';         /** trans='N' */
        public final static char Trans= 'T';            /** trans='T' */
        public final static char ConjTrans= 'C';        /** trans='C' */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    public final static class JOBV {
        private JOBV() {}
        public final static char NoCompute = 'N';       /** eigenvectors are not computed */
        public final static char Compute= 'V';          /** eigenvectors are computed*/
    }
    
    public final static class JOB {
        private JOB() {}
        public final static char All = 'A';             /** all M columns of U are returned in array U */
        public final static char firstInU = 'S';        /** the first min(m,n) columns of U (the left singular
                                                         vectors) are returned in the array U;*/
        public final static char Overwritten = 'O';     /** the first min(m,n) columns of U (the left singular
                                                         vectors) are overwritten on the array A; */
        public final static char NoCompute = 'N';       /** no columns of U (no left singular vectors) are
                                                         computed.*/
    }
    
    public final static class ITYPE {
        private ITYPE(){}
        public final static int first = 1;
        public final static int second = 2;
        public final static int third = 3;
    }
    
    
    
    /* Eigenvector and SVD */
    public static native int dgeev(int matrix_layout, char jobvl, char jobvr,
                                   int n, double[] a, int lda, double[] wr, double[] wi,
                                   double[] vl, int ldvl, double[] vr, int ldvr,
                                   double[] work, int lwork, int[] info);
    
    /**inform java virtual machine that function is defined externally*/
    
    public static native void dsyev(int matrix_layout, char jobz, char uplo, int n,
                                    double[] a, int lda, double[] w, double[] work,
                                    int lwork, int[] info);
 
    private static final long serialVersionUID = 1;
    
}







