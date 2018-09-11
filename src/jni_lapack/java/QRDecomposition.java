package JAMAJniLite;

import JAMAJniLite.*;
import java.lang.*;
/** QR Decomposition.
 <P>
    For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
    orthogonal dormqrmatrix Q and an n-by-n upper triangular matrix R so that
    A = Q*R.
 <P>
    The QR decompostion always exists, even if the matrix does not have
    full rank, so the constructor will never fail.  The primary use of the
    QR decomposition is in the least squares solution of nonsquare systems
    of simultaneous linear equations.  This will fail if isFullRank()
    returns false.
 */

public class QRDecomposition implements java.io.Serializable {
 static {
    /* load library (which will contain lapack functions.)*/
    System.loadLibrary("lapack_QRDecomposition");
    System.loadLibrary("blas_lite"); //dtrsm
 }
     
    /* ------------------------
     *  Class variables
     *  ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] QR;
    private double[] tau;

    private int[] info = new int[]{0};
    double[] work = new double[1000];
    int lwork = 50;
    int[] iwork = new int[24];
    
    /** Row and column dimensions, and min(m, n).*/
    private int m, n, l;

    /* ------------------------
     Constructor
     * ------------------------ */
    public QRDecomposition (Matrix Arg) {
	    m = Arg.getRowDimension();
	    n = Arg.getColumnDimension();
        l = Math.min(m, n);

        QR = Arg.getColumnPackedCopy();
	    int lda = m;
        tau = new double[l];
	    int matrix_layout = QRDecomposition.LAYOUT.ColMajor;
        dgeqrf(matrix_layout, m, n, QR, lda, tau, work, lwork, info);
        
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */

    /** Is the matrix full column rank?
     @return     true if R, and hence A, has full rank.
     */
    
    public boolean isFullRank () {
        for(int j = 0; j < l; j++){
            if (QR[j * m + j] == 0)
                return false;
        }
        return true;
    }
    
    /** no get H */
    
    /** Return the upper triangular factor*/
    public Matrix getR () {
        Matrix X = new Matrix(l,n);
	    double[][] R = X.getArray();
	    for (int i = 0; i < l; i++) {
		    for (int j = 0; j < n; j++) {
                if(i <= j){
                    R[i][j] = QR[i+j*m];
                } else {
                    R[i][j] = 0.0;
                }
		    }
	    }
	    return X;
    }


   /** Generate and return the (economy-sized) orthogonal factor*/
    public Matrix getQ () {
        int lda = m;
        int k = l;
        int matrix_layout = QRDecomposition.LAYOUT.ColMajor;
        double[] Q = new double[m * n];
        for(int i = 0; i < (m * n); i++){
            Q[i] = QR[i];
        }
        dorgqr(matrix_layout, m, n, k, Q, lda, tau, work, lwork, info);
        return new Matrix(Q, m);
        /** here n has to be less than m; otherwise dorgqr will fail ?*/
    }

    /** Least squares solution of A*X = B
     @param B    A Matrix with as many rows as A and any number of columns.
     @return     X that minimizes the two norm of Q*R*X-B.
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is rank deficient.
     */

    public Matrix solve (Matrix B) {
        if (B.getRowDimension() != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!this.isFullRank()) {
            throw new RuntimeException("Matrix is rank deficient.");
        }
        
        // Copy right hand side
        int nx = B.getColumnDimension();
        double[] b = B.getColumnPackedCopy();
        
        
        // Compute Y = transpose(Q)*B
        int matrix_layout = QRDecomposition.LAYOUT.ColMajor;
        char side = QRDecomposition.SIDE.Left;
        char trans = QRDecomposition.TRANSPOSE.Trans;
        int k = l, lda = m, ldb = m;
        dormqr(matrix_layout, side, trans, m, nx, k, QR, lda, tau, b, ldb, work, lwork, info);
        
        // Solve R*X = Y;
        char uplo = QRDecomposition.UPLO.Upper;
        char transa = QRDecomposition.TRANSPOSE.NoTrans;
        char diag = QRDecomposition.DIAG.NoUnit;
        double one = 1.0;
        Matrix.dtrsm(Matrix.LAYOUT.ColMajor, Matrix.SIDE.Left, Matrix.UPLO.Upper, Matrix.TRANSPOSE.NoTrans, Matrix.DIAG.NonUnit, m, n, one, QR, b);
        
        Matrix C = new Matrix(b, B.getRowDimension());
        return C;
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
    
    public final static class SIDE {
        private SIDE() {}
        public final static char Left = 'L';            /** apply Q or Q**T from the Left */
        public final static char Right= 'R';            /** apply Q or Q**T from the Right */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    public final static class DIAG {
        private DIAG() {}
        public final static char Unit = 'U';            /** assumed to be unit */
        public final static char NoUnit= 'N';           /**  not assumed to be unit*/
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
    
    
    /* QR */
    public static native int dgeqrf(int matrix_layout, int m, int n, double[] a,
                                    int lda, double[] tau, double[] work, int lwork,
                                    int[] info);
    
    public static native int dorgqr(int matrix_layout, int m, int n, int k,
                                    double[] a, int lda, double[] tau, double[] work,
                                    int lwork, int[] info);
    
    public static native int dormqr(int matrix_layout, char side, char trans,
                                    int m, int n, int k, double[] a, int lda, double[] tau,
                                    double[] c, int ldc, double[] work,
                                    int lwork, int[] info);


    private static final long serialVersionUID = 1;
    /**inform java virtual machine that function is defined externally*/
    
}







