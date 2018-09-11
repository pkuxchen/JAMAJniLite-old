package JAMAJniLite;

import java.lang.*;

public class LUDecomposition implements java.io.Serializable {
 static {
    /* load library (which will contain blas functions.)*/
    System.loadLibrary("lapack_LUDecomposition");
 }
   
    /* ------------------------
     * Class variables
     * ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] LU;

    /** Row and column dimensions, min(m, n) and INFO.*/
    private int m, n, l, info;
    private int[] inf = new int[]{0};
    
    private double pivsign;

    /** Internal storage of pivot vector.*/
    private int[] ipiv, piv;

    /* ------------------------
     * Constructor
     * ------------------------ */
    public LUDecomposition (Matrix Arg) {
        
        LU = Arg.getColumnPackedCopy();
        m = Arg.getRowDimension();
        n = Arg.getColumnDimension();
        l = Math.min(m, n);
        ipiv = new int[l];
        piv = new int [l];
        pivsign = 1.0;
        int temp;
        
        for (int i = 0; i < l; i++) {
            ipiv[i] = i;
            piv[i] = i;
        }
        int matrix_layout = LUDecomposition.LAYOUT.ColMajor;
        info = dgetrf(matrix_layout, m, n, LU, m, ipiv, inf);

        for (int i = 0; i < l; i++){
            temp = piv[ipiv[i] - 1];
            piv[ipiv[i] - 1] = piv[i];
            piv[i] = temp;
            if (piv[i] != i ){pivsign = -pivsign;}
        }
	}

    /* ------------------------
     * Public Methods
     * ------------------------ */

    /** Is the matrix nonsingular?*/
    public boolean isNonsingular () {
        if(info == 0 ){
            return true;
        } else if (info > 0){
            return false;
        } else {
            return false;
        }
        /* should I add check m == n? */
    }

   /** Return lower triangular factor */
    public Matrix getL () {
        int nrow = m, ncol = ((m >= n) ? n: m);
        Matrix X = new Matrix(nrow, ncol);
        double[][] L = X.getArray();
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                if (i > j) {
                    L[i][j] = LU[i + j * nrow];
                } else if (i == j) {
                    L[i][j] = 1.0;
                } else {
                    L[i][j] = 0.0;
                }
            }
        }
        return X;
    }

    /** Return upper triangular factor*/
    
    public Matrix getU () {
        int nrow = Math.min(m, n), ncol = n;
        Matrix X = new Matrix(nrow, ncol);
        double[][] U = X.getArray();
        
        for (int i = 0; i < nrow; i++){
            for (int j = 0; j < ncol; j++) {
                U[i][j] = ((i <= j) ? LU[i + j * m] : 0.0);
            }
        }
        return X;
    }
        
    /** Return pivot permutation vector */
    public int[] getPivot () {
        int[] p = new int[l];
        for (int i = 0; i < l; i++) {
            p[i] = piv[i];
        }
        return p;
    }
    
    /** Return pivot permutation vector as a one-dimensional double array
     @return     (double) piv
     */
    
    public double[] getDoublePivot () {
        double[] vals = new double[l];
        for (int i = 0; i < l; i++) {
            vals[i] = (double) piv[i];
        }
        return vals;
    }

    /** Determinant
     @return     det(A)
     @exception  IllegalArgumentException  Matrix must be square
     */

    public double det( ) {
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        
        double d = (double) 1.0;
        int temp = 0;
        
        for (int j = 0; j < n; j++){
            d *= LU[j * n + j];
        }
        return (d * pivsign);
    }
    
    
    /** Solve A * X = B */
    public Matrix solve (Matrix B) {
        if (B.getRowDimension() != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        if (!this.isNonsingular()) {
            throw new RuntimeException("Matrix is singular.");
        }
        
        int matrix_layout = LUDecomposition.LAYOUT.ColMajor;
        char trans = LUDecomposition.TRANSPOSE.NoTrans;
        int nrhs = B.getColumnDimension();
        int ldb = m;
        int lda = m;
        double[] b = B.getColumnPackedCopy();
        dgetrs(matrix_layout, trans, m, nrhs, LU, lda, ipiv, b, ldb, inf);
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
    
    /* LU */
    public static native int dgetrf(int matrix_layout, int m, int n, double[] a,
                                    int lda, int[] ipiv, int[] info);
    
    public static native void dgetrs(int matrix_layout, char trans, int n,
                                     int nrhs, double[] a, int lda, int[] ipiv,
                                     double[] b, int ldb, int[] info);
    
    private static final long serialVersionUID = 1;
}
    
    








