package JAMAJniLite;

import JAMAJniLite.*;
import java.lang.*;

public class CholeskyDecomposition implements java.io.Serializable {
	static {
        /* load library (which will contain blas functions.)*/
         System.loadLibrary("lapack_CholeskyDecomposition");
	}

    /* ------------------------
     *  Class variables
     *  ------------------------ */
    
    /** Array for internal storage of decomposition.
     @serial internal array storage.
     */
    private double[] a;

    private int[] info = new int[]{0};

    /** Row and column dimension (square matrix).*/
	private int n;

    /** Symmetric and positive definite flag.*/
	private boolean isspd;

    /* ------------------------
     * Constructor
     * ------------------------ */
    
    /** Cholesky algorithm for symmetric and positive definite matrix.
     Structure to access L and isspd flag.
     @param  Arg   Square, symmetric matrix.
     */
    
	public CholeskyDecomposition (Matrix Arg) {
        
        double[][] A = Arg.getArray();
        a = Arg.getColumnPackedCopy();
        n = Arg.getRowDimension();
        isspd = (Arg.getColumnDimension() == n);
        
        for (int j = 0; j < n; j++){
            for (int k = 0; k < j; k++){
                isspd = isspd & (A[k][j] == A[j][k]);
            }
        }
        
        
        int matrix_layout = CholeskyDecomposition.LAYOUT.ColMajor;
        char uplo = CholeskyDecomposition.UPLO.Lower;
        int lda = n;
        dpotrf(matrix_layout, uplo, n, a, lda, info);
        isspd = isspd & (info[0] == 0);
    }

   /* ------------------------
    * Public Methods
    * ------------------------ */

    /** Is the matrix symmetric and positive definite?*/
	public boolean isSPD () {
		return isspd;
	}

	/** Return triangular factor.*/
	public Matrix getL () {
		Matrix X = new Matrix(n,n);
		double[][] L = X.getArray();
        for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				L[i][j] = this.a[j*n+i];   /* check whether the switch of j and i can be avoid of not*/
			}
		}
		return X;
	}
        
	/** Solve A*X = B
     @param  B   A Matrix with as many rows as A and any number of columns.
     @return     X so that L*L'*X = B
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is not symmetric positive definite.
     */
    
	public Matrix solve (Matrix B) {
		if (B.getRowDimension() != n) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
		if (!isspd) {
			throw new RuntimeException("Matrix is not symmetric positive definite.");
		}
        /*here the isspd only check square matrix...*/
        
        int matrix_layout = CholeskyDecomposition.LAYOUT.ColMajor;
		char uplo = CholeskyDecomposition.UPLO.Lower;
		int nrhs = B.getColumnDimension();
		int ldb = nrhs;
		int lda = n;
		double[] b = B.getColumnPackedCopy();
		dpotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb, info);
		Matrix C = new Matrix(b,B.getRowDimension());
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
        public final static char Lower = 'L';            /** Lower triangular matrix*/
    }

    public final static class JOBV {
        private JOBV() {}
        public final static char NoCompute = 'N';       /** eigenvectors are not computed */
        public final static char Compute= 'V';          /** eigenvectors are computed*/
    }

    public final static class JOB {
        private JOB() {}
        public final static char All = 'A';             /** all M columns of U are returned in array U */
        public final static char firstInU = 'S';        /** the first min(m,n) columns of U (the left singular vectors) are returned in the array U;*/
        public final static char Overwritten = 'O';     /** the first min(m,n) columns of U (the left singularvectors) are overwritten on the array A; */
        public final static char NoCompute = 'N';       /** no columns of U (no left singular vectors) are computed.*/
    }

    public final static class ITYPE {
        private ITYPE(){}
        public final static int first = 1;
        public final static int second = 2;
        public final static int third = 3;
    }
    
    private static final long serialVersionUID = 1;
    
    /**inform java virtual machine that function is defined externally*/
    public static native void dpotrf(int matrix_layout, char uplo, int n,
                                     double[] a, int lda, int[] info);
    
    public static native int dpotri(int matrix_layout, char uplo, int n,
                                    double[] a, int lda, int[] info);
    
    public static native int dpotrs(int matrix_layout, char uplo, int n,
                                    int nrhs, double[] a, int lda, double[] b,
                                    int ldb, int[] info);
    
}




