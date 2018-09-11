import JAMAJniLite.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class JAMAJniLiteExamplesLAPACK {
	private JAMAJniLiteExamplesLAPACK() {}
	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
		System.out.println("###   Testfile for jni_lapack   ### \n \n");
		System.out.println("###   Parameter Preparation   ###");
		int M=3, N=3, K=3;
		double[] A = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
		double[] B = new double[] {4, -2, -6, -2, 10, 9, -6, 9, 14};
		double[] D = new double[9];
		double[] D2 = new double[9];
		int[] ipiv = new int[] {0, 0, 0};
		double[] tau = new double[]{0, 0, 0};
		int[] jpvt = new int[] {1, 2, 3};

		double[] VL = new double[M * N];
		double[] VR = new double[M * N];
		double[] WR = new double[N];
		double[] WI = new double[N];
		double[] S = new double[N];
		double[] U = new double[M * M];
		double[] VT = new double[N * N];
		double[] superb = new double[N - 1];
        int[] info = new int[]{0};
        double[] work = new double[1000];
        int lwork = 3;
        int[] iwork = new int[24];
        //
        // Print the parameters
        //
		printMatrix("Matrix A", LUDecomposition.LAYOUT.RowMajor, A, M, K);
		printMatrix("Matrix B", LUDecomposition.LAYOUT.RowMajor, B, K, N);
        //
        // Set Option
        //
		int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
		char trans = LUDecomposition.TRANSPOSE.NoTrans;
		char Uplo = LUDecomposition.UPLO.Upper;
		char Jobvl = LUDecomposition.JOBV.Compute;
		char Jobvr = LUDecomposition.JOBV.Compute;
		char Jobu = LUDecomposition.JOB.All;
		char Jobvt = LUDecomposition.JOB.All;
		char Jobz = LUDecomposition.JOB.Overwritten;
		int itype = LUDecomposition.ITYPE.first;
        //
        //
		/* ---- LU ---- */
        //
		System.out.println("\n \n");
		System.out.println("###   LU   ###");
        //
		/* dgetrf */
        //
		System.out.println();
		System.out.println("The LU decomposition A=PLU (L has unit diagnoal elements):");
		D = A.clone();
        int m = M;//m = the number of rows of the matrix a
        int n = N;//n = the number of columns of the matrix a
        double[] a = D;//a double precision, dimension(m,n)
        int lda = N;//lda = leading dimension of array a
        //ipiv = pivot indices       
        LUDecomposition.dgetrf( matrix_layout, m, n, a, lda, ipiv, info);
        printMatrix("  The LU matrix of A:",matrix_layout, a, m, n);
        printIntArray("  The permutation vector ipiv: ", ipiv, n);
        //
        /* dgetrs */
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();
        a = D;
        LUDecomposition.dgetrf(matrix_layout, m, n, a, lda, ipiv, info);

        //trans specifies the form of the system of equations
        n = M;//n = order of matrix a
        int nrhs = N;//nrhs = number of columns of matrix b 
        //a is computed by dgetrf
        lda = n;//lda = leading dimension of a
        //ipiv = pivot indices from dgetrf
        double[] b = D2;//b double precision, dimension(n,nrhs)
        int ldb = nrhs;//ldb = leading dimension of b
        LUDecomposition.dgetrs(matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb, info);
        printMatrix("  The solution  matrix of BX=A by LU is:", matrix_layout, b, n, nrhs);
        
        //
        /* ---- Cholesky ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Cholesky   ###");
        //
        /* dpotrf */
        //
        System.out.println();
        D = B.clone();

        char uplo = Uplo;//uplo specifies which part of a is stored
        n = M;//n = order of matrix a
        a = D;//a is symmetric matrix
        lda = M;//lda = leading dimension of a
        CholeskyDecomposition.dpotrf(matrix_layout, uplo, n, a, lda, info);
        printMatrix("  The Cholesky factorization of matrix B is :", matrix_layout, a, n, n);
        
        //
        /*  dpotrs */
        //
        System.out.println();
        D = B.clone();
        D2 = A.clone();

        //uplo specifies which part of a is stored
        n = M;//n = order of matrix a
        a = D;//a is symmetric matrix
        lda = n;//lda = leading dimension of a
        CholeskyDecomposition.dpotrf(matrix_layout, uplo, n, a, lda, info);

        //uplo specifies which part of a is stored
        n = M;//n = order of matrix a
        nrhs = N;//nrhs = number of columns of matrix b
        a = D;//a is computed by dpotrf
        lda = n;//lda = leading dimension of matrix a
        b = D2;//b = right hand matrix, dimension(n,nrhs)
        ldb = nrhs;//ldb = leading dimension of matrix b
        CholeskyDecomposition.dpotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb, info);
        printMatrix("  The solution  matrix of BX=A by Cholesky is :", matrix_layout, b, n, nrhs);
        //
        /* ---- QR ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   QR decomposition   ###");
        //
        /* dgeqrf */
        //
        System.out.println();
        System.out.println("A=QR:");
        D = A.clone();

        m = M;//m = number of rows of matrix a
        n = N;//n = number of columns of matrix a
        a = D;//a double precision, dimension(m,n)
        lda = n;//lda = leading dimension of matrix a
        //tau = scalar factors of the elementary reflectors
        QRDecomposition.dgeqrf(matrix_layout, m, n, a, lda, tau, work, lwork, info);
        printMatrix("  The rewritten matrix A is :", matrix_layout, a, m, n);
        //
        /* dorgqr */
        //
        System.out.println();
        m = M;//m = number of rows of matrix q
        n = N;//n = number of columns of matrix q
        int k = M;//k = number of elementary reflectors whose product defines the matrix q
        a = D;//a is returned by dgeqrf
        lda = M;//lda = leading dimension of a
        //tau is returned by dgeqrf
        QRDecomposition.dorgqr(matrix_layout, m, n, k, a, lda, tau, work, lwork, info);
        printMatrix("  The matrix Q is :", matrix_layout, a, m, n);
        //
        /* ---- Eigenvector ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   Eigenvector   ###");
        //
        /* dgeev */
        //
        System.out.println(" \nGet Eigenvectors of A:");
        D = A.clone();

        char jobvl = Jobvl;//jobvl specifies whether left eigenvectors of a are computed
        char jobvr = Jobvr;//jobvr specifies whether right eigenvectors of a are computed
        n = N;//n = order of matrix a
        a = D;//a double precision
        lda = N;//lda = leading dimension of a
        double[] wr = WR;//wr = real part of eigenvalues
        double[] wi = WI;//wi = image part of eigenvalues
        double[] vl = VL;//vl = left eigenvectors
        int ldvl = N;//ldvl = leading dimension of vl
        double[] vr = VR;//vr = right eigenvectors
        int ldvr = N;//ldvr = leading dimension of vr
        lwork = 30;
        EigenvalueDecomposition.dgeev(matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
        lwork = 3;
        printMatrix("  The overwritten matrix A :", matrix_layout, a, n, n);
        printMatrix("\n  Left eigenvectors (in columns) :", matrix_layout, vl, n, ldvl);
        printMatrix("\n  Right eigenvectors (in columns):", matrix_layout, vr, n, ldvr);
        //
        /* ---- SVD ---- */
        //
        System.out.println("\n \n");
        System.out.println("###   SVD   ###");
        //
        /* dgesvd */
        //
        System.out.println();
        System.out.println(" \nSVD of A: A = U * SIGMA * transpose(V):");
        D = A.clone();

        char jobu = Jobu;//jobu specifies options for computing all or part of the matrix u
        char jobvt = Jobvt;//jobvt specifies options for computing all or part of the matrix vt
        m = M;//m = number of rows of input matrix a
        n = N;//n = number of columns of input matrix a
        a = D;//a double precision
        lda = n;//lda = leading dimension of matri a
        double[] s = S;//s = singular values of a
        double[] u = U;//u = U
        int ldu = M;//ldu = leading dimension of u
        double[] vt = VT;//vt = VT;
        int ldvt = N;//ldvt = leading dimension of vt
        lwork = 50;
        SingularValueDecomposition.dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, m, n);
        printMatrix("\n  The singular values of A :", matrix_layout, s, 1, m);
        printMatrix("\n  Left singular vectors (in columns) :", matrix_layout, u, m, n);
        printMatrix("\n  Right singular vectors (in rows):", matrix_layout, vt, n, n);
        //
        /* dgesdd */
        System.out.println();
        System.out.println(" \nSVD of A calculated by dgesdd:");
        D = A.clone();

        jobu = Jobu;//jobu specifies options for computing all or part of the matrix u 
        m = M;//m = number of rows of input matrix a
        n = N;//n = number of columns of input matrix a
        a = D;//a double precision
        lda = n;//lda = leading dimension of matrix a
        s = S;//s = singular values of a
        u = U;//u = U
        ldu = M;//ldu = leading dimension of u
        vt = VT;//vt = VT;
        ldvt = N;//ldvt = leading dimension of vt
        SingularValueDecomposition.dgesdd(matrix_layout, jobu, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
        printMatrix("  The overwritten matrix A :", matrix_layout, a, m, n);
        printMatrix("\n  The singular values of A :", matrix_layout, s, 1, m);
        printMatrix("\n  Left singular vectors (in columns) :", matrix_layout, u, m, n);
        printMatrix("\n  Right singular vectors (in rows):", matrix_layout, vt, n, n);       
        lwork = 3;
    }

    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
    	System.out.println(prompt);
    	if (layout == LUDecomposition.LAYOUT.ColMajor) {
    		for (int i=0; i<I; i++) {
    			for (int j=0; j<J; j++)
    				System.out.print("\t" + string(X[j*I+i]));
    			System.out.println();
    		}
    	}
    	else if (layout == LUDecomposition.LAYOUT.RowMajor){
    		for (int i=0; i<I; i++) {
    			for (int j=0; j<J; j++)
    				System.out.print("\t" + string(X[i*J+j]));
    			System.out.println();
    		}
    	}
    	else{System.out.println("** Illegal layout setting");}
    }

    private static void printIntArray(String prompt, int[] X, int L) {
    	System.out.println(prompt);
    	for (int i=0; i<L; i++) {
    		System.out.print("\t" + string(X[i]));
    	}
    	System.out.println();
    }

    /** Shorter string for real number. */
    private static String string(double re) {
    	String s="";
    	if (re == (long)re)
    		s += (long)re;
    	else
    		s += re;
    	return s;
    }
}
