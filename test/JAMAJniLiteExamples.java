import JAMAJni.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class JAMAJniLiteExamples {
	private JAMAJniLiteExamples() {}
	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Exaples for JAMAJni   ### \n \n");
        System.out.println("###   Parameter Preparation   ###");
        
        int M=3, N=3;
        double alpha = 2;
        int matrix_layout = Matrix.LAYOUT.RowMajor;
        int[] pivot;

        double[][] a = new double[][] {{12,-51,4},{6,167,-68},{-4,24,-41}};
        double[][] b = new double[][] {{4,-2,-6},{-2,10,9},{-6,9,14}};
        double[][] c;
        //
        //Construct Matrix
        //
        Matrix A = new Matrix(a);
        Matrix B = new Matrix(b);

        Matrix C = new Matrix(M,N);
        Matrix D = new Matrix(M,N);
        Matrix X = new Matrix(M,N);
        Matrix Q = new Matrix(M,N);
        Matrix R = new Matrix(M,N);
        Matrix L = new Matrix(M,N);
        Matrix U = new Matrix(M,N);
        Matrix V = new Matrix(M,N);
        Matrix S = new Matrix(M,N);
		//
        // Print the parameters
        //
        printMatrix("Matrix A", matrix_layout, A.getArray(), M, N);
        printMatrix("Matrix B", matrix_layout, B.getArray(), M, N);
        System.out.println("M =" + M);
        System.out.println("N =" + N);
        System.out.println("alpha =" + alpha);
        //
        /* ---- Basic matrix algebra ---- */
        //
        System.out.println("\n\n###   Basic matrix algebra   ###");
        System.out.println();
        //
        //Matrix Addition
        //
        System.out.println("\n##  Addition: C = A + B  ##");
        C = A.plus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //MINUS
        //
        System.out.println("\n##  Subtraction: C = A - B  ##");
        C = A.minus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES
        //
        System.out.println("\n##  Multiplication: C = A * B  ##");
        C = A.times(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //SCALAR
        //
        System.out.println("\n##  Scalar multiplication: C = alpha * A  ##");
        C = A.times(alpha);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        /* ---- Factorizations ---- */
        //
        System.out.println("\n\n###   Factorizations & Solving Linear Systems   ###");
        System.out.println();
        //
        //LU Decomposition
        //
        System.out.println("\n##  LU Decomposition: A = L*U  ##");
        LUDecomposition LU = A.lu();
        L = LU.getL();
        U = LU.getU();
        pivot = LU.getPivot();
        X = LU.solve(B);
        printMatrix("L = ", matrix_layout, L.getArray(), M, N);
        printMatrix("U = ", matrix_layout, U.getArray(), M, N);
        printIntArray("The permutation vector pivot is : ", pivot, pivot.length);
        printMatrix("The solution matrix of AX=B by LU is", matrix_layout, X.getArray(), M, N);
        //
        //Cholesky Decomposition
        //
        System.out.println("\n##  Cholesky Decomposition of Matrix B  ##");
        CholeskyDecomposition Chol = B.chol(); 
        L = Chol.getL();
        X = Chol.solve(A);
        printMatrix("L = ", matrix_layout, L.getArray(), M, N);
        printMatrix("The solution matrix of BX=A by Cholesky is", matrix_layout, X.getArray(), M, N);
        //
        //QR Decomposition
        //
        System.out.println("\n##  QR Decomposition: A = Q*R  ##");
        QRDecomposition QR = A.qr();
        Q = QR.getQ();
        R = QR.getR();
        X = QR.solve(B);
        printMatrix("Q = ", matrix_layout, Q.getArray(), M, N);
        printMatrix("R = ", matrix_layout, R.getArray(), M, N);
        printMatrix("The solution matrix of AX=B by QR is", matrix_layout, X.getArray(), M, N);
        //
        //SVD Decomposition
        //
        System.out.println("\n##  SVD Decomposition: A = U*S*V  ##");
        SingularValueDecomposition SVD = A.svd();
        U = SVD.getU();
        V = SVD.getV();
        S = SVD.getS();
        printMatrix("U = ", matrix_layout, U.getArray(), M, N);
        printMatrix("V = ", matrix_layout, V.getArray(), M, N);
        printMatrix("S = ", matrix_layout, S.getArray(), M, N);
        //
        //Eigenvalue Decomposition
        //
        System.out.println("##  Eigenvalue Decomposition  ##");
        //
        //Asymmetric
        //
        System.out.println("\n#  Asymmetric: A*V = V*D  #");
        EigenvalueDecomposition Eig = A.eig();
        V = Eig.getV();
        D = Eig.getD();
        printMatrix("V = ", matrix_layout, V.getArray(), M, N);
        printMatrix("D = ", matrix_layout, D.getArray(), M, N);
        //
        //Symmetric
        //
        System.out.println("\n#  Symmetric: B = V*D*V'  #");
        Eig = B.eig();
        V = Eig.getV();
        D = Eig.getD();
        printMatrix("V = ", matrix_layout, V.getArray(), M, N);
        printMatrix("D = ", matrix_layout, D.getArray(), M, N);
    }
    //
    /* Print the matrix X */
    //
    private static void printMatrix(String prompt, int layout, double[][] X, int I, int J) {
        System.out.println(prompt);
        if (layout == Matrix.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j][i]));
                System.out.println();
            }
        }
        else if (layout == Matrix.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i][j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
    }
    //
    /* Print the array X */
    //
    private static void printIntArray(String prompt, int[] X, int L) {
        System.out.println(prompt);
        for (int i=0; i<L; i++) {
            System.out.print("\t" + string(X[i]));
        }
        System.out.println();
    }
    //
    /* Shorter string for real number */
    //
    private static String string(double re) {
        String s="";
        if (re == (long)re)
            s += (long)re;
        else
            s += re;
        return s;
    }
}
