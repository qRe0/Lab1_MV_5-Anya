import java.util.Random;
import java.text.DecimalFormat;

public class Main {
    static final int N = 5; // Размерность матрицы
    static final double detCalc = 11180.; // Определитель матрицы A 5x5, вычисленный в WolframAlpha

    // Метод Гаусса
    static double[] solveWithGauss(double[][] A, double[] b, int N) {
        double[] x = new double[N];

        // Прямой ход
        for (int i = 0; i < N; i++) {
            // Выбор главного элемента
            int maxRow = i;
            for (int k = i + 1; k < N; k++) {
                if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i]))
                    maxRow = k;
            }

            // Перестановка строк
            if (maxRow != i) {
                double[] temp = A[i];
                A[i] = A[maxRow];
                A[maxRow] = temp;

                double tempB = b[i];
                b[i] = b[maxRow];
                b[maxRow] = tempB;
            }

            // Обнуление элементов ниже главного
            for (int k = i + 1; k < N; k++) {
                double factor = A[k][i] / A[i][i];
                b[k] -= factor * b[i];
                for (int j = i; j < N; j++) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }

        // Обратный ход
        for (int i = N - 1; i >= 0; i--) {
            x[i] = b[i] / A[i][i];
            for (int j = i + 1; j < N; j++) {
                x[i] -= A[i][j] * x[j] / A[i][i];
            }
        }

        return x;
    }

    // Метод Холецкого
    public static double[] solveWithCholesky(double[][] A, double[] b, int N) {
        // Разложение на нижнюю L и верхнюю U матрицы
        double[][] L = new double[N][N];
        double[][] U = new double[N][N];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++) {
                double s = 0;
                for (int k = 0; k < j; k++) {
                    s += L[i][k] * U[k][j];
                }
                L[i][j] = A[i][j] - s;
            }

            for (int j = i; j < N; j++) {
                double s = 0;
                for (int k = 0; k < i; k++) {
                    s += L[i][k] * U[k][j];
                }
                U[i][j] = (A[i][j] - s) / L[i][i];
            }
        }

        // Решение Ly = b
        // L - нижняя треугольная матрица
        double[] y = new double[N];
        y[0] = b[0] / L[0][0];

        for (int i = 1; i < N; i++) {
            double s = 0;
            for (int j = 0; j < i; j++) {
                s += L[i][j] * y[j];
            }
            y[i] = (b[i] - s) / L[i][i];
        }

        // Решение Ux = y
        // U - верхняя треугольная матрица
        double[] x = new double[N];
        x[N - 1] = y[N - 1] / U[N - 1][N - 1];

        for (int i = N - 2; i >= 0; i--) {
            double s = 0;
            for (int j = i + 1; j < N; j++) {
                s += U[i][j] * x[j];
            }
            x[i] = (y[i] - s) / U[i][i];
        }

        return x;
    }

    // Вычисление определителя
    static double computeDeterminant(double[][] A, int N) {
        double determinant = 1.0;

        for (int i = 0; i < N; i++) {
            determinant *= A[i][i];
        }

        return determinant;
    }

    public static void printMatrix(double[][] matrix, double[] b, int N, int maxRows, int maxCols) {
        DecimalFormat df = new DecimalFormat("0.0");

        for (int i = 0; i < Math.min(N, maxRows); i++) {
            for (int j = 0; j < Math.min(N, maxCols); j++) {
                System.out.print(String.format("%5s ", df.format(matrix[i][j])));
            }

            if (N > maxCols) {
                System.out.print(" ...");
            }

            System.out.println(" | " + df.format(b[i]));
        }

        if (N > maxRows) {
            for (int i = 0; i < maxCols; i++) {
                System.out.print(String.format("%5s ", "..."));
            }

            if (N > maxCols) {
                System.out.print(" ...");
            }

            System.out.println(" | ...");
        }

        System.out.println();
        DecimalFormat df8 = new DecimalFormat("0.00000000");
    }


    // Вывод вектора
    static void printVector(double[] vector, int N) {
        for (int i = 0; i < N; i++) {
            System.out.printf("| %f |%n", vector[i]);
        }
    }

    static double computeMaxNorm(double[][] A, double[] x, double[] b, int N) {
        double maxNorm = 0.0;

        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            double residual = Math.abs(sum - b[i]);
            if (residual > maxNorm) {
                maxNorm = residual;
            }
        }

        return maxNorm;
    }

    public static void main(String[] args) {
        Random rand = new Random();

        double[][] A = new double[N][N];
        double[] b = new double[N];

        // Заполнение матрицы
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    A[i][j] = 5 * Math.sqrt(i + 1);
                } else {
                    A[i][j] = Math.sqrt((i + 1)) + Math.sqrt((j + 1));
                }
            }
        }

        // Заполнение вектора случайными значениями
        for (int i = 0; i < N; i++) {
            b[i] = 4.5*Math.sqrt(i + 1);
        }

        System.out.println("Matrix:");
        printMatrix(A, b, N, 10, 10);

        // Решение системы методом Гаусса
        long startTime = System.nanoTime();
        double[] xGauss = solveWithGauss(A, b, N);
        System.out.println("\nSolve with Gauss:");
        printVector(xGauss, N);
        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1e9;
        System.out.println("Work time of Gauss: " + duration + " seconds\n");

        startTime = System.nanoTime();
        double[] xLU = solveWithCholesky(A, b, N);
        System.out.println("\nSolve with Cholesky:");
        printVector(xLU, N);
        endTime = System.nanoTime();
        duration = (endTime - startTime) / 1e9;
        System.out.println("Work time of Cholesky: " + duration + " seconds\n");

        // Вычисление максимум-нормы невязки для метода Гаусса
        double maxNormGauss = computeMaxNorm(A, xGauss, b, N);
        System.out.println("Max Norm Residual (Gauss): " + maxNormGauss + "\n");

        // Вычисление максимум-нормы невязки для метода Холецкого
        double maxNormLU = computeMaxNorm(A, xLU, b, N);
        System.out.println("Max Norm Residual (Cholesky): " + maxNormLU + "\n");

        // Вычисление определителя
        double determinantA = computeDeterminant(A, N);
        System.out.println("Det A: " + determinantA);
        System.out.println("DeFacto Det A: " + detCalc);
        System.out.println("Absolute Error: " + Math.abs(determinantA - detCalc));
        System.out.println("Relative Error: " + Math.abs(determinantA - detCalc) / Math.abs(determinantA) * 100 + "%");

        // Освобождение памяти
        A = null;
        b = null;
        xGauss = null;
    }
}
