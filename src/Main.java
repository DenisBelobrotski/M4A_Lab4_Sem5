public class Main {
    public static void main(String... args) {
        double intervalBottom = 0.5;
        double intervalUpper = 1;
        int nodeNum = 10;
        double step = (intervalUpper - intervalBottom) / nodeNum;
        double[] nodes = new double[nodeNum + 1];
        double[][] resultCauchy0;
        double[][] resultCauchy1;
        double[][] resultCauchy2;
        double[] m = {-0.5, -1};

        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = intervalBottom + i * step;
        }

        System.out.println("Узлы:");
        printVector(nodes);
        System.out.println();

        resultCauchy0 = nonHomogeneousCauchyProblem(nodes, step, 0, 0);
        System.out.println("Значения функции u0 в заданных узлах:");
        printVector(resultCauchy0[0]);
        System.out.println("Значения производной функции u0 в заданных узлах:");
        printVector(resultCauchy0[1]);
        System.out.println();

        resultCauchy1 = homogeneousCauchyProblem(nodes, step, 1, 0);
        System.out.println("Значения функции u1 в заданных узлах:");
        printVector(resultCauchy1[0]);
        System.out.println("Значения производной функции u1 в заданных узлах:");
        printVector(resultCauchy1[1]);
        System.out.println();

        resultCauchy2 = homogeneousCauchyProblem(nodes, step, 0, 1);
        System.out.println("Значения функции u2 в заданных узлах:");
        printVector(resultCauchy2[0]);
        System.out.println("Значения производной функции u2 в заданных узлах:");
        printVector(resultCauchy2[1]);
        System.out.println();

        System.out.println("Коэффициенты:");
        printVector(calcCoefs(resultCauchy0, resultCauchy1, resultCauchy2, m, nodeNum));
    }

    private static void printVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            System.out.println(vector[i] + " ");
        }
    }

    private static double calcFunction(double x) {
        return x * x;
    }

    private static double calcP1(double x) {
        return 1 / x;
    }

    private static double calcP2(double x) {
        return -2;
    }

    private static double[][] nonHomogeneousCauchyProblem(double[] nodes, double step, double initialY1, double initialY2) {
        double[][] results = new double[2][nodes.length];
        results[0][0] = initialY1;
        results[1][0] = initialY2;
        for (int i = 0; i < nodes.length - 1; i++) {
            results[0][i + 1] = results[0][i] + step * (results[1][i] + (step / 2) *
                    (calcFunction(nodes[i]) - calcP1(nodes[i]) * results[1][i] - calcP2(nodes[i]) * results[0][i]));
            results[1][i + 1] = results[1][i] + step * (calcFunction(nodes[i] + step / 2) -
                    calcP1(nodes[i] + step / 2) * (results[1][i] + (step / 2) * (calcFunction(nodes[i]) -
                            calcP1(nodes[i]) * results[1][i] - calcP2(nodes[i]) * results[0][i])) -
                    calcP2(nodes[i] + step / 2) * (results[0][i] + (step / 2) * results[1][i]));
        }
        return results;
    }

    private static double[][] homogeneousCauchyProblem(double[] nodes, double step, double initialY1, double initialY2) {
        double[][] results = new double[2][nodes.length];
        results[0][0] = initialY1;
        results[1][0] = initialY2;
        for (int i = 0; i < nodes.length - 1; i++) {
            results[0][i + 1] = results[0][i] + step * (results[1][i] + (step / 2) *
                    (-calcP1(nodes[i]) * results[1][i] - calcP2(nodes[i]) * results[0][i]));
            results[1][i + 1] = results[1][i] + step * (-calcP1(nodes[i] + step / 2) * (results[1][i] + (step / 2) *
                    (-calcP1(nodes[i]) * results[1][i] - calcP2(nodes[i]) * results[0][i])) -
                    calcP2(nodes[i] + step / 2) * (results[0][i] + (step / 2) * results[1][i]));
        }
        return results;
    }

    private static double[] calcCoefs(double[][] result0, double[][] result1, double[][] result2, double[] m, int nodeNum) {
        double[][] mtrA = new double[2][2];
        double[] vectB = new double[2];

        mtrA[0][0] = result1[1][0];
        mtrA[0][1] = result2[1][0];
        mtrA[1][0] = result1[1][nodeNum];
        mtrA[1][1] = result2[1][nodeNum];

        vectB[0] = m[0] - result0[1][0];
        vectB[1] = m[1] - result0[1][nodeNum];

        return gauss(mtrA, vectB);
    }

    private static double[] gauss(double[][] mtrA, double[] vectB) {
        double[][] a = copyMatrix(mtrA);
        double[] b = copyVector(vectB);
        double[] x = new double[b.length];
        double[] xIndexes = new double[b.length];
        double max;
        int maxK;
        int index;

        for (int i = 0; i < x.length; i++) {
            xIndexes[i] = (new Integer(i)).doubleValue();
        }

        for (int k = 0; k < x.length; k++) {
            max = a[k][k];
            maxK = k;
            for (int i = k + 1; i < x.length; i++) {
                if (Math.abs(max) < Math.abs(a[k][i])) {
                    max = a[k][i];
                    maxK = i;
                }
            }
            if (maxK != k) {
                a = swapMatrixColumns(a, k, maxK);
                xIndexes = swapVectorElements(xIndexes, k, maxK);
            }
            for (int j = k; j < x.length; j++) {
                a[k][j] /= max;
            }
            b[k] /= max;
            for (int i = k + 1; i < x.length; i++) {
                for (int j = k + 1; j < x.length; j++) {
                    a[i][j] -= a[i][k] * a[k][j];
                }
                b[i] -= a[i][k] * b[k];
                a[i][k] = 0.0;
            }
        }

        for (int i = x.length - 1; i >= 0; i--) {
            index = (new Double(xIndexes[i])).intValue();
            x[index] = b[i];
            for (int j = i + 1; j < x.length; j++) {
                x[index] -= a[i][j] * x[(new Double(xIndexes[j])).intValue()];
            }
        }

        return x;
    }

    private static double[] copyVector(double[] vector) {
        double[] newVect = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            newVect[i] = vector[i];
        }
        return newVect;
    }

    private static double[] swapVectorElements(double[] vector, int from, int to) {
        double tmp = vector[from];
        vector[from] = vector[to];
        vector[to] = tmp;
        return vector;
    }

    private static double[][] copyMatrix(double[][] mtr) {
        double[][] newMtr = new double[mtr.length][mtr[0].length];
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[i].length; j++) {
                newMtr[i][j] = mtr[i][j];
            }
        }
        return newMtr;
    }

    private static double[][] swapMatrixColumns(double[][] mtr, int from, int to) {
        double tmp;
        for (int i = 0; i < mtr.length; i++) {
            tmp = mtr[i][from];
            mtr[i][from] = mtr[i][to];
            mtr[i][to] = tmp;
        }
        return mtr;
    }

}
