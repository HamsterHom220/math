//#define EPS pow(10,-10)
#include <vector>
#include <iostream>

#ifndef MATRICES_IMPORTED
#include "matrices.cpp"
#define MATRICES_IMPORTED 1
#endif

using namespace std;


ColumnVector* AugmentedMatrix::solve_jacobi(double e) {
    if (lc == rows && ((SquareMatrix)*left_part()).diag_dominant()) {
        SquareMatrix L(rows);
        SquareMatrix D(rows);
        SquareMatrix U(rows);
        for (int i = rows - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                L.arr[i][j] = arr[i][j];
            }
        }
        for (int i = 0; i < rows; i++) {
            D.arr[i][i] = arr[i][i];
        }
        for (int i = 0; i < rows; i++) {
            for (int j = lc - 1; j > 0; j--) {
                U.arr[i][j] = arr[i][j];
            }
        }
        Matrix D_inv = *(D.inverse());
        IdentityMatrix I(rows);
        Matrix alpha = I - (D_inv * (*left_part()));
        ColumnVector beta = D_inv * (*right_part());

        //cout << "alpha:" << endl;
        //cout << alpha;
        //cout << "beta:" << endl;
        //cout << beta;

        // current solution
        ColumnVector x = beta;
        double error = e + 1;

        // improving the accuracy
        int step = 0;
        while (error >= e) {
            ColumnVector prev = x;
            //cout << "x(" << step++ << "):" << endl;
            //cout << prev;

            x = x + D_inv * ((*right_part()) - (*left_part()) * x);
            error = ((ColumnVector)(prev - x)).norm();
            //cout << "e: " << error << endl;
        }
        ColumnVector* res = new ColumnVector(x);
        //cout << "x(" << step++ << "):" << endl;
        //cout << *res;
        return res;
    }
    else {
        cout << "The method is not applicable!" << endl;
        return nullptr;
    }
}
ColumnVector* AugmentedMatrix::solve_seidel(double e) {
    if (lc == rows && ((SquareMatrix)*left_part()).diag_dominant()) {
        SquareMatrix L(rows);
        SquareMatrix D(rows);
        SquareMatrix U(rows);
        for (int i = rows - 1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                L.arr[i][j] = arr[i][j];
            }
        }
        for (int i = 0; i < rows; i++) {
            D.arr[i][i] = arr[i][i];
        }
        for (int i = 0; i < rows; i++) {
            for (int j = lc - 1; j > i; j--) {
                U.arr[i][j] = arr[i][j];
            }
        }
        Matrix D_inv = *(D.inverse());
        IdentityMatrix I(rows);
        Matrix alpha = I - (D_inv * (*left_part()));
        ColumnVector beta = D_inv * (*right_part());


        Matrix L_inv = *(L.inverse());

        SquareMatrix B(rows);
        SquareMatrix C(rows);

        for (int i = rows - 1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                B.arr[i][j] = alpha.arr[i][j];
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = lc - 1; j > i; j--) {
                C.arr[i][j] = alpha.arr[i][j];
            }
        }

        Matrix IB = I - B;
        Matrix IB_1 = *((SquareMatrix)IB).inverse();

        //cout << "beta:" << endl;
        //cout << beta;
        //cout << "alpha:" << endl;
        //cout << alpha;
        //cout << "B:" << endl;
        //cout << B;
        //cout << "C:" << endl;
        //cout << C;
        //cout << "I-B:" << endl;
        //cout << IB;
        //cout << "(I-B)_-1:" << endl;
        //cout << IB_1;

        // current solution
        ColumnVector x = beta;
        double error = e + 1;

        // improving the accuracy
        int step = 0;
        while (error >= e) {
            ColumnVector prev = x;
            //cout << "x(" << step++ << "):" << endl;
            //cout << prev;

            x = L_inv * (*right_part() - U * x);
            error = ((ColumnVector)(prev - x)).norm();
            //cout << "e: " << error << endl;
        }
        ColumnVector* res = new ColumnVector(x);
        //cout << "x(" << step++ << "):" << endl;
        //cout << *res;
        return res;
    }
    else {
        cout << "The method is not applicable!" << endl;
        return nullptr;
    }
}





