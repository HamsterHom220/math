//#define EPS pow(10,-10)
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

#define MATRICES_IMPORTED 1

using namespace std;

class Matrix
{
protected:
    int rows;
    int cols;
    bool valid;
public:
    vector<vector<double>> arr;
    Matrix(int n, int m) {
        rows = n;
        cols = m;
        valid = true;

        for (int i = 0; i < rows; i++) {
            arr.push_back(vector<double>(cols));
            for (int j = 0; j < cols; j++) {
                arr[i][j] = 0;
            }
        }
    }

    int get_cols() const {
        return cols;
    }
    int get_rows() const {
        return rows;
    }

    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);

    Matrix& operator= (Matrix const& A)= default;

    Matrix operator + (Matrix A){
        Matrix res(rows,cols);
        if ((A.rows != rows)||(A.cols != cols)){
            res.valid = false;
            return res;
        }
        for (int i=0;i<rows;i++){
            for (int j=0;j<cols;j++){
                res.arr[i][j] = arr[i][j] + A.arr[i][j];
                /*if (abs(res.arr[i][j]) < EPS) {
                    res.arr[i][j] = 0;
                }*/
            }
        }
        return res;
    }
    Matrix operator - (Matrix A){
        Matrix res(rows,cols);
        if ((A.rows != rows)||(A.cols != cols)){
            res.valid = false;
            return res;
        }
        for (int i=0;i<rows;i++){
            for (int j=0;j<cols;j++){
                res.arr[i][j] = arr[i][j] - A.arr[i][j];
                /*if (abs(res.arr[i][j]) < EPS) {
                    res.arr[i][j] = 0;
                }*/
            }
        }
        return res;
    }
    Matrix operator * (Matrix A){
        Matrix res(rows,A.cols);
        if (cols != A.rows){
            res.valid = false;
            return res;
        }
        for (int i=0;i<rows;i++){
            for (int j=0;j<A.cols;j++){
                res.arr[i][j] = 0;
                for (int k=0;k<cols;k++){
                    res.arr[i][j] += arr[i][k]*A.arr[k][j];
                    /*if (abs(res.arr[i][j]) < EPS) {
                        res.arr[i][j] = 0;
                    }*/
                }
            }
        }
        return res;
    }
    Matrix T(){
        Matrix res(cols,rows);
        for (int i=0;i<cols;i++){
            for (int j=0;j<rows;j++){
                res.arr[i][j] = arr[j][i];
            }
        }
        return res;
    }
    Matrix operator -() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                arr[i][j] = -arr[i][j];
            }
        }
        return *this;
    }
};

ostream& operator << (ostream& stream, Matrix& c)
{
    if (!c.valid){
        stream << "Invalid matrix: "<< c << endl;
    }
    else{
        for (int i = 0; i < c.rows; i++)
        {
            for (int j = 0; j < c.cols-1; j++)
            {
                stream << std::fixed << std::setprecision(2) << c.arr[i][j]<<" ";
            }
            stream << std::fixed << std::setprecision(2) << c.arr[i][c.cols - 1] << endl;
        }
    }
    return stream;
}
istream& operator >> (istream& stream, Matrix& c)
{
    for (int i = 0; i < c.rows; i++)
    {
        for (int j = 0; j < c.cols; j++)
        {
            stream >> c.arr[i][j];
        }
    }
    return stream;
}

class SquareMatrix: public Matrix{
public:
    SquareMatrix(int n) : Matrix(n, n) {}
    SquareMatrix(Matrix const& matrix): Matrix(matrix) {
        if (cols!=rows){
            valid = false;
        }
    }

    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    SquareMatrix& operator= (SquareMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;

    double det();
    Matrix* inverse();
    bool diag_dominant() {
        for (int i = 0; i < rows; i++) {
            double sum = 0;
            for (int j = 0; j < cols; j++) {
                if (j != i) {
                    sum += arr[i][j];
                }
            }
            if (arr[i][i] < sum) {
                return false;
            }
        }
        return true;
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {}
    ColumnVector(Matrix const& matrix) : Matrix(matrix) {
        if (cols != 1) {
            valid = false;
        }
    }
    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    ColumnVector& operator= (ColumnVector const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;

    double dot(ColumnVector& other) {
        double res = 0.0;
        if (rows != other.rows) {
            cout << "Error in dot product: incompatible dimensions" << endl;
        }
        else {
            for (int i = 0; i < rows; i++) {
                res += arr[i][0] * other.arr[i][0];
            }
            /*if (abs(res) < EPS) {
                res = 0;
            }*/
        }
        return res;
    }
    double norm() {
        double res = sqrt(abs(dot(*this)));
        /*if (abs(res) < EPS) {
            res = 0;
        }*/
        return res;
    }
};

class AugmentedMatrix : public Matrix {
private:
    int lc;
public:
    AugmentedMatrix(int n, int m, int left_cols) : Matrix(n, m), lc(left_cols) {}
    AugmentedMatrix(Matrix const& matrix, int left_cols) : Matrix(matrix), lc(left_cols) {}
    AugmentedMatrix(Matrix& A, ColumnVector& b) : Matrix(A.get_rows(), A.get_cols() + 1) {
        if (A.get_rows() != b.get_rows()) {
            valid = false;
        }
        else {
            lc = A.get_cols();
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < lc; j++) {
                    arr[i][j] = A.arr[i][j];
                    // if (abs(arr[i][j]) < EPS) {
                    //     arr[i][j] = 0;
                    // }
                }
            }
            for (int i = 0; i < rows; i++) {
                arr[i][cols - 1] = b.arr[i][0];
                // if (abs(arr[i][lc]) < EPS) {
                //     arr[i][lc] = 0;
                // }
            }
        }
    }
    AugmentedMatrix(Matrix A, Matrix B) : Matrix(A.get_rows(), A.get_cols() + B.get_cols()) {
        if (A.get_rows() != B.get_rows()) {
            valid = false;
        }
        else {
            lc = A.get_cols();
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < lc; j++) {
                    arr[i][j] = A.arr[i][j];
                    // if (abs(arr[i][j]) < EPS) {
                    //     arr[i][j] = 0;
                    // }
                }
            }
            for (int i = 0; i < rows; i++) {
                for (int j = lc; j < cols; j++) {
                    arr[i][j] = B.arr[i][j - lc];
                    // if (abs(arr[i][j]) < EPS) {
                    //     arr[i][j] = 0;
                    // }
                }
            }
        }
    }

    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    AugmentedMatrix& operator=(AugmentedMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;

    Matrix* left_part() {
        Matrix* L = new Matrix(rows, lc);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < lc; j++) {
                L->arr[i][j] = arr[i][j];
            }
        }
        return L;
    }
    Matrix* right_part() {
        Matrix* R = new Matrix(rows, cols - lc);
        for (int i = 0; i < rows; i++) {
            for (int j = lc; j < cols; j++) {
                R->arr[i][j - lc] = arr[i][j];
            }
        }
        return R;
    }

    // Eliminates the given augmented matrix to RREF
    void eliminate();

    // Solves the system represented by the given augmented matrix by elimination to RREF
    ColumnVector* solve() {
        //cout << "step #0:" << endl;
        //cout << *left_part();
        //cout << *right_part();
        eliminate();
        return (ColumnVector*)right_part();
    }
    // Solves the system represented by the given augmented matrix by Jacobi method
    ColumnVector* solve_jacobi(double e);
    // Solves the system represented by the given augmented matrix by Seidel method
    ColumnVector* solve_seidel(double e);
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            arr[i][i] = 1;
        }
    }
    IdentityMatrix(SquareMatrix const& matrix) : SquareMatrix(matrix) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if ((i == j) && (arr[i][j] != 1) || (i != j) && (arr[i][j] != 0)) {
                    valid = false;
                    return;
                }
            }
        }
    }
    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    IdentityMatrix& operator= (IdentityMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;
};

class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix(int n, int row, int col, Matrix& A) : SquareMatrix(n) {
        double fact;
        if (A.arr[col][col] == 0) {
            fact = 0;
        }
        else {
            fact = -A.arr[row][col] / A.arr[col][col];
        }
        for (int i = 0; i < n; i++) {
            arr[i][i] = 1;
        }
        arr[row][col] = fact;
    }
    EliminationMatrix(SquareMatrix const& matrix) : SquareMatrix(matrix) {
        bool fact = false;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if ((i == j) && (arr[i][j] != 1)) {
                    valid = false;
                    return;
                }
                if (arr[i][j] != 0) {
                    if (fact) {
                        valid = false;
                        return;
                    }
                    fact = true;
                }

            }
        }
    }
    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    EliminationMatrix& operator= (EliminationMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;
};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n, int row1, int row2) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            arr[i][i] = 1;
        }
        vector<double> t = arr[row1];
        arr[row1] = arr[row2];
        arr[row2] = t;
    }
    PermutationMatrix(SquareMatrix const& matrix) : SquareMatrix(matrix) {
        bool p = false;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if (arr[i][j] != 0 && arr[i][j] != 1) {
                    valid = false;
                    return;
                }
                if ((i == j) && (arr[i][j] != 1) || (i != j) && (arr[i][j] != 0)) {
                    if (p) {
                        valid = false;
                        return;
                    }
                    p = true;
                }
            }
        }
    }
    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    PermutationMatrix& operator= (PermutationMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;
};

double SquareMatrix::det() {
    // swapping rows changes the sign of det
    // other operations change nothing
    // eventually res*=det(REF)
    double res = 1.0;

    // Reducing the given matrix to REF
    int step = 1;
    SquareMatrix ref = *this;
    for (int k = 0; k < rows - 1; k++) {
        // Pivot selection using partial pivoting (by absolutely maximum element)
        double mx = ref.arr[k][k]; // pivot
        int mx_ind = k; // pivot row
        for (int i = k; i < rows; i++) {
            double cur = ref.arr[i][k];
            if (abs(cur) > abs(mx)) {
                mx = cur;
                mx_ind = i;
            }
        }
        //Swap the current row and the pivot row so that the
        //pivot element becomes the largest absolute value in its column
        if (mx_ind != k) {
            ref = PermutationMatrix(rows, k, mx_ind) * ref;
            res *= -1;
            //cout << "step #" << step++ << ": permutation" << endl;
            //cout << ref;
        }

        // Eliminating the elements under the pivot
        for (int i = k + 1; i < rows; i++) {
            if (ref.arr[i][k] != 0) {
                ref = EliminationMatrix(rows, i, k, ref) * ref;
                //for (int j = k; j < cols; j++) {
                // if (abs(ref.arr[i][j]) < EPS) {
                //     ref.arr[i][j] = 0;
                // }
                //}
                //cout << "step #" << step++ << ": elimination" << endl;
                //cout << ref;
            }
        }
    }

    // Calculating the determinant of REF:
    // The determinant of an upper(or lower) triangular matrix is the product
    // of the main diagonal entries.
    for (int i = 0; i < rows; i++) {
        res *= ref.arr[i][i];
    }
    // if (abs(res) < EPS) {
    //     res = 0.0;
    // }
    //cout << "result:"<<endl;
    //cout<<res;
    return res;
}

void AugmentedMatrix::eliminate() {
    //cout << "Direct way:" << endl;
    // Reducing the given matrix to REF
    int step = 1;
    for (int k = 0; k < lc; k++) {
        // Pivot selection using partial pivoting (by absolutely maximum element)
        double mx = arr[k][k]; // pivot
        int mx_ind = k; // pivot row
        for (int i = k + 1; i < rows; i++) {
            double cur = arr[i][k];
            if (abs(cur) > abs(mx)) {
                mx = cur;
                mx_ind = i;
            }
        }

        //Swap the current row and the pivot row so that the
        //pivot element becomes the largest absolute value in its column
        if (mx_ind != k) {
            *this = AugmentedMatrix(PermutationMatrix(rows, k, mx_ind) * (*this), lc);
            //cout << "step #" << step++ << ": permutation" << endl;
            //cout << *left_part();
            //cout << *right_part();
            //cout << *this;
        }

        // Eliminating the elements under the pivot
        for (int i = k + 1; i < rows; i++) {
            if (arr[i][k] != 0) {
                *this = AugmentedMatrix(EliminationMatrix(rows, i, k, *this) * (*this), lc);
                //for (int j = k; j < cols; j++) {
                // if (abs(arr[i][j]) < EPS) {
                //     arr[i][j] = 0;
                // }
                //}
                //cout << "step #" << step++ << ": elimination" << endl;
                //cout << *left_part();
                //cout << *right_part();
                //cout << *this;
            }
        }
    }
    //cout << "Way back:" << endl;
    // Reducing REF to RREF
    for (int k = lc - 1; k >= 0; k--) {
        // Pivot selection using partial pivoting (by absolutely maximum element)
        double mx = arr[rows - 2][lc - 1]; // pivot
        int mx_ind = k; // pivot row
        for (int i = rows - 2; i >= 0; i--) {
            double cur = arr[i][k];
            if (abs(cur) > abs(mx)) {
                mx = cur;
                mx_ind = i;
            }
        }

        // Eliminating the elements above the pivot
        for (int i = k - 1; i >= 0; i--) {
            if (arr[i][k] != 0) {
                *this = AugmentedMatrix(EliminationMatrix(rows, i, k, *this) * (*this), lc);
                //for (int j = k; j < cols; j++) {
                // if (abs(arr[i][j]) < EPS) {
                //     arr[i][j] = 0;
                // }
                //}
                //cout << "step #" << step++ << ": elimination" << endl;
                //cout << *left_part();
                //cout << *right_part();
                //cout << *this;
            }
        }
    }

    // Diagonal normalization
    //cout << "Diagonal normalization:" << endl;
    for (int i = 0; i < rows; i++) {
        if (arr[i][i] != 0) {
            double div = arr[i][i];
            for (int j = 0; j < cols; j++) {
                arr[i][j] /= div;
                // if (abs(arr[i][j]) < EPS) {
                //     arr[i][j] = 0;
                // }
            }
        }
    }
    //cout << *left_part();
    //cout << *right_part();
    //cout << *this;
}
Matrix* SquareMatrix::inverse() {
    IdentityMatrix I(rows);
    AugmentedMatrix aug((Matrix)(*this), (SquareMatrix)I);
    //cout << "step #0: Augmented Matrix" << endl;
    //cout << aug;
    aug.eliminate();
    //cout << "result:" << endl;
    //cout << *aug.right_part();
    return aug.right_part();
}

// for polynomial regression
class VandermondeMatrix : public Matrix {
public:
    VandermondeMatrix(int pts_num, int deg) : Matrix(pts_num, deg + 1) {}
    VandermondeMatrix(Matrix const& matrix) : Matrix(matrix) {}
    VandermondeMatrix(ColumnVector& x, int deg) : Matrix(x.get_rows(), deg + 1) {
        for (int r = 0; r < rows; r++) {
            for (int d = 0; d < cols; d++) {
                arr[r][d] = pow(x.arr[r][0], d);
            }
        }
    }
    friend ostream& operator << (ostream& stream, Matrix& c);
    friend istream& operator >> (istream& stream, Matrix& c);
    VandermondeMatrix& operator= (VandermondeMatrix const& A) = default;
    using Matrix::operator+;
    using Matrix::operator-;
    using Matrix::operator*;
    using Matrix::T;
    ColumnVector* fit(ColumnVector& y) {
        //cout << "A:" << endl;
        //cout << *this;
        Matrix t = (*this).T(); // F_T
        Matrix left = t * (*this); // F_T * F
        //cout << "A_T*A:" << endl;
        //cout << left;
        Matrix* left_inv = ((SquareMatrix)left).inverse(); // (F_T * F)^-1
        //cout << "(A_T*A)^-1:" << endl;
        //cout << *left_inv;
        Matrix right = t * y; // F_T * y
        //cout << "A_T*b:" << endl;
        //cout << right;
        ColumnVector* ans = new ColumnVector((*left_inv) * right); // alpha (coefs of the polynom)
        //cout << "x~:" << endl;
        //cout << *ans;
        return ans;
    }
};
