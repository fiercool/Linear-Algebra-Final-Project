#include <iostream>
#include "apmatrix.h"
#include "apvector.h"
#include <cmath>
#include <vector>

using namespace std;

void EnterMatrix();
void DisplayMatrix(apmatrix<double> Matrix);

apmatrix<double> LoadMatrix(int n);
apmatrix<double> TransposeMatrix(apmatrix<double> Matrix);
apmatrix<double> AddMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2);
apmatrix<double> MultiplyMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2);
apmatrix<double> RowReduceMatrix(apmatrix<double> Matrix);
apmatrix<double> InvertMatrix(apmatrix<double> Matrix);
apmatrix<double> OLS(apmatrix<double> M, apmatrix<double> y); 
apmatrix<double> QRAlgorithm(apmatrix<double> A); 
apmatrix<double> CholeskyDecomposition(apmatrix<double> A);

apvector<double> CrossProduct(const apvector<double>& u, const apvector<double>& v);


void DistancePointToLine();
void DistanceLineToLine();
void DistancePointToPlane();



apmatrix<double> Matrices[6];  
int enternumber = 0;

int main() {
    int choice, selection, whichmatrix;
    int keepgoing = 0;
    apmatrix<double> Temp, Temp1, Result;
    apmatrix<double> x, y, xmatrix;

    cout << "Welcome to Tanish's Magical Matrix Program!\n\n";

    while (keepgoing == 0) {
        cout << "What would you like to do? \n\n";
        cout << "1. Enter a Matrix\n";
        cout << "2. Display a Matrix\n";
        cout << "3. Transpose a Matrix\n";
        cout << "4. Add two Matrices\n";
        cout << "5. Multiply two Matrices\n";
        cout << "6. Row Reduce a Matrix\n";
        cout << "7. Invert a Matrix\n";
        cout << "8. OLS Regression\n";
        cout << "9. QR Algorithm\n";
        cout << "10. Estimate solution to a system\n";
        cout << "11. CholeskyDecomposition\n";
        cout << "12. Shortest Distance between a point and a line, line/line, point/plane.\n";

        cout << "What do you choose? (type a number 1-12) ==> ";
        cin >> choice;

        switch (choice) {
            case 1:
                EnterMatrix();
                break;
            case 2:
                cout << "\nWhich matrix would you like to display? (choose Matrix 1-6)==> ";
                cin >> whichmatrix;
                Temp = LoadMatrix(whichmatrix);
                DisplayMatrix(Temp);
                break;
            case 3:
                cout << "\nWhich matrix would you like to Transpose? (choose Matrix 1-6)==> ";
                cin >> whichmatrix;
                Temp = LoadMatrix(whichmatrix);
                Temp1 = TransposeMatrix(Temp);
                DisplayMatrix(Temp1);
                break; 
            case 4: 
                int whichmatrix1, whichmatrix2; 
                cout << "\nWhat is the first matrix you want to add? (choose Matrix 1-6)==> ";
                cin >> whichmatrix1;
                Temp = LoadMatrix(whichmatrix1);
                cout << "\nWhat is the second matrix you want to add (Choose Matrix 1-6)==> ";
                cin >> whichmatrix2;
                Temp1 = LoadMatrix(whichmatrix2);
                Result = AddMatrices(Temp, Temp1);
                DisplayMatrix(Result);
                break;
            case 5: 
                
                cout << "\nWhat is the first matrix you want to Multiply? (choose Matrix 1-6)==> ";
                cin >> whichmatrix1;
                Temp = LoadMatrix(whichmatrix1);
                cout << "\nWhat is the second matrix you want to Multiply? (Choose Matrix 1-6)==> ";
                cin >> whichmatrix2;
                Temp1 = LoadMatrix(whichmatrix2);
                Result = MultiplyMatrices(Temp, Temp1);
                DisplayMatrix(Result);
                break;

            case 6:
                cout << "What matrix would you like to row reduce to RREF (choose Matrix 1-6)==> ";
                cin >> whichmatrix1;
                Temp = LoadMatrix(whichmatrix1);
                Result = RowReduceMatrix(Temp);
                DisplayMatrix(Result);
                break;
            
            case 7: 
                cout << "What matrix would you like to find the inverse of (choose Matrix 1-6)==> ";
                cin >> whichmatrix1;
                Temp = LoadMatrix(whichmatrix1);
                Result = InvertMatrix(Temp);
                DisplayMatrix(Result);
                break;
            case 8: {
                int degree, datapoints;
                cout << "\nWhat is the degree of your polynomial? ==> ";
                cin >> degree;

                cout << "\nHow many data points do you want to enter? ==> ";
                cin >> datapoints;

                apmatrix<double> coeffMatrix; 
                apmatrix<double> yMatrix;

                coeffMatrix.resize(datapoints, degree + 1);
                yMatrix.resize(datapoints, 1);

                for (int i = 0; i < datapoints; i++) {
                    double x, y;
                    cout << "Enter x value for point " << i + 1 << ": ";
                    cin >> x;
                    cout << "Enter y value for point " << i + 1 << ": ";
                    cin >> y;

                    yMatrix[i][0] = y;
                    double x_pow = 1.0;
                    for (int j = 0; j <= degree; j++) {
                        coeffMatrix[i][j] = x_pow;
                        x_pow *= x;
                    }
                }

                apmatrix<double> result = OLS(coeffMatrix, yMatrix);

                cout << "Your equation is: y = ";
                for (int i = 0; i <= degree; i++) {
                    double coeff = result[i][0];

                    if (fabs(coeff) < 0.000001) {
                        coeff = 0.0;
                    }

                    if (i == 0) {
                        cout << coeff;
                    } else {
                        cout << " + " << coeff << "x^" << i;
                    }
                }
                cout << endl;

                break;
            }

            case 9: {
                cout << "Which matrix do you want to find eigenvalues of (choose Matrix 1-6)? ";
                cin >> whichmatrix;
                Temp = LoadMatrix(whichmatrix);
                QRAlgorithm(Temp);
                break;
            }

            case 10: { 
                int whichmatrix1, whichmatrix2;
                cout << "\nWhat is the matrix of coefficients (choose Matrix 1-6)? ==> ";
                cin >> whichmatrix1;
                apmatrix<double> A = LoadMatrix(whichmatrix1);

                cout << "\nWhat is the matrix of constants (choose Matrix 1-6)? ==> ";
                cin >> whichmatrix2;
                apmatrix<double> b =LoadMatrix(whichmatrix2);

                if (A.numcols() != b.numrows() || A.numrows() != b.numcols()) {
                    cout << "Coefficient matrix must be nxm and b must be nx1 \n";
                    break;
                }

                apmatrix<double> solution = OLS(A, b);
                cout << "Estimated solution:\n";
                DisplayMatrix(solution);
                break;
            }
            case 11: {
                cout << "Which matrix do you want to perform Cholesky Decomposition on (choose Matrix 1-6)? ";
                cin >> whichmatrix;
                Temp = LoadMatrix(whichmatrix);
                if (Temp.numrows() != Temp.numcols()) {
                    cout << "Cholesky Decomposition requires a square matrix.\n";
                    break;
                }
                apmatrix<double> L = CholeskyDecomposition(Temp);
                cout << "Lower triangular matrix L:\n";
                DisplayMatrix(L);
                break;
            }
                
            case 12: {
            cout << "What would you like to do? \n\n";
            cout << "1. Find the Shortest Distance between a point and a line\n";
            cout << "2. Find the Shortest Distance between a line and a line\n";
            cout << "3. Find the Shortest Distance between a point and a plane\n";
            cout << "What do you choose? (type a number 1-3) ==> ";
            cin >> selection;

            switch (selection) {
                case 1:
                    DistancePointToLine();
                    break;
                case 2:
                    DistanceLineToLine();
                    break;
                case 3:
                    DistancePointToPlane();
                    break;
                default:
                    cout << "Invalid selection.\n";
            }

            break;   

        }  

        default:
            cout << "\nNot an Option";
            continue;
    }

        cout << "\n\nDo you want to keep going? (0=yes, 1=no) ==> ";
        cin >> keepgoing;
    }

    cout << "\n\nSee you next time.";
}



void EnterMatrix() {
    int rows, columns, i, j;
    apmatrix<double> Temp;

    cout << "How many rows would you like your matrix to have? ==> ";
    cin >> rows;
    cout << "How many columns would you like your matrix to have? ==> ";
    cin >> columns;
    Temp.resize(rows, columns);

    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            cout << "What is your row " << i + 1 << " column " << j + 1 << " entry? ==> ";
            cin >> Temp[i][j];
        }
    }

    if (enternumber < 6) {
        Matrices[enternumber] = Temp;
        cout << "\nYour matrix has been stored as matrix " << enternumber + 1 << ". Here is your matrix.\n";
        DisplayMatrix(Temp);
        enternumber++;
    } else {
        cout << "\nYou've already entered 6 matrices. Cannot store more.\n";
    }
    
}

void DisplayMatrix(apmatrix<double> Matrix) {
    int rows = Matrix.numrows();
    int cols = Matrix.numcols();

    if (rows == 0 || cols == 0) {
        cout << "[Empty Matrix]\n";
        return;
    }


    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val = Matrix[i][j];
            if (fabs(val) < .000001){
                val = 0.0;
            } // -0 fix
            cout << val << " ";
        }
        cout << endl;
    }
}

apmatrix<double> LoadMatrix(int n) {
    apmatrix<double> Temp1;
    if (n >= 1 && n <= 6) {
        Temp1 = Matrices[n - 1];
    } else {
        cout << "Invalid matrix number. Returning an empty matrix.\n";
        Temp1.resize(0, 0);
    }
    return Temp1;
}

apmatrix<double> TransposeMatrix(apmatrix<double> Matrix) {
    int rows = Matrix.numrows();
    int cols = Matrix.numcols();
    apmatrix<double> Result(cols, rows);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Result[j][i] = Matrix[i][j];
        }
    }

    return Result;
}

apmatrix<double> AddMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2) {
    int rows1 = Temp1.numrows();
    int cols1 = Temp1.numcols();
    int rows2 = Temp2.numrows();
    int cols2 = Temp2.numcols();
    apmatrix<double> Result;

    if (rows1 == rows2 && cols1 == cols2) {
        Result.resize(rows1, cols1);
        for (int i = 0; i < rows1; i++) {
            for (int j = 0; j < cols1; j++) {
                Result[i][j] = Temp1[i][j] + Temp2[i][j];
            }
        }
    } else {
        cout << "Matrices are not the same size and cannot be added.\n";
        Result.resize(0, 0);  
    }

    return Result;
}

apmatrix<double> MultiplyMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2){
    int rows1 = Temp1.numrows();
    int cols1 = Temp1.numcols();
    int rows2 = Temp2.numrows();
    int cols2 = Temp2.numcols();
    apmatrix<double> Result;

    if (cols1 == rows2){ 
        Result.resize(rows1, cols2);
        for (int i = 0; i < rows1; i++){
            for (int j = 0; j < cols2; j++){
                Result[i][j] = 0;
                for (int k = 0; k < cols1; k++){
                    Result[i][j] += Temp1[i][k] * Temp2[k][j];
                }
            }
        }
    } else { 
        cout << "Matrices are not able to be multiplied. \n";
        Result.resize(0,0);
    }

    return Result;
}



apmatrix<double> RowReduceMatrix(apmatrix<double> Matrix) {
    int rows = Matrix.numrows();
    int cols = Matrix.numcols();
    apmatrix<double> Result = Matrix;

    int pivot_col = 0;
    for (int pivot_row = 0; pivot_row < rows && pivot_col < cols;) {
        //find row with highest abs value in pivot column
        int max_row = pivot_row;
        double max_val = fabs(Result[pivot_row][pivot_col]);

        for (int k = pivot_row + 1; k < rows; k++) {
            double val = fabs(Result[k][pivot_col]);
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }

        //if abs of highest pivot = 0, skip column
        if (max_val < 0.0000001) {
            pivot_col++;
            continue;
        }

        //switch current row with the row with the largest pivot
        if (max_row != pivot_row) {
            for (int j = 0; j < cols; j++) {
                double temp = Result[pivot_row][j];
                Result[pivot_row][j] = Result[max_row][j];
                Result[max_row][j] = temp;
            }
        }

        //normalize pivot
        double pivot = Result[pivot_row][pivot_col];
        for (int j = 0; j < cols; j++) {
            Result[pivot_row][j] /= pivot;
        }

        //subtract other rows by pivot * scalar (eliminate below)
        for (int i = pivot_row + 1; i < rows; i++) {
            double scalar = Result[i][pivot_col];
            for (int j = 0; j < cols; j++) {
                Result[i][j] -= scalar * Result[pivot_row][j];
            }
        }

        pivot_row++;
        pivot_col++;
    }

    //  backward elimination (REF to RREF)
    for (int i = rows - 1; i >= 0; --i) {
        // find leading 1 in row i
        int lead_col = -1;
        for (int j = 0; j < cols; ++j) {
            if (fabs(Result[i][j] - 1.0) < 0.0000001) {
                lead_col = j;
                break;
            }
        }
        if (lead_col == -1) continue;

        // eliminate above pivot
        for (int k = i - 1; k >= 0; --k) {
            double scale = Result[k][lead_col];
            for (int j = 0; j < cols; ++j) {
                Result[k][j] -= scale * Result[i][j];
            }
        }
    }

    return Result;
}

apmatrix<double> InvertMatrix(apmatrix<double> Matrix) {
    int n = Matrix.numrows();
    apmatrix<double> Aug;
    Aug.resize(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Aug[i][j] = Matrix[i][j];
        }
        for (int j = 0; j < n; j++) {
            if (i == j) {
                Aug[i][n + j] = 1;
            } else {
                Aug[i][n + j] = 0;
            }
        }
    }

    apmatrix<double> R = RowReduceMatrix(Aug);
    apmatrix<double> Inv;
    Inv.resize(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Inv[i][j] = R[i][n + j];
        }
    }
    return Inv;
}

apmatrix<double> OLS(apmatrix<double> M, apmatrix<double> y) {
    //(M^T M)^{-1} M^T y
    apmatrix<double> Mt = TransposeMatrix(M);
    apmatrix<double> MtM = MultiplyMatrices(Mt, M);
    apmatrix<double> MtM_inv = InvertMatrix(MtM);
    apmatrix<double> Mt_y = MultiplyMatrices(Mt, y);
    apmatrix<double> Z = MultiplyMatrices(MtM_inv, Mt_y);
    return Z;
}

apmatrix<double> QRAlgorithm(apmatrix<double> A) {
    int n = A.numrows();
    int maxIter = 100;
    double tol = 0.00001;

    apmatrix<double> Ak = A;

    for (int iter = 0; iter < maxIter; ++iter) {
        //gram
        apmatrix<double> Q(n, n);
        apmatrix<double> R(n, n);

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i)
                Q[i][j] = Ak[i][j];

            for (int k = 0; k < j; ++k) {
                double dot = 0.0;
                for (int i = 0; i < n; ++i)
                    dot += Q[i][k] * Ak[i][j];

                R[k][j] = dot;
                for (int i = 0; i < n; ++i)
                    Q[i][j] -= dot * Q[i][k];
            }

            double norm = 0.0;
            for (int i = 0; i < n; ++i)
                norm += Q[i][j] * Q[i][j];

            R[j][j] = sqrt(norm);
            for (int i = 0; i < n; ++i)
                Q[i][j] /= R[j][j];
        }

        Ak = MultiplyMatrices(R, Q);
    }


    cout << "Estimated eigenvalues:\n";
    for (int i = 0; i < n; ++i) {
        cout << Ak[i][i] << endl;
    }
    return Ak;
}

apmatrix<double> CholeskyDecomposition(apmatrix<double> A) {
    int n = A.numrows();
    apmatrix<double> L(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                //diagonal L[i][i]
                L[i][j] = sqrt(A[i][i] - sum);
            } else {
                // off-diagonal L[i][j]
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }

    return L;
}


void DistancePointToLine() {
    //honestly not sure if this is supposed to be in 3 dimentions, all the examples of questions in the textbook were 3d so i figured id just do it like that.
   
    double Ax, Ay, Az, vx, vy, vz, Px, Py, Pz;

    cout << "What is the x-coordinate of point A on the line? ==> ";
    cin >> Ax;
    cout << "What is the y-coordinate of point A on the line? ==> ";
    cin >> Ay;
    cout << "What is the z-coordinate of point A on the line? ==> ";
    cin >> Az;

    cout << "What is the x-component of the direction vector v of the line? ==> ";
    cin >> vx;
    cout << "What is the y-component of the direction vector v of the line? ==> ";
    cin >> vy;
    cout << "What is the z-component of the direction vector v of the line? ==> ";
    cin >> vz;

    cout << "What is the x-coordinate of point P? ==> ";
    cin >> Px;
    cout << "What is the y-coordinate of point P? ==> ";
    cin >> Py;
    cout << "What is the z-coordinate of point P? ==> ";
    cin >> Pz;

    //AP = P - A
    double APx = Px - Ax;
    double APy = Py - Ay;
    double APz = Pz - Az;


    //proj_v(AP) = [(AP ⋅ v) / (v ⋅ v)] * v
    double dot_AP_v = APx * vx + APy * vy + APz * vz;
    double v_dot_v = vx * vx + vy * vy + vz * vz;

    
    double scalar = dot_AP_v / v_dot_v;
    double projx = scalar * vx;
    double projy = scalar * vy;
    double projz = scalar * vz;

    //  AP - proj_v(AP)
    double perp_x = APx - projx;
    double perp_y = APy - projy;
    double perp_z = APz - projz;

    // mag of perp
    double distance = sqrt(perp_x * perp_x + perp_y * perp_y + perp_z * perp_z);

    cout << "The shortest distance from point P to the line is: " << distance << endl;
}

void DistanceLineToLine() {
    double A1x, A1y, A1z, v1x, v1y, v1z, A2x, A2y, A2z, v2x, v2y, v2z;

    cout << "Enter the x-coordinate of point A1 on line 1: ==> ";
    cin >> A1x;
    cout << "Enter the y-coordinate of point A1 on line 1: ==> ";
    cin >> A1y;
    cout << "Enter the z-coordinate of point A1 on line 1: ==> ";
    cin >> A1z;

    cout << "Enter the x-component of direction vector v1: ==> ";
    cin >> v1x;
    cout << "Enter the y-component of direction vector v1: ==> ";
    cin >> v1y;
    cout << "Enter the z-component of direction vector v1: ==> ";
    cin >> v1z;

    cout << "Enter the x-coordinate of point A2 on line 2: ==> ";
    cin >> A2x;
    cout << "Enter the y-coordinate of point A2 on line 2: ==> ";
    cin >> A2y;
    cout << "Enter the z-coordinate of point A2 on line 2: ==> ";
    cin >> A2z;

    cout << "Enter the x-component of direction vector v2: ==> ";
    cin >> v2x;
    cout << "Enter the y-component of direction vector v2: ==> ";
    cin >> v2y;
    cout << "Enter the z-component of direction vector v2: ==> ";
    cin >> v2z;


    double Cx = v1y * v2z - v1z * v2y;
    double Cy = v1z * v2x - v1x * v2z;
    double Cz = v1x * v2y - v1y * v2x;

    double crossMag = sqrt(Cx * Cx + Cy * Cy + Cz * Cz);

  -
    if (fabs(crossMag) < 0.00001) {
        
        cout << "Lines are parallel or intersect. \n";
        return;
    }

    double diffx = A2x - A1x;
    double diffy = A2y - A1y;
    double diffz = A2z - A1z;

    double dot = diffx * Cx + diffy * Cy + diffz * Cz;

    double distance = fabs(dot) / crossMag;

    cout << "The shortest distance between the two lines is: " << distance << endl;
}



void DistancePointToPlane() {
    double A, B, C, D, Px, Py, Pz;

    cout << "What is the A value in the equation Ax + By + Cz + D = 0? ==> ";
    cin >> A;
    cout << "What is the B value in the equation Ax + By + Cz + D = 0? ==> ";
    cin >> B;
    cout << "What is the C value in the equation Ax + By + Cz + D = 0? ==> ";
    cin >> C;
    cout << "What is the D value in the equation Ax + By + Cz + D = 0? ==> ";
    cin >> D;

    double norm_n_sq = A*A + B*B + C*C;
    double Qx, Qy, Qz;

    if (fabs(C) > 0.000001) {
        Qx = 0.0;
        Qy = 0.0;
        Qz = D / C;
    } else if (fabs(B) > 0.000001) {
        Qx = 0.0;
        Qz = 0.0;
        Qy = D / B;
    } else {
        Qy = 0.0;
        Qz = 0.0;
        Qx = D / A;
    }
    
    cout << "Enter the x-coordinate of point P: ==> ";
    cin >> Px;
    cout << "Enter the y-coordinate of point P: ==> ";
    cin >> Py;
    cout << "Enter the z-coordinate of point P: ==> ";
    cin >> Pz;

   
    double dot_PQ_n = Px * A + Py * B + Pz * C + D;
    double norm_n = sqrt(norm_n_sq);
    double distance = fabs(dot_PQ_n) / norm_n;
    cout << "The shortest distance from point P to the plane is: " << distance << endl;
}