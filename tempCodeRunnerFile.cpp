#include <iostream>
#include "apmatrix.h"
#include <cmath>

using namespace std;

//define functions
void EnterMatrix();
void DisplayMatrix(apmatrix<double> Matrix);
apmatrix<double> LoadMatrix(int n);
apmatrix<double> TransposeMatrix(apmatrix<double> Matrix);
apmatrix<double> AddMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2);
apmatrix<double> MultiplyMatrices(apmatrix<double> Temp1, apmatrix<double> Temp2);
apmatrix<double> RowReduceMatrix(apmatrix<double> Matrix);


apmatrix<double> Matrices[6];  
int enternumber = 0;

int main() {
    int choice, whichmatrix;
    int keepgoing = 0;
    apmatrix<double> Temp, Temp1, Result;

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
        cout << "What do you choose? (type a number 1-8) ==> ";
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
                cout << "
What matrix would you like to row reduce to RREF (choose Matrix 1-6)==> ";
                cin >> whichmatrix1;
                Temp = LoadMatrix(whichmatrix1);
                Result = RowReduceMatrix(Temp);
                DisplayMatrix(Result);
                break;

            default:
                cout << "
Not an Option";
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
            cout << Matrix[i][j] << " ";
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
    for (int pivot_row = 0; pivot_row < rows && pivot_col < cols; pivot_row++, pivot_col++) {
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
        }

        //if abs of highest pivot = 0 (w/ floating point), skip column
        if (max_val < 0.0000001){ 
            pivot_col++;
            continue;
        }
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
            Result[pivot_row][j] = Result[pivot_row][j] / pivot;
        }

        //subtract other rows by pivot * scalar 
        for (int i = pivot_row + 1; i < rows; i++) {
            double scalar = Result[i][pivot_col];
            for (int j = 0; j < cols; j++) {
                Result[i][j] -= scalar * Result[pivot_row][j];
            }
        }
    }

    // added this block to perform back-substitution and get RREF
    for (int pivot_row = rows - 1; pivot_row >= 0; pivot_row--) {
        int pivot_col = -1;
        for (int j = 0; j < cols; j++) {
            if (fabs(Result[pivot_row][j]) > 0.0000001) {
                pivot_col = j;
                break;
            }
        }
        if (pivot_col == -1) continue;

        for (int i = pivot_row - 1; i >= 0; i--) {
            double scalar = Result[i][pivot_col];
            for (int j = 0; j < cols; j++) {
                Result[i][j] -= scalar * Result[pivot_row][j];
            }
        }
    }

    return Result;
}
