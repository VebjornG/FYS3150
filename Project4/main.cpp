#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>

using namespace std;

void randomSpinMatrix(int **spin_matrix, int L);

int main(int argc, char *argv[])
{
    //int n_spins = 2;
    int L = 2;

    //int **spin_matrix;
    int **spin_matrix = new int *[L]; // allokerte denne her i stede for inni randomSpinMatrix. Da funket det.

    randomSpinMatrix(spin_matrix, L);

    cout << endl;

    //cout << spin_matrix[1][1] << endl;

    for(int i = 0; i < L; i++) {
        for(int j = 0; j < L; j++) {
            //spin_matrix[i][j] = 0;
            cout << setw(3) << spin_matrix[i][j] << " ";
            //cout << i << ", " << j;
            //cout << ": " << spin_matrix[i][j] << endl;

        }
        cout << endl;
    }

    return 0;
}

void randomSpinMatrix(int **spin_matrix, int L) {
    int N = L;
    clock_t time;
    time = clock();
    long int seed = time;
    //long int seed = 248527847;
    srand(seed);

    // Initialize spin_matrix by dynamic memory allocation

    //spin_matrix = new int *[N]; // denne var allokert her, men ble flyttet ut fra funksjonen.
    for (int i = 0; i < N; i++) {
        spin_matrix[i] = new int[N];
    }

    // Initialize all spin_matrix elements to zero
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            spin_matrix[i][j] = 0;
            //cout << setw(3) << spin_matrix[i][j] << " ";
        }
       //cout << endl;
    }

    cout << endl;

    // Random numbers in matrix
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            // (condition) ? (if_true) : (if_false)
            int r = (rand()/((double)RAND_MAX)) > 0.5 ? 1 : -1;
            spin_matrix[i][j] = r;
            cout << setw(3) << spin_matrix[i][j] << " ";
        }
        cout << endl;
    }
}

