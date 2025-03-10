#include <stdio.h>

#define NR 88
#define NZ 88

void fill_matrix(int A[NR][NZ]) {
    // Initialize with zeros
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            A[i][j] = 0;
        }
    }
    
    // Initialization
    A[0][0] = 1;
    A[NR-1][NZ-1] = NR * NZ;
    
    for (int i = 0; i < NR - 1; i++) {
        A[i+1][0] = A[i][0] + i + 2;
        A[NR-i-2][NZ-1] = A[NR-i-1][NZ-1] - (i + 2);
    }
    
    for (int i = 0; i < NR - 1; i++) { 
        for (int j = 1; j < NZ - i; j++) { 
            A[i][j] = A[i][j-1] + (i + j);
            A[NR-i-1][NZ-j-1] = A[NR-i-1][NZ-j] - (i + j);
        }
    }
}

void save_matrix_to_file(int A[NR][NZ], const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    
    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            fprintf(file, "%d ", A[i][j]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
}

int main() {
    int A[NR][NZ];
    fill_matrix(A);
    save_matrix_to_file(A, "matrix.txt");
    return 0;
}
