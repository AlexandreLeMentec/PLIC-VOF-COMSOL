// EXAMPLE EXTERNAL FUNCTION 

#include <math.h>
#include <stdlib.h>
#include <string.h>
 
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

static const char *error = NULL;

void fill_matrix(int Nr, int Nz, int A[Nr][Nz]) { // generation de la matrice (il faut que Nr = Nz en l'instant T)
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nz; j++) {
            A[i][j] = 0;
        }
    }
    
    A[0][0] = 1;
    A[Nr-1][Nz-1] = Nr * Nz;
    
    for (int i = 0; i < Nr - 1; i++) {
        A[i+1][0] = A[i][0] + i + 2;
        A[Nr-i-2][Nz-1] = A[Nr-i-1][Nz-1] - (i + 2);
    }
    
    for (int i = 0; i < Nr - 1; i++) { 
        for (int j = 1; j < Nz - i; j++) { 
            A[i][j] = A[i][j-1] + (i + j);
            A[Nr-i-1][Nz-j-1] = A[Nr-i-1][Nz-j] - (i + j);
        }
    }
    return 1;
}

void fill_dr(int Ntot, int Nr, int Nz, int A[Nr][Nz], double dr[Ntot], const double **inReal){
    int ilist[Nz];
    double hlist[Nz];

    for (int i = 0; i < Nr - 1; i++) {
        for (int j = 0; i < Nz - 1; i++) {
            ilist[j] = A[i][j];
            hlist[j] = inReal[A[i][j]][0];
        }
        int hmin = 1000000.0;
        for (int j = 0; i < Nz - 1; i++) {
            if (hlist[j]<hmin){
                hmin = hlist[j];
            }
        }
        for (int j = 0; i < Nz - 1; i++) {
            dr[ilist[j]] = hmin;
        }
    }
    return 1;
}

int fill_dz(int Ntot, int Nr, int Nz, int A[Nr][Nz], double dz[Ntot], const double **inReal){
    int ilist[Nr];
    double hlist[Nr];

    for (int i = 0; i < Nz - 1; i++) {
        for (int j = 0; i < Nr - 1; i++) {
            ilist[j] = A[i][j];
            hlist[j] = inReal[A[i][j]][0];
        }
        int hmin = 1000000.0;
        for (int j = 0; i < Nr - 1; i++) {
            if (hlist[j]<hmin){
                hmin = hlist[j];
            }
        }
        for (int j = 0; i < Nr - 1; i++) {
            dz[ilist[j]] = hmin;
        }
    }
    return 1;
}

EXPORT int init(const char *str) {
    return 1;
}
   
EXPORT const char * getLastError() {
    return error;
}

EXPORT int eval(const char *func,
                                int nArgs,
                                const double **inReal,
                                const double **inImag,
                                int blockSize,
                                double *outReal,
                                double *outImag) {
    int i;
    int Nr = floor(inReal[0][1]);
    int Nz = floor(inReal[0][2]);
    int reconst_mat[Nr][Nz];
    fill_matrix(Nr,Nz,reconst_mat);

    if (strcmp("dr", func) == 0) {
        if (nArgs != 3) {
        error = "Three argument expected";
        return 0;
        }
        double dr_list[blockSize];
        fill_dr(blockSize,Nr,Nz,reconst_mat,dr_list,inReal);
        for (i = 0; i < blockSize; i++) {
            outReal[i] = dr_list[i];
        }
        return 1;
    }
    else if (strcmp("dz", func) == 0){
        if (nArgs != 3) {
        error = "Three argument expected";
        return 0;
        }
        double dz_list[blockSize];
        fill_dz(blockSize,Nr,Nz,reconst_mat,dz_list,inReal);
        for (i = 0; i < blockSize; i++) {
        outReal[i] = dz_list[i];
        }
        return 1;
    }
    else {
        error = "Unknown function";
        return 0;
    }
}




