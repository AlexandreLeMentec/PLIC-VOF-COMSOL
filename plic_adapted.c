#include <math.h>
#include <stdlib.h>
#include <string.h>
 
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

typedef struct {
    double x;
    double z;
} tuple;

 
static const char *error = NULL;

EXPORT int init(const char *str) {
    return 1;
  }
  
EXPORT const char * getLastError() {
    return error;
  }

EXPORT int eval(const char *func, int nArgs, const double **inReal,
                const double **inImag,int blockSize,double *outReal,
                double *outImag ){
    int i, j;

    if (strcmp("plic_vof", func) == 0) {
        if (nArgs != 7) {
          error = "Seven argument expected";
          return 0;
        }
        for (i = 0; i < blockSize; i++) {
            // Appel variables 
            // float C, float nx, float nz, float dx, float dz, float x, float z

            double C = inReal[0][i];
            double nx = inReal[1][i];
            double nz = inReal[2][i];
            double dx = inReal[3][i];
            double dz = inReal[4][i];
            double x = inReal[5][i];
            double z = inReal[6][i];

            // Calcul variables
            double nuo = x/ dx;
            double mx = fabs(nx)*dx/(fabs(nz)*dz + fabs(nx)*dx);
            double mz = fabs(nz)*dz/(fabs(nz)*dz + fabs(nx)*dx);

            //
            double nx_cond = 0.5 + copysign(0.5,nx);                      // condition nx > = 0   
            double mxz_cond = 0.5 + copysign(0.5,(mz-mx));                // condition mx < = mz
            double phi_tr = mxz_cond*nx_cond*(mx/(2.*mz)*(nuo+1./3.)/(nuo+1./2.)) + (1.-mxz_cond)*nx_cond*(mz/(2.*mx)*(nuo+mz/(3.*mx))/(nuo+1./2.)) + (1.-nx_cond)*mxz_cond*(mx/(2.*mz)*(nuo+2./3.)/(nuo+1./2)) + (1-nx_cond)*(1-mxz_cond)*(mz/(2*mx)*(nuo+1.-mz/(3*mx))/(nuo+1/2));
            double phi_crit = nx_cond*(2.*mx/(3.*mz)*(powf(nuo,3.))/(nuo+1./2.)) + (1.-nx_cond)*(2.*mx/(3.*mz)*(powf((nuo+1.),3.))/(nuo+1./2.)) ;
            //

            double M = 3*powf(mx,2.)*mz*(nuo+0.5);
            double phio = asin(1. - 2.*C/phi_crit);
            double alpha = 0.        ;
            // double phimax 

            // #################################### CONDITIONNER ALPHA ##################################
            
            if (C < phi_tr)             
            {
                if (C < phi_crit || nx < 0.)                    // equation (1 et 6) et (4 et 9)
                {
                double alpha_a = (-1 + 2*cos(phio/3 + acos(-1)/6 )) ;
                double alpha_b = (-1 + 2*sin(phio/3)) ;
                double alpha = powf(M*phi_crit/2,1/3) * (nx_cond*alpha_a + (1. - nx_cond)*alpha_b)  ;   
                }
                else                 // equation 2 et 7
                {
                double alpha_c = powf(2.*C - phi_crit + 2.*sqrt(C*(C - phi_crit)),1/3) + powf(2.*C - phi_crit - 2*sqrt(C*C - phi_crit),1./3.);
                double alpha = ( alpha_c + copysign(1.0,(mz-mx))*pow(phi_crit,1/3))*pow(M/2,1/3) ;
                }
            }
            else
            {
                if((mz-mx)>0.)
                {
                    // equation 3 et 5
                    double alpha = mz * C + mx/2. * (nuo + 1/3 + nx_cond*1.0/3.0)/(nuo+1/2) ; // Ajout d'une condition sur le 1/3 -> 2/3
                }
                else
                {
                    if(nx<0.)
                    {
                        // equation 10
                        double alpha = 1/2 * (mz - 2*mx*(nuo+1) - sqrt(4*powf(mx,2)*pow(nuo+1,2)- 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)));
                    }
                    else
                    {
                        // equation 8 
                        double alpha = 1/2 * (mz - 2*mx*nuo + sqrt(4*powf(mx,2)*powf(nuo,2) + 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)));
                    }
                }
            }

            // #################################### calcul de la surface d'interface ##################################
            tuple candidates[4] = {
                {(alpha - mz) / mx, 1.0},  // sigz1
                {1.0, (alpha - mx) / mz},  // sigx1
                {alpha / mx, 0.0},         // sigz0
                {0.0, alpha / mz}          // sigx0
                };
            tuple results[2];
            int count = 0;
            for (int i = 0; i < 4; i++) {
                int mask = (candidates[i].x >= 0.0f && candidates[i].x <= 1.0f) &
                        (candidates[i].z >= 0.0f && candidates[i].z <= 1.0f);
                results[count] = candidates[i];  
                count += mask;
                if (count == 2) break;
            }
            // results[0] => {sigx1, sigz1}, results[0].x = sigx1
            
            tuple results_coord[2];
            for (int i = 0; i < 2; i++) {
                results_coord[i].x = dx * results[0].z + x     ;
                results_coord[i].z = dz * (1- results[0].x)+z  ;
            }  
            outReal[i] = sqrt((results_coord[1].x - results_coord[0].x)*(results_coord[1].x - results_coord[0].x) + (results_coord[1].z - results_coord[0].z)*(results_coord[1].z - results_coord[0].z));
        }
        return 1;
    }
    else {
        error = "Unknown function";
        return 0;
    }
  }