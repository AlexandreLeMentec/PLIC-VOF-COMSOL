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

            float nuo = x/ dx;
            float mx = fabs(nx)*dx/(fabs(nz)*dz + fabs(nx)*dx);
            float mz = fabs(nz)*dz/(fabs(nz)*dz + fabs(nx)*dx);

            
            float alpha = 0.0;
            
            // float phimax 

            // #################################### CONDITIONNER ALPHA ##################################
            if (nx*nz != 0.0){
                //
                float nx_cond = 0.5 + copysign(0.5,nx);                      // condition nx > = 0   
                float mxz_cond = 0.5 + copysign(0.5,(mz-mx));                // condition mx < = mz
                float phi_tr = mxz_cond*nx_cond*(mx/(2.*mz)*(nuo+1.0/3.0)/(nuo+0.5)) + (1.0-mxz_cond)*nx_cond*(mz/(2.0*mx)*(nuo+mz/(3.0*mx))/(nuo+0.5)) + (1.0-nx_cond)*mxz_cond*(mx/(2.0*mz)*(nuo+2.0/3.0)/(nuo+0.5)) + (1.0-nx_cond)*(1.0-mxz_cond)*(mz/(2.0*mx)*(nuo+1.-mz/(3.0*mx))/(nuo+0.5));
                float phi_crit = nx_cond*(2.0*mx/(3.0*mz)*(powf(nuo,3.0))/(nuo+0.5)) + (1.0-nx_cond)*(2.0*mx/(3.0*mz)*(powf((nuo+1.0),3.0))/(nuo+0.5)) ;
                //

                float M = 3.0*powf(mx,2.0)*mz*(nuo+0.5);
                float phio = asin(1.0 - 2.0*C/phi_crit);
                if (C < phi_tr)             
                {
                    if (C < phi_crit || nx < 0.0)                    // equation (1 et 6) et (4 et 9)
                    {
                    float alpha_a = (-1.0 + 2.0*cos(phio/3.0 + acos(-1.0)/6.0 )) ;
                    float alpha_b = (-1.0 + 2.0*sin(phio/3.0)) ;
                    alpha = copysign(1.0,nx)*powf(M*phi_crit*0.5,1.0/3.0) * (nx_cond*alpha_a + (1.0 - nx_cond)*alpha_b)  ;   

                    }
                    else                 // equation 2 et 7
                    {
                    float alpha_c = powf(2.*C - phi_crit + 2.*sqrt(C*(C - phi_crit)),1.0/3.0) + powf(2.*C - phi_crit - 2*sqrt(C*C - phi_crit),1.0/3.0);
                    alpha = ( alpha_c + copysign(1.0,(mz-mx))*powf(phi_crit,1.0/3.0))*powf(M/2.0,1.0/3.0) ;
                    }
                }
                else
                {
                    if((mz-mx)>=0.0)
                    {
                        // equation 3 et 5
                        alpha = mz * C + mx/2.0 * (nuo + 1.0/3.0 + nx_cond*1.0/3.0)/(nuo+0.5) ; // Ajout d'une condition sur le 1/3 -> 2/3
                    }
                    else
                    {
                        if(nx<0.0)
                        {
                            // equation 10
                            alpha = 0.5 * (mz + 2.0*mx*(nuo+1.0) - powf(4.0*powf(mx,2.0)*powf(nuo+1.0,2.0)- 8 * powf(mx,2.0)*(nuo+0.5)*C - 1.0/3.0*powf(mz,2.0),0.5));
                        }
                        else
                        {
                            // equation 8 
                            alpha = 0.5 * (mz - 2.0*mx*nuo + powf(4.0*powf(mx,2.0)*powf(nuo,2.0) + 8.0 * powf(mx,2.0)*(nuo+0.5)*C - 1.0/3.0*powf(mz,2.0),0.5));
                        }
                    }
                }
            }
            else {
                alpha = C;
                if (C < 0.99 && C > 0.01){
                    outReal[i] = dx*nx+dz*nz;
                }
                else{
                    outReal[i] = 0.0;     
                }
             return 1;   
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
            for (int j = 0; j < 4; j++) {
                int mask = (candidates[j].x >= 0.0f && candidates[j].x <= 1.0f) &
                        (candidates[j].z >= 0.0f && candidates[j].z <= 1.0f);
                results[count] = candidates[j];  
                count += mask;
                if (count == 2) break;
            }
            // results[0] => {sigx1, sigz1}, results[0].x = sigx1
            
            tuple results_coord[2];
            if (count == 2 && C < 0.9 && C > 0.1){
                for (int j = 0; j < 2; j++) {
                    results_coord[j].x = dx * results[j].x + x     ;
                    results_coord[j].z = dz * (1.0- results[j].z)+z*dz  ;
                }  
                float Si = sqrt((results_coord[1].x - results_coord[0].x)*(results_coord[1].x - results_coord[0].x) + (results_coord[1].z - results_coord[0].z)*(results_coord[1].z - results_coord[0].z));
                outReal[i] = Si;
            }
            else {
                outReal[i] = 0.0;
            }
        }
        return 1;
    }
    else {
        error = "Unknown function";
        return 0;
    }
  }