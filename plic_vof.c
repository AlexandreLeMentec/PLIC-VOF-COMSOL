#include <math.h>

typedef struct {
    float x;
    float z;
} tuple;

float plic_cyl(float C, float nx, float nz, float dx, float dz, float x, float z)
    {
        float nuo = x/ dx;
        float mx = fabs(nx)*dx/(fabs(nz)*dz + fabs(nx)*dx);
        float mz = fabs(nz)*dz/(fabs(nz)*dz + fabs(nx)*dx);

        //
        float nx_cond = 0.5 + copysign(0.5,nx);                      // condition nx > = 0   
        float mxz_cond = 0.5 + copysign(0.5,(mz-mx));                // condition mx < = mz
        float phi_tr = mxz_cond*nx_cond*(mx/(2.*mz)*(nuo+1./3.)/(nuo+1./2.)) + (1.-mxz_cond)*nx_cond*(mz/(2.*mx)*(nuo+mz/(3.*mx))/(nuo+1./2.)) + (1.-nx_cond)*mxz_cond*(mx/(2.*mz)*(nuo+2./3.)/(nuo+1./2)) + (1-nx_cond)*(1-mxz_cond)*(mz/(2*mx)*(nuo+1.-mz/(3*mx))/(nuo+1/2));
        float phi_crit = nx_cond*(2.*mx/(3.*mz)*(powf(nuo,3.))/(nuo+1./2.)) + (1.-nx_cond)*(2.*mx/(3.*mz)*(powf((nuo+1.),3.))/(nuo+1./2.)) ;
        //

        float M = 3*powf(mx,2.)*mz*(nuo+0.5);
        float phio = asin(1. - 2.*C/phi_crit);
        float alpha = 0.        ;
        // float phimax 

        // #################################### CONDITIONNER ALPHA ##################################
        
        if (C < phi_tr)             
        {
            if (C < phi_crit || nx < 0.)                    // equation (1 et 6) et (4 et 9)
            {
            float alpha_a = (-1 + 2*cos(phio/3 + acos(-1)/6 )) ;
            float alpha_b = (-1 + 2*sin(phio/3)) ;
            float alpha = powf(M*phi_crit/2,1/3) * (nx_cond*alpha_a + (1. - nx_cond)*alpha_b)  ;   
            }
            else                 // equation 2 et 7
            {
            float alpha_c = powf(2.*C - phi_crit + 2.*sqrt(C*(C - phi_crit)),1/3) + powf(2.*C - phi_crit - 2*sqrt(C*C - phi_crit),1./3.);
            float alpha = ( alpha_c + copysign(1.0,(mz-mx))*pow(phi_crit,1/3))*pow(M/2,1/3) ;
            }
        }
        else
        {
            if((mz-mx)>0.)
            {
                // equation 3 et 5
                float alpha = mz * C + mx/2. * (nuo + 1/3 + nx_cond*1.0/3.0)/(nuo+1/2) ; // Ajout d'une condition sur le 1/3 -> 2/3
            }
            else
            {
                if(nx<0.)
                {
                    // equation 10
                    float alpha = 1/2 * (mz - 2*mx*(nuo+1) - sqrt(4*powf(mx,2)*pow(nuo+1,2)- 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)));
                }
                else
                {
                    // equation 8 
                    float alpha = 1/2 * (mz - 2*mx*nuo + sqrt(4*powf(mx,2)*powf(nuo,2) + 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)));
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
        float Si = sqrt((results_coord[1].x - results_coord[0].x)*(results_coord[1].x - results_coord[0].x) + (results_coord[1].z - results_coord[0].z)*(results_coord[1].z - results_coord[0].z));
        return Si;
    }



    
