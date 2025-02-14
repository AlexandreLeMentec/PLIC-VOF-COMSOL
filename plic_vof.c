#include <math.h>

float plic_cyl(float C, float nx, float nz, float dx, float dz, float x)
    {
        float nuo = x/ dx;
        float mx = fabs(nx)*dx/(fabs(nz)*dz + fabs(nx)*dx); // TODO: Checker si c'est dans le bon sens 
        float mz = fabs(nz)*dz/(fabs(nz)*dz + fabs(nx)*dx);

        //
        float nx_cond = 0.5 + copysign(0.5,nx);                      // condition nx > = 0   
        float mxz_cond = 0.5 + copysign(0.5,(mz-mx));                // condition mx < = mz
        float phi_tr = mxz_cond*nx_cond*(mx/(2.*mz)*(nuo+1./3.)/(nuo+1./2.)) + (1.-mxz_cond)*nx_cond*(mz/(2.*mx)*(nuo+mz/(3.*mx))/(nuo+1./2.)) + (1.-nx_cond)*mxz_cond*(mx/(2.*mz)*(nuo+2./3.)/(nuo+1./2)) + (1-nx_cond)*(1-mxz_cond)*(mz/(2*mx)*(nuo+1.-mz/(3*mx))/(nuo+1/2));
        float phi_crit = nx_cond*(2.*mx/(3.*mz)*(powf(nuo,3.))/(nuo+1./2.)) + (1.-nx_cond)*(2.*mx/(3.*mz)*(powf((nuo+1.),3.))/(nuo+1./2.)) ;
        //

        float M = 3*powf(mx,2.)*mz*(nuo+0.5);
        float phio = asin(1. - 2.*C/phi_crit);
        // float phimax 


        // equation (1 et 6) et (4 et 9)
        float alpha_a = (-1 + 2*cos(phio/3 + acos(-1)/6 )) ;
        float alpha_b = (-1 + 2*sin(phio/3)) ;
        float alpha = powf(M*phi_crit/2,1/3) * (nx_cond*alpha_a + (1. - nx_cond)*alpha_b)  ;          

        // equation 2 et 7
        float alpha_c = powf(2.*C - phi_crit + 2.*sqrt(C*(C - phi_crit)),1/3) + powf(2.*C - phi_crit - 2*sqrt(C*C - phi_crit),1./3.);
        float alpha = ( alpha_c + mxz_cond*pow(phi_crit,1/3))*pow(M/2,1/3) ;          // CONDITION à vérifier
        
        // equation 3 et 5
        float alpha = mz * C + mx/2. * (nuo + 2/3)/(nuo+1/2) ;

        // equation 8 
        float alpha = 1/2 * (mz - 2*mx*nuo + sqrt(4*powf(mx,2)*powf(nuo,2) + 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)))    ;

        // equation 10
        float alpha = 1/2 * (mz - 2*mx*(nuo+1) - sqrt(4*powf(mx,2)*pow(nuo+1,2)- 8 * powf(mx,2)*(nuo+1/2)*C - 1/3*powf(mz,2)))   ;

        




    }



    
