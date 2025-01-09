/*************************************************************************************************
UDF for Cross rheology model
*************************************************************************************************/ 
/* Additional info
This rheology model differs from the Cross model of Fluent as follows:
 - It features a high-shear-rate Newtonian viscosity mu_inf, which is equal to the solvent, here H2O.

Requires
#include "udf.h"
*/ 

#include "udf.h"


struct FourParameterModel {
    real K, n, mu_0, mu_inf;
};


/* Compute viscosity from Cross model for current shear rate */
real ViscosityCross(struct FourParameterModel Cross, real SR)
{
	real eta = Cross.mu_inf+(Cross.mu_0-Cross.mu_inf)/(1.0+pow(Cross.K*SR, 1.0-Cross.n));
	
	return eta;
}

DEFINE_PROPERTY(rheology_cross,c,t)
{
#if !RP_HOST
	real SR, eta;

	/* Rheological properties of the fluid */
	struct FourParameterModel Cross;

	/* Does appaerently not work, host to node issue
	Cross.K = RP_Get_Real("cross/k/scaled");
	Cross.n = RP_Get_Real("cross/n/scaled");
	Cross.mu_0 = RP_Get_Real("cross/mu_0/scaled");
	Cross.mu_inf = RP_Get_Real("cross/mu_inf/scaled");
	*/
	
	/* Cross coefficient PAC2 
	Cross.K = 0.0109;
	Cross.n = 0.414;
	Cross.mu_0 = 0.0721;
	Cross.mu_inf = 0.00102;	
	*/
	
	/* Cross coefficient PAC4 */
	Cross.K = 0.0261;
	Cross.n = 0.392;
	Cross.mu_0 = 0.21;
	Cross.mu_inf = 0.00102;	
	
		
	SR = C_STRAIN_RATE_MAG(c, t);
	
	eta = ViscosityCross(Cross, SR);
	
	return eta;
#endif
}
