/*************************************************************************************************
UDF for Cross rheology model
*************************************************************************************************/ 
/* Additional info
This rheology model differs from the Cross model of fluent as follows:
 - It features a high-shear-rate Newtonian viscosity mu_inf, which is equal to the solvent, here H2O.
 - It has the exponent = n_Cr, the Fluent version has the exponent = 1 - n_Cr

Requires
#include "udf.h"
*/ 

#include "ParticleTrajectory.h"




DEFINE_PROPERTY(rheology_cross,c,t)
{
	real SR, eta;

	/* Rheological properties of the fluid */
	struct FourParameterModel Cross;
	Cross.lambda = RP_Get_Real("cross/lambda");
	Cross.n = RP_Get_Real("cross/n");
	Cross.mu_0 = RP_Get_Real("cross/mu_0");
	Cross.mu_inf = RP_Get_Real("cross/mu_inf");
	
	SR = C_STRAIN_RATE_MAG(c, t);
	
	eta = ViscosityCross(Cross, SR);
	
	return eta;
}
