/*************************************************************************************************
UDF for Cross rheology model for parallel computations/hpc
*************************************************************************************************/ 
/* Additional info
This rheology model differs from the Cross model of Fluent as follows:
 - It features a high-shear-rate Newtonian viscosity mu_inf, which is equal to the solvent, here H2O.

Requires
#include "udf.h"


A DEFINE_ADJUST function that is called at the beginning of each iteration gets the RP variables
on the host only and distributes them to the nodes. This is because the macro RP_Get_... cannot be
called within DEFINE_PROPERTY functions or any other UDF that passes cells and threads, since they
are called multiple times within a single process and (potentially) uneven between processors.
Moreover, RP_Get_... functions are valid for the host only, since it is the only process that can
access session variables.

Note that the structure variable Cross is defined globally in order to be accessed by the
DEFINE_PROPERTY function with the values passed in the DEFINE_ADJUST.

Also note that the DEFINE_ADJUST is called only during runtime. This means that when you initialize
the flow the viscosity will not be updated (you need to run one iteration).

*/ 


#include "udf.h"


/* Definition of structure FourParameterModel */
struct FourParameterModel {
    real K, n, mu_0, mu_inf;
};


/* Compute viscosity from Cross model for current shear rate */
real ViscosityCross(struct FourParameterModel Cross, real SR)
{
	real eta = Cross.mu_inf+(Cross.mu_0-Cross.mu_inf)/(1.0+pow(Cross.K*SR, 1.0-Cross.n));
	
	return eta;
}


/* Compute viscosity based on Cross model coefficients and current strain rate */
DEFINE_PROPERTY(rheology_cross,c,t)
{
#if !RP_HOST
	real SR, eta;

	/* Rheological properties of the fluid */
	struct FourParameterModel Cross;
	
	/* Cross coefficient PAC2 */
	Cross.K = 0.0109;
	Cross.n = 0.414;
	Cross.mu_0 = 0.0721;
	Cross.mu_inf = 0.00102;
	
	/* Cross coefficient PAC4 */
	/*
	Cross.K = 0.0261;
	Cross.n = 0.392;
	Cross.mu_0 = 0.21;
	Cross.mu_inf = 0.00102;	
	*/
	
	SR = C_STRAIN_RATE_MAG(c, t);
	
	eta = ViscosityCross(Cross, SR);
	
	return eta;
#endif
}
