/*************************************************************************************************
VLES turbulence model of Johansen et al. (2004)
*************************************************************************************************/ 
/* How to
- Requires udf.h and model.h
- Use realizable k-eps
- Enable 1 user defined scalar (but switch of calculation of the scalar by the solver)
- Enable 1 user defined memory

Requires
#include "udf.h"
#include "mem.h"

Changes to version provided by Jan Erik Olsen/SINTEF
- C_UDMI(c,t,0) --> C_UDMI(c,t,TKE)
- cxboolean VLES_model_on   = TRUE; --> ... = RP_Get_Boolean("models/vles");

*/

#include "model.h"

/* !!! This may have to go into an define adjust function similar to the rheology model coefficients */
cxboolean VLES_model_on   = RP_Get_Boolean("models/vles"); /* cxboolean VLES_model_on   = TRUE; */
real C3                   = 1.22;

/* VLES turbulent viscosity */
DEFINE_TURBULENT_VISCOSITY(VLES_turb_visc, c, t)
{
	real mu_t, mu_t_standard, mu_t_filter, filtersize;
	filtersize = pow(C_VOLUME(c,t), 0.33333); 
	mu_t = 0.0;  

	if ( (C_D(c,t) > 1.0e-6) && (C_K(c,t) > 1.0e-6) )
	{
		mu_t_standard = 0.09 * C_R(c,t) * SQR(C_K(c,t)) / C_D(c,t);
		mu_t_filter   = mu_t_standard * C3 * filtersize * C_D(c,t) / pow(C_K(c,t), 1.5);
		
		if ((mu_t_standard > mu_t_filter) && (VLES_model_on))
		{
		  mu_t = mu_t_filter; 
		}
		else
		{
		  mu_t = mu_t_standard;
		}
	}
	
	return mu_t;
}

/*Calculates gradient of filtersize if the filtersize is based on mesh size */
DEFINE_ADJUST(grid_gradient, domain)
{
	Thread *t;
	cell_t c;
	real s;
	cxboolean t_store_uds;
	
	thread_loop_c(t,domain)
	{
		begin_c_loop(c,t)
		{
			C_UDSI(c,t,0) = pow(C_VOLUME(c,t),1./3.);
			t_store_uds    = ( NNULLP(T_STORAGE_R(t,SV_UDSI_G(0))) ); /* t_store_uds checks if the gradient of C_UDSI(c,t, 0) exists */
			
			if (t_store_uds)
				{
					C_UDMI(c,t,TKE) = ( C_U(c,t)*C_UDSI_G(c,t,0)[0] + C_V(c,t)*C_UDSI_G(c,t,0)[1] + C_W(c,t)*C_UDSI_G(c,t,0)[2] ) * 0.75 * pow(C_D(c,t),2./3.) * pow(C_UDSI(c,t,0),-1./3.);
				}
				else
				{
					C_UDMI(c,t,TKE) = 0.0;
				}
		}
		end_c_loop(c,t)
	}
}


/* TKE source term in fluid/phase zone */ 
DEFINE_SOURCE(k_grid_gradient,c,t,dS,eqn)
{
	real source = 0.0;
	  
	if (VLES_model_on) source  = -C_UDMI(c,t,TKE);
	{
		dS[eqn] = 0.0;
	}
	
	return source;
}

