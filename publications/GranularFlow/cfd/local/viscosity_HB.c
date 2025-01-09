
#include "udf.h"
#include "mem.h"

/* Initialisation */
DEFINE_INIT(initialisation, mix_domain)
{
	/* DEFINE_INIT always refers to the mixture domain. */
	
	int phase_domain_index;
	cell_t cell;
	Thread *cell_thread;
	Domain *subdomain;
	real xc[ND_ND];
		
	/* loop over all subdomains (phases) in the superdomain (mixture) */
	sub_domain_loop(subdomain, mix_domain, phase_domain_index)
	{
		/* loop if secondary phase */
		if (DOMAIN_ID(subdomain) == 1)
		{
			/* loop over all cell threads in the secondary phase domain */
			thread_loop_c(cell_thread,subdomain)
			{
				/* loop over all cells in secondary phase cell threads */
				begin_c_loop_all(cell,cell_thread)
				{
					C_UDMI(cell, cell_thread, 1) = 0.0;
				}
				end_c_loop_all(cell,cell_thread)
			}
		}
	}
}



DEFINE_PROPERTY(viscosity_HB,c,t)
{

/* eta = [tau0 + k*(S_dot^n - (tau0/eta0)^n)]/S_dot ; S_dot > tau0/eta0
eta = eta0; ; S_dot < tau0/eta0
 */

real angle, cohesion, eta,eta0,k,n,tau0,S_dot,S_crit;
real NV_VEC(cell_position);

/* Get cell centroid location:
cell_position[0] = x
cell_position[1] = y
cell_position[2] = z in 3D only */
C_CENTROID(cell_position,c,t);


/* Coefficients */
angle = 30.0; /* Angle of repose [°] */
cohesion = 1.0; /* Cohesion coefficient [Pa] */
tau0 = cohesion + C_P(c,t)*tan(angle*M_PI/180.0); /* Yield stress [Pa] */
k = 1.0; /* Consistency index */
n = 1.0; /* Flow index  */
eta0 = 1.0E6; /* Limiting dynamic viscosity for no flow */

/* Strain rates */
S_dot = C_STRAIN_RATE_MAG(c,t);
S_crit = MAX(tau0, 1.0E-10) / MAX(eta0, 1.0E-10);

/* H-B equation */
if (S_dot <= S_crit)
{
	eta = eta0;
/* 	C_UDMI(c, t, 1) = 1.0; /* Mark cell in order to freeze velocity in cell */
}
else
{
	eta = (tau0 + k*(pow(S_dot, n) - pow(S_crit, n)))/S_dot;
	/* eta = (tau0 + k*(pow(S_dot, n) - pow(S_crit, n)))/S_dot; */
/* 	C_UDMI(c, t, 1) = 0.0; /* Mark cell in order to freeze velocity in cell */
}


/*
CX_Message("Sdot is %10.3e\n", S_dot);
CX_Message("Viscosity is %10.3e\n", eta);
*/

return eta;
}

DEFINE_ADJUST(freeze_velocities,d)
{
	Domain *mix_domain = Get_Domain(1); /* Get mixture domain */
	int phase_domain_index;
	cell_t cell;
	Thread *cell_thread;
	Domain *subdomain;
	real xc[ND_ND];
	
	/* loop over all subdomains (phases) in the superdomain (mixture) */
	sub_domain_loop(subdomain, mix_domain, phase_domain_index)
	{
		/* loop if secondary phase */
		if (DOMAIN_ID(subdomain) == 1)
		{
			/* loop over all cell threads in the secondary phase domain */
			thread_loop_c (cell_thread,subdomain)
			{
				/* loop over all cells in secondary phase cell threads */
				begin_c_loop_all (cell,cell_thread)
				{
 					if (C_UDMI(cell, cell_thread, 1) == 1.0)
					{
						/* Freeze velocity in cell */
						C_U(cell,cell_thread) = 0.0;
						C_V(cell,cell_thread) = 0.0;
						C_W(cell,cell_thread) = 0.0;
					}
				}
				end_c_loop_all (cell,cell_thread)
			}
		}
	}
}