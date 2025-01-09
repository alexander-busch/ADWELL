/*************************************************************************************************
UDF for specifying a steady-state velocity profile at 2D channel inlet
Fluent R17.2 UDF manual, 8.1.3. (page 496)
*************************************************************************************************/
/* Additional info
Parabolic velocity profile

Dimensional version adjusts for potentially changing sediment bed with height h_bottom

Requires
#include "udf.h"
*/ 

#include "ParticleTrajectory.h"


/* Parabolic velocity profile

Dimensional version adjusts for potentially changing sediment bed with height h_bottom

Requires
#include "udf.h"
*/ 
DEFINE_PROFILE(par_vel_pro, thread, position)
{
	real x[ND_ND]; /* Position vector */
	real y, u1, u2;
	face_t f;
	real m = 1.0/7.0;
	real uTmax = (m+1.0)*(m+2.0)*u_bulk/2.0;
	
	
	/* Loop all inlet faces */
	begin_f_loop(f, thread)
	{
		/* Get vector x to current face center  */
		F_CENTROID(x, f, thread);
		
		if (rp_turb)  /* Turbulent velocity profile */
		{
			if (x[1] < h_0)
			{
				u1 = 0.0;
			}
			else if (x[1] < (h_0+(h-h_0)/2.0))
			{
				u1 = uTmax*pow(2.0*((x[1]-h_0)/(h-h_0)),m);
			}
			else
			{
				u1 = uTmax*pow(2.0*(1.0-(x[1]-h_0)/(h-h_0)),m);
			}
		}
		else /* Laminar velocity profile */
		{	
			
			/* Non-dimensional y coordinate and velocity u = f(x2) */
			y = 2.0*(x[1] - 0.5*h) / h; 
			u1 = 3.0/2.0*u_bulk*(1.0 - y*y);
			
			/* Dimensional y coordinate and velocity u = f(x2) */
			y = x[1];
			u2 = 2.0*u_bulk-2.0*u_bulk/pow(((h_0-h)/2.0),2.0)*pow((y-(h_0 +(h- h_0)/2.0)),2.0);
		}
		
		/* Assign velocity */
		F_PROFILE(f, thread, position) = u1;
	}
	end_f_loop(f, thread)
}


/* Experimental velocity profile

Based on polynomials fitted toe xperimental data

Requires
#include "udf.h"
*/ 
DEFINE_PROFILE(exp_vel_pro, thread, position)
{
	real x[ND_ND]; /* Position vector */
	real y, u;
	real factor = 1.0;
	face_t f;
	real p1, p2, p3, p4, p5, p6;
	
	/* Determine coefficients of polynomial based on current fluid */	
	int fluid = RP_Get_Integer("fluid");
	if (fluid == 1) /* H2O */
	{
		p1 = -5.699e+08;
		p2 = -4.222e+07 ;
		p3 = -1.055e+06;
		p4 = -1.161e+04;
		p5 = -57.62;
		p6 = 0.01091;
	}
	else if (fluid == 2) /* PAC2 */
	{
		p1 = 0.0;
		p2 = 0.0;
		p3 = 0.0;
		p4 = -780.4;
		p5 = -16.33;
		p6 = 0.002025;
	}
	else /* PAC4 */
	{
		p1 = 0.0;
		p2 = 0.0;
		p3 = 0.0;
		p4 = -1202.0;
		p5 = -24.86;
		p6 = 0.02218;
		
		/* Scale 0.085 velocity profile to 0.048 */
		if (u_bulk == 0.048)
		{
			factor = 0.048/0.085;
		}
		else
		{
			factor = 1.0;
		}
		
	}
	
	/* Loop all inlet faces */
	begin_f_loop(f, thread)
	{
		/* Get vector x to current face center  */
		F_CENTROID(x, f, thread);
		
		y = x[1]-h;

		u = factor * (p1*pow(y,5.0) + p2*pow(y,4.0) + p3*pow(y,3.0) + p4*pow(y,2.0) + p5*pow(y,1.0) + p6);
		
		/* Assign velocity */
		F_PROFILE(f, thread, position) = u;
	}
	end_f_loop(f, thread)
}
