/* Granular dense flow correction header file */




/*************************************************************************************************
Include
*************************************************************************************************/ 

#include <string.h>
#include <stdio.h>

#include "udf.h"
#include "unsteady.h"
#include "mem.h"
#include "sg.h"



/*************************************************************************************************
Definitions
*************************************************************************************************/ 

// General
#define g (9.81)

// Fluid
#define rho_f (1000.0) // (RP_Get_Real("rho_f")) // Fluid density [kg/m³]
/* #define lambda (RP_Get_Real("pac_cross_lambda")); //
#define n (RP_Get_Real("pac_cross_n")); //
#define mu_0 (RP_Get_Real("pac_cross_mu0")); //
#define mu_inf (eta_h2o); //
 */
 
// Solid
#define rho_s (2560.0) // (RP_Get_Real("d_s")) // Solid density [kg/m³]
#define angle_rep (30.0 * M_PI / 180.0) // Angle of repose
#define alpha_pack (0.63) // Maximum packing density
#define alpha_jamm (0.96825*alpha_pack) // Jamming density of FLOW3D
#define alpha_fric (0.55) // Friction density, particles start to touch each other
#define c_0 (1.0) // Cohesion
#define d_p (0.001) // (RP_Get_Real("d_s")) // Solid particle diameter [m]
#define v_cr (1.41*c_0*pow(d_p*g*(rho_s-rho_f)/rho_f, 0.5)) // Treshold velocity (Bagnold 1941)

// Mesh
#define dx (0.0005) // Grid cell length in x
#define gradVOF_cr ((alpha_pack-alpha_fric)/dx) // Critical magnitude of solid volume fraction gradient




/*************************************************************************************************
Structures
*************************************************************************************************/ 
/* struct gradient {
	real dx, dy;
#if RP_3D
	real dz;
#endif
};
*/

/*************************************************************************************************
Functions
*************************************************************************************************/ 
//extern struct gradient get_gradVOF(cell_t c, Thread *t);

// extern real get_constant(cell_t c, Thread *t, real v_mag)


// Get Volume fraction gradients from Fluent solver for current cell
real comp_constant(cell_t c, Thread *t, real v_mag)
{
	real constant, angle;
	real gradVOF[ND_ND], e_y[ND_ND];
	real alpha_s = C_VOF(c, t);

	// Check if gradient exists, does not exist in first iteration/after loading
	if (NULL != THREAD_STORAGE(t, SV_VOF_G))
	{
		gradVOF[0] = C_VOF_G(c, t)[0];
		gradVOF[1] = C_VOF_G(c, t)[1];
#if RP_3D
		gradVOF[2] = C_VOF_G(c, t)[2];
#endif
		//Message0("VOF gradient exists\n ");
	}

	// Check if reconstruction gradient exists, does not exist in first iteration/after loading
	if (NULL != THREAD_STORAGE(t, SV_VOF_RG))
	{
		gradVOF[0] = C_VOF_RG(c, t)[0];
		gradVOF[1] = C_VOF_RG(c, t)[1];
#if RP_3D
		gradVOF[2] = C_VOF_RG(c, t)[2];
#endif
		//Message0("VOF reconstruction gradient exists\n ");
	}
	else
	{
		gradVOF[0] = 0.0;
		gradVOF[1] = 0.0;
#if RP_3D
		gradVOF[2] = 0.0;
#endif
	}


	// Define unit vector in opposite direction of gravity
	ND_SET(e_y[0], e_y[1], e_y[2], 0.0, 1.0, 0.0);

	// Compute slope of bed based on volume fraction gradient and unit vector
	angle = acos(NV_DOT(gradVOF, e_y) / (NV_MAG(gradVOF)));

	// Check if cell has dense, slow solids
	if ((alpha_s >= alpha_fric) && (v_mag <= v_cr))
	{
		if (NV_MAG(gradVOF) < gradVOF_cr)
		{
			constant = -alpha_s*1.0*pow(g / 4.0 / d_p, 0.5)*(alpha_s - alpha_fric) / (alpha_pack - alpha_fric)*rho_s;
			constant = constant*((v_cr - v_mag) / v_cr);
		}
		else
		{
			if (fabs(angle) <= angle_rep)
			{
				constant = 0.0;
			}
			else
			{
				constant = -alpha_s*1.0*pow(g / 4.0 / d_p, 0.5)*(alpha_s - alpha_fric) / (alpha_pack - alpha_fric)*rho_s;
				constant = constant*((v_cr - v_mag) / v_cr);
			}
		}
	}
	else
	{
		constant = 0.0;
	}

	return constant;
}

