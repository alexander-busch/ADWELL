#include "GranularDenseFlowCorrection.h"


/* Structure test
struct gradient get_gradVOF(cell_t c, Thread *t)
{

	struct gradient gradVOF;

	// Check if gradient exists, does not exist in first iteration/after loading
	if (NULL != THREAD_STORAGE(t, SV_VOF_G))
	{
		gradVOF.x = C_VOF_G(c, t)[0];
		gradVOF.y = C_VOF_G(c, t)[1];
#if RP_3D
		gradVOF.z = C_VOF_G(c, t)[2];
#endif
		//Message0("VOF gradient exists\n ");
	}

	if (NULL != THREAD_STORAGE(t, SV_VOF_RG))
	{
		gradVOF.x = C_VOF_RG(c, t)[0];
		gradVOF.y = C_VOF_RG(c, t)[1];
#if RP_3D
		gradVOF.z = C_VOF_RG(c, t)[2];
#endif
		//Message0("VOF reconstruction gradient exists\n ");
	}
	else
	{
		gradVOF.x = 0.0;
		gradVOF.y = 0.0;
#if RP_3D
		gradVOF.z = 0.0;
#endif
	}

	return gradVOF;
}
*/


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
