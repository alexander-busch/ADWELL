#include "GranularDenseFlowCorrection.h"


DEFINE_ON_DEMAND(MarkVOFzones)
{
	//Marks solid volume fraction areas
	//Needs one UDM!

	Domain *d = Get_Domain(3);
	Thread *t;
	cell_t c;

	real con, con2, source, v_mag, angle;
	real v[ND_ND], gradVOF[ND_ND], e_y[ND_ND];
	real alpha_s;

	// Initialize UDM with zero
	thread_loop_c(t, d) //Loop over all cell threads in solid domain
	{
		begin_c_loop_all(c, t) // Loop over all cells in solid phase cell threads
		{
			C_UDMI(c, t, 0) = 0.0;
		}
		end_c_loop_all(c, t)
	}

	// Assign integers to various solid phase zones
	thread_loop_c(t, d) //Loop over all cell threads in solid domain
	{
		begin_c_loop_all(c, t) // Loop over all cells in solid phase cell threads
		{
			alpha_s = C_VOF(c, t);

			// Check if gradient exists, does not exist in first iteration/after loading
			if (NULL != THREAD_STORAGE(t, SV_VOF_RG))
			{
				gradVOF[0] = C_VOF_RG(c, t)[0];
				gradVOF[1] = C_VOF_RG(c, t)[1];
#if RP_3D
				gradVOF[3] = C_VOF_RG(c, t)[2];
#endif
			}
			else
			{
				gradVOF[0] = 0.0;
				gradVOF[1] = 0.0;
#if RP_3D
				gradVOF[3] = 0.0;
#endif
			}

			ND_SET(e_y[0], e_y[1], e_y[2], 0.0, 1.0, 0.0);
			angle = acos(NV_DOT(gradVOF, e_y) / (NV_MAG(gradVOF)));

			// Get solid velocities and compute magnitude
			v[0] = C_U(c, t);
			v[1] = C_V(c, t);
#if RP_3D
			v[2] = C_W(c, t);
#endif
			v_mag = NV_MAG(v);

			if ((alpha_s >= alpha_fric) && (v_mag <= v_cr))
			{
				real dvofdx = gradVOF[0];
				real dvofdy = gradVOF[1];
				real gradVOF_mag = NV_MAG(gradVOF);
				real gradVOF_tr = gradVOF_cr;
				
				if (NV_MAG(gradVOF)<gradVOF_cr)
				{
					C_UDMI(c, t, 0) = 4.0;
				}
				else
				{
					if (fabs(angle) <= angle_rep)
					{
						C_UDMI(c, t, 0) = 3.0;
					}
					else
					{
						C_UDMI(c, t, 0) = 2.0;
					}
				}
			}
			else
			{
				C_UDMI(c, t, 0) = 0.0;
			}
		}
		end_c_loop_all(c, t)
	}
} // End







DEFINE_ON_DEMAND(showgrad)
{
	Domain *domain;
	Thread *t; domain = Get_Domain(3);
	if (!Data_Valid_P()) return;
	Message0(" >>> entering show-grad: \n ");

	thread_loop_c(t, domain)
	{
		Material *m = THREAD_MATERIAL(t);
		int nspe = MIXTURE_NSPECIES(m);
		int nspm = nspe - 1;
		Message0("::::\n ");
		Message0(":::: Reconstruction Gradients :::: \n ");
		Message0("::::\n ");
		if (NNULLP(THREAD_STORAGE(t, SV_P_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of P is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_U_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of U is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_V_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of V is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_W_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of W is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_T_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of T is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_H_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of H is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_K_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of K is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_D_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of D is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_O_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of O is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_NUT_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of NUT is available \n ");
		}
		if (nspe && NNULLP(THREAD_STORAGE(t, SV_Y_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of Species is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_VOF_RG)))
		{
			Message0("....show-grad:Reconstruction Gradient of volume fraction is available \n ");
		}
		/********************************************************************/
		/********************************************************************/
		/********************************************************************/
		/********************************************************************/
		Message0("::::\n ");
		Message0(":::: Gradients :::: \n ");
		Message0("::::\n ");
		if (NNULLP(THREAD_STORAGE(t, SV_P_G)))
		{
			Message0("....show-grad:Gradient of P is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_U_G)))
		{
			Message0("....show-grad:Gradient of U is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_V_G)))
		{
			Message0("....show-grad:Gradient of V is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_W_G)))
		{
			Message0("....show-grad:Gradient of W is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_T_G)))
		{
			Message0("....show-grad:Gradient of T is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_H_G)))
		{
			Message0("....show-grad:Gradient of H is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_K_G)))
		{
			Message0("....show-grad:Gradient of K is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_D_G)))
		{
			Message0("....show-grad:Gradient of D is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_O_G)))
		{
			Message0("....show-grad:Gradient of O is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_NUT_G)))
		{
			Message0("....show-grad:Gradient of NUT is available \n ");
		}
		if (nspe && NNULLP(THREAD_STORAGE(t, SV_Y_G)))
		{
			Message0("....show-grad:Gradient of Species is available \n ");
		}
		if (NNULLP(THREAD_STORAGE(t, SV_VOF_G)))
		{
			Message0("....show-grad:Gradient of volume fraction is available \n ");
		}
	}
}

