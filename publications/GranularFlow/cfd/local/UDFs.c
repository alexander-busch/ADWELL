/*******************************************************************
UDF for for modifying solid velocity in dense and slow granular flows 
based on Flow3D model. Additional solid momentum source term to model
jamming of granular material.
*******************************************************************/
//#include "GranularDenseFlowCorrection.h"




#include "udf.h"
#include "mem.h"


/*************************************************************************************************
Definitions
*************************************************************************************************/ 

// General
#define g (9.81)

// Fluid
#define rho_f (RP_Get_Real("rho_f")) // Fluid density [kg/m³]
/* #define lambda (RP_Get_Real("pac_cross_lambda")); //
#define n (RP_Get_Real("pac_cross_n")); //
#define mu_0 (RP_Get_Real("pac_cross_mu0")); //
#define mu_inf (eta_h2o); //
 */
 
// Solid
#define rho_s (RP_Get_Real("rho_s")) // Solid density [kg/m³]
#define angle_intfric (RP_Get_Real("angle_intfric") * M_PI / 180.0) // Angle of repose
#define cohesion (1.0) // Cohesion coefficient [Pa]
#define alpha_pack (RP_Get_Real("alpha_pack"))// Maximum packing density
#define alpha_jamm (0.96825*alpha_pack) // Jamming density of FLOW3D
#define alpha_fric (RP_Get_Real("alpha_fric")) // Friction density, particles start to touch each other
#define c_0 (1.0) // Cohesion
#define d_s (RP_Get_Real("d_s")) // Solid particle diameter [m]
#define v_cr (1.41*c_0*pow(d_s*g*(rho_s-rho_f)/rho_f, 0.5)) // Treshold velocity (Bagnold 1941)

// Mesh
#define dx (0.001) // Grid cell length in x
#define gradVOF_cr ((alpha_pack-alpha_fric)/(6*d_s)) // Critical magnitude of solid volume fraction gradient





// Compute magnitude of deviatoric part of strain rate tensor = magnitude of shear rate tensor
real comp_shear_rate_mag(cell_t c, Thread *t)
{
	// Velocity gradient components
	real dudx = C_DUDX(c, t);
	real dudy = C_DUDY(c, t);
	real dvdx = C_DVDX(c, t);
	real dvdy = C_DVDY(c, t);
#if RP_3D
	real dudz = C_DUDZ(c, t);
	real dvdz = C_DVDZ(c, t);
	real dwdz = C_DWDZ(c, t);
	real dwdx = C_DWDX(c, t);
	real dwdy = C_DWDY(c, t);
#endif

	// Trace of strain rate tensor
#if RP_2D
	real trace = dudx + dvdy;
#else
	real trace = dudx + dvdy + dwdz;
#endif

	// Diagonal components of deviatoric part of strain rate tensor
	real S11 = dudx - trace/3;
	real S22 = dvdy - trace/3;
#if RP_3D
	real S33 = dwdz - trace/3;
#endif
	
	// Magnitude of shear rate tensor based on definition of C_STRAIN_RATE_MAG(c,t) in mem.h
#if RP_2D
	real ShearRateTensor_Mag = \
	sqrt(2.0*(pow(S11, 2.0) + pow(S22, 2.0)) + pow(dudy+dvdx, 2.0));
#else
	real ShearRateTensor_Mag = \
		sqrt(2.0*(pow(S11, 2.0) + pow(S22, 2.0) + pow(S33, 2.0)) + pow(dudz + dwdx, 2.0) + pow(dudy + dvdx, 2.0) + pow(dwdy + dvdz, 2.0));
#endif
	
	return ShearRateTensor_Mag;
}


// Compute frictional pressure
real comp_p_fric(cell_t c, Thread *t)
{
	real alpha_s, p_fric, Fr, exp_nom, exp_denom;
	alpha_s = C_VOF(c, t); //solid volume fraction

	// Cheng et al. (2017) coefficients
	Fr = 0.05;
	exp_nom = 3.0;
	exp_denom = 5.0;

	// Jackson and Johnson (1987) coefficients
	Fr = 0.1;
	exp_nom = 2.0;
	exp_denom = 5.0;

	if (alpha_s > alpha_fric)
	{
		// Syamlal et al. (1993)
		// p_fric = alpha_s * pow(10.0, 25.0) * pow((alpha_s - alpha_fric), 10.0);

		// Johnson and Jackson (1987)
		p_fric = (Fr*alpha_s) * (pow((alpha_s - alpha_fric), exp_nom) / pow((alpha_pack - alpha_s), exp_denom));

		// Schneiderbauer et al. (2012) / Chialvo et al. (2012)
		// real strain_rate_mag = comp_shear_rate_mag(c, t); // Magnitude of deviatoric part of strain rate tensor
		// p_fric = 4.0*rho_s*pow(0.2*d_s*strain_rate_mag / (alpha_pack - alpha_s), 2.0);
	}
	else
	{
		p_fric = 0.0;
	}
	return p_fric;
}




// Constant of FLOW3D source term
real comp_constant(cell_t c, Thread *t, real v_mag)
{
	real constant, angle;
	real gradVOF[ND_ND], gravity[ND_ND];
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
	ND_SET(gravity[0], gravity[1], gravity[2], 0.0, -g, 0.0);

	// Compute slope of bed based on volume fraction gradient and unit vector
	angle = acos(NV_DOT(gradVOF, gravity) / (NV_MAG(gradVOF)*NV_MAG(gravity)));

	// Check if cell has dense, slow solids
	if ((alpha_s >= alpha_fric) && (v_mag <= v_cr))
	{
		// Check if cell is in quasi-static bed by evaluating gradient of solid volume fraction
		if (NV_MAG(gradVOF) < gradVOF_cr)
		{
			// Add FLOW3D source term
			constant = -alpha_s*1.0*pow(g / 4.0 / d_s, 0.5)*(alpha_s - alpha_fric) / (alpha_pack - alpha_fric)*rho_s;
			constant = constant*((v_cr - v_mag) / v_cr);
		}
		else
		{
			// Check if cell is satisfying angle of repose
			if (fabs(angle) <= angle_intfric)
			{
				// Mark cell to set velocity to zero
				C_UDMI(c, t,0) = 1.0;
			}
			else
			{
				// Add FLOW3D source term
				constant = -alpha_s*1.0*pow(g / 4.0 / d_s, 0.5)*(alpha_s - alpha_fric) / (alpha_pack - alpha_fric)*rho_s;
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



// Solid x-velocity source term
DEFINE_SOURCE(sol_x_mom_src,c,t,dS,eqn)
{
	real v_mag = sqrt(C_VMAG2(c, t));
	real constant = comp_constant(c, t, v_mag);
	real source = constant*C_U(c, t);
	dS[eqn] = -source*fabs(C_U(c, t)) / v_mag;
	return source;
}

/* Solid y-velocity source term */
DEFINE_SOURCE(sol_y_mom_src,c,t,dS,eqn)
{
	real v_mag = sqrt(C_VMAG2(c, t));
	real constant = comp_constant(c, t, v_mag);
	real source = constant*C_V(c, t);
	dS[eqn] = -source*fabs(C_V(c, t)) / v_mag;
	return source;
}


// Solid granular temperature source term
DEFINE_SOURCE(sol_gran_temp_src,c,t,dS,eqn)
{
	real v_mag = sqrt(C_VMAG2(c, t));
	real constant = comp_constant(c, t, v_mag);
	real source = constant*C_GT(c, t);
	dS[eqn] = -source*fabs(C_GT(c, t)) / v_mag;
	return source;
}


// Set velocity to zero for fully packed cells
DEFINE_ADJUST(freeze,d)
{
	Domain *domain;
	Thread *t; 
	domain = Get_Domain(3); // 1 = mixture, 2 = fluid, 3 = solid
	cell_t c;
	
	thread_loop_c(t, domain) //Loop over all cell threads in solid domain
	{
		begin_c_loop_all(c, t) // Loop over all cells in solid phase cell threads
		{
			if ((C_VOF(c, t) >= 0.99*alpha_pack) || (C_UDMI(c,t,0) == 1.0))
			{
				C_U(c, t) = 0.0;
				C_V(c, t) = 0.0;
				C_GT(c, t) = 0.0;
#if RP_3D
				C_W(c, t) = 0.0;
#endif
			}
			else
			{

			}
		}
		end_c_loop_all(c, t)
	} // End of thread loop
} // End of function


// Limit granular temperature
DEFINE_ADJUST(limit_GT, d)
{
	Domain *domain;
	Thread *t;
	domain = Get_Domain(3); // 1 = mixture, 2 = fluid, 3 = solid
	cell_t c;
	real GT_lim = 0.00001;

	thread_loop_c(t, domain) //Loop over all cell threads in solid domain
	{
		begin_c_loop_all(c, t) // Loop over all cells in solid phase cell threads
		{
			if (C_GT(c, t) <= GT_lim)
			{
				C_GT(c, t) = GT_lim;
			}
		}
		end_c_loop_all(c, t)
	} // End of thread loop
} // End of function



DEFINE_PROPERTY(pressure_fric, c, t)
{
	real p_fric = comp_p_fric(c, t);
	return p_fric;
}


// Radial distribution function
DEFINE_PROPERTY(rad_distr_func, c, t)
{
	// Solid volume fraction
	real alpha_s = C_VOF(c, t);
	
	// Carnahan and Starling (1969)
	real g_0_CS = 1.0/(1.0-alpha_s) + 3.0*alpha_s/(2.0*pow(1.0 - alpha_s,2.0)) + pow(alpha_s,2.0)/(2.0*pow(1-alpha_s,3.0));
	
	// Iddir and Arastoopour(2005)
	real g_0_IA = 1.0 / (1.0 - alpha_s) / alpha_pack;

	// Ding and Gidaspow (1990), Laux (1997)
	real g_0 = 0.6 / (1.0 - pow(alpha_s / alpha_pack, 1.0 / 3.0));

	// Schneiderbauer & Pirker (2012)
	// real g_0 = fmin(g_0_CS, g_0_IA);
	

	
	return g_0;
}


DEFINE_PROPERTY(viscosity_fric,c,t)
{
	// Get mixture level thread
	Thread *t_mix = THREAD_SUPER_THREAD(t);
	
	real eta, eta_lim;
	real NV_VEC(cell_position);
	
	// Cell centroid location
	C_CENTROID(cell_position, c, t);

	// Limiting dynamic viscosity for no flow 
	eta_lim = 1.0E5;
	
	// Various pressures
	real p_oper = RP_Get_Real("operating-pressure");
	real p_stat = C_P(c,t_mix); p_stat = C_P(c,t);
	real p_gran = C_GP(c, t);
	real p_tot = p_stat + p_gran;
	real p_abs = p_stat + p_oper; p_abs = p_tot + p_oper;
	real p_fric = comp_p_fric(c, t);

	// Average normal stress formulation for frictional viscosity
	real normalstress = fmax(p_gran,p_fric);
	
	// Other factors for frictional viscosity
	real alpha_s = C_VOF(c, t); //solid volume fraction
	real sin_angle = sin(angle_intfric);
	real strain_rate_mag = C_STRAIN_RATE_MAG(c, t); // Magnitude of strain rate tensor
	strain_rate_mag = comp_shear_rate_mag(c, t); // Magnitude of deviatoric part of strain rate tensor

	if (alpha_s > alpha_fric)
	{
		// Frictional viscosity
		eta = normalstress*sin(angle_intfric) / (pow(2.0, 0.5)*strain_rate_mag);

		// Limit frictional viscosity
		if ((eta > eta_lim)) // || (strain_rate_mag <=0.1)
		{
			eta = eta_lim;
			//C_UDMI(c, t, 0) = 1.0;
		}
	}
	else
	{
		eta = 0.0;
	}

	return eta;
}

