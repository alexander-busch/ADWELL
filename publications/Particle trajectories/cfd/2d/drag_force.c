/*************************************************************************************************
UDF for computing DPM particle drag force for shear-thinning fluids

in the form required by Fluent: (18 Cd Re/24)
- where in general both c_D and Re are functions of a viscosity that varies with shear rate
- and the effective shear rate as seen by the particle is computed as the vector sum of Fluents
	shear rate and a shear rate induced by the particles relative velocity.

Two DPM drag macros DEFINE_DPM_DRAG are provided:
- DEFINE_DPM_DRAG(drag_force_st,Re,p) for shear-thinning fluids as described above and based on
	drag laws as defined in section "Functions" below.
- DEFINE_DPM_DRAG(drag_force_n,Re,p) for Newtonian fluids based on Schiller & Naumann (1933)
  
  
See Fluent R17.2 UDF manual, 2.5.3. (page 183)
*************************************************************************************************/ 
/* Additional info
DEFINE_DPM_DRAG is used to specify the drag between particles and fluid as a dimensionless group:
	(18 Cd Re/24)
	as it appears in the drag force per unit particle mass (eq. (2.19), page 183).

Drag coefficient of Fluent ("Spherical drag law") based on
	R. Clift, J. R. Grace and M. E. Weber "Bubbles, Drops, and Particles" (1978).

Requires
	#include "udf.h"
*/ 


#include "ParticleTrajectory.h"



/*************************************************************************************************
Functions (Various drag laws) 
*************************************************************************************************/
/* Additional info

*/ 
/* Coefficient of drag - Schiller-Naumann (1933) */
real c_D_SchillerNaumann(real Re)
{
	real c_D;
	if (Re > 1000.0)
	{
		c_D = 0.44;
	}
	else
	{
		c_D = (24.0 / Re)*(1.0 + 0.15*pow(Re, 0.687));
	}
	return (c_D);
}


/* Viscoelastic correction of drag_force, Acharya (1976)
- corrects Acharya drag law for viscoelastic effects
- uses particle Reynolds number based on  v_r / d_p */
real c_D_correction_Acharya(struct FourParameterModel Carreau, struct PowerLaw PL, real Re, real SR, real c_D)
{
	/* Time scale lambda */
	real lambda;
	int fluid = RP_Get_Integer("fluid");
	if (fluid == 1) /* H2O */
	{
		lambda = 0.0;
	}
	else if (fluid == 2) /* PAC2 */
	{
		/* Carreau time constant */
		lambda = Carreau.lambda;
	}
	else /* PAC4 */
	{
		/* PAC4 K_FNSC and n_FNSC from Busch et al. (2017) */
		real lambda = pow(0.084 / 2.0 / PL.K, 1.0 / (0.964 - PL.n));

		/* Carreau time constant */
		lambda = Carreau.lambda;
	}

	/* Weissenberg number */
	real Wi = lambda * SR;

	/* Correction of drag coefficent */
	c_D = c_D*(1 - 0.18* pow(Re*Wi, 0.19));

	return (c_D);
}


/* Coefficient of drag - Acharya (1976)
	- determines large Re drag coefficient based on PL model
	- the functions determines PL coefficients from the Carreau model for a given shear rate */
real c_D_Acharya(struct FourParameterModel Carreau, real v_r, real d_p, real Re, real SR)
{
	real f1, f2, f3;
	real c_D;
	struct PowerLaw PL;
	PL	= Carreau2PL(Carreau, SR);
	
	/* Drag law coefficients */
	f1 = pow(3.0,(3.0*PL.n-3.0)/2.0) \
		* (33.0*pow(PL.n,5.0)-64.0*pow(PL.n,4.0)-11.0*pow(PL.n,3.0)+97.0*pow(PL.n,2.0)+16.0*PL.n) \
		/ (4.0*pow(PL.n,2.0)*(PL.n+1.0)*(PL.n+2.0)*(2.0*PL.n+1.0));	
	f2 = 10.5*PL.n - 3.5;
	f3 = 0.32*PL.n + 0.13;
	
	/* Drag law coefficient correction of Kawase and Ulbrecht (1981) */
	f1 = pow(3.0,(3.0*PL.n-3.0)/2.0) \
		* (-22.0*pow(PL.n,2.0)+29.0*PL.n+2.0) \
		/ (PL.n*(PL.n+2.0)*(2.0*PL.n+1.0));	
		
	/* Drag law coefficient correction of Kawase and Moo-Young (1986) */
	f1 = pow(3.0,(3.0*PL.n-3.0)/2.0) \
		* (-7.0*pow(PL.n,2.0)+4.0*PL.n+26.0) \
		/ (5.0*(PL.n+2.0));
	
	/* Particle Reynolds number */
	Re  = rho_f*pow(v_r,(2.0-PL.n))*pow(d_p,PL.n)/PL.K;

	/* Drag coefficient */
	if (Re < 1.0)
	{
		c_D = 24.0/Re*f1;
	}
	else
	{
		c_D = 24.0/Re*f1 + f2/(pow(Re,f3));
	}

	/* Viscoelastic correction of drag_force */
	c_D = c_D_correction_Acharya(Carreau, PL, Re, SR, c_D);

	return (c_D);
}


/* Coefficient of drag - Chhabra & Uhlherr (1980)
	- corrects Schiller & Naumann drag coefficient for nN effects based on Carreau model
	- uses Newtonian zero-shear viscosity based particle Reynolds number */
real c_D_ChhabraUhlherr(struct FourParameterModel Carreau, real Re0, real SR)
{
	/* Carreau number */
	real Ca = Carreau.lambda * SR; /* SR_p; */
	
	/* Schiller-Naumann coefficient of drag */	
	real c_D = c_D_SchillerNaumann(Re0);
	
	/* nN Schiller-Naumann correction factor */
	real correction = 1.0 + 0.65*(Carreau.n-1.0)*pow(Ca,0.2);
		
	/* nN corrected Schiller-Naumann coefficient of drag */
	return (c_D*correction);
}


/* Coefficient of drag - Fluents spherical drag law
	- suggested valid range 20 < Re < 260
	- code is taken from https://www.cfd-online.com/Forums/fluent-udf/100432-how-add-velocity-dpm.html
	- is also available as macro in Fluent "drag_force = SphereDragCoeff(Re);"
	- "drag_force = SphereDragCoeff(p->Re);" makes use of Reynolds number provided to DPM drag law */
real c_D_Fluent_spherical(real Re)
{
	real drag_force;
	/* Code snippet from CFD online */
	if (Re < 0.01)
	{
		drag_force=18.0;
	}
	else if (Re < 20.0)
	{
		real w = log10(Re);
		drag_force = 18.0 + 2.367*pow(Re,0.82-0.05*w) ;
	}
	else
	{
		drag_force = 18.0 + 3.483*pow(Re,0.6305) ;
	}

	/* Fluent macro */
	drag_force = SphereDragCoeff(Re);

	/* Coefficient of drag */
	real c_D = 24.0*drag_force/18.0/Re;
	return (c_D);
}



/*************************************************************************************************
UDF particle drag force for non-Newtonian liquids based on several nN drag law models
*************************************************************************************************/ 
/* Additional info
Computational procedure:
1. Relative particle velocity and magnitude of relative particle velocity
2. Effective (Newtonian) shear rate as seen by the particle
3. True viscosity as seen by the particle based on Cross or Carreau constitutive equation
4. Reynolds number, either using efffective viscosity from 3. or depending on drag law
5. Drag law, one of the above defined functions
6. Return (18 Cd Re/24) to Fluent solver
*/ 
DEFINE_DPM_DRAG(drag_force_st,Re,p)
{	
	/* Tracked_Particle *p cell index of the cell that the particle is currently in */
	cell_t c_p = P_CELL(p);
	
	/* Tracked particle *p pointer to the thread the particle currently is in */
	Thread *t_p = P_CELL_THREAD(p);

	/* Tracked particle position vector */
	real x[ND_ND];

	/* Particle diameter */
	real d_p = P_DIAM(p);
	
	/* Fluid density */
	real rho_l = C_R(c_p,t_p);
	
	/* Velocity and shear rate vectors */
	real u[ND_ND], v[ND_ND], v_r[ND_ND], a[ND_ND], e_m[ND_ND], e_p[ND_ND], SR_vector[ND_ND];
	
	/* Looping variable and spatial dimensions */
	int i, idim = ND_ND;
	
	/* Models to be used */
	int model_SR = RP_Get_Integer("model/sr");
	int model_SR_p = RP_Get_Integer("model/sr_p");
	char *model_R = RP_Get_String("model/rheology");
	char *model_c_D = RP_Get_String("model/drag");

	/* Rheological properties of the fluid */
	struct FourParameterModel Cross;
	Cross.lambda = RP_Get_Real("cross/lambda");
	Cross.n = RP_Get_Real("cross/n");
	Cross.mu_0 = RP_Get_Real("cross/mu_0");
	Cross.mu_inf = RP_Get_Real("cross/mu_inf");
	struct FourParameterModel Carreau;
	Carreau.lambda = RP_Get_Real("carreau/lambda");
	Carreau.n = RP_Get_Real("carreau/n");
	Carreau.mu_0 = RP_Get_Real("carreau/mu_0");
	Carreau.mu_inf = RP_Get_Real("carreau/mu_inf");
	
	/* Shear rates & viscosity */
	real SR, SR_p, SR_f = C_STRAIN_RATE_MAG(c_p, t_p);
	real eta, eta_Fluent, Re_Fluent, Re0;
	
	/* Drag law */
	real w, c_D, c_D_Fluent, drag_force;

	
	
	/* ----------------------------------------------------------------------------------------------------------- */ 
	/* Relative velocity and magnitude of relative velocity */
	/* ----------------------------------------------------------------------------------------------------------- */ 
	
	/* Fluid velocity components of background shear flow */
	u[0] = C_U(c_p,t_p);
	u[1] = C_V(c_p,t_p);
	#if RP_3D
	u[2] = C_W(c_p,t_p);
	#endif
	
	/* Loop all spatial dimensions */ 
	for (i=0; i < idim; i++)
	{
		/* Particle position components */
		x[i] = P_POS(p)[i];
		real x1 = x[0];
		real x2 = x[1];

		/* Particle velocity components */
		v[i] = P_VEL(p)[i];
		
		/* Particle relative velocity components */
		v_r[i] = u[i]-v[i];
	}
	
	/* Velocity magnitudes */
	real u_mag = NV_MAG(u);
	real v_r_mag = NV_MAG(v_r);
	
	
	/* ----------------------------------------------------------------------------------------------------------- */ 
	/* Effective (Newtonian) shear rate as seen by the particle */
	/* ----------------------------------------------------------------------------------------------------------- */ 
	
	/* Particle relative-velocity induced  shear rate */
	if (model_SR_p == 1) /* e.g. Acharya (1976) */
	{
		SR_p = v_r_mag/d_p;
	}
	else if (model_SR_p == 2) /* e.g. Chhabra & Uhlherr (1980) */
	{
		SR_p = v_r_mag/(d_p/2.0);
	}
	else if (model_SR_p == 3) /* Renaud et al. (2004) */
	{
		SR_p = v_r_mag/(d_p/sqrt(6.0));
	}
	else if (model_SR_p == 4) /* tbd */
	{}
	else /*  */
	{
		
	}
		
	/* Total shear based on current shear rate model */
	if (model_SR == 1) /* STJ, Literature */
	{
		SR = pow( pow(SR_f,2.0) + pow(SR_p,2.0), 0.5);
	}
	else if (model_SR == 2) /* Addition of shear rate vectors */
	{
		NV_VS_VS(SR_vector, =, u, *, SR_f/u_mag, +, v_r, *, SR_p/v_r_mag);
		SR = NV_MAG(SR_vector );
	}
	else if (model_SR == 3) /* Substraction of shear rate vectors */
	{
		NV_VS_VS(SR_vector, =, u, *, SR_f/u_mag, -, v_r, *, SR_p/v_r_mag);
		SR = NV_MAG(SR_vector );
	}
	else if (model_SR == 4) /*  */
	{
		for (i = 0; i < idim; i++)
		{
			/* Unit vector components of fluid velocity directions */
			e_m[i] = u[i] / u_mag;

			/* Unit vector components of particle relative velocity directions */
			e_p[i] = v_r[i] / v_r_mag;
		}
	
		#if RP_2D
		real factor = abs(e_m[0]*e_p[1]-e_m[1]*e_p[0]);
		#endif
		
		#if RP_3D
		NV_CROSS(a,e_m,e_p);
		real factor = NV_MAG(a);
		#endif
		SR = pow(pow(SR_f, 2.0) + pow(factor*SR_p, 2.0), 0.5);
	}
	else /*  */
	{

	}
	
	
	/* ----------------------------------------------------------------------------------------------------------- */	
	/* True viscosity as seen by the particle */
	/* ----------------------------------------------------------------------------------------------------------- */
	
	int fluid = RP_Get_Integer("fluid");

	/* Viscosity as determined by Fluent */
	eta_Fluent = C_MU_L(c_p, t_p);

	if (fluid == 1) /* H2O */
	{
		eta = eta_Fluent;
	}
	else /* Shear-thinning case, evaluate effective viscosity with computed total shear rate */
	{
		if (strcmp("Cross",model_R) == 0) /* Rheological model: CROSS */
		{
			eta = ViscosityCross(Cross, SR);
		}
		else if (strcmp("Carreau",model_R) == 0) /* Rheological model: CARREAU */
		{
			eta = ViscosityCarreau(Carreau, SR);
		}
		else /*  */
		{
		
		}
	}
	
	/* Turbulent viscosity */
	/* if (rp_turb)
	{
		
	}
	*/

	
	/* ----------------------------------------------------------------------------------------------------------- */	
	/* Reynolds number */
	/* ----------------------------------------------------------------------------------------------------------- */
	
	/* Fluent default particle Reynolds number based on Fluent macro */
	Re_Fluent = RE_NUMBER(rho_l,v_r_mag,d_p,eta_Fluent);
	
	/* Fluent default particle Reynolds number based on effective viscosity, should be identical to previous */
	Re_Fluent = rho_l*v_r_mag*d_p / eta_Fluent;

	/* True particle Reynolds number based on true effective viscosity, should be identical to previous */
	Re = rho_l*v_r_mag*d_p / eta;
	
	/* True particle Reynolds number based on Fluent macro */
	Re = RE_NUMBER(rho_l, v_r_mag, d_p, eta);
	
	if (strcmp("CU",model_c_D) == 0)
	{
		/* Particle Reynolds number based on Carreau zero-shear viscosity as required for Chhabra & Uhlherr (1980) drag law */
		Re0 = rho_l*v_r_mag*d_p / Carreau.mu_0;
	}
	
	/* ----------------------------------------------------------------------------------------------------------- */	
	/* Drag law */
	/* ----------------------------------------------------------------------------------------------------------- */
		
	if (strcmp("FS",model_c_D) == 0) /* Fluents spherical drag law */
	{
		c_D_Fluent = c_D_Fluent_spherical(Re_Fluent);
		c_D = c_D_Fluent_spherical(Re);
	}
	else if (strcmp("SN",model_c_D) == 0) /* Schiller & Naumann (1935) */
	{
		c_D_Fluent = c_D_SchillerNaumann(Re_Fluent);
		c_D = c_D_SchillerNaumann(Re);
	}
	else if (strcmp("A",model_c_D) == 0) /* Acharya (1976) and correction for viscoelastic effects */
	{
		c_D_Fluent = c_D_Acharya(Carreau, v_r_mag, d_p, Re, SR_f);
		c_D = c_D_Acharya(Carreau, v_r_mag, d_p, Re, SR);
	}
	else if (strcmp("CU",model_c_D) == 0) /* Chhabra & Uhlherr (1980) */
	{
		c_D_Fluent = c_D_ChhabraUhlherr(Carreau, Re0, SR_f);
		c_D = c_D_ChhabraUhlherr(Carreau, Re0, SR);
	}
	else /*  */
	{
	}

	/* Drag force as required by Fluent including potential non-Newtonian viscosity correction (first term)*/
	drag_force = eta/eta_Fluent * 18.0*c_D*Re/24.0;

	/* ----------------------------------------------------------------------------------------------------------- */
	/* Output in TUI */
	/* ----------------------------------------------------------------------------------------------------------- */	
	
	Message0("\n particle_drag_coefficient.c \
	\t x %f \
	\t y %f \
	\t SR_f %f \
	\t SR_p %f \
	\t SR %f \
	\t eta_F %f \
	\t eta %f \
	\t Re_Fluent %f \
	\t Re %f \
	\t c_D_Fluent %f \
	\t c_D %f \
	\
	\t \n", x[0], x[1], SR_f, SR_p, SR, eta_Fluent, eta, Re_Fluent, Re, c_D_Fluent, c_D);
	
	return (drag_force);
}


/*************************************************************************************************
UDF particle drag force for Newtonian fluids based on Schiller & Naumann (1933)
*************************************************************************************************/ 

DEFINE_DPM_DRAG(drag_force_n,Re,p)
{	
	
	/* Schiller & Naumann (1935) */
	real c_D;
	if (Re > 1000.0)
	{
		c_D = 0.44;
	}
	else
	{
		c_D = (24.0 / Re)*(1.0 + 0.15*pow(Re, 0.687));
	}
	
	/* Drag force as required by Fluent including potential non-Newtonian viscosity correction (first term)*/
	real drag_force = 18.0*c_D*Re/24.0;
	
	return (drag_force);
}
