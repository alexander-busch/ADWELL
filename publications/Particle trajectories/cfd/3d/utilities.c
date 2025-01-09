#include "ParticleTrajectory.h"



/* Compute PL coefficients from Cross model for current shear rate */
struct PowerLaw Cross2PL(real lambda, real n, real mu_0, real mu_inf, real SR)
{
	struct PowerLaw PL;
	
    PL.K = (pow(SR,n*(mu_inf*pow(lambda,2.0*n)*pow(SR,2.0*n) \
		+2.0*mu_0*pow(lambda,n)*pow(SR,n)+mu_0)/(mu_inf*pow(lambda,2.0*n) \
		*pow(SR,2.0*n)+pow(SR,n)*(mu_0+mu_inf)*pow(lambda,n)+mu_0)) \
		*mu_inf*pow(lambda,n) \
		+ pow(SR,pow(lambda,n)*pow(SR,n)*n*(mu_0-mu_inf) \
		/(mu_inf*pow(lambda,2.0*n)*pow(SR,2.0*n)+pow(SR,n) \
		*(mu_0+mu_inf)*pow(lambda,n)+mu_0))*mu_0) \
		/(1.0+pow(lambda,n)*pow(SR,n));
		
    PL.n = (mu_inf*pow(lambda,2.0*n)*pow(SR,2.0*n)-((-n-1.0) \
		*mu_inf+mu_0*(n-1.0))*pow(SR,n)*pow(lambda,n)+mu_0) \
		/(mu_inf*pow(lambda,2.0*n)*pow(SR,2.0*n)+pow(SR,n) \
		*(mu_0+mu_inf)*pow(lambda,n)+mu_0);

	return PL;
}


/* Compute PL coefficients from Carreau model for current shear rate */
struct PowerLaw Carreau2PL(struct FourParameterModel Carreau, real SR)
{
	struct PowerLaw PL;
	
	/* PL coefficients based on current shear rate */
	if (SR < (1.0/Carreau.lambda)) /* Newtonian zero-shear viscosity region */
	{
		PL.K = Carreau.mu_0;
		PL.n = 1.0;
	}
	else /* Shear-thinning region */
	{
		PL.K = Carreau.mu_0*pow(Carreau.lambda,(Carreau.n-1.0));
		PL.n = Carreau.n;
	}

	return PL;
}


/* Compute viscosity from Cross model for current shear rate */
real ViscosityCross(struct FourParameterModel Cross, real SR)
{
	real eta = Cross.mu_inf+(Cross.mu_0-Cross.mu_inf)/(1.0+pow(Cross.lambda*SR, Cross.n));
	
	return eta;
}


/* Compute viscosity from Carreau model for current shear rate */
real ViscosityCarreau(struct FourParameterModel Carreau, real SR)
{
	real eta = Carreau.mu_inf+(Carreau.mu_0-Carreau.mu_inf)*pow(1.0+pow(Carreau.lambda*SR, 2.0), (Carreau.n-1.0)/2.0);
	
	return eta;
}