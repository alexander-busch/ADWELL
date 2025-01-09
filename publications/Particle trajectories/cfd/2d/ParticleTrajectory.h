/* Particle Trajectory header file */




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

/* Case geometry */
#define h (RP_Get_Real("h_channel")) /* Inlet height [m] */
#define h_0 (0) /* Initial bed height [m] */

/* Physical constants */
#define gravity (9.81) /* gravitational acceleration */

/* Fluid properties */
#define u_bulk (RP_Get_Real("u_bulk")) /* Bulk velocity [m/s] */
#define eta_h2o (0.001002) /* Dynamic (shear) viscosity of water */
#define rho_f (1000.0) /* density water */

/* Rheological properties of the fluid (Cross model) */
/* #define lambda (RP_Get_Real("pac_cross_lambda")); */
/* #define n (RP_Get_Real("pac_cross_n")); */
/* #define mu_0 (RP_Get_Real("pac_cross_mu0")); */
/* #define mu_inf (eta_h2o); */

/* Solid properties */
#define d_s (0.0003) /* particle diameter */
#define rho_s (2650.0) /* solid density */
#define s (rho_s/rho_f) /* Relative density */



/*************************************************************************************************
Forward declarations
*************************************************************************************************/ 

/* struct PowerLaw;
struct FourParameterModel; 
*/


/*************************************************************************************************
Structures
*************************************************************************************************/ 

struct PowerLaw {
  real K, n;
};

struct FourParameterModel {
    real lambda, n, mu_0, mu_inf;
};

/*Declares a PowerLaw structure
typedef struct PowerLaw Struct_PowerLaw;
typedef struct PowerLaw Struct_Carreau;
*/

/*************************************************************************************************
Functions
*************************************************************************************************/ 
/*
Cross2PL(real, real, real, real, real);
Carreau2PL(struct FourParameterModel, real);
*/

extern struct PowerLaw Cross2PL(real lambda, real n, real mu_0, real mu_inf, real SR);
extern struct PowerLaw Carreau2PL(struct FourParameterModel Carreau, real SR);

extern real ViscosityCross(struct FourParameterModel Cross, real SR);
extern real ViscosityCarreau(struct FourParameterModel Carreau, real SR);
