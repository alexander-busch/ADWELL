/******************************************************************
UDF for customizing the default Syamlal drag law in ANSYS Fluent
*******************************************************************/

/* ANSYS Fluent 17.2 UDF Manual 2.4.3.3. */

/* The following UDF, named custom_drag, can be used to customize the default Syamlal drag law in
ANSYS Fluent. The default drag law uses 0.8 (for void < 0.85) and 2.65 (void>0.85) for bfac. This results
in a minimum fluid velocity of 25 cm/s. The UDF modifies the drag law to result in a minimum fluid
velocity of 8 cm/s, using 0.28 and 9.07 for the bfac parameters. */

#include "udf.h"
#define pi 4.*atan(1.)
#define diam2 3.e-4

DEFINE_EXCHANGE_PROPERTY(custom_drag,cell,mix_thread,s_col,f_col)
{
Thread *thread_g, *thread_s;
real x_vel_g, x_vel_s, y_vel_g, y_vel_s, abs_v, slip_x, slip_y,
rho_g, rho_s, mu_g, reyp, afac,
bfac, void_g, vfac, fdrgs, taup, k_g_s;
/* find the threads for the gas (primary) */
/* and solids (secondary phases) */
thread_g = THREAD_SUB_THREAD(mix_thread, s_col);/* gas phase */
thread_s = THREAD_SUB_THREAD(mix_thread, f_col);/* solid phase*/
/* find phase velocities and properties*/
x_vel_g = C_U(cell, thread_g);
y_vel_g = C_V(cell, thread_g);
x_vel_s = C_U(cell, thread_s);
y_vel_s = C_V(cell, thread_s);
slip_x = x_vel_g - x_vel_s;
slip_y = y_vel_g - y_vel_s;
rho_g = C_R(cell, thread_g); rho_s = C_R(cell, thread_s);
mu_g = C_MU_L(cell, thread_g);
/*compute slip*/
abs_v = sqrt(slip_x*slip_x + slip_y*slip_y);
/*compute Reynolds number*/
reyp = rho_g*abs_v*diam2/mu_g;
/* compute particle relaxation time */
taup = rho_s*diam2*diam2/18./mu_g;
void_g = C_VOF(cell, thread_g);/* gas vol frac*/
/*compute drag and return drag coeff, k_g_s*/
afac = pow(void_g,4.14);
if(void_g<=0.85)
bfac = 0.281632*pow(void_g, 1.28);
else
bfac = pow(void_g, 9.076960);
vfac = 0.5*(afac-0.06*reyp+sqrt(0.0036*reyp*reyp+0.12*reyp*(2.*bfacafac)+
afac*afac));
fdrgs = void_g*(pow((0.63*sqrt(reyp)/
vfac+4.8*sqrt(vfac)/vfac),2))/24.0;
k_g_s = (1.-void_g)*rho_s*fdrgs/taup;
return k_g_s;
}