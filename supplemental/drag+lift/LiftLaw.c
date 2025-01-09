/******************************************************************
UDF for implementing the lift model of Tomiyama et al. (2002)
*******************************************************************/

/* ANSYS Fluent 17.2 UDF Manual 2.4.3.4. */

/*  */

#include "udf.h"
#include "flow.h"
#define Eo_l1 4.0
#define Eo_l2 10.0
DEFINE_EXCHANGE_PROPERTY(custom_lift,c,t,i,j)
{
/* i -- liquid-phase; j -- vapor-phase */
Thread **pt = THREAD_SUB_THREADS(t);
real v_x=0., v_y=0., v_z=0.;
real vel, Rev, Eo, d2, T_sfc, sigma;
real lift_coeff, lift_co, wk_co;;
real diam = C_PHASE_DIAMETER(c,pt[j]);
real rho_v = C_R(c,pt[j]);
real rho_l = C_R(c,pt[i]);
real mu_l = C_MU_L(c,pt[i]);
real gravity = NV_MAG(M_gravity);
Property *prop =
DOMAIN_COMPLEX_PROP_PROPERTY(DOMAIN_INTERACTION(root_domain),
COMPLEX_PROP_sfc_tension_coeff,
i,j);
T_sfc = (sg_temperature && NNULLP(THREAD_STORAGE(pt[i],SV_T)))? C_T(c,pt[i]) : T_REF;
if(prop == NULL || PROPERTY_METHOD(prop,0) == PROP_METHOD_NONE)
Error("Lift-Tomiyama: Please set value for surface tension !");
sigma = generic_property(c,t,prop,(Property_ID)0,T_sfc);
if(sigma <= 0.) Error("Lift-Tomiyama: Please set nonzero value for surface tension !");
/* calculate bubble Reynolds Number */
v_x = C_U(c,pt[j]) - C_U(c,pt[i]);
v_y = C_V(c,pt[j]) - C_V(c,pt[i]);
#if RP_3D
v_z = C_W(c,pt[j]) - C_W(c,pt[i]);
#endif
vel = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
Rev = RE_NUMBER(rho_l,vel,diam,mu_l);
d2 = diam*diam;
Eo = gravity*(rho_l-rho_v)*d2/sigma;
if (Eo <= Eo_l1)
wk_co = 0.0;
else if (Eo < Eo_l2)
wk_co = -0.096*Eo + 0.384;
else
wk_co = -0.576;
lift_co = 0.288*tanh(0.121*Rev);
lift_coeff = lift_co + wk_co;
return lift_coeff;
}