/***********************************************************************
UDF for computing particle drag coefficient (18 Cd Re/24)
curve as suggested by R. Clift, J. R. Grace and M. E. Weber
"Bubbles, Drops, and Particles" (1978)
************************************************************************/

/* ANSYS Fluent 17.2 UDF Manual 2.5.3. */

/* DEFINE_DPM_DRAG is used to specify the drag between particles and fluid as a dimensionless
group (18 Cd Re/24) as it appears in the drag force per unit particle mass (eq. (2.19), page 183). */



#include "udf.h"
DEFINE_DPM_DRAG(particle_drag_force,Re,p)
{
real w, drag_force;
if (Re < 0.01)
{
drag_force=18.0;
return (drag_force);
}
else if (Re < 20.0)
{
w = log10(Re);
drag_force = 18.0 + 2.367*pow(Re,0.82-0.05*w) ;
return (drag_force);
}
else
/* Note: suggested valid range 20 < Re < 260 */
{
drag_force = 18.0 + 3.483*pow(Re,0.6305) ;
return (drag_force);
}
}