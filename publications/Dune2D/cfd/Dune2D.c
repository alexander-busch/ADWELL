
/* Dune2D */

/* Fluent */
#include "udf.h"
#include "unsteady.h"
#include "mem.h"
#include "dynamesh_tools.h"
#include "sg.h"

/* Zone Identifiers */
#define BC_movingbed_ID (8) /* ID of MovingBed boundary */

/* Operating point */
#define FluidBulkVel (0.81) /* Fluid bulk velocity */
#define InitialBedHeight (0.005) /* Initial bed height */

/* Case geometry */
#define ChannelHeight (0.04) /* inlet height - TODO compute from mesh and make dynamic for inflow profile */
#define A_Pipe (M_PI*ChannelHeight*ChannelHeight/4.0) /* Pipe cross-sectional area */
#define Q_f (FluidBulkVel*A_Pipe) /* Volumetric fluid flow rate */
#define delta_z (1.0) /* z-spacing - Unit length as 2D */ 

/* Physical constants */
#define gravity (9.81) /* gravitational acceleration */

/* Fluid properties */
#define rho_f (1000.0) /* density water */
#define sigma_t (0.9) /* turbulent Prandt-Schmidt number */

/* Solid properties */
#define d_s (0.0012	) /* particle diameter */
#define k_s (2.5*d_s) /* sediment roughness height */
#define rho_s (2500.0) /* solid density */
#define phi (30.0 * M_PI / 180.0) /* Angle of repose */
#define alpha_s_bed (0.55) /* Solids fraction in bed load layer */
#define alpha_f_bed (1-alpha_s_bed) /* Fluids fraction in bed load layer, "Bed porosity", "Void fraction" */
#define s (rho_s/rho_f) /* Relative density */

/* Bed load model */
#define GravitationalStress (rho_f * (s - 1.0) * gravity * d_s) /* Denominator of Shields number */
#define C_q_0 (12.0 / ( 1.0 - alpha_f_bed)) /* Nominator is 12 in Tron, why? */
#define C_grad_h (1.5)

/* UDM enumerators */
enum UDMnames { WSS, BLTR, ECT, EST, SMF, SSRC, FSRC, NBH, DHDX, MBWF, SSV, N_REQ_UDM };
/* 		Description
WSS		Wall shear stress on "moving bed"
BLTR	Bed load transport rate
ECT		Exner convective term
EST		Exner source term
SMF		Net solid volumetric flux from "moving bed" wall into wall-adjacent grid cell
SSRC	Solid mass source in wall-adjacent grid cell
FSRC	Fluid mass source in wall-adjacent grid cell
NBH		New bed height
DHDX	Geometrical gradient of "moving bed" face
MBWF	Cell flag <--> Cell is moving bed wall-adjacent
SSV		Settling velocity of solids in "moving bed" wall-adjacent grid cell
*/


/* Global variables */
real delta_x; /* x-spacing */
real U;  /* Current fluid bulk velocity */

/* Function declaration */
real CompFaceGradient(face_t, Thread *);
real CompDeltaX(face_t, Thread *);
real CompSetVel(real, real, real);
void UpdateUDM(enum UDMnames, face_t, Thread *, real);
void WriteUDMtoTextFile(real, Thread *);
real CompShearStress_mag(face_t, Thread *);
real CompBedLoadTranRate(real, real, face_t, Thread *);
real CompDepositionRate(face_t, Thread *);
real CompEntrainmentRate(face_t, Thread *);




/* ******************************************************************************************** */
/* Various functions */
/* ******************************************************************************************** */

/* Compute face gradient dh/dx for current face */
real CompFaceGradient(face_t face, Thread *mix_face_thread)
{
	real x_node[2], y_node[2], dx, dy, dhdx;
	int n;
	Node *node_p;

	/* Get x- and y-coordinates of face nodes */
	f_node_loop(face, mix_face_thread, n)
	{
		node_p = F_NODE(face, mix_face_thread, n);
		x_node[n] = NODE_X(node_p); /* x-coordinates of nodes */
		y_node[n] = NODE_Y(node_p); /* y-coordinates of nodes */
	}

	/* Compute dx and dy based on which node comes first in array */
	if (x_node[0] < x_node[1])
	{
		dx = x_node[1] - x_node[0];
		dy = y_node[1] - y_node[0];
	}
	else
	{
		dx = x_node[0] - x_node[1];
		dy = y_node[0] - y_node[1];
	}


	dhdx = (y_node[1] - y_node[0]) / (x_node[1] - x_node[0]);
	dhdx = (y_node[0] - y_node[1]) / (x_node[0] - x_node[1]);


	/* Compute gradient */
	dhdx = dy / dx;

	/* ????????????????????????????????????????????????????????????????????????? */

	/* TODO Check if by using the following macros the face gradient is obtained more conveniently */
	real A[ND_ND];		/* area normal vector, boundary face area normals always point out of the domain */
	real ds;			/* distance between the cell centroid and the face centroid */
	real es[ND_ND]; 	/* unit normal vector in the direction from centroid of cell c0 to the face centroid */
	real A_by_es;		/* Projected area */
	real dr0[ND_ND];	/* value vector that connects the centroid of cell c0 to the face centroid */

	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);


	/* Mixture thread values */
	BOUNDARY_FACE_GEOMETRY(face, mix_face_thread, A, ds, es, A_by_es, dr0);

	/* Fluid thread values */
	BOUNDARY_FACE_GEOMETRY(face, flu_face_thread, A, ds, es, A_by_es, dr0);

	/* Solid thread values */
	BOUNDARY_FACE_GEOMETRY(face, sol_face_thread, A, ds, es, A_by_es, dr0);

	/* ????????????????????????????????????????????????????????????????????????? */


	return dhdx;
}


/* Compute delta x for current face */
real CompDeltaX(face_t face, Thread *mix_face_thread)
{
	real x_node[ND_ND], dx;
	int n;
	Node   *node_p;

	/* Loop all face nodes and get node coordinates */
	f_node_loop(face, mix_face_thread, n)
	{
		node_p = F_NODE(face, mix_face_thread, n);
		x_node[n] = NODE_X(node_p); /* Get x-coordinates of nodes */
	}

	/* Compute delta_x */
	dx = fabs(x_node[1] - x_node[0]);

	return dx;
}


/* Re-compute solids settling velocity */
real CompSetVel(real v_set, real mu_eff, real gamma_dot_global)
{
	real Re_p, c_D, gamma_dot_eff, mu_eff2;

	/* Effective strain rate */
	/* gamma_dot_eff = sqrt(pow(gamma_dot_global, 2.0) + pow(2.0*v_set / d_s, 2.0)); */

	/* Effective viscosity based on effective strain rate */
	/* mu_eff2 = (0.0026 - 0.001) / (1.0 + pow(0.008*gamma_dot_eff, 0.37)) + 0.001; */
	
	/* Particle Reynolds number */
	/* Re_p = rho_f*v_set*d_s / mu_eff2; */
	
	/* Particle Reynolds number */
	Re_p = rho_f*v_set*d_s / mu_eff;

	/* Coefficient of drag (Schiller-Naumann) */
	c_D = (24.0 / Re_p)*(1.0 + 0.15*pow(Re_p, 0.687));

	/* Settling velocity */
	v_set = sqrt(4.0*d_s*gravity*(s - 1.0) / (3.0*c_D));

	return v_set;
}


/* Update F_UDMI and C_UDMI for current Enumerator value from current face */
void UpdateUDM(enum UDMnames Enumerator, face_t face, Thread *mix_face_thread, real new_value)
{
	/* Wall adjacent mixture cell index */
	cell_t mix_c0 = F_C0(face, mix_face_thread);

	/* Wall adjacent mixture cell Thread */
	Thread *mix_ct0 = THREAD_T0(mix_face_thread);

	/* Save to face UDM */
	F_UDMI(face, mix_face_thread, Enumerator) = new_value;

	/* Save to cell UDM */
	C_UDMI(mix_c0, mix_ct0, Enumerator) = new_value;
}


/* Write all UDM values to text file*/
void WriteUDMtoTextFile(real flow_time, Thread *mix_face_thread)
{
	/* Local variables */
	FILE *fp;
	face_t face;
	real face_pos[ND_ND];
	char *line = "--------------------------------------------------------"
		"--------------------------------------------------------"
		"--------------------------------------------------------"
		"--------------------------------------------------------"
		"--------------------------------------------------------"
		"\n";
	char *heading = "Time\t\t"
		"Face#\t"
		"x_face\t\t"
		"y_face = h_b\t"
		"tau_w\t\t"
		"q_b\t\t"
		"Convective\t"
		"Source\t\t"
		"E + D\t\t"
		"mass_src_s\t"
		"mass_src_f\t"
		"h_b_n\t \t"
		"DHDX\t\t"
		"vel_set_s \n";
	char *definitions = "%f\t"
		"%u\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"%f\t"
		"\n";
	char filename[100];
	sprintf(filename, "./udm-values/UDM_values_t=%f.txt", flow_time);

	fp = fopen(filename, "a");
	fprintf(fp, line);
	fprintf(fp, "t = %f\n", flow_time);
	fprintf(fp, line);
	fprintf(fp, heading);
	fprintf(fp, line);

	begin_f_loop(face, mix_face_thread)
	{
		/* Get face coordinates of current face = Initial face positions */
		F_CENTROID(face_pos, face, mix_face_thread);

		fprintf(fp, definitions, \
			flow_time, \
			face, \
			face_pos[0], \
			face_pos[1], \
			F_UDMI(face, mix_face_thread, WSS), \
			F_UDMI(face, mix_face_thread, BLTR), \
			F_UDMI(face, mix_face_thread, ECT), \
			F_UDMI(face, mix_face_thread, EST), \
			F_UDMI(face, mix_face_thread, SMF), \
			F_UDMI(face, mix_face_thread, SSRC), \
			F_UDMI(face, mix_face_thread, FSRC), \
			F_UDMI(face, mix_face_thread, NBH), \
			F_UDMI(face, mix_face_thread, DHDX), \
			F_UDMI(face, mix_face_thread, SSV));
	}
	end_f_loop(face, mix_face_thread)

		fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fclose(fp);
}



/* ******************************************************************************************** */
/* Boundary and initial conditions */
/* ******************************************************************************************** */

/* Parabolic velocity profile at inlet */
DEFINE_PROFILE(parabolic_velocity_profile, phase_face_thread, position)
{

	/* This UDF defines a custom velocity profile for the inlet boundary zone and is hooked to 	the appropriate velocity phase in Fluent in the relevant boundary condition dialog box.
	Appropriate phase variables will be passed to the function by the solver at run time. See UDF manual page 299 */

	/* Get mixture domain */
	Domain *mix_domain = Get_Domain(1);

	/* Get mixture face thread of "moving bed" */
	Thread *mix_face_thread = Lookup_Thread(mix_domain, BC_movingbed_ID);

	/* Local variables  */
	real	face_pos[ND_ND]; /* Face position vector */
	Node	*node_p;
	face_t	face = 298; /* First face */
	real	x_face, y_face, x_node, y_node; /* Face and node x- & y-coordinates */
	real	u, h, H, A_Bed, h_bed_inlet;

	real	current_time = CURRENT_TIME;
	int n;

	/* Get current bed height at inlet  */
	F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates of first cell */
	x_face = face_pos[0]; /* Get x-component of current face */

	f_node_loop(face, mix_face_thread, n)
	{
		node_p = F_NODE(face, mix_face_thread, n);
		x_node = NODE_X(node_p); /* x-coordinate of node */
		y_node = NODE_Y(node_p); /* y-coordinate of node */

		/* Node is upstream <--> inlet node */
		if (x_node < x_face) { h_bed_inlet = y_node; }
	}

	H = ChannelHeight;
	h = h_bed_inlet;

	/* Current cross-sectional area of bed based on pipe cross-sectional area */
	A_Bed = (pow((H / 2.0), 2.0)*acos(1.0 - 2.0*h / H) - (H / 2.0 - h_bed_inlet)*pow((H*h - h*h), 0.5));

	/* Current fluid bulk velocity based on existing cross-sectional flow area */
	U = Q_f / (A_Pipe - A_Bed);

	/* Loop all boundary faces and compute corresponding x-velocity */
	begin_f_loop(face, phase_face_thread)
	{
		/* Get face coordinates */
		F_CENTROID(face_pos, face, phase_face_thread);
		y_face = face_pos[1]; /* Get y-component of current face */

		/* x-velocity in current cell */
		u = 2.0*U - 2.0*U / pow(((h - H) / 2.0), 2.0)*pow((y_face - (h + (H - h) / 2.0)), 2.0);

		/* Assign to current face */
		F_PROFILE(face, phase_face_thread, position) = u;
	}
	end_f_loop(face, phase_face_thread)

}

/* Solid mass source term for "moving bed"-adjacent grid cell */
DEFINE_SOURCE(sol_mass_SRC, cell, sol_cell_thread, dS, eqn)
{
	/* This macro will loop through all cells of the domain and will get passed the pointer to the cell thread of the solid phase as it is hooked to the solid mass transport equation? */

	/* Solid phase cell thread --> Mixture cell Thread */
	Thread *mix_cell_thread = THREAD_SUPER_THREAD(sol_cell_thread);

	real source;

	/* Check if the current cell is "moving bed"-adjacent */
	if (C_UDMI(cell, mix_cell_thread, MBWF) == 1.0)
	{
		/* Compute solid mass source with units [kg/m³/s] */
		source = C_UDMI(cell, mix_cell_thread, SMF)*rho_s*delta_x*delta_z / C_VOLUME(cell, mix_cell_thread); /* Solids fraction in bed load layer times total source */

		/* Update cell UDM */
		C_UDMI(cell, mix_cell_thread, SSRC) = source;
	}
	else
	{
		/* Cell is not "moving bed"-adjacent <--> Source term is zero */
		source = 0.0;
		dS[eqn] = 0.0;
	}
	return source;
}

/* Fluid mass source term for "moving bed"-adjacent grid cell */
DEFINE_SOURCE(flu_mass_SRC, cell, flu_cell_thread, dS, eqn)
{
	/* This macro will loop through all cells of the domain and will get passed the pointer to the cell thread of the fluid phase as it is hooked to the fluid mass transport equation? */

	/* Fluid phase cell thread --> Mixture cell Thread */
	Thread *mix_cell_thread = THREAD_SUPER_THREAD(flu_cell_thread);

	real source;

	/* Check if the current cell is "moving bed"-adjacent */
	if (C_UDMI(cell, mix_cell_thread, MBWF) == 1.0)
	{
		/* Compute fluid mass source with units [kg/m³/s] */
		source = alpha_f_bed / alpha_s_bed*C_UDMI(cell, mix_cell_thread, SMF)*rho_f*delta_x*delta_z / C_VOLUME(cell, mix_cell_thread);

		/* Update cell UDM */
		C_UDMI(cell, mix_cell_thread, FSRC) = source;
	}

	else
	{
		/* Cell is not "moving bed"-adjacent <--> Source term is zero */
		source = 0.0;
		dS[eqn] = 0.0;
	}
	return source;
}

/* Initialisation */
DEFINE_INIT(initialisation, mix_domain)
{
	/* DEFINE_INIT always refers to the mixture domain.
	So far, only the UDMs are initialized, this should be phase-independent, therefore a loop in the mixture domain should be sufficent. */


	/* Mixture face thread of "moving bed" wall */
	Thread *mix_face_thread = Lookup_Thread(mix_domain, BC_movingbed_ID);

	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);

	/* Local variables */
	cell_t cell, mix_c0;
	face_t face;
	Node   *node_p;
	Thread *mix_cell_thread, *mix_ct0;

	real face_pos[ND_ND], x_node[2], mu_eff, v_set;
	real flow_time = CURRENT_TIME;

	/* Check if sufficient amopunt of UDMs are available */
	if (N_UDM < N_REQ_UDM)
	{
		Message("Number of UDM's : %d\n", N_UDM);
		Message("Number of used UDM's : %d\n", N_REQ_UDM);
		Error("ERROR : Number of UDM's too small");
	}

	/* Loop over all cell threads in the mixture domain */
	thread_loop_c(mix_cell_thread, mix_domain)
	{
		/* Loop over all cells in current cell thread and sets UDMs to zero */
		begin_c_loop_all(cell, mix_cell_thread)
		{
			/* Loop through enumerations of UDMs */
			for (int i = WSS; i<N_REQ_UDM; i++)
			{
				C_UDMI(cell, mix_cell_thread, i) = 0.0;
			}
		}
		end_c_loop_all(cell, mix_cell_thread)
	}

	/* Loop over all faces of "moving bed" wall */
	begin_f_loop(face, mix_face_thread)
	{
		/* Get face coordinates of current face = Initial face positions */
		F_CENTROID(face_pos, face, mix_face_thread);

		/* Loop through enumerations of UDMs */
		for (int i = WSS; i < N_REQ_UDM; i++)
		{
			if (i == NBH)
			{
				F_UDMI(face, mix_face_thread, i) = face_pos[1];
			}
			else if (i == SSV)
			{
				F_UDMI(face, mix_face_thread, i) = v_set;
			}
			else if (i == MBWF)
			{
				F_UDMI(face, mix_face_thread, i) = 1.0;
			}
			else
			{
				F_UDMI(face, mix_face_thread, i) = 0.0;
			}
		}

		/* Wall adjacent mixture cell index & Thread */
		mix_c0 = F_C0(face, mix_face_thread); /* Get cell of current face */
		mix_ct0 = THREAD_T0(mix_face_thread); /* Get cell thread of current face thread*/

		/* Wall adjacent fluid cell index & Thread */
		cell_t flu_c0 = F_C0(face, flu_face_thread);
		Thread *flu_ct0 = THREAD_T0(flu_face_thread);

		/* Effective viscosity in wall-adjacent grid cell */
		mu_eff = C_MU_EFF(flu_c0, flu_ct0);

		/* Compute Stokes settling velocity */
		v_set = pow(d_s, 2.0)*(rho_s - rho_f)*gravity / (18.0*mu_eff);

		/* Set wall-adjacent cell UDMs */
		C_UDMI(mix_c0, mix_ct0, NBH) = face_pos[1]; /* bed height = face y-component */
		C_UDMI(mix_c0, mix_ct0, MBWF) = 1.0; /* Set flag <--> This cell is a "moving bed"-adjacent cell */
		C_UDMI(mix_c0, mix_ct0, SSV) = v_set; /* Stokes settling velocity */

		/* Compute mesh delta_x from first face */
		if (face == 298) { delta_x = CompDeltaX(face, mix_face_thread); }
	}
	end_f_loop(face, mix_face_thread)

		/* Write current values of all UDMs to text file */
		WriteUDMtoTextFile(flow_time, mix_face_thread);
}

/* Rename UDMs */
DEFINE_ON_DEMAND(rename_UDMs)
{
	Set_User_Memory_Name(0, "UDM0: Shear stress acting on bed");
	Set_User_Memory_Name(1, "UDM1: Bed load transport rate");
	Set_User_Memory_Name(2, "UDM2: Exner convective term");
	Set_User_Memory_Name(3, "UDM3: Exner source term");
	Set_User_Memory_Name(4, "UDM4: Net solid volumetric flux from bed load layer into wall-adjacent grid cell");
	Set_User_Memory_Name(5, "UDM5: Solid mass source in wall-adjacent grid cell");
	Set_User_Memory_Name(6, "UDM6: Fluid mass source in wall-adjacent grid cell");
	Set_User_Memory_Name(7, "UDM7: Bed height");
	Set_User_Memory_Name(8, "UDM8: Geometrical gradient of bed");
	Set_User_Memory_Name(9, "UDM9: Cell flag <--> Cell is moving bed wall-adjacent");
	Set_User_Memory_Name(10, "UDM10: Settling velocity of solids in bed load layer wall-adjacent grid cell");
}

/* ******************************************************************************************** */
/* Exner equation */
/* ******************************************************************************************** */

/* Compute shear stress for current face */
real CompShearStress_mag(face_t face, Thread *flu_face_thread)
{
	real ShearForce_mag, ShearStress_mag, area;
	real FC[ND_ND], A[ND_ND];

	/* Get face area vector */
	F_AREA(A, face, flu_face_thread);

	/* Get face center position */
	F_CENTROID(FC, face, flu_face_thread);

	/* Get face area */
	area = NV_MAG(A);

	/* Compute shear force magnitude */
	ShearForce_mag = NV_MAG(F_STORAGE_R_N3V(face, flu_face_thread, SV_WALL_SHEAR));

	/* Compute shear stress magnitude based on face area */
	ShearStress_mag = ShearForce_mag / area;

	return ShearStress_mag;
}

/* Compute volumetric bed load transport rate for current face */
real CompBedLoadTranRate(real teta, real teta_cr, face_t face, Thread *mix_face_thread)
{
	cell_t flu_c0;
	Thread *flu_ct0;
	real q_x0, q_x, u_x;

	/* Get x-component of fluid velocity of wall-adjacent grid cell */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	flu_c0 = F_C0(face, flu_face_thread);
	flu_ct0 = THREAD_T0(flu_face_thread);
	u_x = C_U(flu_c0, flu_ct0);

	/* Bed load transport rate for horizontal bed */
	q_x0 = C_q_0 * sqrt(gravity * (s - 1.0) * pow(d_s, 3.0) * teta) * (teta - teta_cr);

	/* Bed load transport rate for inclined bed*/
	q_x = q_x0 * (u_x / fabs(u_x) - C_grad_h * F_UDMI(face, mix_face_thread, DHDX)); /* TODO: What happens if dhdx is larger than 2/3? */

	return q_x;
}

/* Compute volumetric bed load transport rate incl. all variables for upstream face */
real CompBedLoadTranRate_us(face_t face_us, Thread *mix_face_thread)
{
	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);

	/* Wall adjacent fluid cell index & Thread */
	cell_t flu_c0 = F_C0(face_us, flu_face_thread);
	Thread *flu_ct0 = THREAD_T0(flu_face_thread);

	/* Local variables */
	real face_pos[ND_ND];
	real dhdx, beta, ShearStress_mag, teta, angle, mu_eff, dstar, teta_cr0, teta_cr, q_x_us;
	real y_face, current_bed_height;

	/* Geometrical gradient of current face */
	dhdx = CompFaceGradient(face_us, mix_face_thread);

	/* Bed slope */
	beta = atan(dhdx);

	/* Wall shear stress magnitude acting on current (sloped) face */
	ShearStress_mag = CompShearStress_mag(face_us, flu_face_thread);

	/* Dimensionless Shields number for sloped face */
	teta = ShearStress_mag / GravitationalStress;

	/* Effective viscosity in wall-adjacent grid cell */
	mu_eff = C_MU_EFF(flu_c0, flu_ct0);

	/* Dimensionless particle diameter */
	dstar = d_s*pow((s - 1.0)*gravity / pow(mu_eff / rho_f, 2.0), 0.3333333333333);

	/* Critical Shields number for horizontal bed */
	teta_cr0 = (0.24 / dstar) + 0.055 * (1.0 - exp(-0.02 * dstar));

	/* Critical Shields number with slope correction for inclined bed */
	angle = phi + beta;
	if (angle < 0.0) angle = 0.0; /* Limitation concept taken from Tron - TODO: Check with Export */
	if (angle > 0.5*M_PI) angle = 0.5 * M_PI; /* Limitation concept taken from Tron - TODO: Check with Export */
	teta_cr = teta_cr0 * sin(angle) / sin(phi);

	/* Get current bed height - TODO check approach from UDM vs from face*/
	F_CENTROID(face_pos, face_us, mix_face_thread); /* Get face coordinate of upstream cell */
	y_face = face_pos[1]; /* Get y-component of upstream face from face coordinate*/
	current_bed_height = F_UDMI(face_us, mix_face_thread, NBH); /* Get y-component of upstream face from UDM */

	/* Bed load transport rate per unit bed width of current face */
	/* Bed load transport only occurs */
	/* - if the bed shear stress is above the critical shear stress <--> teta > teta_c */
	/* - and if a bed exists <--> current_bed_height > 0 */
	if ((teta > teta_cr) && (current_bed_height > 0.0))
	{
		q_x_us = CompBedLoadTranRate(teta, teta_cr, face_us, mix_face_thread);
	}
	else
	{
		q_x_us = 0.0;
	}

	return q_x_us;
}

/* Compute deposition flux for current face */
real CompDepositionRate(face_t face, Thread *mix_face_thread)
{
	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);

	/* Wall adjacent mixture cell index & Thread */
	cell_t mix_c0 = F_C0(face, mix_face_thread);
	Thread *mix_ct0 = THREAD_T0(mix_face_thread);

	/* Wall adjacent fluid cell index & Thread */
	cell_t flu_c0 = F_C0(face, flu_face_thread);
	Thread *flu_ct0 = THREAD_T0(flu_face_thread);

	/* Wall adjacent solid cell index & Thread */
	cell_t sol_c0 = F_C0(face, sol_face_thread);
	Thread *sol_ct0 = THREAD_T0(sol_face_thread);

	/* Get flow current time */
	real current_time = CURRENT_TIME;

	/* Get required cell values */
	real V = C_V(sol_c0, sol_ct0); /* Solid v-velocity, TODO compare with v_set */
	real alpha_f = C_VOF(flu_c0, flu_ct0); /* Fluid volume fraction */
	real alpha_s = C_VOF(sol_c0, sol_ct0); /* Solid volume fraction */
	real mu_eff = C_MU_EFF(flu_c0, flu_ct0); /* Effective viscosity */
	real gamma_dot = C_STRAIN_RATE_MAG(flu_c0, flu_ct0); /* Strain rate magnitude of background fluid */

	/* Local variables */
	real v_set, v_set_iter, Re_p, n, D;
	real delta_v_set = 1.0;

	/* Solids settling velocity guess from previous time step */
	v_set = C_UDMI(mix_c0, mix_ct0, SSV);

	/* Iteratively compute solids settling velocity */
	while (fabs(delta_v_set) > 0.01)
	{
		/* Re-compute solids settling velocity */
		v_set_iter = CompSetVel(v_set, mu_eff, gamma_dot);

		/* Difference of original and recomputed settling velocity */
		delta_v_set = v_set - v_set_iter;

		/* Overwrite settling velocity with re-computed settling velocity */
		v_set = v_set_iter;
	}

	/* Save current solids settling velocity to UDMs for next time step */
	UpdateUDM(SSV, face, mix_face_thread, v_set);

	/* Particle Reynolds number */
	Re_p = rho_f*v_set*d_s / mu_eff;

	/* Hindered settling exponent n of Zaki & Richardson */
	n = (0.27*pow(Re_p, 0.9) + 5.1) / (0.1*pow(Re_p, 0.9) + 1.0);

	/* Deposition rate incl. hindered settling */
	D = -alpha_s*(pow(alpha_f, n)*v_set);

	return D;
}

/* Compute entrainment flux for current face */
real CompEntrainmentRate(face_t face, Thread *mix_face_thread)
{
	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);

	/* Wall adjacent mixture cell index & Thread */
	cell_t mix_c0 = F_C0(face, mix_face_thread);
	Thread *mix_ct0 = THREAD_T0(mix_face_thread);

	/* Wall adjacent fluid cell index & Thread */
	cell_t flu_c0 = F_C0(face, flu_face_thread);
	Thread *flu_ct0 = THREAD_T0(flu_face_thread);

	/* Wall adjacent solid cell index & Thread */
	cell_t sol_c0 = F_C0(face, sol_face_thread);
	Thread *sol_ct0 = THREAD_T0(sol_face_thread);

	/* Local variables */
	real mu_lam, mu_turb, mu_eff, d_eff, d_turb;
	real x_cell[ND_ND], face_pos[ND_ND];
	real y1, y2, dalphady, E;

	/* Get required cell values */
	mu_lam = C_MU_L(flu_c0, flu_ct0); /* Laminar viscosity */
	if (rp_turb)
	{
		mu_turb = C_MU_T(flu_c0, flu_ct0); /* Turbulent viscosity */
		mu_eff = C_MU_EFF(flu_c0, flu_ct0); /* Effective viscosity */
	}
	real alpha_s = C_VOF(sol_c0, sol_ct0); /* Solid volume fraction */

	/* Get solid volume fraction gradient - TODO Why does it not work? */
	if (NULL != THREAD_STORAGE(sol_ct0, SV_VOF_G)) /* Checks if gradient exists, does not in first iteration */
	{
		dalphady = C_VOF_G(sol_c0, sol_ct0)[1]; /* Get solid volume fraction gradient at cell center */
	}
	else
	{
		dalphady = 0.0;
	}

	/* Laminar diffusivity */
	d_eff = mu_lam;

	if (rp_turb)
	{
		d_turb = mu_turb / sigma_t; /* Turbulent diffusivity */
		d_eff += d_turb;
	}

	/* Compute solid fraction gradient, from wall-adjacent grid cell center to" moving bed" */
	y1 = F_UDMI(face, mix_face_thread, NBH); /* Get y-component of face coordinate from UDMI */
	F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates */
	y1 = face_pos[1]; /* Get y-component of face coordinate from Face */
	C_CENTROID(x_cell, mix_c0, mix_ct0); /* Get cell center coordinates */
	y2 = x_cell[1]; /* Get y-component of cell coordinate from Face */
	dalphady = (alpha_s - alpha_s_bed) / (y2 - y1); /* Compute solid volume fraction gradient */

	/* Entrainment rate */
	E = -d_turb*dalphady / rho_s;

	return E;
}

/* Solve Exner equation with first-order upwind scheme */
DEFINE_EXECUTE_AT_END(solve_exner)
{
	/* Mixture domain and Thread */
	Domain *mix_domain = Get_Domain(1);
	Thread *mix_face_thread = Lookup_Thread(mix_domain, BC_movingbed_ID);

	/* Fluid and solid wall face threads */
	Thread *flu_face_thread = THREAD_SUB_THREAD(mix_face_thread, 0);
	Thread *sol_face_thread = THREAD_SUB_THREAD(mix_face_thread, 1);

	face_t face, face_us;

	real flow_time = CURRENT_TIME;
	real delta_t = CURRENT_TIMESTEP;
	real face_pos[ND_ND], A[ND_ND];
	real cellvolume, facearea;

	real dhdx, alpha, beta, teta, angle, mu_eff, dstar, teta_cr0, teta_cr;
	real tau_wall, ShearStress_mag, ShearStress_x, q_x, q_x_us, q_x_down, current_bed_height, new_bed_height, delta_h;
	real D, E, convectiveterm, sourceterm;

	/* NOTE: Fluent does not access the faces in the movingbed face loop in the correct order due to an odd face numbering concept.
	Check numbering with script "Ini_CheckCellNumbering_Standard.c" ,
	then implement and test correct stepwise treatment in "Ini_CheckCellNumbering_Stepwise.c".
	Current face IDs
	-------------------------------------------
	|     |     |     |     |     |     |     |
	| 298 |  0  |  1  | ... | 296 | 297 | 299 |
	|	  |     |     |     |     |     |     |
	-------------------------------------------
	*/

	/* Compute new bed height */
	begin_f_loop(face, mix_face_thread) /* Loop all faces of" moving bed" wall */
	{
		/* Wall adjacent mixture cell index & Thread */
		cell_t mix_c0 = F_C0(face, mix_face_thread);
		Thread *mix_ct0 = THREAD_T0(mix_face_thread);

		/* Wall adjacent fluid cell index & Thread */
		cell_t flu_c0 = F_C0(face, flu_face_thread);
		Thread *flu_ct0 = THREAD_T0(flu_face_thread);

		/* Wall adjacent solid cell index & Thread */
		cell_t sol_c0 = F_C0(face, sol_face_thread);
		Thread *sol_ct0 = THREAD_T0(sol_face_thread);


		/* Geometrical gradient of current face */
		dhdx = CompFaceGradient(face, mix_face_thread);

		/* Bed slope */
		beta = atan(dhdx);

		/* Wall shear stress magnitude acting on current (sloped) face */
		ShearStress_mag = CompShearStress_mag(face, flu_face_thread);

		/* Save Wall shear stress magnitude to UDMs */
		UpdateUDM(WSS, face, mix_face_thread, ShearStress_mag);

		/* Dimensionless Shields number for sloped face */
		teta = ShearStress_mag / GravitationalStress;

		/* Effective viscosity in wall-adjacent grid cell */
		mu_eff = C_MU_EFF(flu_c0, flu_ct0);

		/* Dimensionless particle diameter */
		dstar = d_s*pow((s - 1.0)*gravity / pow(mu_eff / rho_f, 2.0), 0.3333333333333);

		/* Critical Shields number for horizontal bed */
		teta_cr0 = (0.24 / dstar) + 0.055 * (1.0 - exp(-0.02 * dstar));

		/* Critical Shields number with slope correction for inclined bed */
		angle = phi + beta;
		if (angle < 0.0) angle = 0.0; /* Limitation concept taken from Tron - TODO: Check with Export */
		if (angle > 0.5*M_PI) angle = 0.5 * M_PI; /* Limitation concept taken from Tron - TODO: Check with Export */
		teta_cr = teta_cr0 * sin(angle) / sin(phi);

		/* Get current bed height - TODO check approach from UDM vs from face*/
		F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates of current cell */
		current_bed_height = face_pos[1]; /* Get y-component of current face from face coordinate*/
		current_bed_height = F_UDMI(face, mix_face_thread, NBH); /* Get y-component of current face from UDM */

		/* Entrainment rate */
		/* Entrainment does not occur */
		/* - if the critical shields number is not exceeded <--> teta <= teta_c */
		/* - or if the sediment bed is not existing  <--> current_bed_height <= 0.0 */
		if ((teta <= teta_cr) || (current_bed_height <= 0.0))
		{
			E = 0.0;
		}
		else
		{
			E = CompEntrainmentRate(face, mix_face_thread);
		}

		/* Deposition rate */
		D = CompDepositionRate(face, mix_face_thread);

		/* Save net volumetric solid flux to UDMs */
		UpdateUDM(SMF, face, mix_face_thread, E + D);

		/* EXNER sourceterm */
		sourceterm = delta_t / alpha_s_bed * (E + D);

		/* Save Exner sourceterm to UDMs */
		UpdateUDM(EST, face, mix_face_thread, sourceterm);




		/* ????????????????????????????????????????????????????????????????????????? */

		/* ToDo - Test face area and volume values */
		cellvolume = C_VOLUME(mix_c0, mix_ct0);
		F_AREA(A, face, mix_face_thread);
		facearea = NV_MAG(A);

		/* ????????????????????????????????????????????????????????????????????????? */

		/* Bed load transport rate per unit bed width of current face */
		/* Bed load transport only occurs */
		/* - if the bed shear stress is above the critical shear stress <--> teta > teta_c */
		/* - and if a bed exists <--> current_bed_height > 0 */
		if ((teta > teta_cr) && (current_bed_height > 0.0))
		{
			q_x = CompBedLoadTranRate(teta, teta_cr, face, mix_face_thread);
		}
		else
		{
			q_x = 0.0;
		}

		/* Save bed load transport rate to UDMs*/
		UpdateUDM(BLTR, face, mix_face_thread, q_x);


		/* Bed load transport rate  of upstream face */
		/* The upstream face depends on the direction of the bed load transport rate. */

		/* Check different face cases and determine upstream bed load transport rate */
		if (face == 0) /* Second face */
		{
			if (q_x > 0.0)
			{
				/* Positive fluid velocity direction --> face 298 */
				face_us = 298;
			}
			else
			{
				/* Negative fluid velocity direction --> face 298 */
				face_us = 1;
			}
			q_x_us = CompBedLoadTranRate_us(face_us, mix_face_thread);
		}
		else if (face == 298) /* First face */
		{
			if (q_x > 0.0)
			{
				/* Positive fluid velocity direction --> Ghost face */

				/* BC-Alternative 1: Upstream transport rate = current transport rate */
				q_x_us = q_x;

				/* BC-Alternative 2: Upstream transport rate = 0 */
				q_x_us = 0.0;

				/* BC-Alternative 3: Relate to wall shear stress of inlet velocity profile */
				/* current_bed_height = F_UDMI(face, mix_face_thread, NBH); */
				/* delta_h = ChannelHeight - current_bed_height; */
				/* tau_wall = mu_eff*8.0*FluidBulkVel / delta_h; */
				/* ... */
			}
			else
			{
				/* Negative fluid velocity direction --> face 0 */
				q_x_us = CompBedLoadTranRate_us(0, mix_face_thread);
			}
		}
		else if (face == 297) /* Third-last face */
		{
			if (q_x > 0.0)
			{
				/* Positive fluid velocity direction --> face 296 */
				face_us = face - 1;
			}
			else
			{
				/* Negative fluid velocity direction --> face 299 */
				face_us = face + 2;
			}
			q_x_us = CompBedLoadTranRate_us(face_us, mix_face_thread);
		}
		else if (face == 299) /* Last face */
		{
			if (q_x > 0.0)
			{
				/* Positive fluid velocity direction --> Ghost face */
				q_x_us = q_x;
			}
			else
			{
				/* Negative fluid velocity direction --> face 297 */
				q_x_us = CompBedLoadTranRate_us(297, mix_face_thread);
			}
		}
		else /* All other faces */
		{
			if (q_x > 0.0)
			{
				/* Positive fluid velocity direction --> face 296 */
				face_us = face - 1;
			}
			else
			{
				/* Negative fluid velocity direction --> face 299 */
				face_us = face + 1;
			}
			q_x_us = CompBedLoadTranRate_us(face_us, mix_face_thread);
		}

		/* Delta x of current face */
		delta_x = CompDeltaX(face, mix_face_thread);

		/* EXNER convective term */
		convectiveterm = delta_t / alpha_s_bed / delta_x*(q_x - q_x_us);

		/* Save convective term to UDMs*/
		UpdateUDM(ECT, face, mix_face_thread, convectiveterm);

		/* Compute new y-position of current face */
		new_bed_height = current_bed_height - convectiveterm - sourceterm;

		/* Save new y-position of current face to UDM */
		if (new_bed_height > 0.0)
		{
			UpdateUDM(NBH, face, mix_face_thread, new_bed_height);
		}
		else
		{
			UpdateUDM(NBH, face, mix_face_thread, 0.0);
		}
	}
	end_f_loop(face, mix_face_thread)

		/* Write current values of all UDMs to text file */
		WriteUDMtoTextFile(flow_time, mix_face_thread);
}



/* ******************************************************************************************** */
/* Update deformable boundary nodes */
/* ? Compute node from face values of field variables http://www.eureka.im/280.html */
/* https://www.cfd-online.com/Forums/fluent-udf/173970-node-data-usage-udf.html */
/* ******************************************************************************************** */
DEFINE_GRID_MOTION(update_node_positions, domain, dt, time, dtime)
{
	/* Mixture domain and thread */
	Thread *mix_face_thread = DT_THREAD((Dynamic_Thread *)dt);

	/* Local variables */
	face_t face, face_ds, face_us, face_min_gradient, face_max_peak;
	cell_t mix_c0;
	Thread *mix_ct0;
	Node   *node_p;
	real x_face, y_face_new, y_face_new_us, y_face_new_ds, dhdx; /* Face coordinates */
	real delta_x_node, x_node, y_node, x_node_us, x_node_ds, y_node_us_new, y_node_ds_new, dy, y_n_new; /* Node coordinates */
	real flow_time = CURRENT_TIME;
	real face_pos[ND_ND];
	int ReposeViolation = 0;
	int n;

	/* Local variables for searching faces */
	real y_face_new_pre, dhdx_pre;
	real y_face_new_max = 0.0;
	real dhdx_min = 0.0;



	/* ************************************************************ */
	/* Loop all faces */
	/* Determine new node y-positions from face position --> N_UDMI */
	/* Determine face gradients, check face angle vs. angle of repose and set violation flag */
	/* ************************************************************ */

	/* Modify first face and first node - Face ID 298 */
	face = 298; /* Set face index of first face */
	face_ds = 0; /* Set face index of downstream face */

	/* Current first face x-coordinate */
	F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates */
	x_face = face_pos[0]; /* Get x-component of first face */

	/* New first and downstream face y-coordinates */
	y_face_new = F_UDMI(face, mix_face_thread, NBH); /* Get new y-coordinate of first face */
	y_face_new_ds = F_UDMI(face_ds, mix_face_thread, NBH); /* Get new y-coordinate of downstream face */

	/* Loop the faces two nodes and compute new node y-positions */
	f_node_loop(face, mix_face_thread, n)
	{
		node_p = F_NODE(face, mix_face_thread, n);
		x_node = NODE_X(node_p); /* x-coordinate of node */
		y_node = NODE_Y(node_p); /* y-coordinate of node */

		/* Check if node is downstream or upstream of face x-position and average accordingly. */
		if (x_node < x_face) /* Node is upstream */
		{
			/* Compute new node position based on current face and downstream face y-position */
			y_node_us_new = y_face_new - 0.5*(y_face_new_ds - y_face_new);

			/* Save new node position to N_UDMI */
			N_UDMI(node_p, 0) = y_node_us_new;

			/* Save upstream node x-coordinate */
			x_node_us = x_node;
		}
		else /* Node is downstream */
		{
			/* Compute new node position based on current face and downstream face y-position */
			y_node_ds_new = y_face_new + 0.5*(y_face_new_ds - y_face_new);

			/* Save new node position to N_UDMI */
			N_UDMI(node_p, 0) = y_node_ds_new;

			/* Save downstream node y-coordinate */
			x_node_ds = x_node;
		}
	}

	/* Compute gradient of first face */
	dhdx = (y_node_ds_new - y_node_us_new) / (x_node_ds - x_node_us);

	/* Set flag if angle of repose is violated */
	if (fabs(atan(dhdx))>phi) { ReposeViolation = 1; }

	/* Loop all boundary faces and all corresponding nodes - This faceloop will only modify downsstream nodes as the upstream node is already modified */
	begin_f_loop(face, mix_face_thread)
	{
		if (face == 298)
		{
			/* First face, do nothing as this has been modified already */
		}
		else
		{
			/* Current face x-coordinate */
			F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates */
			x_face = face_pos[0]; /* Get x-component of current face */

			/* New current and downstream face x-coordinates */
			y_face_new = F_UDMI(face, mix_face_thread, NBH); /* Get new y-coordinate of current face */

			/* New y-component of downstream face */
			if (face == 297)
			{
				face_ds = face + 2; /* Second last face --> Downstream face 299 */
			}
			else if (face == 299) /* Last face */
			{
				face_ds = face; /* Last face --> Downstream Ghost face */
			}
			else /* All other faces */
			{
				face_ds = face + 1;
			}

			/* Get new y-coordinate of downstream face */
			y_face_new_ds = F_UDMI(face_ds, mix_face_thread, NBH);

			/* Loop the faces two nodes and compute new node y-positions */
			f_node_loop(face, mix_face_thread, n)
			{
				node_p = F_NODE(face, mix_face_thread, n);
				x_node = NODE_X(node_p); /* x-coordinate of node */
				y_node = NODE_Y(node_p); /* y-coordinate of node */

				/* Check if node is downstream or upstream of face x-position and average accordingly. */
				if (x_node < x_face) /* Node is upstream */
				{
					/* Node has been modified in previous loop, get node position from UDM */
					y_node_us_new = N_UDMI(node_p, 0);

					/* Save upstream node x-coordinate */
					x_node_us = x_node;
				}
				else /* Node is downstream */
				{
					/* Compute new node position based on current face and downstream face y-position */
					y_node_ds_new = y_face_new + 0.5*(y_face_new_ds - y_face_new);

					/* Save new node position to N_UDMI */
					N_UDMI(node_p, 0) = y_node_ds_new;

					/* Save downstream node y-coordinate */
					x_node_ds = x_node;
				}
			}

			/* Compute gradient */
			delta_x_node = x_node_ds - x_node_us;
			dhdx = (y_node_ds_new - y_node_us_new) / delta_x_node;

			/* Compute gradient of first face */
			dhdx = (y_node_ds_new - y_node_us_new) / (x_node_ds - x_node_us);

			/* Set flag if angle of repose is violated */
			if (fabs(atan(dhdx))>phi) { ReposeViolation = 1; }


			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

			/* TEST - Search for largest negative gradient */
			if (dhdx < dhdx_min)
			{
				dhdx_min = dhdx;
				face_min_gradient = face;
			}
			dhdx_pre = dhdx;

			/* TEST - Search for highest peak */
			if (y_face_new > y_face_new_max)
			{
				y_face_new_max = y_face_new;
				face_max_peak = face;
			}
			y_face_new_pre = y_face_new;

			/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
		}
	}
	end_f_loop(face, mix_face_thread)


		/* ************************************************************ */
		/* Loop while ReposeViolation == TRUE */
		/* Set ReposeViolation == FALSE */
		/* Loop all faces */
		/* Determine face gradients  --> F_UDMI */
		/* Check if angle of repose is violated at current cell */
		/* Check if angle of repose is larer or smaller than zero */
		/* Compute required delta_y to satisfy angle of repose */
		/* Modify nodes accordingly --> N_UDMI */
		/* Set ReposeViolation == TRUE */
		/* ************************************************************ */


		while (ReposeViolation == 1)
		{
			/* Reset repose violation flag */
			ReposeViolation = 0;

			/* Loop all" moving bed" faces */
			begin_f_loop(face, mix_face_thread)
			{
				/* Current face x-coordinate */
				F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates */
				x_face = face_pos[0]; /* Get x-component of current face */

				/* Loop the faces two nodes and get all new node coordinates */
				f_node_loop(face, mix_face_thread, n)
				{
					node_p = F_NODE(face, mix_face_thread, n);
					x_node = NODE_X(node_p); /* x-coordinate of node */
					y_node = NODE_Y(node_p); /* y-coordinate of node */

					/* Check if node is downstream or upstream of face x-position */
					if (x_node < x_face) /* Node is upstream */
					{
						/* Get node position from UDM */
						y_node_us_new = N_UDMI(node_p, 0);

						/* Save upstream node x-coordinate */
						x_node_us = x_node;
					}
					else /* Node is downstream */
					{
						/* Get node position from UDM */
						y_node_ds_new = N_UDMI(node_p, 0);

						/* Save upstream node x-coordinate */
						x_node_ds = x_node;
					}
				}

				/* Compute gradient */
				delta_x_node = x_node_ds - x_node_us;
				dhdx = (y_node_ds_new - y_node_us_new) / delta_x_node;

				/* Check for violation of angle of repose of current face */
				if (fabs(atan(dhdx)) > phi + 0.0001)
				{
					/* Set violation flag */
					ReposeViolation = 1;

					/* Compute required dy of nodes such that angle of repose is met */
					dy = delta_x_node / 2.0 * (fabs(dhdx) - tan(phi));

					/* Check for uphill or downhill slope */
					if (atan(dhdx) > 0.0) /* Slope is positive, i.e. uphill */
					{
						f_node_loop(face, mix_face_thread, n)
						{
							node_p = F_NODE(face, mix_face_thread, n);
							x_node = NODE_X(node_p); /* x-coordinate of node */

							/* Get new node y-coordinate from UDM */
							y_n_new = N_UDMI(node_p, 0);

							if (x_node < x_face) /* Node is upstream */
							{
								N_UDMI(node_p, 0) = y_n_new + dy;
							}
							else /* Node is downstream */
							{
								N_UDMI(node_p, 0) = y_n_new - dy;
							}
						}
					}
					else /* Slope is negative, i.e. downhill */
					{
						f_node_loop(face, mix_face_thread, n)
						{
							node_p = F_NODE(face, mix_face_thread, n);
							x_node = NODE_X(node_p); /* x-coordinate of node */

							y_n_new = N_UDMI(node_p, 0);

							if (x_node < x_face) /* Node is upstream */
							{
								N_UDMI(node_p, 0) = y_n_new - dy;
							}
							else /* Node is downstream */
							{
								N_UDMI(node_p, 0) = y_n_new + dy;
							}
						}
					}
				}
			}
			end_f_loop(face, mix_face_thread)
		}


	/* ************************************************************ */
	/* Modify mesh */
	/* ************************************************************ */

	/* Set deforming flag on adjacent cell zone --> cells adjacent to the deforming wall will also be deformed, in order to avoid skewness. */
	SET_DEFORMING_THREAD_FLAG(mix_face_thread->t0);

	/* Loop all boundary faces and all corresponding nodes, get new node y-coordinates from UDM and update mesh */
	begin_f_loop(face, mix_face_thread)
	{
		f_node_loop(face, mix_face_thread, n)
		{
			/* Get current node */
			node_p = F_NODE(face, mix_face_thread, n);

			/* Get new node y-position from N_UDMI */
			y_n_new = N_UDMI(node_p, 0);

			/* Assign new node y-coordinate to node */
			/* - if the new node y-position is not below the bottom channel wall <--> y_n_new  >= 0.0 */
			/* - if the new node y-position is not above the top channel wall <--> y_n_new  < ChannelHeight */
			/* - if the current node has not been updated previously */
			if (y_n_new >= 0.0 && y_n_new  < ChannelHeight && NODE_POS_NEED_UPDATE(node_p))
			{
				/* Set flag to indicate that the current node's position has been updated, so that it will not be updated during a future pass through the loop */
				NODE_POS_UPDATED(node_p);

				/* Overwrite the current nodes y-coordinate with the new y-coordinate */
				NODE_Y(node_p) = y_n_new;
			}
		}
		Update_Face_Metrics(face, mix_face_thread);

		/* Get updated face position */
		F_CENTROID(face_pos, face, mix_face_thread); /* Get face coordinates of current cell */

		/* Update all "new bed height" UDMs */
		UpdateUDM(NBH, face, mix_face_thread, face_pos[1]);

		/* Geometrical gradient of updated face */
		dhdx = CompFaceGradient(face, mix_face_thread);

		/* Save updated gradient to UDMs */
		UpdateUDM(DHDX, face, mix_face_thread, dhdx);
	}
	end_f_loop(face, mix_face_thread)
}
