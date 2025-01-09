/* Parabolic velocity profile at inlet */ 
DEFINE_PROFILE(parabolic_velocity_profile, phase_face_thread, position)
{
#if !RP_HOST
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
			if (x_node < x_face) {h_bed_inlet = y_node;}
		}

	H = ChannelHeight;
	h = h_bed_inlet;
	
	/* Current cross-sectional area of bed based on pipe cross-sectional area */
	A_Bed = (pow((H / 2.0), 2.0)*acos(1.0 - 2.0*h / H) - (H / 2.0 - h_bed_inlet)*pow((H*h - h*h), 0.5));
	
	/* Current fluid bulk velocity based on existing cross-sectional flow area */
	U = Q_f / (A_Pipe - A_Bed);
	
	/* Loop all boundary faces and compute corresponding x-velocity */
	begin_f_loop(face,phase_face_thread)
	{
		/* Get face coordinates */
		F_CENTROID(face_pos,face,phase_face_thread);
		y_face = face_pos[1]; /* Get y-component of current face */
		
		/* x-velocity in current cell */
		u =2.0*U-2.0*U/pow(((h-H)/2.0),2.0)*pow((y_face-(h+(H-h)/2.0)),2.0);
		
		/* Assign to current face */
		F_PROFILE(face,phase_face_thread,position) = u;
	}
	end_f_loop(face,phase_face_thread)
#endif
}