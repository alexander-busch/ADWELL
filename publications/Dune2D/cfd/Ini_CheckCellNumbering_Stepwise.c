/* Ini_CheckFluentCellNumbering_Stepwise.c


/file/read-case mesh.cas OK
/define/user-defined/user-defined-memory 1
/define/user-defined/compiled-functions unload "libudf"
!rmdir /s /q libudf
!if exist "MovingBedCellValues.txt" del "MovingBedCellValues.txt"
/define/user-defined/compiled-functions compile "libudf" y "Ini_CheckCellNumbering_Stepwise.c" , ,
/define/user-defined/compiled-functions load "libudf"
/define/user-defined/function-hooks/initialization "initialisation::libudf" ,
/solve/initialize/hyb-initialization OK
/display/contour mixture udm-0 , ,


/file/read-case-data case_initialized.cas OK
/define/user-defined/compiled-functions unload "libudf"
!rmdir /s /q libudf
!if exist "MovingBedCellValues.txt" del "MovingBedCellValues.txt"
/define/user-defined/compiled-functions compile "libudf" y "Ini_CheckCellNumbering_Stepwise.c" , ,
/define/user-defined/compiled-functions load "libudf"
/define/user-defined/function-hooks/initialization "initialisation::libudf" ,
/solve/initialize/hyb-initialization OK
/display/contour mixture udm-0 , ,


*/



/* Fluent */
#include "udf.h"
#include "mem.h"



/* ******************************************************************************************** */ 
/* Initialisation */ 
/* ******************************************************************************************** */

DEFINE_INIT(initialisation, domain)
{
	cell_t c,c0;
	Thread *t,*tf,*ct0;
	face_t f;
	/*real NV_VEC(xf);*/
	real xf[ND_ND];
	int i=0;
	FILE *fp;
	tf = Lookup_Thread(domain, 8);
   
   	/* Create a text file */
	fp = fopen("MovingBedCellValues.txt","a");
	fprintf(fp,"x \t cf# \t i \t udm-f \t udm-c \n");
   
   	
	/* First cell */
	c0=F_C0(198,tf);
	ct0=THREAD_T0(tf);
	F_CENTROID(xf, 198, tf); /* get face coordinates = Initial node positions */
		
	F_UDMI(198,tf,0) = i;
	C_UDMI(c0,ct0,0) = i;
		
	fprintf(fp,"%f \t %u \t %u \t %g \t %g \n",xf[0],198,i,F_UDMI(198,tf,0),C_UDMI(c0,ct0,0));  
		
	i=i+1;   
   
	begin_f_loop (f, tf)
	{
		if (f==198) /* f==0 || f==198 */
		{		

		}
		
		/* Do nothing */
		else
		{
		
		c0=F_C0(f,tf);
		ct0=THREAD_T0(tf);
		
		F_CENTROID(xf, f, tf); /* get face coordinates = Initial node positions */
		
		F_UDMI(f,tf,0) = i;
		C_UDMI(c0,ct0,0) = i;
		
		fprintf(fp,"%f \t %u \t %u \t %g \t %g \n",xf[0],f,i,F_UDMI(f,tf,0),C_UDMI(c0,ct0,0));  
		
		i=i+1;
		}

    }
	end_f_loop (f, tf)
	
	
	/* Last cell 
	c0=F_C0(0,tf);
	ct0=THREAD_T0(tf);
	F_CENTROID(xf, 0, tf); /* get face coordinates = Initial node positions 
		
	F_UDMI(0,tf,0) = i;
	C_UDMI(c0,ct0,0) = i;
		
	fprintf(fp,"%f \t %u \t %u \t %g \t %g \n",xf[0],0,i,F_UDMI(0,tf,0),C_UDMI(c0,ct0,0));  
	*/	

	fclose(fp);
}

