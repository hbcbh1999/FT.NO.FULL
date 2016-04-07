/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*				testseslib.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/
#include <cdecs.h>
#include <vmalloc.h>

LOCAL	boolean	file_is_readable(const char*);
LOCAL	void	usage(void);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    FORTRAN	void	FORTRAN_NAME(s2get)(int*,int*,double*,int*,int*,int*);
    FORTRAN	void	FORTRAN_NAME(s2geti)(int*,int*,double*,int*,int*,int*);
    FORTRAN	void	FORTRAN_NAME(s2eos)(int*,double*,double*,double*,
                                            double*,double*);
    FORTRAN	void	FORTRAN_NAME(s2eosi)(int*,double*,double*,double*,
                                             double*,double*);
    FORTRAN	void	FORTRAN_NAME(s2hugi)(int*,double*,double*,double*,int*);
    FORTRAN	void	FORTRAN_NAME(s2shki)(int*,double*,double*,double*,int*);
    FORTRAN	void	FORTRAN_NAME(s2abti)(int*,double*,double*,double*,int*);
#if defined(cray)
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,_fcd);
#else /* defined(cray) */
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,const char*,int);
#endif /* defined(cray) */

typedef struct	_S2DIR  {int lcmx, nrs, lcfw[10];} S2DIR;
    EXPORT	S2DIR	FORTRAN_NAME(s2dir);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

int C_MAIN_PROGRAM(int argc,char **argv)
{
	IMPORT	S2DIR	FORTRAN_NAME(s2dir);
	static const int  DEFAULT_TABLE_DIM = 20000;/*TOLERANCE*/
	static const char *fmt = "%lf %lf";
	boolean    inverse = NO;
	static const char *libname = "../../../databases/gas/sesdata/m101394.bin";
	double	   *tbls;
	double	   rho;
	double	   v0[3], v[6];
	double	   p[3], e[3], T[3];
	char	   s[Gets_BUF_SIZE];
	int	   mfl;
	int	   dim_tbls = DEFAULT_TABLE_DIM;
	int	   i, ids2;
	int	   ir = 1;
	int	   lcnt = 1;
	int	   lu = 10;
	int	   ifl;

	for (i = 0; i < argc; i++)
	{
	    if (strncmp(argv[i],"-i",2)==0)
	    	inverse = YES;
	    if (strncmp(argv[i],"-seslib",7)==0)
	    {
	    	libname = argv[++i];
	    }
	    if (strncmp(argv[i],"-tsize",5)==0)
	    {
	    	if (((sscanf(argv[++i],"%d",&dim_tbls) != 1)) ||
	    		(dim_tbls < 0))
	    	{
	    	    (void) fprintf(stderr,"Invalid size for table array\n");
	    	    usage();
	    	    return 0;
	    	}
	    }
	}

	if (sscanf(argv[argc-1],"%d",&ids2) != 1)
	{
	    (void) fprintf(stderr,"Can't read ids2\n");
	    usage();
	    return 0;
	}

	if (file_is_readable(libname) == NO)
	{
	    (void) fprintf(stderr,"Can't read %s\n",libname);
	    usage();
	    return 0;
	}

	uni_array(&tbls,dim_tbls,FLOAT);
	if (tbls == NULL)
	{
	    (void) fprintf(stderr,"Can't allocate table array\n");
	    usage();
	    return 0;
	}

	FORTRAN_NAME(s2dir).lcmx = dim_tbls;
	FORTRAN_NAME(s2dir).nrs = 1;
	for (i = 0; i < 10; i++)
	    FORTRAN_NAME(s2dir).lcfw[i] = 0;

#if defined(cray)
	SFORTRAN_NAME(oplib_uf)(&lu,_cptofcd(libname,strlen(libname)));
#else /* defined(cray) */
	SFORTRAN_NAME(oplib_uf)(&lu,libname,(int)strlen(libname));
#endif /* defined(cray) */

	if (inverse == YES)
	    FORTRAN_NAME(s2geti)(&ir,&ids2,tbls,&lcnt,&lu,&ifl);
	else
	    FORTRAN_NAME(s2get)(&ir,&ids2,tbls,&lcnt,&lu,&ifl);

	if (ifl == 0)
	{
	    (void) fprintf(stderr,"Material %d not found\n",ids2);
	    usage();
	    return 0;
	}
	if (ifl < 0)
	{
	    dim_tbls -= ifl;
	    free(tbls);
	    uni_array(&tbls,dim_tbls,FLOAT);
	    if (tbls == NULL)
	    {
	    	(void) fprintf(stderr,"Can't allocate table array\n");
	    	usage();
	    	return 0;
	    }
	    FORTRAN_NAME(s2dir).lcmx = dim_tbls;
	    FORTRAN_NAME(s2dir).nrs = 1;
	    for (i = 0; i < 10; i++)
	    	FORTRAN_NAME(s2dir).lcfw[i] = 0;
	    if (inverse == YES)
	    	FORTRAN_NAME(s2geti)(&ir,&ids2,tbls,&lcnt,&lu,&ifl);
	    else
	    	FORTRAN_NAME(s2get)(&ir,&ids2,tbls,&lcnt,&lu,&ifl);
	}

	if (inverse == YES)
	{
	    while (1)
	    {
		(void) fprintf(stderr,"Enter a density (gram/cc) and ");
		(void) fprintf(stderr,"internal energy (mj/kg): ");
		(void) Gets(s);
		if ((s[0] == '\0') || (sscanf(s,fmt,&rho,e) != 2))
		    break;
	    	FORTRAN_NAME(s2eosi)(&ir,tbls,&rho,e,p,T);
	    	(void) printf("p = %g GPa, ",             p[0]);
		(void) printf("dp/drho = %g GPa/(g/cc), ",p[1]);
	    	(void) printf("dp/de = %g GPa/(mj/kg)\n", p[2]);
	    	(void) printf("T = %g (K), ",                 T[0]);
		(void) printf("dT/drho = %g (K)/(grams/cc), ",T[1]);
	    	(void) printf("dT/de = %g K/(mj/kg)\n",       T[2]);
		v0[0] = rho;
		v0[1] = p[0];
		v0[2] = e[0];
		(void) fprintf(stderr,"Enter a density (gram/cc) on ");
		(void) fprintf(stderr,"the Hugoiniot: ");
	        fscan_float(stdin,v);
	    	FORTRAN_NAME(s2hugi)(&ir,tbls,v0,v,&mfl);
		switch (mfl)
		{
		case 1:
		    (void) fprintf(stderr,"State on Hugoniot\n");
		    (void) fprintf(stderr,"\trho = %g gram/cc\n",v[0]);
		    (void) fprintf(stderr,"\tp = %g GPa\n",v[1]);
		    (void) fprintf(stderr,"\te = %g mj/kg\n",v[2]);
		    (void) fprintf(stderr,"\tT = %g K\n",v[3]);
		    (void) fprintf(stderr,"\tShock velocity = %g km/sec\n",
		    	           v[4]);
		    (void) fprintf(stderr,
				   "\tParticle velocity = %g km/sec\n",v[5]);
		    break;
		case 0:
		    (void) fprintf(stderr,"Shock velocity is imaginary\n");
		    break;
		case -1:
		    (void) fprintf(stderr,"No Solution");
		    break;

		}
		(void) fprintf(stderr,"Enter an energy (mj/kg) on ");
		(void) fprintf(stderr,"the Hugoiniot: ");
	        fscan_float(stdin,v+2);
	    	FORTRAN_NAME(s2shki)(&ir,tbls,v0,v,&mfl);
		switch (mfl)
		{
		case 1:
		    (void) fprintf(stderr,"State on Hugoniot\n");
		    (void) fprintf(stderr,"\trho = %g gram/cc\n",v[0]);
		    (void) fprintf(stderr,"\tp = %g GPa\n",v[1]);
		    (void) fprintf(stderr,"\te = %g mj/kg\n",v[2]);
		    (void) fprintf(stderr,"\tT = %g K\n",v[3]);
		    (void) fprintf(stderr,"\tShock velocity = %g km/sec\n",
		    	           v[4]);
		    (void) fprintf(stderr,"\tParticle velocity = %g km/sec\n",
				   v[5]);
		    break;
		case 0:
		    (void) fprintf(stderr,"Shock velocity is imaginary\n");
		    break;
		case -1:
		    (void) fprintf(stderr,"No Solution");
		    break;

		}
		(void) fprintf(stderr,"Enter a density (gram/cc) on ");
		(void) fprintf(stderr,"the adiabat: ");
	        fscan_float(stdin,v);
	    	FORTRAN_NAME(s2abti)(&ir,tbls,v0,v,&mfl);
		switch (mfl)
		{
		case 1:
		    (void) fprintf(stderr,"State on Adiabat\n");
		    (void) fprintf(stderr,"\trho = %g gram/cc\n",v[0]);
		    (void) fprintf(stderr,"\tp = %g GPa\n",v[1]);
		    (void) fprintf(stderr,"\te = %g mj/kg\n",v[2]);
		    (void) fprintf(stderr,"\tT = %g K\n",v[3]);
		    (void) fprintf(stderr,"\tSound speed = %g km/sec\n",
		    	           v[4]);
		    (void) fprintf(stderr,"\tParticle velocity = %g km/sec\n",
				   v[5]);
		    break;
		case 0:
		    (void) fprintf(stderr,"Sound speed is imaginary\n");
		    break;
		case -1:
		    (void) fprintf(stderr,"Integration did not converge");
		    (void) fprintf(stderr,"State on Adiabat\n");
		    (void) fprintf(stderr,"\trho = %g gram/cc\n",v[0]);
		    (void) fprintf(stderr,"\tp = %g GPa\n",v[1]);
		    (void) fprintf(stderr,"\te = %g mj/kg\n",v[2]);
		    (void) fprintf(stderr,"\tT = %g K\n",v[3]);
		    (void) fprintf(stderr,"\tSound speed = %g km/sec\n",
		    	           v[4]);
		    (void) fprintf(stderr,"\tParticle velocity = %g km/sec\n",
				   v[5]);
		    break;

		}
	    }
	}
	else
	{
	    (void) printf("tbls[0] = %g, tbls[1] = %g\n",tbls[0],tbls[1]);
	    (void) printf("tbls[2] = %g, tbls[3] = %g\n",tbls[2],tbls[3]);
	    while (1)
	    {
	        (void) fprintf(stderr,"Enter a density (g/cc) and ");
		(void) fprintf(stderr,"temperature (K): ");
		(void) Gets(s);
		if ((s[0] == '\0') || (sscanf(s,fmt,&rho,T) != 2))
		    break;
	    	FORTRAN_NAME(s2eos)(&ir,tbls,&rho,T,p,e);
	    	(void) printf("p = %g GPa, ",             p[0]);
		(void) printf("dp/drho = %g GPa/(g/cc), ",p[1]);
	    	(void) printf("dp/dT = %g GPa/K\n",       p[2]);
	    	(void) printf("e = %g mj/kg, ",                   e[0]);
		(void) printf("de/drho = %g (mj/kg)/(grams/cc), ",e[1]);
	    	(void) printf("de/dT = %g (mj/kg)/K\n",           e[2]);
		(void) printf("sound speed squared = %g (km2/sec2)\n",
			p[1] + (p[0]/(rho*rho) - e[1])*p[2]/e[2]);
	    }
	}
	return 0;
}

LOCAL	void	usage(void)
{
	(void) fprintf(stderr,
		"Usage: testseslib [-i] [-tsize n] [-seslib filename] sesid\n");
}		/*end usage*/

LOCAL	boolean	file_is_readable(
	const char	*fname)
{
	FILE	*fp;
	if ((fp = fopen(fname,"r")) == NULL)
	    return NO;
	(void) fclose(fp);
	return YES;
}		/*end file_is_readable*/
