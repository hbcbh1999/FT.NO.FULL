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


#include <cdecs.h>
#include <sys/types.h>
#include <time.h>

/*
*			mkseslib
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Usage:	mkseslib [-update] [-preserve] filename
*/

typedef	enum {NEW_SESAME_LIBRARY=0,UPDATE_SESAME_LIBRARY=1} SesLibType;

LOCAL	boolean	file_is_readable(const char*);
LOCAL	void	usage(void);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
#   if defined(cray)
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,_fcd);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,_fcd);
    FORTRAN	void	FORTRAN_NAME(update)(double*,int*,_fcd,
				             int*,int*,int*,int*,int*,int*);
#else /* defined(cray) */
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,const char*,int);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,const char*,int);
    FORTRAN	void	FORTRAN_NAME(update)(double*,int*,const char*,
				             int*,int*,int*,int*,int*,int*,int);
#endif /* defined(cray) */
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

EXPORT	int	C_MAIN_PROGRAM(int argc, char **argv)
{
	char		cmd[256];
	char		fmt[Gets_BUF_SIZE];
	char            *sesdat;
	char		seslib[Gets_BUF_SIZE], nseslib[Gets_BUF_SIZE];
	static const char	*scratchlib = "seslib-scratch";
	double		date;
	boolean		svolib = NO;
	int		lp = 6, lup = 7, lib = 8, lm = 9, lnew = 10;
	int		inew = NEW_SESAME_LIBRARY;
	int		ncpw;
	int		i;
	time_t		clock;
	struct tm	*ltm;

	/*Set date*/
	clock = time(NULL);
	ltm = localtime(&clock);
	date = 10000.0*ltm->tm_mday + 100.0*(ltm->tm_mon+1) + ltm->tm_year;

	if (argc < 2)
	{
	    usage();
	    return 0;
	}
	for (i = 0; i < argc; i++)
	{
	    if (strcmp(argv[i],"-update")==0)
	    	inew = UPDATE_SESAME_LIBRARY;
	    if (strcmp(argv[i],"-preserve")==0)
	    {
	    	inew = UPDATE_SESAME_LIBRARY;
	    	svolib = YES;
	    }
	}
	sesdat = argv[argc-1];
	if (file_is_readable(sesdat) == NO)
	{
	    (void) fprintf(stderr,"ERROR,  can't open %s\n",sesdat);
	    usage();
	    return 0;
	}

	ncpw = sizeof(double)/sizeof(char);
	(void) sprintf(fmt,"(a%d)",ncpw);
	(void) sprintf(seslib,"%s.bin",basename(sesdat));

#if defined(cray)
	SFORTRAN_NAME(oplib_f)(&lup,_cptofcd(sesdat,strlen(sesdat)));
	SFORTRAN_NAME(oplib_uf)(&lm,_cptofcd(scratchlib,strlen(scratchlib)));
#else /* defined(cray) */
	SFORTRAN_NAME(oplib_f)(&lup,sesdat,(int)strlen(sesdat));
	SFORTRAN_NAME(oplib_uf)(&lm,scratchlib,(int)strlen(scratchlib));
#endif /* defined(cray) */

	if (inew == UPDATE_SESAME_LIBRARY)
	{
	    (void) sprintf(nseslib,"n%s",seslib);
#if defined(cray)
	    SFORTRAN_NAME(oplib_uf)(&lib,_cptofcd(seslib,strlen(seslib)));
	    SFORTRAN_NAME(oplib_uf)(&lnew,_cptofcd(nseslib,strlen(nseslib)));
#else /* defined(cray) */
	    SFORTRAN_NAME(oplib_uf)(&lib,seslib,(int)strlen(seslib));
	    SFORTRAN_NAME(oplib_uf)(&lnew,nseslib,(int)strlen(nseslib));
#endif /* defined(cray) */
	}
	else
	{
#if defined(cray)
	    SFORTRAN_NAME(oplib_uf)(&lnew,_cptofcd(seslib,strlen(seslib)));
#else /* defined(cray) */
	    SFORTRAN_NAME(oplib_uf)(&lnew,seslib,(int)strlen(seslib));
#endif /* defined(cray) */
	}

#if defined(cray)
	FORTRAN_NAME(update)(&date,&ncpw,_cptofcd(fmt,strlen(fmt)),
		             &inew,&lup,&lib,&lm,&lnew,&lp);
#else /* defined(cray) */
	FORTRAN_NAME(update)(&date,&ncpw,fmt,&inew,&lup,&lib,&lm,
			     &lnew,&lp,(int)strlen(fmt));
#endif /* defined(cray) */

	if ((inew == UPDATE_SESAME_LIBRARY) && (svolib == NO))
	{
	    (void) sprintf(cmd,"mv -f %s %s",nseslib,seslib);
	    system(cmd);
	}
	(void) sprintf(cmd,"/bin/rm -f %s",scratchlib);
	system(cmd);
}

LOCAL	void	usage(void)
{
	(void) fprintf(stderr,"Usage: ");
	(void) fprintf(stderr,"mkseslib [-update] [-preserve] filename\n");
}

LOCAL	boolean	file_is_readable(
	const char	*fname)
{
	FILE	*fp;
	if ((fp = fopen(fname,"r")) == NULL)
	    return NO;
	(void) fclose(fp);
	return YES;
}		/*end file_is_readable*/
