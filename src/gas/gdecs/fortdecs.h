ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
c Front Tracking is a numerical method for the solution of partial differential equations 
c whose solutions have discontinuities.  
c 
c 
c Copyright (C) 1999 by The University at Stony Brook. 
c 
c  
c This library is free software; you can redistribute it and/or
c modify it under the terms of the GNU Lesser General Public
c License as published by the Free Software Foundation; either
c version 2.1 of the License, or (at your option) any later version.
c 
c This library is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c Lesser General Public License for more details.
c 
c You should have received a copy of the GNU Lesser General Public
c License along with this library; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c
c			fortdecs.h
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
c
#if !defined(_FORTDECS_H)
#define _FORTDECS_H

#if defined(cray)

#   define amax0	max
#   define amax1	max
#   define amin0	min
#   define amin1	min
#   define alog10	log10
#   define alog	log

#endif /* defined(cray) */

#   if defined(cray)

#   	define REAL    real*16
#   	define INTEGER integer*16
#   	define COMPLEX complex*32
#   	define NWPC    5

#   else /* defined(cray) */

#   	define REAL real*8
#   	define INTEGER integer*8
#   	if defined(_AIX)
#   	    define COMPLEX double complex
#   	else /* defined(_AIX) */
#   	    define COMPLEX complex*16
#   	endif /* defined(_AIX) */
#   	define NWPC    10

#   endif /* defined(cray) */

#   define Cmplx dcmplx

#   if !defined(cray)

#   	define aint		dint
#   	define anint		dnint
#   	define nini		idnint
#   	define abs		dabs
#   	define cabs		cdabs
#   	define amod		dmod
#   	define sign		dsign
#   	define dim		ddim
#   	define amax0		dmax0
#   	define amax1		dmax1
#   	define amin0		dmin0
#   	define amin1		dmin1
#   	define Real		dble
#   	define Aimag		dimag
#   	define conjg		dconjg
#   	define sqrt		dsqrt
#   	define csqrt		cdsqrt
#   	define exp		dexp
#   	define cexp		cdexp
#   	define alog		dlog
#   	define clog		cdlog
#   	define alog10		dlog10
#   	define sin		dsin
#   	define csin		cdsin
#   	define cos		dcos
#   	define ccos		cdcos
#   	define tan		dtan
#   	define ctan		cdtan
#   	define asin		dasin
#   	define acos		dacos
#   	define atan		datan
#   	define atan2		datan2
#   	define sinh		dsinh
#   	define cosh		dcosh
#   	define tanh		dtanh

#   endif /* !defined(cray) */

#endif /* !defined(_FORTDECS_H) */
