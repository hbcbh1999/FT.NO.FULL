c
c			gflw.F
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
c
#include <gdecs/fortdecs.h>
c-------------------------------------------------------------------------------
c			   artvsc()
c
c Computes the coefficients of artificial viscosity used in the
c Split Lax-Wendroff finite difference step
c-------------------------------------------------------------------------------
      subroutine artvsc(v0,c0,v1,c1,vsc0,vsc1,vsc2,chi1,chi3,cms2,lam0,
     1  vsccf,spcoef,ivs)
      implicit none
      REAL v0(1), c0(1), v1(1), c1(1), vsc0(1), vsc1(1), vsc2(1)
      REAL chi1(1), chi3(1), cms2(1), lam0(1), vsccf(1), spcoef(1)
      REAL tol
      integer ivs, i
      integer ivsm1
      ivsm1 = ivs - 1
c Compute eigenvalues and maximum wave speed
      do 100 i = 1,ivs
          chi1(i) = v0(i) - c0(i)
          chi3(i) = v0(i) + c0(i)
 100  continue
      do 200 i = 1,ivsm1
	  tol =  0.00025*(abs(v0(i))+c0(i)+abs(v0(i+1))+c0(i+1))
	  if (c1(i) .le. tol) then
		vsc0(i) = vsccf(i)*abs(v0(i+1)-v0(i))
		vsc1(i) = 0.0
		vsc2(i) = 0.0
	  else
          	cms2(i) = 2.0*c1(i)*c1(i)
          	lam0(i) = v1(i) - c1(i)
          	vsc0(i) = vsccf(i)*abs(chi1(i+1) - chi1(i))
          	vsc1(i) = vsccf(i)*abs(v0(i+1) - v0(i))
          	vsc2(i) = vsccf(i)*abs(chi3(i+1) - chi3(i))
          	vsc2(i) = (vsc0(i) - 2.0*vsc1(i) + vsc2(i))/cms2(i)
          	vsc1(i) = (vsc1(i) - vsc0(i))/c1(i) -
     1                   (v1(i) + lam0(i))*vsc2(i)
          	vsc0(i) = vsc0(i) - lam0(i)*(vsc1(i) + lam0(i)*vsc2(i))
	endif
 200  continue
      return
      end
c-------------------------------------------------------------------------------
c			lwvec1()
c
c	Performs the first half of the split Lax-Wendroff method 
c-------------------------------------------------------------------------------
      subroutine lwvec1(r0,e0,m10,m20,p0,v0,r1,e1,m11,m21,v1,
     1	gecor,gm1cor,k1,k2,vsc1,vsc2,lmf,lm,
     2  gdt,isgrav,ivs)
      implicit none
      REAL r0(1), e0(1), m10(1), m20(1), p0(1), v0(1)
      REAL r1(1), e1(1), m11(1), m21(1), v1(1)
      REAL gecor(1), gm1cor(1), k1(1), k2(1)
      REAL vsc1(1), vsc2(1), lmf(1)
      REAL lm, gdt(1)
      integer isgrav, ivs, i
      if (isgrav .ne. 0) then
	  do 100 i=1,ivs
	      gecor(i) = gdt(i)*m11(i)
	      gm1cor(i) = gdt(i)*r1(i)
 100      continue
      endif
      do 200 i=1,ivs
	  lmf(i) = 0.5*(vsc2(i) + lm)
          k1(i) = 0.5*(1.0 + vsc1(i)) + lmf(i)*v0(i) 
          k2(i) = 0.5*(1.0 - vsc1(i)) - lmf(i)*v0(i+1) 
          r1(i) = k2(i)*r0(i+1) + k1(i)*r0(i)
          e1(i) = k2(i)*e0(i+1) + k1(i)*e0(i) -
     1                   lmf(i)*(v0(i+1)*p0(i+1) - v0(i)*p0(i))
          m11(i) = k2(i)*m10(i+1) + k1(i)*m10(i) -
     1                   lmf(i)*(p0(i+1) - p0(i))
          m21(i) = k2(i)*m20(i+1) + k1(i)*m20(i)
 200  continue
      if (isgrav .ne. 0) then
	  do 300 i=1,ivs
	      e1(i) = e1(i) + gecor(i)
	      m11(i) = m11(i) + gm1cor(i)
 300      continue
      do 400 i=1,ivs
          v1(i) = m11(i)/r1(i)
 400  continue
      endif
      return
      end
c-------------------------------------------------------------------------------
c			lwvec2()
c
c	Performs the second half of the split Lax-Wendroff method 
c-------------------------------------------------------------------------------
      subroutine lwvec2(r0,e0,m10,m20,r1,e1,m11,m21,
     1		p1,v1,rans,eans,m1ans,m2ans,k1,k2,k3,vsc0,
     2          ivs,lm,gdt,isgrav)
      implicit none
      REAL r0(1), e0(1), m10(1), m20(1)
      REAL r1(1), e1(1), m11(1), m21(1), p1(1), v1(1)
      REAL rans(1), eans(1), m1ans(1), m2ans(1)
      REAL k1(1), k2(1), k3(1), vsc0(1)
      REAL lm, gdt(1)
      integer isgrav, ivs, i
      REAL lmd2
      lmd2 = 0.5*lm
      do 100 i=1,ivs
          k1(i) = lmd2*vsc0(i) 
          k3(i) = lmd2*vsc0(i+1) 
          k2(i) = 1.0 - k1(i) - k3(i)
          rans(i) = k1(i)*r0(i) + k2(i)*r0(i+1) + k3(i)*r0(i+2)
     1               - lm*(m11(i+1) - m11(i))
          eans(i) = k1(i)*e0(i) + k2(i)*e0(i+1) + k3(i)*e0(i+2)
     1               - lm*(v1(i+1)*(e1(i+1) + p1(i+1)) - 
     2                        v1(i)*(e1(i) + p1(i)))
          m1ans(i) = k1(i)*m10(i) + k2(i)*m10(i+1) + k3(i)*m10(i+2)
     1               - lm*(v1(i+1)*m11(i+1) + p1(i+1) -
     2                        v1(i)*m11(i) - p1(i))
          m2ans(i) = k1(i)*m20(i) + k2(i)*m20(i+1) + k3(i)*m20(i+2)
     1               - lm*(v1(i+1)*m21(i+1) - v1(i)*m21(i))
 100  continue
      if (isgrav .ne. 0) then
	   do 200 i=1,ivs
              eans(i) = eans(i) + gdt(i)*(m11(i+1) + m11(i))
              m1ans(i) = m1ans(i) + gdt(i)*(r1(i+1) + r1(i))
 200      continue
	endif
      return
      end
