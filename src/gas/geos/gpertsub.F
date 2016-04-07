c
c				gpertsub.F
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
c
c	Contains subroutines for use in sinusoidal perturbation problems.
c
c	See I. G. Currie, Fundamental Mechanics of Fluids, Chapter 6.,
c	or see Lamb's Hydrodynamics for the incompressible analysis;
c	for the compressible analysis, use Crocco's equation in place
c	of Bernoulli's equation and the wave equation in place of
c	Laplace's equation.
c
#include <gdecs/fortdecs.h>
      subroutine khstate(x, y, t, rho, yintfc, ybdry, amp, kx, gspeed, 
     &    pspeed, phi, Vx, Vy, Drho, Dp, press, vel, gamma)
      implicit none
      REAL x, y, Vx, Vy, Drho, Dp
      REAL t, yintfc, ybdry, amp, kx, gspeed
      REAL pspeed, phi, rho, press, vel, gamma
      common /misc/ i
      COMPLEX i
      REAL Aimag
      COMPLEX cexp, ccosh, csinh, Cmplx, csqrt
      REAL csq, denfac
      COMPLEX expon, omega, factor, sqmach, difvel, iky
c
c				khstate():
c
c	Calculates the state variables rho, e, mx, and my for the
c	position x,y at time t.
c
      omega = kx*Cmplx(pspeed, gspeed)
      expon = (-i)*kx*amp*cexp(i*(kx*x-omega*t-phi))
      csq = gamma*press/rho
      difvel = Cmplx(pspeed, gspeed)-vel
      sqmach = difvel**2/csq
      iky = kx*csqrt(1.-sqmach)
      factor = (kx/iky)*ccosh(iky*(y-ybdry))/csinh(iky*(yintfc-ybdry))
      denfac = Real(factor*expon*sqmach)
      Drho = rho*denfac
      Vx = Real(factor*expon*difvel)
      factor = csinh(iky*(y-ybdry))/csinh(iky*(yintfc-ybdry))
      Vy = Aimag(factor*expon*difvel)
      Dp = press*denfac
      end
#if DONT_COMPILE
C  Old version of khstate,  saved for comparison of gravity
      subroutine khstate(rho, e, mx, my, x, y, t, yintfc, ybdry, 
     &   amp, kx, gspeed, pspeed, phi, dens, press, vel, gamma, grav)
      implicit none
      REAL rho, e, mx, my, x, y
      REAL t, yintfc, ybdry, amp, kx, gspeed
      REAL pspeed, phi, dens, press, vel, gamma
      REAL grav
      common /misc/ i
      COMPLEX i
      REAL Aimag
      COMPLEX cexp, ccosh, csinh, Cmplx, csqrt
      REAL csq, denfac
      COMPLEX expon, omega, factor, sqmach, difvel, iky
      omega = kx*Cmplx(pspeed, gspeed)
      expon = (-i)*kx*amp*cexp(i*(kx*x-omega*t-phi))
      csq = gamma*press/dens
      dens = dens*(1.-(gamma-1.)*grav*(y-yintfc)/csq)**(1./(gamma-1.))
      csq = gamma*press/dens
      difvel = Cmplx(pspeed, gspeed)-vel
      sqmach = difvel**2/csq
      iky = kx*csqrt(1.-sqmach)
      factor = (kx/iky)*ccosh(iky*(y-ybdry))/csinh(iky*(yintfc-ybdry))
      denfac = Real(factor*expon*sqmach)+1.
      rho = dens*denfac
      mx = rho*(vel+Real(factor*expon*difvel))
      factor = csinh(iky*(y-ybdry))/csinh(iky*(yintfc-ybdry))
      my = rho*Aimag(factor*expon*difvel)
      e = .5*(mx*mx+my*my)/rho+press*denfac**gamma/(gamma-1.)
      end
#endif /*DONT_COMPILE*/
      subroutine velint(vpot, x, t, yintfc, ybdry, amp, kx, 
     &   gspeed, pspeed, phi, dens, press, vel, gamma, grav)
      implicit none
      REAL vpot, x, t, yintfc, ybdry, amp
      REAL kx, gspeed, pspeed, phi, dens, press
      REAL vel, gamma, grav
      common /misc/ i
      COMPLEX i
      REAL Aimag
      COMPLEX cexp, ccosh, csinh, Cmplx, csqrt
      REAL y, csq
      COMPLEX expon, omega, factor, sqmach, difvel, iky
c
c				velint():
c
c	Calculates the line integral velint of the velocity from the point
c	on the interface at 0 to the point on the interface at x.
c	This integeral is given by the difference between values of the
c	velocity potential at these two points.
c
      omega = kx*Cmplx(pspeed, gspeed)
      expon = (-i)*kx*amp*cexp((-i)*(omega*t+phi))*(cexp(i*kx*x)-1.)
      y = yintfc+(1./kx)*Real(expon)
      csq = gamma*press/dens
      dens = dens*(1.-(gamma-1.)*grav*(y-yintfc)/csq)**(1./(gamma-1.))
      csq = gamma*press/dens
      difvel = Cmplx(pspeed, gspeed)-vel
      sqmach = difvel**2/csq
      iky = kx*csqrt(1.-sqmach)
      factor = (1./iky)*ccosh(iky*(y-ybdry))/csinh(iky*(yintfc-ybdry))
      vpot = vel*x+Aimag(difvel*factor*expon)
      end
      subroutine speed(gspeed, incspd, pspeed, kx, rhoa, rhob, 
     &   csqra, csqrb, dista, distb, vela, velb, grvty, sigma)
      implicit none
      REAL gspeed, incspd, pspeed, kx, rhoa, rhob
      REAL csqra, csqrb, dista, distb, vela, velb
      REAL grvty, sigma
      common /misc/ i
      COMPLEX i
      REAL cabs, tanh, sqrt, Aimag
      COMPLEX delta, Cmplx
      integer n
      REAL hterm, b, a, psp0, gsp0, stterm
      REAL rtterm, rhobot, rhoaot, eps
      COMPLEX spd, del
      REAL temp, temp1
c
c				speed():
c
c	Calculates the propagation speed of a sinusoidal perturbation of the
c	interface between two fluids. The "growth speed" is the imaginary part
c,
c	and the "propagation speed" is the real part, of omega / k.  Here the
c	two fluids are labeled "a" (for above) and "b" (for below); the
c	variables dista and distb are distances from the unperturbed
c	interface to the Neumann boundary for the corresponding fluid.
c
      a = 0.
      b = 1.
      i = Cmplx(a, b)
      rhoaot = rhoa/tanh(kx*dista)
      rhobot = rhob/tanh(kx*distb)
      hterm = rhoaot*rhobot*(vela-velb)**2/(rhoaot+rhobot)
      rtterm = (rhoa-rhob)*grvty/kx
      stterm = (-kx)*sigma
      psp0 = 0.
      gsp0 = sqrt((hterm+rtterm+stterm)/(rhoaot+rhobot))
      incspd = gsp0
c
c		Newton's method:
c
      spd = Cmplx(psp0, gsp0)
      eps = .00001*cabs(spd)
      n = 0
      del = eps+1.
   1  if (cabs(del) .le. eps) goto  6
         n = n+1
         if (n .le. 10) goto 5
            write (6,  2) 
   2        format ('ERROR:  Newton''s method failed')
            write (6,  3) psp0, gsp0
   3        format ('spd0 = ', 1pe20.6, 1pe20.6)
            temp1 = Real(spd)
            temp = Aimag(spd)
            write (6,  4) temp1, temp
   4        format ('spd = ', 1pe20.6, 1pe20.6)
            stop 
   5     del = delta(spd, kx, rhoa, rhob, csqra, csqrb, dista, distb, 
     &      vela, velb, rtterm+stterm)
         spd = spd+del
c####
c		write(,"del = ":c,Real(del):e(20,6),Aimag(del):e(20,6))
c		write(,"spd = ":c,Real(spd):e(20,6),Aimag(spd):e(20,6))
c		write(,:s(2))
c####
         goto  1
   6  pspeed = Real(spd)
      gspeed = Aimag(spd)
c####
c	write(,"n = ":c,n:"I3")
c	write(,"gspeed = ":c,gspeed:e(20,6)," pspeed = ":c,pspeed:e(20,6))
c	write(,"gsp0   = ":c,gsp0:e(20,6)," psp0   = ":c,psp0:e(20,6))
c	write(,:s(3))
c####
      end
      COMPLEX function delta(speed, kx, rhoa, rhob, csqra, csqrb
     &   , dista, distb, vela, velb, ans)
      implicit none
      REAL kx, rhoa, rhob, csqra, csqrb, dista
      REAL distb, vela, velb, ans
      COMPLEX speed
      common /misc/ i
      COMPLEX i
      COMPLEX ctanh, csqrt
      COMPLEX root, sqdifv, valueb, valuea, sqmach, difvel
      COMPLEX lnderb, tnh, lndera
      difvel = speed-vela
      sqdifv = difvel**2
      sqmach = sqdifv/csqra
      root = csqrt(1.-sqmach)
      tnh = ctanh(kx*root*dista)
      valuea = rhoa*sqdifv/(root*tnh)
      lndera = (sqmach/(1.-sqmach)+2.+kx*dista*(1./tnh-tnh)*sqmach/root)
     &   /difvel
      difvel = speed-velb
      sqdifv = difvel**2
      sqmach = sqdifv/csqrb
      root = csqrt(1.-sqmach)
      tnh = ctanh(kx*root*distb)
      valueb = rhob*sqdifv/(root*tnh)
      lnderb = (sqmach/(1.-sqmach)+2.+kx*distb*(1./tnh-tnh)*sqmach/root)
     &   /difvel
c####
c	write(,"valuea = ":c,Real(valuea):e(20,6),Aimag(valuea):e(20,6))
c	write(,"lndera = ":c,Real(lndera):e(20,6),Aimag(lndera):e(20,6))
c	write(,"valueb = ":c,Real(valueb):e(20,6),Aimag(valueb):e(20,6))
c	write(,"lnderb = ":c,Real(lnderb):e(20,6),Aimag(lnderb):e(20,6))
c####
      delta = (-(valuea+valueb-ans))/(lndera*valuea+lnderb*valueb)
      return
      end
      COMPLEX function ccosh(x)
      implicit none
      COMPLEX x
      common /misc/ i
      COMPLEX i
      COMPLEX cexp
      COMPLEX expon
      expon = cexp(x)
      ccosh = .5*(expon+1./expon)
      return
      end
      COMPLEX function csinh(x)
      implicit none
      COMPLEX x
      common /misc/ i
      COMPLEX i
      COMPLEX cexp
      COMPLEX expon
      expon = cexp(x)
      csinh = .5*(expon-1./expon)
      return
      end
      COMPLEX function ctanh(x)
      implicit none
      COMPLEX x
      common /misc/ i
      COMPLEX i
      COMPLEX cexp
      COMPLEX expon, exponi
      expon = cexp(x)
      exponi = 1./expon
      ctanh = (expon-exponi)/(expon+exponi)
      return
      end
