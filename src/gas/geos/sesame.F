c
c		sesame.F
c
#if defined(SESAME_CODE)
#include <gdecs/fortdecs.h>
#define DONT_COMPILE 0
#if defined(float)
#	define	I2FLOAT	dfloat
#	define ERRTEMP -1.d99
#	define ZERO 0.d0
#	define TENTH 1.0d-01
#	define FIFTH 2.0d-01
#	define HALF 0.5d0
#	define ONE 1.0d0
#	define ONEPOINTTWO 1.2d0
#	define TWO 2.0d0
#	define TWO0ONE 201.d0
#	define ONENINENINE 199.d0
#	define TWOHUNDRED 200.d0
#	define THREE0ONE 301.d0
#	define THREE0THREE 303.d0
#	define FIVE0TWO 502.d0
#	define FIVE0THREE 503.d0
#	define FIVE0FOUR 504.d0
#	define FIVE0FIVE 505.d0
#	define ONETHOUSAND 1000.d0
#	define EPS3 1.0d-03
#	define EPS6 1.0d-06
#	define FIVEEPS6 5.d-06
#	define EPS12 1.0d-12
#	define BIG10 1.0d+10
#else /*defined(float)*/
#	define	I2FLOAT	float
#	define ERRTEMP -1.e30
#	define ZERO 0.0
#	define TENTH 0.1
#	define FIFTH 0.2
#	define HALF 0.5
#	define ONE 1.0
#	define ONEPOINTTWO 1.2
#	define TWO 2.0
#	define TWO0ONE 201.0
#	define ONENINENINE 199.0
#	define TWOHUNDRED 200.0
#	define THREE0ONE 301.0
#	define THREE0THREE 303.0
#	define FIVE0TWO 502.0
#	define FIVE0THREE 503.0
#	define FIVE0FOUR 504.0
#	define FIVE0FIVE 505.0
#	define ONETHOUSAND 1000.0
#	define EPS3 1.0e-03
#	define EPS6 1.0e-06
#	define FIVEEPS6 5.e-06
#	define EPS12 1.0e-12
#	define BIG10 1.0e+10
#endif /*defined(float)*/
cses2i
cses2i----------------------------------------------------------
cses2i
cses2i                 file ses2i - 5 oct 1977
cses2i
cses2i  this file is a fortran listing and documentation of
cses2i  subprograms needed to implement sesame tables in a hydro
cses2i  code or other computer program which uses equation of
cses2i  state data.  the routines are classified into several
cses2i  types, as shown below.
cses2i
cses2i  these routines are written in ans fortran so that users
cses2i  should not have to adapt them to their local compilers.
cses2i  however, there may be certain exceptions to this rule.
cses2i     1. input and output may be system dependent, so it
cses2i  is confined to the two routines listed under a.
cses2i
cses2i  a. generalized input/output - these routines call system
cses2i     i/o.  the user is expected to check and modify them,
cses2i     if necessary.
cses2i        outbuf - binary write or end file mark
cses2i        inbuf  - binary read and end file test
cses2i     outbuf is used to create or update the eos library.
cses2i     inbuf is used to read the library and should be
cses2i     compatible with outbuf.
cses2i
cses2i  b. library preprocessor - subroutine update either
cses2i     creates or adds data tables to a binary library file
cses2i     from a data tape written in bcd/card image format.
cses2i     although this routine is system independent, the
cses2i     user must provide the library file, a scratch file
cses2i     and a file on which the new library (created or updated)
cses2i     will be returned.  if the user specifies the creation
cses2i     of a new library, the old library file is not used.
cses2i     if the user specifies an update to the existing library,
cses2i     that library remains intact after the update.  it is the
cses2i     responsibility of the user to save the new library
cses2i     after execution.
cses2i     externals: inbuf, outbuf, rdcard, s2bldm, retmat,
cses2i                bdindx, s2srtm, s2bldl, s2updt, skipmf,
cses2i                copmat, copmts, s2mod, asort2, upsum,
cses2i                compak, uptpgm, combin, tbfch.
cses2i
cses2i  c. packaged eos routines - these routines are intended
cses2i     for hydro codes which use density and temperature as
cses2i     the independent variables and pressure and internal
cses2i     energy as the dependent variables.
cses2i        s2get - get data from library
cses2i        s2eos - compute an eos point
cses2i     externals: inbuf, tbfch, la301a, isrch, ratfn2.
cses2i
cses2i  d. packaged sin code routines - these routines are to be
cses2i     incorporated into a hydro code, such as sin, which
cses2i     uses the inverted form of an equation of state.
cses2i     density and internal energy are the independent
cses2i     variables, and pressure and temperature are the
cses2i     dependent variables.
cses2i        s2geti - get data from library
cses2i        s2eosi - compute an eos point
cses2i     externals: inbuf, tbfch, la302a, inv301, isrch,
cses2i                isrchk, ratfn1, ratfn2.
cses2i
cses2i  e. test problem - subroutine tprb1 accesses data
cses2i     for several materials specified by the user and
cses2i     computes and prints hugoniot curves.
cses2i     externals: s2geti, s2shki.
cses2i
cses2i  f. auxiliary packaged routines (optional)
cses2i        s2hugi - compute point on a hugoniot
cses2i        s2shki - compute point on a hugoniot
cses2i        s2abti - compute point on an adiabat
cses2i
cses2i  g. search and interpolation - basic routines to operate
cses2i     on data strings read in from library.  not usually
cses2i     called directly - la301a, la302a, wdw301, inv301,
cses2i     isrch, isrchk, ratfn1, ratfn2, la401a, la302b,
cses2i     la302s.
cses2i
cses2i  h. library manipulation - basic routines to get data
cses2i     from a library.  not usually called directly -
cses2i     tbfch, copmat, rdcard, skipmf, copmts, libnit.
cses2i
cses2i  i. library checkout - cklib is a subroutine used to check
cses2i     a library for self consistency and data range errors.
cses2i     externals: inbuf, skipmf.
cses2i
cses2i----------------------------------------------------------
cses2i
#if defined(MAKE_SES_LIB)
      subroutine outbuf(lu,z,n)
      implicit none
      integer lu
      REAL    z(1)
      integer n
c---------------------------------------------------------------
c
c  subroutine:  outbuf(lu,z,n)
c
c  purpose:     write data or end file mark onto disc or tape
c
c  arguments:   lu (input) - number of tape, disc, etc.
c               z  (input) - name of local array
c               n  (input) - number of words to be written
c                 for n.le.0, writes an end file mark
c
c  remarks:     this routine is just an binary write or
c               end file mark.  the user may wish to modify it
c               and make corresponding changes to inbuf.
c
c  remarks:     system binary write and end file.
c
c  programmer:  g. t. rood and g. i. kerley, t-4.
c
c  date:        7 december 1976
c
c---------------------------------------------------------------
c-----local variables
      integer i
      character*5 fm
      data fm/'$fm'/
c-----end local variables
      if (n .gt. 0) then
        write(lu)(z(i),i=1,n)
      else
	write(lu) fm
      end if
      return
      end
#endif /*defined(MAKE_SES_LIB)*/
      subroutine inbuf(lu,z,n,im)
      implicit none
      integer lu
      REAL    z(1)
      integer n
      integer im
c---------------------------------------------------------------
c
c  subroutine:  inbuf(lu,z,n,im)
c
c  purpose:     read from disc or tape into local array
c
c  arguments:   lu  (input) - number of tape, disc, etc.
c               z  (in/out) - name of local array
c               n   (input) - number of words to be read
c               im (output) - end file flag
c                 im = 0, if end file mark encountered
c                 im = 1, if no end file mark
c
c  remarks:     this routine is just an binary read and test
c               for end of file.  it assumes the data records
c               were written with outbuf and should be made
c               compatible with that routine.  note that the
c               two branch end file test is not standard.
c
c  externals:   system binary read, end file test.
c
c  programmer:  g. t. rood and g. i. kerley, t-4.
c
c  date:        6 december 1976
c
c---------------------------------------------------------------
c-----local variables
      integer i
      character*5 fm, test
      data fm/'$fm'/
c-----end local variables
      read(lu,end=1) test
      if (test .eq. fm) then
	im = 0
	return
      end if
      backspace lu
      read(lu,end=1)(z(i),i=1,n)
      im = 1
      return
 1    im = 0
      return
      end
      subroutine s2get(ir,ids2,tbls,lcnt,lu,ifl)
      implicit none
      integer ir, ids2
      REAL    tbls(1)
      integer lcnt, lu, ifl
c---------------------------------------------------------------
c
c  subroutine:  s2get(ir,ids2,tbls,lcnt,lu,ifl)
c
c  purpose:     load eos tables from library
c
c  arguments:   ir    (input) - material region number
c               ids2  (input) - sesame 2 id number
c               tbls (in/out) - array for storage of tables
c               lcnt (in/out) - current word in array, updated
c                               by this routine
c               lu    (input) - logical unit number for library
c               ifl  (output) - error flag
c                 ifl = 1, if data is loaded with no problems
c                 ifl = 0, if data cannot be located
c                 ifl = -number of extra words needed, if
c                       storage allocated is insufficient
c
c               user also supplies a common block having the
c               form:  common/s2dir/lcmx,nrs,lcfw
c                 lcmx = length of array tbls
c                 nrs  = number of regions in problem,
c                        including non-sesame regions
c                 lcfw = array used for directory to tbls
c               this common block also used by routine s2eos.
c
c  remarks:     loads data into array tbls, beginning with word
c               tbls(lcnt) and computing new value of lcnt.
c
c  externals:   inbuf, tbfch.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        4 october 1977
c
c---------------------------------------------------------------
c-----local variables
      integer id, i, loc, itest, nmax, ieof, nwds
c-----end local variables
c-----common block variables
      integer lcmx, nrs, lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      if(nrs .lt. 1) then
	ifl = 0
	return
      endif
      id = iabs(ids2)
      do i=1,nrs
        loc = lcfw(i)
        if(loc .eq. 0) go to 1
        itest = tbls(loc)
        if(id .eq. itest) then
          lcfw(ir) = loc
          ifl = 1
          return
	endif
 1    enddo
      nmax = lcmx-lcnt+1
      ifl = nmax-5
      if(ifl .lt. 0) return
      rewind lu
      itest = -1
      do while (id .ne. itest)
	ieof = 1
	do while (ieof .eq. 1)
          call inbuf(lu,tbls(lcnt),1,ieof)
	enddo
        call inbuf(lu,tbls(lcnt),5,ieof)
        if(ieof .eq. 0) then
	  ifl = 0
	  return
        endif
        itest = tbls(lcnt)
      enddo
      nwds = tbls(lcnt+4)
      nwds = nwds+nwds+5
      ifl = nmax-nwds
      if(ifl .lt. 0) return
      backspace lu
      call inbuf(lu,tbls(lcnt),nwds,ieof)
      if(ieof .eq. 0) then
	ifl = 0
	return
      endif
      backspace lu
      loc = lcnt+nwds
      nmax = nmax-nwds
      call tbfch(TWO0ONE,tbls(loc),tbls(lcnt),lu,nmax,ifl)
      if(ifl .lt. 1) return
      tbls(lcnt+1) = tbls(loc+2)
      nmax = lcmx-lcnt-1
      call tbfch(THREE0ONE,tbls(lcnt+2),tbls(lcnt),lu,nmax,ifl)
      if(ifl .lt. 1) return
      loc = lcnt
      tbls(lcnt) = id
      lcnt = lcnt+2+ifl
 3    lcfw(ir) = loc
      ifl = 1
      return
      end
      subroutine s2eos(ir,tbls,r,t,p,e)
      implicit none
      integer ir
      REAL    r, t, tbls(1), p(3), e(3)
c---------------------------------------------------------------
c
c  subroutine:  s2eos(ir,tbls,r,t,p,e)
c
c  purpose:     compute eos from sesame 2 tables
c
c  arguments:   ir   (input) - material region number
c               tbls (input) - array for storage of tables
c               r    (input) - density in mg/m**3
c               t    (input) - temperature in kelvins
c               p,e (output) - pressure, internal energy vectors
c                 p(1),e(1) = pressure(gpa), energy(mj/kg)
c                 p(2),e(2) = density derivatives
c                 p(3),e(3) = temperature derivatives
c
c  remarks:
c
c  externals:   la301a.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        23 december 1976
c
c---------------------------------------------------------------
c-----local variables
      integer loc
c-----end local variables
c-----common block variables
      integer lcmx, nrs, lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      loc = lcfw(ir)+2
      call la301a(tbls(loc),r,t,p,e)
      return
      end
      subroutine s2geti(ir,ids2,tbls,lcnt,lu,ifl)
      implicit none
      integer ir, ids2
      REAL    tbls(1)
      integer lcnt, lu, ifl
c---------------------------------------------------------------
c
c  subroutine:  s2geti(ir,ids2,tbls,lcnt,lu,ifl)
c
c  purpose:     load inverted eos tables from library
c
c  arguments:   ir    (input) - material region number
c               ids2  (input) - sesame 2 id number
c               tbls (in/out) - array for storage of tables
c               lcnt (in/out) - current word in array, updated
c                               by this routine
c               lu    (input) - logical unit number for library
c               ifl  (output) - error flag
c                 ifl = 1, if data is loaded with no problems
c                 ifl = 0, if data cannot be located
c                 ifl = -number of extra words needed, if
c                       storage allocated is insufficient
c
c               user also supplies a common block having the
c               form:  common/s2dir/lcmx,nrs,lcfw
c                 lcmx = length of array tbls
c                 nrs  = number of regions in problem,
c                        including non-sesame regions
c                 lcfw = array used for directory to tbls
c               this common block also used by routine s2eosi.
c
c  remarks:     loads data into array tbls, beginning with word
c               tbls(lcnt) and computing new value of lcnt.
c
c  externals:   inbuf, tbfch, inv301.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        4 october 1977
c
c---------------------------------------------------------------
c-----local variables
      integer loc, id, i, itest, nmax, ieof, nwds
      integer nr, nt, lds
      REAL    r0
c-----end local variables
c-----common block variables
      integer lcmx, nrs, lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      if(nrs .lt. 1) then
        ifl = 0
        return
      endif
      id = -iabs(ids2)
      do i=1,nrs
        loc = lcfw(i)
        if(loc .eq. 0) go to 1
        itest = tbls(loc)
        if(id .eq. itest) then
          lcfw(ir) = loc
          ifl = 1
          return
	endif
 1    enddo
      nmax = lcmx-lcnt+1
      ifl = nmax-5
      if(ifl .lt. 0) return
      rewind lu
      itest = -1
      do while (id .ne. itest)
	ieof = 1
	do while (ieof .eq. 1)
         call inbuf(lu,tbls(lcnt),1,ieof)
	enddo
        call inbuf(lu,tbls(lcnt),5,ieof)
        if(ieof .eq. 0) then
          ifl = 0
          return
        endif
        itest = -tbls(lcnt)
      enddo
      nwds = tbls(lcnt+4)
      nwds = nwds+nwds+5
      ifl = nmax-nwds
      if(ifl .lt. 0) return
      backspace lu
      call inbuf(lu,tbls(lcnt),nwds,ieof)
      if(ieof .eq. 0) then
        ifl = 0
        return
      endif
      backspace lu
      loc = lcnt+nwds
      nmax = nmax-nwds
      call tbfch(TWO0ONE,tbls(loc),tbls(lcnt),lu,nmax,ifl)
      if(ifl .lt. 1) return
      tbls(lcnt+1) = tbls(loc+2)
      nmax = lcmx-lcnt-1
      call tbfch(THREE0ONE,tbls(lcnt+2),tbls(lcnt),lu,nmax,ifl)
      if(ifl .lt. 1) return
      loc = lcnt
      nr = tbls(loc+2)
      nt = tbls(loc+3)
      nwds = ifl+nr+3*nt
      ifl = nmax-nwds
      if(ifl .lt. 0) return
      r0 = tbls(loc+1)
      call inv301(tbls(loc+2),r0,lds)
      tbls(lcnt) = id
      lcnt = lcnt+2+lds
      lcfw(ir) = loc
      ifl = 1
      return
      end
      subroutine s2eosi(ir,tbls,r,e,p,t)
      implicit none
      integer ir
      REAL    tbls(1), r, e, p(3), t(3)
c---------------------------------------------------------------
c
c  subroutine:  s2eosi(ir,tbls,r,e,p,t)
c
c  purpose:     compute inverted eos from sesame 2 tables
c
c  arguments:   ir   (input) - material region number
c               tbls (input) - array for storage of tables
c               r    (input) - density in mg/m**3
c               e    (input) - internal energy in mj/kg
c               p,t (output) - pressure, temperature vectors
c                 p(1),t(1) = pressure(gpa), temperature(kelvin)
c                 p(2),t(2) = density derivatives
c                 p(3),t(3) = energy derivatives
c
c  remarks:
c
c  externals:   la302a.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        4 january 1977
c
c---------------------------------------------------------------
c-----local variables
      integer loc
c-----end local variables
c-----common block variables
      integer lcmx, nrs, lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      loc = lcfw(ir)+2
      call la302a(tbls(loc),r,e,p,t)
      return
      end
      subroutine s2hugi(ir,tbls,v0,v,mfl)
      implicit none
      integer ir
      REAL    tbls(1), v0(3), v(6)
      integer mfl
c---------------------------------------------------------------
c
c  subroutine:  s2hugi(ir,tbls,v0,v,mfl)
c
c  purpose:     compute a point on the hugoniot
c
c  arguments:   ir   (input) - material region number
c               tbls (input) - array for storage of tables
c               v0   (input) - initial state vector
c                 v0(1) = density, mg/m**3
c                 v0(2) = pressure, gpa
c                 v0(3) = internal energy, mj/kg
c               v   (in/out) - final state vector
c                 v(1) = density (input)
c                 v(2) = pressure
c                 v(3) = internal energy
c                 v(4) = temperature, kelvins
c                 v(5) = shock velocity, km/sec
c                 v(6) = particle velocity, km/sec
c               mfl (output) - error flag
c                 mfl = 1,  if calculation is okay
c                 mfl = 0,  if shock velocity is imaginary
c                 mfl = -1, if no solution can be found
c
c  remarks:
c
c  externals:   s2eosi.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        12 january 1977
c
c---------------------------------------------------------------
c-----local variables
      REAL p(3), t(3)
      REAL a, b, e, etest
c-----end local variables
      b = HALF*(ONE/v0(1)-ONE/v(1))
      a = v0(3)+b*v0(2)
      v(3) = v0(3)
      etest = BIG10
      do while (etest .gt. ZERO)
        e = v(3)
        call s2eosi(ir,tbls,v(1),e,p,t)
        if(p(1) .gt. BIG10) then
          mfl = -1
          return
        endif
        v(2) = p(1)
        v(4) = t(1)
        v(3) = a+b*v(2)
        etest = abs(v(3)-e)-EPS6*(abs(v(3))+EPS3)
      enddo
      v(6) = TWO*b*(v(2)-v0(2))
      if(v(6) .le. ZERO) then
        mfl = 0
        v(5) = ZERO
        v(6) = ZERO
      else
        mfl = 1
        v(6) = sqrt(v(6))
        v(5) = HALF*v(6)/v0(1)/b
      endif
      return
      end
      subroutine s2shki(ir,tbls,v0,v,mfl)
      implicit none
      integer ir
      REAL    tbls(1), v0(3), v(6)
      integer mfl
c---------------------------------------------------------------
c
c  subroutine:  s2shki(ir,tbls,v0,v,mfl)
c
c  purpose:     compute a point on the hugoniot
c
c  arguments:   ir   (input) - material region number
c               tbls (input) - array for storage of tables
c               v0   (input) - initial state vector
c                 v0(1) = density, mg/m**3
c                 v0(2) = pressure, gpa
c                 v0(3) = internal energy, mj/kg
c               v   (in/out) - final state vector
c                 v(1) = density
c                 v(2) = pressure
c                 v(3) = internal energy (input)
c                 v(4) = temperature, kelvins
c                 v(5) = shock velocity, km/sec
c                 v(6) = particle velocity, km/sec
c               mfl (output) - error flag
c                 mfl = 1, if calculation is okay
c                 mfl = 0, if no solution can be found
c
c  remarks:
c
c  externals:   s2eosi.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        23 march 1977
c
c---------------------------------------------------------------
c-----local variables
      integer n
      REAL p(3), r(3), f(3)
      REAL b, x
c-----end local variables
      mfl = 0
      b = HALF/(v0(1)*(v(3)-v0(3)))
      r(2) = v0(1)
      f(2) = v0(1)
      do while (f(2) .gt. ZERO)
        r(1) = r(2)
        f(1) = f(2)
        r(2) = ONEPOINTTWO*r(2)
        call s2eosi(ir,tbls,r(2),v(3),p,v(4))
        f(2) = r(2)-b*(r(2)-v0(1))*(p(1)+v0(2))
      enddo
      r(3) = r(2)
      do n=1,250
        v(1) = r(3)
        x = f(2)/(f(2)-f(1))
        if(x*(ONE-x) .lt. TENTH) x=HALF
        r(3) = r(2)+x*(r(1)-r(2))
        if(abs(ONE-r(3)/v(1)) .lt. EPS6) then
          v(2) = p(1)
          v(6) = sqrt((p(1)-v0(2))*(ONE/v0(1)-ONE/v(1)))
          v(5) = v(6)/(ONE-v0(1)/v(1))
          mfl = 1
          return
        endif
        call s2eosi(ir,tbls,r(3),v(3),p,v(4))
        f(3) = r(3)-b*(r(3)-v0(1))*(p(1)+v0(2))
        if(f(1)*f(3) .gt. ZERO) then
          r(1) = r(3)
          f(1) = f(3)
        else
          r(2) = r(3)
          f(2) = f(3)
	endif
      enddo
      return
      end
      subroutine s2abti(ir,tbls,v0,v,mfl)
      implicit none
      integer ir
      REAL    tbls(1), v0(3), v(6)
      integer mfl
c---------------------------------------------------------------
c
c  subroutine:  s2abti(ir,tbls,v0,v,mfl)
c
c  purpose:     compute a point on the adiabat
c
c  arguments:   ir   (input) - material region number
c               tbls (input) - array for storage of tables
c               v0   (input) - initial state vector
c                 v0(1) = density, mg/m**3
c                 v0(2) = pressure, gpa
c                 v0(3) = internal energy, mj/kg
c               v   (in/out) - final state vector
c                 v(1) = density (input)
c                 v(2) = pressure
c                 v(3) = internal energy
c                 v(4) = temperature, kelvins
c                 v(5) = sound velocity, km/sec
c                 v(6) = particle velocity, km/sec
c               mfl (output) - error flag
c                 mfl = 1, if calculation is okay
c                 mfl = 0, if sound velocity is imaginary
c                 mfl = -1, if the integration did not converge
c
c  remarks:
c
c  externals:   s2eosi.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        7 february 1977
c
c---------------------------------------------------------------
c-----local variables
      integer i, n, m
      REAL    p(3), t(3)
      REAL    esave
      REAL    a, b, e, e1, etest
      REAL    c, r, r0, dr, up
c-----end local variables
      n = 5
      dr = FIFTH*(v(1)-v0(1))
      if(abs(dr) .lt. EPS6) then
        call s2eosi(ir,tbls,v(1),v0(3),v(2),t)
        c = v(3)+v(2)*v(4)/v(1)**2
        mfl = 1
        if(c .le. ZERO) then
          c = ZERO
          mfl = 0
	endif
        v(5) = sqrt(c)
        v(3) = v0(3)
        v(4) = t(1)
        v(6) = ZERO
        return
      endif
      do m=1,10
        r0 = v0(1)
        v(2) = v0(2)
        v(3) = v0(3)
        c = ZERO
        v(6) = ZERO
        mfl = 1
        do i=1,n
          r = r0+dr
          b = HALF*(ONE/r0-ONE/r)
          a = v(3)+b*v(2)
          e = a+b*dr*c
	  etest = BIG10
	  do while (etest .gt. ZERO)
            e1 = e
            call s2eosi(ir,tbls,r,e1,p,t)
            e = a+b*p(1)
            etest = abs(e-e1)-EPS6*(abs(e)+EPS3)
	  enddo
          v(5) = c
          c = (p(1)-v(2))/dr
          if(c .le. ZERO) then
            mfl = 0
            c = ZERO
	  endif
          up = v(6)-sqrt(c/r/r0)*dr
          r0 = r
          v(2) = p(1)
          v(3) = e
          v(4) = t(1)
          v(6) = up
        enddo
        c = c+HALF*(c-v(5))
        if(n .eq. 5) then
          esave = v(3)
          n = n+n
          dr = HALF*dr
        else
          etest = abs(v(3)-esave)-FIVEEPS6*(abs(v(3))+EPS3)
          if(etest .lt. ZERO) then
            if(c .gt. ZERO) then
              v(5) = sqrt(c)
	    else
              mfl = 0
              v(5) = ZERO
            endif
            return
          endif
        endif
      enddo
      mfl = -1
      if(c .gt. ZERO) then
        v(5) = sqrt(c)
      else
        v(5) = ZERO
      endif
      return
      end
      subroutine la301a(tbls,r,t,p,e)
      implicit none
      REAL tbls(*), r, t, p(3), e(3)
c---------------------------------------------------------------
c
c  subroutine:  la301a(tbls,r,t,p,e)
c
c  purpose:     search/interpolate for pressure and internal
c               energy as functions of density and temperature.
c
c  arguments:   tbls (input) - data string of type 301
c               r    (input) - density
c               t    (input) - temperature
c               p   (output) - pressure vector of length 3
c               e   (output) - energy vector of length 3
c                 p(1),e(1) = pressure and internal energy
c                 p(2),e(2) = density derivatives
c                 p(3),e(3) = temperature derivatives
c
c  remarks:
c
c  externals:   isrch, ratfn2.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        19 november 1976
c
c---------------------------------------------------------------
c-----local variables
c-----end local variables
c-----external variables
      integer isrch
c-----end external variables
c-----common block variables
      integer locx, i, nr, locy, j, nt, locz, nz
      common/rtblk2/locx,i,nr,locy,j,nt,locz,nz
c-----end common block variables
      nz = 1
      nr = tbls(1)
      nt = tbls(2)
      locx = 3
      locy = locx+nr
      i = isrch(r,tbls(locx),nr)
      if(i*(nr-i) .le. 0) then
        if(i .lt. 1) i=1
        if(i .gt. nr-1) i=nr-1
      endif
      j = isrch(t,tbls(locy),nt)
      if(j*(nt-j) .le. 0) then
        if(j .lt. 1) j=1
        if(j .gt. nt-1) j=nt-1
      endif
      locz = locy+nt
      call ratfn2(r,t,tbls,p)
      locz = locz+nr*nt
      call ratfn2(r,t,tbls,e)
      return
      end
      subroutine la302a(tbls,r,e,p,t)
      implicit none
      REAL tbls(*), r, e, p(3), t(3)
c---------------------------------------------------------------
c
c  subroutine:  la302a(tbls,r,e,p,t)
c
c  purpose:     search/interpolate for pressure and temperature
c               as functions of density and internal energy.
c
c  arguments:   tbls (input) - data string of type 302
c               r    (input) - density
c               e    (input) - energy
c               p   (output) - pressure vector of length 3
c               t   (output) - temperature vector of length 3
c                 p(1),t(1) = pressure and temperature
c                 p(2),t(2) = density derivatives
c                 p(3),t(3) = energy derivatives
c
c  remarks:
c
c  externals:   isrch, ratfn1, ratfn2.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        3 january 1977
c
c---------------------------------------------------------------
c-----local variables
      REAL et, decdr
c-----end local variables
c-----external variables
      integer isrch
c-----end external variables
c-----common block variables
      integer locx, kx, locy, ky, ix, n
      common/rtblk1/locx,kx,locy,ky,ix,n
      integer locr, ir, nr, loce, ie, ne, locz, nz
      common/rtblk2/locr,ir,nr,loce,ie,ne,locz,nz
c-----end common block variables
      nr = tbls(1)
      ne = tbls(2)
      locx = 3
      kx = 1
      locy = locx+nr+ne
      ky = 1
      ix = isrch(r,tbls(locx),nr)
      if(ix*(nr-ix) .le. 0) then
        if(ix .lt. 1) ix=1
        if(ix .gt. nr-1) ix=nr-1
      endif
      n = nr
      call ratfn1(r,tbls,p)
      et = e-p(1)
      decdr = p(2)
      if(et .lt. ZERO) et=ZERO
      locr = locx
      ir = ix
      loce = locx+nr
      ie = isrch(et,tbls(loce),ne)
      if(ie*(ne-ie) .le. 0) then
        if(ie .lt. 1) ie=1
        if(ie .gt. ne-1) ie=ne-1
      endif
      locz = locy+nr
      nz = 1
      call ratfn2(r,et,tbls,p)
      p(2) = p(2)-decdr*p(3)
      locz = locz+nr*ne
      call ratfn2(r,et,tbls,t)
      t(2) = t(2)-decdr*t(3)
      return
      end
      subroutine wdw301(dstr,rmin,rmax,tmin,tmax,loc)
      implicit none
      REAL dstr(*), rmin, rmax, tmin, tmax
      integer loc
c---------------------------------------------------------------
c
c  subroutine:  wdw301(dstr,rmin,rmax,tmin,tmax,loc)
c
c  purpose:     contract data string of type 301 by deleting
c               data outside specified ranges.
c
c  arguments:   dstr (input) - data string to be contracted
c               rmin (input) - lower density limit
c               rmax (input) - upper density limit
c               tmin (input) - lower temperature limit
c               tmax (input) - upper temperature limit
c               loc (output) - length of new data string
c
c  remarks:
c
c  externals:   isrch
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        4 october 1977
c
c---------------------------------------------------------------
c-----local variables
      integer nr, nt, imin, imax, jmin, jmax
      integer i, j
c-----end local variables
c-----external variables
      integer isrch
c-----end external variables
      nr = dstr(1)
      nt = dstr(2)
      imin = isrch(rmin,dstr(4),nr-2)+1
      imax = isrch(rmax,dstr(4),nr-2)+1
      if(rmax .gt. dstr(imax+2)) imax=imax+1
      if(imax .lt. imin+3) imax=imin+3
      do while (imax .gt. nr)
        imax = imax-1
        imin = imin-1
      enddo
      jmin = isrch(tmin,dstr(nr+4),nt-2)+1
      jmax = isrch(tmax,dstr(nr+4),nt-2)+1
      if(tmax .gt. dstr(nr+jmax+2)) jmax=jmax+1
      if(jmin .eq. 1) jmin=2
      if(jmax .lt. jmin+2) jmax=jmin+2
      do while (jmax .gt. nt)
        jmax = jmax-1
        jmin = jmin-1
      enddo
      loc = 2
      do i=imin,imax
        loc = loc+1
        dstr(loc) = dstr(i+2)
      enddo
      loc = loc+1
      dstr(loc) = dstr(3+nr)
      do j=jmin,jmax
        loc = loc+1
        dstr(loc) = dstr(nr+j+2)
      enddo
      do i=imin,imax
        loc = loc+1
        dstr(loc) = dstr(nr+nt+i+2)
      enddo
      do j=jmin,jmax
        do i=imin,imax
          loc = loc+1
          dstr(loc) = dstr(nt+i+2+nr*j)
	enddo
      enddo
      do i=imin,imax
        loc = loc+1
        dstr(loc) = dstr(nr+nt+nr*nt+i+2)
      enddo
      do j=jmin,jmax
        do i=imin,imax
          loc = loc+1
          dstr(loc) = dstr(nt+nr*(nt+j)+i+2)
	enddo
      enddo
      nr = imax-imin+1
      nt = jmax-jmin+2
      dstr(1) = nr
      dstr(2) = nt
      return
      end
      subroutine inv301(dstr,r0,lds)
      implicit none
      REAL     dstr(*), r0
      integer lds
c---------------------------------------------------------------
c
c  subroutine:  inv301(dstr,r0,lds)
c
c  purpose:     invert data string of type 301 to type 302.
c
c  arguments:   dstr (input) - data string to be inverted
c               r0   (input) - approximate density of solid
c               lds (output) - length of new data string
c
c  remarks:
c               this routine overwrites locations following the
c               data string.  it expands the string by nr words,
c               where nr is the number of densities.  it also
c               uses 3*nt words as temporary storage, where nt
c               is the number of temperatures.
c
c  externals:   isrchk, ratfn1.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        4 october 1977
c
c---------------------------------------------------------------
c-----local variables
      integer loct, lcec, locp, loce, locn, imax
      integer i, j, jj
      REAL z(2), q, et
c-----end local variables
c-----external variables
      integer isrchk
c-----end external variables
c-----common block variables
      integer locx, nr, locy, ky, jx, nt
      common/rtblk1/locx,nr,locy,ky,jx,nt
c-----end common block variables
      nr = dstr(1)
      nt = dstr(2)
      loct = 3+nr
      lcec = loct+nt
      locp = lcec+nr
      loce = locp+nr*nt
      locn = loce+nr*nt
      imax = 2*nr*nt
      do i=1,imax
        dstr(locn-i) = dstr(locn-i-nr)
      enddo
      do i=1,nr
        jj = loce+i-1
        q = EPS12*abs(dstr(jj))
        dstr(lcec+i-1) = dstr(jj)
        dstr(jj) = ZERO
        do j=2,nt
          jj = jj+nr
          dstr(jj) = dstr(jj)-dstr(lcec+i-1)
          if(dstr(jj)-dstr(jj-nr) .lt. q) dstr(jj)=dstr(jj-nr)+q
	enddo
      enddo
      i = isrchk(r0,dstr(4),nr-2,1)+1
      do j=1,nt
        dstr(locn+j-1) = dstr(loct+j-1)
        dstr(loct+j-1) = dstr(loce+i-1+nr*(j-1))
      enddo
      do i=1,nr
        locx = loce+i-1
        do j=1,nt
          et = dstr(loct+j-1)
          jx = isrchk(et,dstr(locx+nr),nt-2,nr)+1
          locy = locp+i-1
          ky = nr
          call ratfn1(et,dstr,z)
          dstr(locn+nt+j-1) = z(1)
          locy = locn
          ky = 1
          call ratfn1(et,dstr,z)
          dstr(locn+nt+nt+j-1) = z(1)
	  enddo
        do j=1,nt
          dstr(locp+i-1+nr*(j-1)) = dstr(locn+nt+j-1)
          dstr(locx+nr*(j-1)) = dstr(locn+nt+nt+j-1)
	enddo
      enddo
      lds = locn-1
      return
      end
      function isrch(x,tbls,n)
      implicit none
      integer isrch
      REAL    x, tbls(1)
      integer n
c---------------------------------------------------------------
c
c  function:    isrch(x,tbls,n)
c
c  purpose:     find index of x in an array tbls.
c
c  arguments:   x    (input) - value to be located
c               tbls (input) - table to be searched
c               n    (input) - number of values to be searched
c               the value of the function = index i, where
c                 tbls(i).le.x.lt.tbls(i+1), or
c                 i=0 or i=n if x is outside range of table.
c
c  remarks:
c               table values must be contiguous and in ascending
c               order.  for more general routine, see isrchk.
c
c  externals:   none.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        19 november 1976
c
c---------------------------------------------------------------
c-----local variables
      integer j, jp
c-----end local variables
      isrch = 0
      j = n+1
      do while (j-isrch .ne. 1)
        jp = HALF*(j+isrch)
        if(x .lt. tbls(jp)) then
          j = jp
        else
          isrch = jp
        endif
      enddo
      return
      end
      function isrchk(x,tbls,n,k)
      implicit none
      integer isrchk
      REAL    x, tbls(1)
      integer n, k
c---------------------------------------------------------------
c
c  function:    isrchk(x,tbls,n,k)
c
c  purpose:     find index of x in an array tbls.  table values
c               need not be contiguous and can be in either
c               ascending or descending order.
c
c  arguments:   x    (input) - value to be located
c               tbls (input) - table to be searched
c               n    (input) - number of values to be searched
c               k    (input) - spacing between values in table
c               the value of the function = index i, where
c                 tbls(1+k*(i-1)).le.x.lt.tbls(1+k*i), or
c                 tbls(1+k*(i-1)).ge.x.gt.tbls(1+k*i), or
c                 i=0 or i=n if x is outside range of table.
c
c  remarks:
c
c  externals:   none.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        19 november 1976
c
c---------------------------------------------------------------
c-----local variables
      integer j, ki, jp
      REAL    s
c-----end local variables
      isrchk = 0
      j = n+1
      ki = 1-k
      s = tbls(ki+k*n)-tbls(1)
      do while (j-isrchk .ne. 1)
        jp = HALF*(j+isrchk)
        if(s*(x-tbls(ki+k*jp)) .lt. ZERO) then
          j = jp
        else
          isrchk = jp
        endif
      enddo
      return
      end
      subroutine ratfn1(x,tbls,y)
      implicit none
      REAL x, tbls(1), y(2)
c---------------------------------------------------------------
c
c  subroutine:  ratfn1(x,tbls,y)
c
c  purpose:     interpolate for a function y(x) and its
c               derivative from tables located in array tbls.
c
c  arguments:   x    (input) - independent variable
c               tbls (input) - array containing tables
c               y   (output) - vector of length 2, where
c                 y(1) = value of function
c                 y(2) = derivative of function
c               the routine also requires a common block,
c               common/rtblk1/locx,kx,locy,ky,i,n
c                 locx = location of x vector
c                 kx   = spacing of x vector
c                 locy = location of y vector
c                 ky   = spacing of y vector
c                 i    = index of x and y vectors
c                 n    = length of x and y vectors
c
c  remarks:     uses rational function method with quadratic
c               estimate of derivatives at the mesh points.
c
c  externals:   none, but a search routine must be called
c               first, to compute index i.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        18 november 1976
c
c---------------------------------------------------------------
c-----local variables
      integer ix, iy
      REAL    q, d, r, s, sp, c2, dm, sm, c1, c3, c4
c-----end local variables
c-----common block variables
      integer locx, kx, locy, ky, i, n
      common/rtblk1/locx,kx,locy,ky,i,n
c-----end common block variables
      ix = locx+kx*(i-1)
      iy = locy+ky*(i-1)
      q = x-tbls(ix)
      d = tbls(ix+kx)-tbls(ix)
      r = d-q
      s = (tbls(iy+ky)-tbls(iy))/d
      if(i .le. 1) then
        sp = (tbls(iy+ky+ky)-tbls(iy+ky))/(tbls(ix+kx+kx)-tbls(ix+kx))
        c2 = (sp-s)/(tbls(ix+kx+kx)-tbls(ix))
        if(s*(s-d*c2) .le. ZERO) c2=s/d
        y(1) = tbls(iy)+q*(s-r*c2)
        y(2) = s+(q-r)*c2
        return
      endif
      dm = tbls(ix)-tbls(ix-kx)
      sm = (tbls(iy)-tbls(iy-ky))/dm
      c1 = (s-sm)/(d+dm)
      if(i .ge. n-1) then
        y(1) = tbls(iy)+q*(s-r*c1)
        y(2) = s+(q-r)*c1
        return
      endif
      if(i .le. 2) then
        if(sm*(sm-dm*c1) .le. ZERO) c1=(s-sm-sm)/d
      endif
      sp = (tbls(iy+ky+ky)-tbls(iy+ky))/(tbls(ix+kx+kx)-tbls(ix+kx))
      c2 = (sp-s)/(tbls(ix+kx+kx)-tbls(ix))
      c3 = abs(c2*r)
      c4 = c3+abs(c1*q)
      if(c4 .gt. ZERO) c3=c3/c4
      c4 = c2+c3*(c1-c2)
      y(1) = tbls(iy)+q*(s-r*c4)
      y(2) = s+(q-r)*c4+d*(c4-c2)*(ONE-c3)
      return
      end
      subroutine ratfn2(x,y,tbls,z)
      implicit none
      REAL x, y, tbls(1), z(3)
c---------------------------------------------------------------
c
c  subroutine:  ratfn2(x,y,tbls,z)
c
c  purpose:     interpolate for a function z(x,y) and its
c               derivatives from tables located in array tbls.
c
c  arguments:   x,y  (input) - independent variables
c               tbls (input) - array containing tables
c               z   (output) - vector of length 3, where
c                 z(1) = value of function
c                 z(2) = x derivative of function
c                 z(3) = y derivative of function
c               the routine also requires a common block,
c               common/rtblk2/locx,ix,nx,locy,iy,ny,locz,nz
c                 locx = location of x vector
c                 ix   = index of x vector
c                 nx   = length of x vector
c                 locy = location of y vector
c                 iy   = index of y vector
c                 ny   = length of y vector
c                 locz = location of z(x,y) array
c                 nz   = spacing of z array
c               a second common block - common/intord/ifn -
c               is used to specify the order of interpolation.
c               set ifn = 1 for the bilinear form.
c
c  remarks:     unless bilinear form is specified, routine
c               uses rational function method with quadratic
c               estimate of derivatives at the mesh points.
c
c  externals:   none, but a search routine must be called
c               first, to compute indices ix and iy.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        24 february 1977
c
c---------------------------------------------------------------
c-----local variables
      integer k, i, iz, kz, ibr, nbr
      REAL    zz(4), zd(4)
      REAL    q, d, r, s, sp, c2, dm, sm, c1, c3, c4
      REAL    dx, qx
c-----end local variables
c-----common block variables
      integer locx, ix, nx, locy, iy, ny, locz, nz
      common/rtblk2/locx,ix,nx,locy,iy,ny,locz,nz
      integer ifn
c     common/intord/ifn
c-----end common block variables
      data ifn/0/
      k = 0
      i = locx+ix-1
      iz = locz+nz*(ix-1+nx*(iy-1))
      kz = nz
      ibr = ix
      nbr = nx-ix
      q = x-tbls(i)
 1    d = tbls(i+1)-tbls(i)
      r = d-q
 2    k = k+1
      s = (tbls(iz+kz)-tbls(iz))/d
      zz(k) = HALF*(tbls(iz)+q*s)
      zd(k) = HALF*s
      if(ifn .eq. 1) go to 6
      if(ibr .gt. 1) go to 3
      sp = (tbls(iz+kz+kz)-tbls(iz+kz))/(tbls(i+2)-tbls(i+1))
      c2 = (sp-s)/(tbls(i+2)-tbls(i))
      if(s*(s-d*c2) .le. ZERO) c2=s/d
      zz(k) = zz(k)-q*r*c2
      zd(k) = zd(k)+(q-r)*c2
      go to (7,8,9,10),k
 3    dm = tbls(i)-tbls(i-1)
      sm = (tbls(iz)-tbls(iz-kz))/dm
      c1 = (s-sm)/(d+dm)
      if(nbr .gt. 1) go to 4
      zz(k) = zz(k)-q*r*c1
      zd(k) = zd(k)+(q-r)*c1
      go to (7,8,9,10),k
 4    if(ibr .gt. 2) go to 5
      if(sm*(sm-dm*c1) .le. ZERO) c1=(s-sm-sm)/d
 5    sp = (tbls(iz+kz+kz)-tbls(iz+kz))/(tbls(i+2)-tbls(i+1))
      c2 = (sp-s)/(tbls(i+2)-tbls(i))
      c3 = abs(c2*r)
      c4 = c3+abs(c1*q)
      if(c4 .gt. ZERO) c3=c3/c4
      c4 = c2+c3*(c1-c2)
      zz(k) = zz(k)-q*r*c4
      zd(k) = zd(k)+(q-r)*c4+d*(c4-c2)*(ONE-c3)
 6    go to (7,8,9,10),k
 7    iz = iz+nx*nz
      go to 2
 8    dx = d
      qx = q
      i = locy+iy-1
      kz = kz*nx
      iz = iz-kz
      ibr = iy
      nbr = ny-iy
      q = y-tbls(i)
      go to 1
 9    iz = iz+nz
      go to 2
 10   zz(2) = (zz(2)-zz(1))/d
      zz(4) = (zz(4)-zz(3))/dx
      z(1) = zz(1)+zz(2)*q+zz(3)+zz(4)*qx
      z(2) = zz(4)+zd(1)+(zd(2)-zd(1))*q/d
      z(3) = zz(2)+zd(3)+(zd(4)-zd(3))*qx/dx
      return
      end
      subroutine tbfch(ctnr,z,xtb,lu,nmax,num)
      implicit none
      REAL ctnr, z(1), xtb(1)
      integer lu, nmax, num
c---------------------------------------------------------------
c
c  subroutine:  tbfch(ctnr,z,xtb,lu,nmax,num)
c
c  purpose:     get a data table from sesame 2 library.
c
c  arguments:   ctnr (input) - catalogue number of data table
c               z   (in/out) - array where table is to be stored
c               xtb  (input) - material index
c               lu   (input) - logical unit number for library
c               nmax (input) - maximum length of array z
c               num (output) - number of words loaded into z
c                 num = 0, if the table cannot be found
c                 if nmax is too small, num = -number of extra
c                 words needed
c
c  remarks:
c               before calling, the pointer on unit lu must be
c               set to the first word of the material file.
c
c  externals:   inbuf.
c
c  programmer:  g. t. rood and g. i. kerley, t-4.
c
c  date:        13 december 1976
c
c---------------------------------------------------------------
c-----local variables
      integer i, j, k, ieof
c-----end local variables
      j = 5
      k = xtb(j)
      do i=1,k
        if(ctnr .eq. xtb(j+i)) then
          num = xtb(j+i+k)
          if(num .gt. nmax) then
           num = nmax-num
           return
          endif
          do k=1,i
            call inbuf(lu,z,1,ieof)
            if(ieof .eq. 0) then
              num = 0
              return
	    endif
          enddo
          call inbuf(lu,z,num,ieof)
          if(ieof .eq. 0) then
            num = 0
            return
          endif
          do k=1,i+1
            backspace lu
          enddo
          return
	endif
      enddo
      num = 0
      return
      end
#if defined(MAKE_SES_LIB)
      subroutine rdcard (lu,tbl,tblid,nwds,nwpc,fmt)
      implicit none
      integer lu
      REAL    tbl(1), tblid
      integer nwds, nwpc
      character fmt*(*)
c----------------------------------------------------------------
c
c  subroutine   rdcard (lu,tbl,tblid,nwds,nwpc,fmt)
c
c  purpose      routine reads bcd card image data from logical
c               unit -lu-, decides whether hollerith or numeric,
c               and if numeric, it decides whether integer or
c               floating point.
c
c  arguments    lu    (in) - logical unit of bcd input file
c               tbl  (out) - scm array to hold data
c               tblid (in) - table id number
c               nwds  (in) - number of words to be read
c               nwpc  (in) - number of words per card image
c               fmt   (in) - bcd format (8a10-cdc, 20a4-ibm)
c
c  remarks      if tblid.le.19999 data is assumed hollerith, and
c               -nwds- is number of hollerith characters to be read.
c               if tblid.ge.20000 data is assumed numeric.  five
c               integer flags are read with each card to indicate
c               type of each corresponding number.  all data
c               numbers are read with e15.0 format, but if the
c               integer flag says it should be integer type, it
c               is converted and stored back in itself.
c
c  programmer   c. forest,  g. rood,   t-4
c
c  date         jan 1977
c
c----------------------------------------------------------------
c-----local variables
      REAL    tmp(5)
      INTEGER imp(5)
      REAL    rkr
      INTEGER irkr
      REAL    rkp
      INTEGER irkp
      integer iflg(5)
      integer ncard, ncds, i, j, k, l, m
      equivalence (tmp,imp)
      equivalence(rkr,irkr),(rkp,irkp)
c-----end local variables
      if(tblid .gt. ONENINENINE) then
c------ read numeric table
        m=-5
        ncard=(nwds+4)/5
        do i=1,ncard
          m=m+5
          read (lu,100) (tmp(k),k=1,5), (iflg(k),k=1,5)
  100     format(5e15.0,5i1)
          do j=1,5
            if(iflg(j) .eq. 2) imp(j)=tmp(j)
            tbl(j+m)=tmp(j)
	  enddo
        enddo
        return
      endif
c------ read hollerith data
      ncds=(nwds+nwpc-1)/nwpc
      j=1-nwpc
      do l=1,ncds
        j=j+nwpc
        k=j+nwpc-1
        if(k .gt. nwds) k=nwds
        call rcrdline(lu,tbl,j,k,fmt)/*JWG*/
c       read (lu,fmt) (tbl(i),i=j,k)/*JWG*/
      enddo
      return
      end
cJWG--
      subroutine rcrdline(lu,tbl,j,k,fmt)
      implicit none
      integer lu
      REAL    tbl(1)
      integer j, k
      character fmt*(*)
c----------------------------------------------------------------
c
c  subroutine:  rcrdline(lu,fmt,tbl,j,k)
c
c  purpose:     Transfer the contents of one character line into
c               the array tbl between indicies j and k
c
c  arguments    lu    (in) - logical unit of bcd input file
c               tbl  (out) - scm array to hold data
c               nwpc  (in) - number of words per card image
c               fmt   (in) - bcd format (8a10-cdc, 20a4-ibm)
c----------------------------------------------------------------
c-----local variables
      integer i
      character*80 cline
      REAL         rline(NWPC)
      equivalence(cline,rline)
c-----end local variables
      read (lu,'(a80)') cline
      do i=j,k
        tbl(i) = rline(i-j+1)
      end do
      return
      end
cJWG--
      subroutine libnit (dir,lb)
      implicit none
      REAL    dir(1)
      integer lb
c----------------------------------------------------------------
c
c  subroutine:  libnit (dir,lb)
c
c  purpose:     to initialize  *sesame ii*  library by reading
c               master library directory from file -lb- into scm.
c
c  arguments:   dir (out) - scm array to contain master library
c                           directory
c               lb   (in) - logical unit of library file
c
c  remarks:     routine rewinds library file, buffers in first
c               three words of directory.  number of words in
c               second record to be buffered in is calculated
c               from first word of first record, then second
c               record is buffered in.
c
c  externals:   inbuf
c
c  programmer:  g. rood,   t-4
c
c  date:        1 sep 1976
c
c               system independent
c
c----------------------------------------------------------------
c-----local variables
      integer j, ieof, i2
c-----end local variables
      j=3
      rewind lb
      call inbuf(lb,dir,j,ieof)
      if(ieof .eq. 0) then
        dir(1)=ZERO
        dir(2)=ZERO
        dir(3)=ZERO
        return
      endif
      i2=dir(1)
      if(i2  .eq.  ZERO) return
      i2=3*i2
      call inbuf(lb,dir(j+1),i2,ieof)
      return
      end
      subroutine skipmf (lu,nf)
      implicit none
      integer lu, nf
c----------------------------------------------------------------
c
c  subroutine:  skipmf (lu,nf)
c
c  purpose:     to skip over -nf- material files.
c
c  arguments:   lu (in) - logical unit number
c               nf (in) - number of files to be skipped over
c
c  remarks:     first checks for end-of-information, then reads
c               one word from each succeeding record until
c               end-of-file is read.  loops until it has done
c               this -nf- times, or read eoi, whichever occurs
c               first.  if initially positioned immediately before
c               an eof, then skipmf simply skips one file mark and
c               returns independent of the value of nf.
c
c  externals:   inbuf
c
c  programmer:  g. rood,   t-4
c
c  date:        21 jan 1977
c
c----------------------------------------------------------------
c-----local variables
      integer i, ieof
      REAL    z(1)
c-----end local variables
      do i=1,nf
        call inbuf(lu,z,1,ieof)
        if(ieof .eq. 0) return
	do while (ieof .ne. 0)
          call inbuf(lu,z,1,ieof)
	enddo
      enddo
      return
      end
      subroutine la401a(tbls,y,j,ifl)
      implicit none
      REAL    tbls(1), y(12)
      integer j, ifl
c---------------------------------------------------------------
c
c  subroutine:  la401a(tbls,y,j,ifl)
c
c  purpose:     calculate quantities on vaporization/sublimation
c               curve as functions of a specified variable.
c
c  arguments:   tbls (input) - data string of type 401
c               y   (in/out) - dimensioned array for variables
c                              and their derivatives
c                 y(1),y(2)   = vapor pressure and derivative
c                 y(3),y(4)   = temperature and derivative
c                 y(5),y(6)   = vapor density and derivative
c                 y(7),y(8)   = liquid density and derivative
c                 y(9),y(10)  = vapor energy and derivative
c                 y(11),y(12) = liquid energy and derivative
c               j    (input) - specifies y(2*j-1) to be the
c                             independent variable
c               ifl (output) - error flag
c                 ifl = 0  if out of range.  values on the
c                          boundary are returned.
c                 ifl.gt.0 if in range
c
c  remarks:     certain derivatives are rigorously zero or
c               infinite at the critical point.  this routine
c               will not reproduce such behavior.
c
c  externals:   isrchk, ratfn1.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        26 may 1977
c
c---------------------------------------------------------------
c-----local variables
      REAL    x
      integer i
c-----end local variables
c-----external variables
      integer isrch, isrchk
c-----end external variables
c-----common block variables
      integer locx, kx, locy, ky, ix, n
      common/rtblk1/locx,kx,locy,ky,ix,n
c-----end common block variables
      n = tbls(1)
      x = y(j+j-1)
      y(j+j) = ONE
      locx = 2+n*(j-1)
      kx = 1
      ky = 1
      ix = isrchk(x,tbls(locx),n,1)
      ifl = ix*(n-ix)
      if(ifl .gt. 0) go to 2
      if(ix .eq. n) go to 1
      ix = 1
      x = tbls(locx)
      go to 2
 1    ix = ix-1
      x = tbls(locx+ix)
 2    locy = 2
      do i=1,6
        if(i .eq. j) go to 3
        call ratfn1(x,tbls,y(i+i-1))
c  if out of range, derivatives are set = ZERO
        if(ifl .eq. 0) y(i+i)=ZERO
        locy = locy+n
 3    enddo
      return
      end
      subroutine la302b(tbls,r,e,p)
      implicit none
      REAL tbls(1), r, e(2), p(2)
c---------------------------------------------------------------
c
c  subroutine:  la302b(tbls,r,e,p)
c
c  purpose:     calculate cold curve from data string 302
c
c  arguments:   tbls (input) - data string of type 302
c               r    (input) - density
c               e,p (output) - doubly-dimensioned arrays
c                 e(1),p(1) = internal energy and pressure
c                 e(2),p(2) = density derivatives
c
c  remarks:
c
c  externals:   isrch, ratfn1.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        25 may 1977
c
c---------------------------------------------------------------
c-----local variables
      integer ne
c-----end local variables
c-----external variables
      integer isrch
c-----end external variables
c-----common block variables
      integer locx, kx, locy, ky, i, nr
      common/rtblk1/locx,kx,locy,ky,i,nr
c-----end common block variables
      nr = tbls(1)
      ne = tbls(2)
      locx = 3
      kx = 1
      locy = locx+nr+ne
      ky = 1
      i = isrch(r,tbls(locx),nr)
      if(i*(nr-i) .le. 0) then
        if(i .lt. 1) i=1
        if(i .gt. nr-1) i=nr-1
      endif
      call ratfn1(r,tbls,e)
      locy = locy+nr
      call ratfn1(r,tbls,p)
      return
      end
      subroutine la302s(tbls,r,e,p,t,ibr)
      implicit none
      REAL    tbls(1), r, e, p(3), t(3)
      integer ibr
c---------------------------------------------------------------
c
c  subroutine:  la302s(tbls,r,e,p,t,ibr)
c
c  purpose:     search/interpolate for pressure and temperature
c               as functions of density and internal energy.
c
c  arguments:   tbls (input) - data string of type 302
c               r    (input) - density
c               e    (input) - energy
c               p   (output) - pressure vector of length 3
c               t   (output) - temperature vector of length 3
c                 p(1),t(1) = pressure and temperature
c                 p(2),t(2) = density derivatives
c                 p(3),t(3) = energy derivatives
c               ibr  (input) - branch parameter
c                 ibr = 0, calculate pressure and temperature
c                 ibr = 1, calculate temperature only
c                 ibr = 2, calculate pressure only.
c
c  remarks:
c               adapted from sesame 2 subroutine la302a.
c
c  externals:   isrch, ratfn1, ratfn2.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        9 september 1977
c
c---------------------------------------------------------------
c-----local variables
      REAL z(2)
      REAL et, decdr
c-----end local variables
c-----external variables
      integer isrch
c-----end external variables
c-----common block variables
      integer locx, kx, locy, ky, ix, n
      common/rtblk1/locx,kx,locy,ky,ix,n
      integer locr, ir, nr, loce, ie, ne, locz, nz
      common/rtblk2/locr,ir,nr,loce,ie,ne,locz,nz
c-----end common block variables
      nr = tbls(1)
      ne = tbls(2)
      locx = 3
      kx = 1
      locy = locx+nr+ne
      ky = 1
      ix = isrch(r,tbls(locx),nr)
      if(ix*(nr-ix) .le. 0) then
        if(ix .lt. 1) ix=1
        if(ix .gt. nr-1) ix=nr-1
      endif
      n = nr
      call ratfn1(r,tbls,z)
      et = e-z(1)
      decdr = z(2)
      if(et .lt. ZERO) et=ZERO
      locr = locx
      ir = ix
      loce = locx+nr
      ie = isrch(et,tbls(loce),ne)
      if(ie*(ne-ie) .le. 0) then
        if(ie .lt. 1) ie=1
        if(ie .gt. ne-1) ie=ne-1
      endif
      locz = locy+nr
      nz = 1
      if(ibr .ne. 1) then
        call ratfn2(r,et,tbls,p)
        p(2) = p(2)-decdr*p(3)
        if(ibr .eq. 2) return
      endif
      locz = locz+nr*ne
      call ratfn2(r,et,tbls,t)
      t(2) = t(2)-decdr*t(3)
      return
      end
      subroutine tprb1(ids2,n,lu,lp)
      implicit none
      integer ids2(1), n, lu, lp
c---------------------------------------------------------------
c
c  subroutine:  tprb1(ids2,n,lu,lp)
c
c  purpose:     calculate and print hugoniots for materials
c               specified by the user.  for illustration
c               and testing of sesame library procedures.
c
c  arguments:   ids2 (input) - list of material numbers
c               n    (input) - number of materials specified
c               lu   (input) - unit number of sesame library
c               lp   (input) - unit number of line printer
c
c  externals:   s2geti, s2shki.
c
c  programmer:  g. i. kerley, t-4.
c
c  date:        14 october 1977
c
c---------------------------------------------------------------
c-----local variables
      integer loc, ir, lcnt, ifl
      integer k, mflag, i, j
      REAL    fk
      REAL    rd(3), r(6,50)
c-----end local variables
c-----common block variables
c  storage for the tables is in common block eosdat
      REAL tbls(12000)
      common/eosdat/tbls
c  common block s2dir is a directory
      integer lcmx, nrs, lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      lcmx = 12000
      nrs = n
      if(nrs .gt. 10) nrs=10
      do ir=1,nrs
        lcfw(ir) = 0
      enddo
c  load eos tables from library, on unit lu, into common
c  block eosdat.  s2geti is used to obtain internal energy
c  as an independent variable.
      lcnt = 1
      do ir=1,nrs
        call s2geti(ir,ids2(ir),tbls,lcnt,lu,ifl)
        if(ifl .lt. 0) go to 4
        if(ifl .gt. 0) go to 3
        write(lp,2)ids2(ir)
 2      format(1x,8hmaterial,i5,21h not found in library)
 3    enddo
      go to 5
 4    nrs = ir-1
c  hugoniot calculation. see routine s2shki for an example
c  of use of the routine s2eosi.
 5    do ir=1,nrs
        loc = lcfw(ir)
        if(loc .eq. 0) go to 10
c  initial state for hugoniot calculation
        rd(1) = tbls(loc+1)
        rd(2) = ZERO
        rd(3) = ZERO
c  hugoniot calculation using s2shki
        k = 0
 6      k = k+1
        fk = k
        r(3,k) = HALF*fk*fk
        call s2shki(ir,tbls,rd,r(1,k),mflag)
c  if calculation fails, temperature is set to negative value
        if(mflag .eq. 0) r(4,k)=ERRTEMP
        if(r(2,k) .ge. ONETHOUSAND) go to 7
        if(k .lt. 50) go to 6
c  print results
 7      write(lp,8)ids2(ir),(r(i,1),i=1,6)
 8      format(//12x,41hshock hugoniot for sesame material number,
     1    i5//7x,1hr,10x,1hp,10x,1he,10x,1ht,10x,2hus,9x,2hup/5x,
     2    5hmg/m3,7x,3hgpa,7x,5hmj/kg,5x,6hkelvin,6x,6hkm/sec,5x,
     3    6hkm/sec//1x,6e11.3)
        write(lp,9)((r(i,j),i=1,6),j=2,k)
 9      format(1x,6e11.3)
 10   enddo
      if(n .eq. nrs) return
      write(lp,11)
 11   format(/1x,40hinsufficient storage for other materials)
      return
      end
      subroutine copmts(l1,l2,n,index)
      implicit none
      integer l1, l2, n, index(1)
c----------------------------------------------------------------------
c
c  subroutine  copmts(l1,l2,n,index)
c
c  purpose     to copy a material file to another unit
c
c  arguments   l1        (input)   copy from unit no.
c              l2        (input)   copy to unit no.
c              n         (input)   no. of material to be copied
c              index     (input)   no. of words in index records
c
c  remarks     none
c
c  externals   copmat
c
c  programmer  j.abdallah,jr.
c
c  date         26 feb 1978
c
c----------------------------------------------------------------------
c
c-----local variables
      integer j
c-----end local variables
c-----common block variables
      REAL buff(20000)
      common/sescom/buff
c-----end common block variables
      rewind l1
      rewind l2
      do j=1,n
        if(j .ne. 1) call skipmf(l1,1)
        call copmat(l1,l2,buff,index(j))
      enddo
      call outbuf(l2,buff,0)
      return
      end
      subroutine update(date,ncpw,fmt,inew,lup,lib,lm,lnew,lp)
      implicit none
      REAL date
      integer ncpw
      character fmt*(*)
      integer inew, lup, lib, lm, lnew, lp
c----------------------------------------------------------------------
c
c  subroutine  update(date,ncpw,fmt,inew,lup,lib,lm,lnew,lp)
c
c  purpose     to create / update a sesame ii library
c
c  arguments   date      (input)   floating pt. date mmddyy
c              ncpw      (input)   no. of hollerith characters/word
c              fmt       (input)   format for reading hollerith card
c              inew      (input)   1=update old library
c                        (input)   0=create new library
c              lup       (input)   coded update file unit no.
c              lib       (input)   current library unit no.
c              lm        (input)   unit no. of scratch file
c              lnew      (input)   unit no, of created/updated
c                                  library
c              lp        (input)   printer unit no.
c
c  remarks     none
c
c  externals   s2bldm,s2srtm,s2bldl,s2updt,copmts
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c------------------------------------------------------------------------
c-----local variables
      REAL     rmats(400)
      integer index(400)
      integer nmats
c-----end local variables
      write(lp,1000)
c . . build material file from update file
      call s2bldm(lup,lm,lnew,lp,inew,date,ncpw,fmt,nmats,rmats,index)
c . . if no materials found go back
      if(nmats .le. 0) go to 100
c . . sort material file according to material and table id
      write(lp,1010)
      call s2srtm(lm,lnew,nmats,rmats,index)
      if(inew .eq. 1) go to 10
c . . build new library here
      call copmts(lnew,lm,nmats,index)
      call s2bldl(lm,lnew,date,ONE)
      write(lp,1050)
      return
c . . update existing library here
 10    write(lp,1030)
      call s2updt(lnew,lib,lm,date,nmats,rmats,index)
      write(lp,1060)
      return
 100  write(lp,1040)
      return
 1000 format(1x,39hcreating material file from update file)
 1010 format(1x,21hsorting material file)
 1020 format(1x,20hcreating new library)
 1030 format(1x,25hupdating existing library)
 1040 format(1x,31herror in processing update file)
 1050 format(1x,17hcreation complete)
 1060 format(1x,15hupdate complete)
      end
      subroutine s2bldm(lup,lm,lscr,lp,inew,date,ncpw,fmt,nmats,
     $                  rmats,index)
      implicit none
      integer lup
      integer lm
      integer lscr, lp, inew
      REAL    date
      integer ncpw
      character fmt*(*)
      integer nmats
      REAL    rmats(1)
      integer index(1)
c----------------------------------------------------------------------
c
c  subroutine  s2bldm(lup,lm,lscr,lp,inew,date,ncpw,fmt,nmats,rmats,index)
c
c  purpose     to build a material file from a coded
c              update file
c
c  arguments   lup       (input)   update file unit no.
c              lm        (input)   output material file unit no.
c              lscr      (input)   scratch file unit no.
c              lp        (input)   printer unit no.
c              inew      (input)   0=create new library
c                                  1=update existing library
c              date      (input)   floating pt. date
c              ncpw      (input)   no of characters per word
c              fmt       (input)   format for reading card image
c              nmats     (output)  number of materials obtained on lm
c                                  0=error
c              rmats     (output)  array of material ids
c              index     (output)  array no of words in material indices
c
c
c  remarks     a material file is a file in sesame format except
c              that the  file directory has been excluded
c
c  externals   retmat,bdindx,inbuf,outbuf,upsum
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      REAL hindex(35), tblids(15)
      REAL rmat, date1, date2, vers
      integer nwtab(15)
      integer i, ntab, iflag
      integer nwi, nw
      integer j, im
c-----end local variables
c-----common block variables
      REAL buff(20000)
      common/sescom/buff
c-----end common block variables
      rewind lm
      rewind lup
      i=0
 10   call retmat(lup,lscr,ncpw,fmt,rmat,date1,date2,vers,
     $           ntab,tblids,nwtab,iflag)
      if(iflag .ne. 0) go to 1000
      if(date1 .ne. ZERO) go to 20
      date1=date
      date2=date
      go to 30
 20   if(date2 .eq. ZERO) date2=date
 30   if((inew .eq. 0) .and. (vers .eq. ZERO)) vers=ONE
      i=i+1
      rmats(i)=rmat
      call bdindx(rmat,date1,date2,vers,ntab,tblids,nwtab,hindex)
      call upsum(lp,i,rmat,date1,ntab,tblids,nwtab)
      nwi=ntab+ntab+5
      index(i)=nwi
      call outbuf(lm,hindex,nwi)
      rewind lscr
      do j=1,ntab
        nw=nwtab(j)
        if(nw .eq. 0) go to 100
        call inbuf(lscr,buff,nw,im)
        call outbuf(lm,buff,nw)
 100  enddo
      call outbuf(lm,buff,0)
      go to 10
 1000 nmats=i
      call outbuf(lm,buff,0)
      return
      end
      subroutine retmat(lup,ls,ncpw,fmt,rmat,date1,date2,vers,
     $                 ntab,tblids,nwtab,iflag)
      implicit none
      integer   lup
      integer   ls
      integer   ncpw
      character fmt*(*)
      REAL      rmat, date1, date2, vers
      integer   ntab
      REAL      tblids(1)
      integer   nwtab(1), iflag
c-----------------------------------------------------------------------
c
c  subroutine  retmat(lup,ls,ncpw,fmt,rmat,date1,date2,vers,
c                     ntab,tblids,nwtab,iflag)
c
c  purpose     to retrieve the data for one material from the
c              update file
c
c  arguments   lup       (input)   update file unit no.
c              ls        (input)   unit to which tables are written
c              ncpw      (input)   no. of chars. / word
c              fmt       (input)  format for reading card image
c              date1     (output)  creation date
c              date2     (output)  date of last update
c              vers      (output)  version no.
c              rmat      (output)  material id
c              ntab      output)  no. of tables
c              tblids    (output)  array of table ids
c              nwtab     (output)  no of words per table
c              iflag     (output)  1=end of file
c
c
c  remarks     data tables are written on ls in the
c              order they are read off the update file
c
c  externals   rdcard,outbuf
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c-----local variables
      integer idel, jflag, ifm
      integer nwdc
      integer idis, nwpc, j
      REAL ratid, tblid
      REAL datea, dateb, uvers
      save ifm, ratid, tblid, nwdc, idis, datea, dateb, uvers
c-----end local variables
c-----common block variables
      REAL buff(20000)
      common/sescom/buff
c-----end common block variables
      data idel/1hd/
      data jflag/0/
      iflag=0
      ntab=0
c . . read marker card
      rewind ls
      if (jflag .ne. 0) go to 110
      jflag=1
 10   read(lup,100) ifm,ratid,tblid,nwdc,idis,datea,dateb,uvers
 100  format(i2,2f6.0,i6,3x,a1,2f9.0,f4.0)
 110  if((ntab .eq. 0) .and. (ifm .eq. 2)) then
        iflag=1
        return
      endif
      if(ifm .eq. 2) return
      if((ntab .ne. 0) .and. (ifm .eq. 0)) return
      date1=datea
      date2=dateb
      vers=uvers
      ntab=ntab+1
      tblids(ntab)=tblid
      if(tblid .lt. TWOHUNDRED) nwdc=(nwdc+ncpw-1)/ncpw
      nwtab(ntab)=nwdc
      nwpc=80/ncpw
      if(nwpc*ncpw .lt. 80) nwpc=nwpc+1
      rmat=ratid
      if(idis .eq. idel) then
        nwtab(ntab)=0
      else
        call rdcard(lup,buff,tblid,nwdc,nwpc,fmt)
        call outbuf(ls,buff,nwdc)
      endif
      go to 10
      end
      subroutine bdindx(rmat,date1,date2,vers,ntab,tblids,nwtab,hindex)
      implicit none
      REAL    rmat, date1, date2, vers
      integer ntab
      REAL    tblids(1)
      integer nwtab(1)
      REAL    hindex(1)
c----------------------------------------------------------------------
c
c  subroutine  bdindx(rmat,date1,date2,vers,ntab,tblids,nwtab,hindex)
c
c  purpose     to build the material index
c
c  arguments   rmat      (input)   material i.d.
c              date1     (input)   creation date
c              date2     (input)   date of last update
c              vers      (input)   version number
c              ntab      (input)   no. of tables
c              tblids    (input)   array of table ids
c              nwtab     (input)   array of no. of words/table
c              hindex    (output)   floating point material index
c
c  remarks     none
c
c  externals   none
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer j
c-----end local variables
      hindex(1)=rmat
      hindex(2)=date1
      hindex(3)=date2
      hindex(4)=vers
      hindex(5)=I2FLOAT(ntab)
      do j=1,ntab
        hindex(5+j)=tblids(j)
        hindex(5+j+ntab)=I2FLOAT(nwtab(j))
      enddo
      return
      end
      subroutine s2srtm(lm,lms,n,rmats,index)
      implicit none
      integer lm
      integer lms
      integer n
      REAL    rmats(1)
      integer index(1)
c----------------------------------------------------------------------
c
c  subroutine  s2srtm(lm,lms,n,rmats,index)
c
c  purpose     to sort material file in increasing material and
c              table ids
c
c  arguments   lm        (input)   unit no. of input material file
c              lms       (input)   unit no. of output sorted file
c              n         (input)   no.of materials on file
c              rmats     (input)   array of material ids
c              index     (input)   no of words in each material index
c
c  remarks     the call to s2srtm may be skipped if file is already
c              sorted
c
c  externals   asort2,inbuf,outbuf,tbfch,compak
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer j, ipt1, nx
      integer nwds
      integer iflg, nt, ipt2, k
      integer nt3, ipt3, k5
      integer ifl, ntw
      REAL    rm, t1
c-----end local variables
c-----common block variables
      REAL buff(20000)
      common/sescom/buff
c-----end common block variables
      rewind lm
      rewind lms
      do j=1,n
        buff(j)=index(j)
      enddo
      call asort2(n,rmats,buff(1))
      ipt1=n+1
      nx=0
      do j=1,n
        rm=rmats(j)
 10     if(nx .ne. 0) go to 12
        rewind lm
        go to 14
 12     call skipmf(lm,1)
 14     nx=nx+1
        nwds=index(nx)
        call inbuf(lm,buff(ipt1),nwds,iflg)
        if(nx .ne. n) go to 15
        nx=0
 15     if(buff(ipt1) .eq. rm) go to 20
        go to 10
 20     nt=buff(ipt1+4)
        ipt2=ipt1+nwds
        if(nt .eq. 1) go to 50
        t1=buff(ipt1+5)
        do k=2,nt
          if(buff(ipt1+4+k) .lt. t1) go to 100
          t1=buff(ipt1+4+k)
        enddo
c . . no sort on tables required - copy them
 50     call outbuf(lms,buff(ipt1),nwds)
        do k=1,nt
          ntw=buff(ipt1+4+k+nt )
          if(ntw .eq. 0) go to 80
          call inbuf(lm,buff(ipt2),ntw,iflg)
          call outbuf(lms,buff(ipt2),ntw)
 80     enddo
        go to 900
c . . table sort required
 100    call compak(buff(ipt1),buff(ipt2))
        nt3=buff(ipt2+4)
        ipt3=ipt2+nt3+nt3+5
        call asort2(nt,buff(ipt1+5),buff(ipt1+5+nt))
        call outbuf(lms,buff(ipt1),nwds)
        backspace lm
        do k=1,nt
          k5=ipt1+4+k
          if(buff(k5+nt) .eq. ZERO) go to 130
          call tbfch(buff(k5),buff(ipt3),buff(ipt2),lm,9765,ifl)
          call outbuf(lms,buff(ipt3),ifl)
 130    enddo
 900    call outbuf(lms,buff,0)
 1000 enddo
      call outbuf(lms,buff,0)
      do j=1,n
        index(j)=buff(j)
      enddo
      return
      end
      subroutine s2bldl(lm,lnew,date,vers)
      implicit none
      integer lm, lnew
      REAL    date, vers
c----------------------------------------------------------------------
c
c  subroutine  s2bldl(lm,lnew,date,vers)
c
c  purpose     to build  a sesame ii library from a material file
c
c
c  arguments   lm        (input)   material file unit no.
c              lnew      (input)   library unit no.
c              date      (input)   floating point date
c              vers      (input)   version number
c
c  remarks     s2bldl builds the master directory from the material
c              file and copies it and the data from the
c              material file to the library unit  to create a
c              sesame ii type file
c
c  externals   skipmf,inbuf,copmat
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer iad, im, in, ia, iptm, iptn, ipta, nm
      integer iflg, nt, nindx, j, nw, nm3, ndir
      integer ipt, nwi, n, il, n3, ipts, nsp
      integer nut
c-----end local variables
c-----common block variables
      REAL buff(12000), bufi(2000)
      INTEGER ibuff(6000)
      common/sescom/buff,ibuff,bufi
c-----end common block variables
      iad=1
      im=3
      in=3000
      ia=0
      iptm=im
      iptn=in
      ipta=ia
      rewind lm
c . . build tables for directory
      nm=0
      go to 2
 1    call skipmf(lm,1)
 2    call inbuf(lm,bufi,5,iflg)
      if(iflg .eq. 0) go to 1000
      nm=nm+1
      iptm=iptm+1
      iptn=iptn+1
      ipta=ipta+1
      ibuff(ipta)=iad
      buff(iptm)=bufi(1)
      nt=bufi(5)
      nindx=5+nt+nt
      backspace lm
      call inbuf(lm,bufi,nindx,iflg)
      buff(iptn)=I2FLOAT(nindx)
      iad=iad+nindx+1
      do j=1,nt
        nw=bufi(5+nt+j)
        iad=iad+nw+1
      enddo
      iad=iad+1
      go to 1
c . . move directory to proper places
 1000 buff(1)=I2FLOAT(nm)
      buff(2)=date
      buff(3)=vers
      rewind lnew
      call outbuf(lnew,buff,3)
      nm3=3*nm
      ndir=4+nm3+2
c . . build second record of directory
c . . move material no. table
      ipt=0
      do j=1,nm
        ipt=ipt+1
        im=im+1
        buff(ipt)=buff(im)
      enddo
c . . move index table
      iptn=in
      do j=1,nm
        ipt=ipt+1
        iptn=iptn+1
        buff(ipt)=buff(iptn)
      enddo
c . . move address table
      do j=1,nm
        ipt=ipt+1
        ia=ia+1
        iad=ibuff(ia)+ndir
        buff(ipt)=I2FLOAT(iad)
      enddo
      call outbuf(lnew,buff,nm3)
      call outbuf(lnew,buff,0)
c . . copy data tables from lm to lnew
      rewind lm
      iptn=in
      do j=1,nm
        iptn=iptn+1
          buff(j)=buff(iptn)
      enddo
      ipt=nm+1
      do j=1,nm
        if(j .ne. 1) call skipmf(lm,1)
        nwi=buff(j)
        call copmat(lm,lnew,buff(ipt),nwi)
      enddo
      call outbuf(lnew,buff,0)
      return
      end
      subroutine s2updt(lm,lib,lscr,date,nm,rmats,index)
      implicit none
      integer lm, lib, lscr
      REAL    date
      integer nm
      REAL    rmats(1)
      integer index(1)
c----------------------------------------------------------------------
c
c
c  subroutine  s2updt(lm,lib,lscr,date,nm,rmats,index)
c
c  purpose     to update the current sesame library and build
c              a new sesame library
c
c
c  arguments   lm        (input)   unit no. of material file containing
c                                  updates - it is replaced by the
c                                  sesame ii library upon exit
c              lib       (input)   unit no of current library
c              lscr      (input)   unit no. of scratch file
c              nm        (input)   no of materials on material file
c              rmats     (input)   array of material ids on material
c                                  file
c              index     (input)   no. of words in each material index
c                                  on material file
c              date      (input)   floating point date
c
c  remarks     none
c
c  externals   inbuf,outbuf,uptpgm,skipmf,s2mod,copmat,s2bldl
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c-----local variables
      integer iupdt(200)
      integer iflg, im, n, il, n3, ipts, nsp, nut, j
      integer idj, nwi
      REAL    vers
c-----end local variables
c-----common block variables
      REAL buff(20000)
      common/sescom/buff
c-----end common block variables
c . . get sesame directory
      rewind lib
      rewind lm
      rewind lscr
      call inbuf(lib,buff,3,iflg)
      im=0
      n=buff(1)
      vers=buff(3)
      il=0
      n3=n+n+n
      ipts=n+1
      nsp=15000-ipts +1
      call inbuf(lib,buff(4),n3,iflg)
      call uptpgm(n,buff(4),nm,rmats,nut,iupdt)
      do j=1,n
        buff(j)=buff(3+n+j)
      enddo
      do j=1,nut
        idj=iupdt(j)
        if(idj) 600,400,500
 400    il=il+1
        im=im+1
        nwi=buff(il)
        if(im .gt. 1) call skipmf(lm,1)
        call skipmf(lib,1)
        call s2mod(lm,lib,lscr,buff(ipts),nsp,index(im),nwi)
        go to 1000
 500    im=im+1
        if(im .gt. 1) call skipmf(lm,1)
        call copmat(lm,lscr,buff(ipts),index(im))
        go to 1000
 600    il=il+1
        call skipmf(lib,1)
        nwi=buff(il)
        call copmat(lib,lscr,buff(ipts),nwi)
 1000 enddo
      call outbuf(lscr,buff,0)
      vers=vers+ONE
      call s2bldl(lscr,lm,date,vers)
      return
      end
      subroutine copmat(l1,l2,store,nwi)
      implicit none
      integer l1, l2
      REAL    store(1)
      integer nwi
c----------------------------------------------------------------------
c
c  subroutine  copmat(l1,l2,store,nwi)
c
c  purpose     to copy a material from one file to another
c
c  arguments   l1        (input)   copy from unit no.
c              l2        (input)   copy to unit no.
c              store     (input)   storage area for data records
c              nwi       (input)   no. of words in material index
c                                  if nwi.le.0 the number of words
c                                  is computed
c
c  remarks     l1 and l2 must point to the index record
c
c  externals   inbuf,outbuf
c
c  programmer   j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer nw, iflg, nt, ipt, j
c-----end local variables
      nw=5
      if(nwi .gt. 0) nw=nwi
      call inbuf(l1,store,nw,iflg)
      nt=store(5)
      if(nwi .gt. 0) go to 10
      nw=5+nt+nt
      backspace l1
      call inbuf(l1,store,nw,iflg)
 10   ipt=nw+1
      call outbuf(l2,store,nw)
c . . write remaining data tables
      do j=1,nt
        nw=store(5+j+nt)
        call inbuf(l1,store(ipt),nw,iflg)
        call outbuf(l2,store(ipt),nw)
      enddo
      call outbuf(l2,store,0)
      return
      end
      subroutine s2mod(lm,lib,lmu,store,lens,nwm,nwl)
      implicit none
      integer lm, lib, lmu
      REAL    store(1)
      integer lens, nwm, nwl
c----------------------------------------------------------------------
c
c  subroutine  s2mod(lm,lib,lmu,store,lens,nwm,nwl)
c
c  purpose     to modify the library material file using the
c              binary update file
c
c  arguments   lm        (input)   unit no. of update file
c              lib       (input)   library version of material unit no.
c              lmu       (input)   unit no of updated material
c              store     (input)   storage area
c              lens      (input)   length of storage area
c              nwm       (input)   no. of words in mat. index on update
c                                  file
c              nwl       (input)   no. of words in index of library
c                                  version
c
c  remarks     lm and lib must point to the index record
c
c  externals   combin,compak,inbuf, outbuf,tbfch
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer im1, im, il, in1, in, istx, lnx, iflg
      integer nt, nwds, k, num
      REAL    tid
c-----end local variables
      im1=1
      im=1+nwm
      il=im+nwm
      in1=il+nwl
      in=in1+nwl+nwm+1
      istx=in+nwl+nwm
      lnx=lens-istx+1
      call inbuf(lm,store(im1),nwm,iflg)
      call inbuf(lib,store(il),nwl,iflg)
      call combin(store(im1),store(il),store(in1))
      call compak(store(im1),store(im))
      call compak(store(in1),store(in))
c . . check for total deletion
      nt=store(in+4)
      if(nt .eq. 0) return
      nwds=nt+nt+5
      call outbuf(lmu,store(in),nwds)
      backspace lm
      backspace lib
      do k=1,nt
        tid=store(in+4+k)
        if(store(im+4) .eq. ZERO) go to 60
c . . search update file for table
        call tbfch(tid,store(istx),store(im),lm,lnx,num)
        if(num .gt. 0) go to 80
c . .   search lib file for table
 60     call tbfch(tid,store(istx),store(il),lib,lnx,num)
 80     call outbuf(lmu,store(istx),num)
      enddo
      call outbuf(lmu,store,0)
      return
      end
      subroutine compak(dexa,dexb)
      implicit none
      REAL dexa(1), dexb(1)
c----------------------------------------------------------------------
c
c  subroutine  compak(dexa,dexb)
c
c  purpose     to compact material index by removing
c              table deletion indicators
c
c  arguments   dexa      (input)   index to be compacted
c              dexb      (output)  compacted index
c
c  remarks     none
c
c  externals   none
c
c  programmer  j.abdallah
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c . . copy first four words of index
c-----local variables
      integer na, nb, j, j5, k5, ipta, ipt, i
c-----end local variables
      dexb(1)=dexa(1)
      dexb(2)=dexa(2)
      dexb(3)=dexa(3)
      dexb(4)=dexa(4)
c . . compute the number of tables in the compacted index
      na=dexa(5)
      nb=0
      do j=1,na
        if(dexa(j+na+5) .eq. ZERO) go to 10
        nb=nb+1
 10   enddo
      dexb(5)=I2FLOAT(nb)
      if(nb .eq. 0) return
c . . move table ids into index
      j5=5
      k5=5
      ipta=5+na
      ipt=5+nb
      i=0
      do j=1,na
        ipta=ipta+1
        j5=j5+1
        if(dexa(ipta) .eq. ZERO) go to 30
        i=i+1
        ipt=ipt+1
        k5=k5+1
        dexb(ipt)=dexa(ipta)
        dexb(k5)=dexa(j5)
        if(i .ge. nb) go to 40
 30   enddo
 40   continue
      return
      end
      subroutine uptpgm(nl,rmatl,nm,rmatm,nut,ipgm)
      implicit none
      integer nl
      REAL    rmatl(1)
      integer nm
      REAL    rmatm(1)
      integer nut, ipgm(1)
c----------------------------------------------------------------------
c
c  subroutine  uptpgm(nl,rmatl,nm,rmatm,nut,ipgm)
c
c  purpose     uptpgm(nl,rmatl,nm,rmatm,nut,ipgm)
c
c  arguments
c              nl        (input)   no. of materials on library
c              rmatl     (input)   material ids on library
c              nm        (input)   no. of materials on update file
c              rmatm    (input)    array of materials on update
c                                  file
c              nut       (output)  no. of materials in updated
c                                  library
c              ipgm      (output)  program array
c                                  -1=use old material as is
c                                   0=update material
c                                   1=use new material as is
c
c  remarks     the ipgm array is used to control processing in
c              s2updt
c
c  externals   none
c
c  programmer  j.abdallah
c
c  date      20 feb 1979
c
c----------------------------------------------------------------------
c-----local variables
      integer il, im
      REAL    rsm, rsl
c-----end local variables
      il=1
      im=1
      nut=0
 10   rsm=rmatm(im)
      rsl=rmatl(il)
      nut=nut+1
      if(rsl-rsm) 100,200,300
 100  il=il+1
      ipgm(nut)=-1
      go to 500
 200  ipgm(nut)=0
      il=il+1
      im=im+1
      go to 500
 300  ipgm(nut)=1
      im=im+1
 500  if((im .le. nm) .and. (il .le. nl)) go to 10
      if((im .gt. nm) .and. (il .gt. nl)) go to 600
      nut=nut+1
      if(im .gt. nm) go to 100
      if(il .gt. nl) go to 300
 600  continue
      return
      end
      subroutine combin(dexm,dexl,dexn)
      implicit none
      REAL dexm(1), dexl(1), dexn(1)
c----------------------------------------------------------------------
c
c  subroutine  combin(dexm,dexl,dexn)
c
c  purpose     to combine the material index of the update file
c              with the index of the library version to form a
c              new material index
c
c  arguments   dexm      (input)   index of material on update file
c              dexl      (input)   index of material on library
c              dexn      (output)  new updated index
c
c  remarks     the resulting index includes deletion indicators
c              (if any)
c
c  externals   none
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c . . compute the total number of tables
c . . in updated material
c-----local variables
      integer ntm, ntl, ntot, k, j, im, il
      integer nwm, nwl
      REAL    tm, tl, vers
c-----end local variables
      ntm=dexm(5)
      ntl=dexl(5)
      ntot=ntl
      do k=1,ntm
        tm=dexm(5+k)
        do j=1,ntl
          if(tm .eq. dexl(j+5)) go to 20
        enddo
        ntot=ntot+1
 20   enddo
c . . build table array
      im=1
      il=1
      do j=1,ntot
        tm=dexm(5+im)
        nwm=dexm(5+im+ntm)
        tl=dexl(5+il)
        nwl=dexl(5+il+ntl)
        if(im .gt. ntm) go to 180
        if(il .gt. ntl) go to 140
        if(tm-tl) 140,160,180
 140    im=im+1
        dexn(5+j)=tm
        dexn(5+j+ntot)=I2FLOAT(nwm)
        go to 500
 160    im=im+1
        il=il+1
        dexn(j+5)=tm
        dexn(j+5+ntot)=I2FLOAT(nwm)
        go to 500
 180    il=il+1
        dexn(5+j)=tl
        dexn(5+j+ntot)=I2FLOAT(nwl)
 500  enddo
c . . move the remainder of the new index to appropriate spots
      dexn(1)=dexl(1)
      dexn(2)=dexl(2)
      dexn(3)=dexm(3)
      vers=dexm(4)
      if(vers .eq. ZERO) vers=dexl(4)+ONE
      dexn(4)=vers
      dexn(5)=I2FLOAT(ntot)
      return
      end
      subroutine asort2(n,a,b)
      implicit none
      integer  n
      REAL     a(1), b(1)
c----------------------------------------------------------------------
c
c  subroutine  asort2(n,a,b)
c
c  purpose     to sort a and b according to increasing a
c
c  arguments   n         (input)   no. of values to be sorted
c              a         (in/out)  array to be sorted
c              b         (in/out)  array corresponding to a
c
c  remarks     none
c
c  externals   none
c
c  programmer  j.abdallah,jr.
c
c  date        20 feb 1979
c
c----------------------------------------------------------------------
c
c-----local variables
      integer i, i1, j
      REAL    asave, bsave
c-----end local variables
      if(n .eq. 1) return
      i=1
 50   i1=i+1
      do j=i1,n
        if(a(j) .gt. a(i)) go to 100
        asave=a(i)
        a(i)=a(j)
        a(j)=asave
        bsave=b(i)
        b(i)=b(j)
        b(j)=bsave
 100  enddo
      i=i+1
      if(i .lt. n) go to 50
      return
      end
      subroutine upsum(lp,i,rmat,date,ntab,tblids,nwtab)
      implicit none
      integer lp,i
      REAL    rmat, date
      integer ntab
      REAL    tblids(1)
      integer nwtab(1)
c----------------------------------------------------------------------
c
c  subroutine  upsum(lp,i,rmat,date,ntab,tblids,nwtab)
c
c  purpose     to print the update summary
c
c  arguments   lp        (input)   printer unit no.
c              i         (input)   relative mat no. in update
c                                  file
c
c
c              rmat      (input)   material id
c              ntab      (input)   no. of tables
c              tblids    (input)   array of table ids
c              nwtab     (input)   no. of words in each table
c
c  remarks     none
c
c  externals   none
c
c  programmer  j.abdallah
c
c  date        27 feb 1979
c
c----------------------------------------------------------------------
c-----local variables
      integer matid, idate, j, itbl, nwds
c-----end local variables
      if(i .gt. 1) go to 10
      write(lp,1000)
 1000 format(1h1,18x,24hsesame ii update summary,//,
     $3x,11hmaterial id,5x,8hcreation,7x,8htable id,
     $5x,23hno.of words(0=deletion))
 10   matid=rmat
      idate=date
      do j=1,ntab
        itbl=tblids(j)
        nwds=nwtab(j)
        if(j .gt. 1) go to 30
        write(lp,1010) i,matid,idate,itbl,nwds
        go to 40
 30     write(lp,1020) itbl,nwds
 40   enddo
 1010 format(1x,i2,2x,i6,9x,i6,7x,i6,9x,i6)
 1020 format(33x,i6,9x,i6)
      return
      end
#endif /*defined(MAKE_SES_LIB)*/
#endif /*defined(SESAME_CODE)*/
