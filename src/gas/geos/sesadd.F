c
c			sesadd.F
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
c
#if defined(SESAME_CODE) && defined(PHASE_CODE)
#include <gdecs/fortdecs.h>
#if defined(float)
#define TWO0ONE 201.d0
#define THREE0ONE 301.d0
#define FOUR0ONE 401.d0
#else /*defined(float)*/
#define TWO0ONE 201.0
#define THREE0ONE 301.0
#define FOUR0ONE 401.0
#endif /*defined(float)*/
c
c     Contains functions similar to SESAME routines, but more specialized.
c
c
      subroutine s4get(ir,ids2,tbls,lcnt,lu,ifl)
      implicit none
      integer ir,ids2
      REAL tbls(1)
      integer lcnt,lu,ifl
c---------------------------------------------------------------
c
c  subroutine:  s4get(ir,ids2,tbls,lcnt,lu,ifl)
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
c               tbls can be declared lcm on the cdc 7600.
c
c  externals:   inbuf, tbfch.
c
c  programmer:  Lisa Osterman NYU
c
c  date:        4 April 1988
c
c---------------------------------------------------------------
c****      lcm tbls
c-----local variables
      integer id, i, loc, itest, nmax, ieof, nwds
c-----end local variables
c-----common block variables
      integer lcmx,nrs,lcfw(10)
      common/s2dir/lcmx,nrs,lcfw
c-----end common block variables
      if(nrs.lt.1) then
        ifl = 0
        return
      endif
      id = iabs(ids2)
      do i=1,nrs
        loc = lcfw(i)
        if(loc.ne.0) then
          itest = tbls(loc)
          if(id.eq.itest) then
            lcfw(ir) = loc
            ifl = 1
            return
          endif
        endif
      enddo
      nmax = lcmx-lcnt+1
      ifl = nmax-5
      if(ifl.lt.0) return
      rewind lu
      ieof = 1
      itest = -1
      do while (id.ne.itest)
        do while (ieof.eq.1)
          call inbuf(lu,tbls(lcnt),1,ieof)
        enddo
        call inbuf(lu,tbls(lcnt),5,ieof)
        if(ieof.eq.0) then
          ifl = 0
          return
        endif
        itest = tbls(lcnt)
      enddo
      nwds = tbls(lcnt+4)
      nwds = nwds+nwds+5
      ifl = nmax-nwds
      if(ifl.lt.0) return
      backspace lu
      call inbuf(lu,tbls(lcnt),nwds,ieof)
      if(ieof.eq.0) then
        ifl = 0
        return
      endif
      backspace lu
      loc = lcnt+nwds
      nmax = nmax-nwds
      call tbfch(TWO0ONE,tbls(loc),tbls(lcnt),lu,nmax,ifl)
      if(ifl.lt.1) return
      tbls(lcnt+1) = tbls(loc+2)
      nmax = lcmx-lcnt-1
      call tbfch(THREE0ONE,tbls(lcnt+2),tbls(lcnt),lu,nmax,ifl)
      if(ifl.lt.1) return
      loc = lcnt
      tbls(lcnt) = id
      lcnt = lcnt+2+ifl
      nmax = nmax - ifl
      tbls(lcnt) = -1
      do while (tbls(lcnt).ne.id)
        call inbuf(lu,tbls(lcnt),5,ieof)
        if(ieof.eq.0) then
          ifl = 0
          return
        endif
      enddo
      itest = tbls(lcnt)
      if(id.ne.itest) then
        ifl = 0
        return
      endif
      nwds = tbls(lcnt+4)
      nwds = nwds+nwds+5
      ifl = nmax-nwds
      if(ifl.lt.0) return
      backspace lu
      call inbuf(lu,tbls(lcnt),nwds,ieof)
      if(ieof.eq.0) then
        ifl = 0
        return
      endif
      backspace lu
      loc = lcnt
      call tbfch(FOUR0ONE,tbls(loc),tbls(lcnt),lu,nmax,ifl)
      if(ifl.lt.1) then
c--------Material does not have a 401 table
	  ifl = 1
	  return
      endif
      lcnt = lcnt +ifl
      lcfw(ir) = loc
      ifl = 1
      return
      end
#endif /*defined(SESAME_CODE) && defined(PHASE_CODE)*/
