c
c				gsestoc.F
c
c	Copyright 1999 by The University at Stony Brook, All rights reserved.
c
#if defined(SESAME_CODE)
#include <gdecs/fortdecs.h>
#define stdin 5
#define stdout 6
#define stderr 0
      subroutine oplib_uf(lu,libnam)
      implicit none
      integer lu
      character libnam*(*)
c-------------------------------------------------------------------------------
c
c	subroutine oplib_uf(libnam,lu)
c	character libnam*(*)
c	integer lu
c
c	opens a file with given unit number lu for UNFORMATED io
c-------------------------------------------------------------------------------
      open(unit=lu,file=libnam,form='unformatted')
      return
      end
      subroutine oplib_f(lu,libnam)
      implicit none
      integer lu
      character libnam*(*)
c-------------------------------------------------------------------------------
c
c	subroutine oplib_f(libnam,lu)
c	character libnam*(*)
c	integer lu
c
c	opens a file with given unit number lu for FORMATED io
c-------------------------------------------------------------------------------
      open(unit=lu,file=libnam)
      return
      end
      subroutine cllib(lu)
      implicit none
      integer lu
c-------------------------------------------------------------------------------
c
c	subroutine cllib(lu)
c	integer lu
c
c	closes the sesame library file with given unit number lu
c-------------------------------------------------------------------------------
      close(unit=lu)
      return
      end
#endif /*defined(SESAME_CODE)*/
