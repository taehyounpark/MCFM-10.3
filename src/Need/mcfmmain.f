!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine mcfmmain(r,er)
          use Integration
          use parseinput
      implicit none
c***********************************************************************
c                                                                      *
c  This is the main program for MCFM                                   *
c                                                                      *
c  The sequence of calls should always be:                             *
c   call mcfm_init          : basic variable initialization, print-out *
c   call mcfm_vegas(warmup) : warm-up the Vegas grid                   *
c   call mcfm_vegas(accum)  : accumulate results                       *
c   call mcfm_exit          : final processing and print-out           *
c                                                                      *
c***********************************************************************

      real(dp):: integ,integ_err,r,er
      logical:: dryrun
      common/dryrun/dryrun

      integer*4 old_cw

      call f_fpu_fix_start(old_cw)

c basic variable initialization, print-out
      call mcfm_init()

      call mcfm_vegas_adaptive(integ,integ_err)

c final processing and print-out
      r=integ
      er=integ_err
      call mcfm_exit(integ,integ_err)
      end

