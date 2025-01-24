!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_mVpmH(p,mu0)
      implicit none
      include 'types.f'
c==== this is the scale setting routine for
c===  mV+mH, its not really a dynmaic scale, but its
c==== more convienent than carrying around junk inthe input.dAT
      include 'kprocess.f'
      include 'mxpart.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kWHbbar) .or.
     &   (kcase==kWH1jet) .or.
     &   (kcase==kWHgaga) .or.
     &   (kcase==kWH__WW) .or.
     &   (kcase==kZH__WW) .or.
     &   (kcase==kZHgaga) .or.
     &   (kcase==kZHbbar) .or.
     &   (kcase==kZH1jet)) then
         mu0=mass2+mass3
      else
        write(6,*)'dynamicscale mV+mH not supported for this process.'
        stop
      endif

      return
      end

