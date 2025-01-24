!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine scaleset_ptj1(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  ptj1, which is the pt of the leading jet in the event
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'Rcut.f'
      integer isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt

      if((kcase==ktwojet) .or.
     &   (kcase==ktwo_ew)) then

c--- first work out whether this point is real radiation or not
        if (abs(p(npart+2,4))  >  1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif

c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets
        call genclust2(p,Rcut,pjet,isub)

c        write(6,*) 'partons:'
c        if (p(5,4)  >=  1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4)  >=  1.e-8_dp) write(6,*) 'pt6',pt(6,p)
c      write(6,*) 'jets:'
c        if (pjet(5,4)  >=  1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4)  >=  1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)

c--- restore old value of jets
        jets=oldjets

c--- asssign scale
        mu0=max(pt(3,pjet),pt(4,pjet),pt(5,pjet))
c        write(6,*) mu0,pt(3,pjet),pt(4,pjet),pt(5,pjet)

      else
        write(6,*) 'dynamicscale ptj1'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

