!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_Msqptj1sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+ptj1^2), where M is the mass of the particle (34)
c--- and ptj1 is the pt of the leading jet in the event

      include 'mxpart.f'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'breit.f'
      include 'Rcut.f'
      integer:: isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt,ptj1,ptj2

      if((kcase==kW_1jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==khjetma)) then

c--- first work out whether this point is real radiatio or not
        if (abs(p(npart+2,4)) > 1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif

c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets
        call genclust2(p,Rcut,pjet,isub)

c        write(6,*) 'partons:'
c        if (p(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,p)
c      write(6,*) 'jets:'
c        if (pjet(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)

c--- restore old value of jets
        jets=oldjets

c--- order according to pt
      ptj1=pt(5,pjet)
      ptj2=pt(6,pjet)
      if (ptj2 > ptj1) ptj1=ptj2

c--- assign scale
        mu0=mass3**2+ptj1**2
        mu0=sqrt(abs(mu0))

      else
        write(6,*) 'dynamicscale sqrt(M^2+ptj1^2)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

