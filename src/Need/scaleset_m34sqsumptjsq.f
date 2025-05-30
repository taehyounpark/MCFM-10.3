!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_m34sqsumptjsq(p,mu0)
          use scaleset_m
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(m(34)^2+sum(ptj^2)), where m(34) is the invariant mass of (3+4)
c--- and the sum is over the pt-squared of each jet in the event
      include 'mxpart.f'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'Rcut.f'
      real(dp):: dot
      integer:: isub,oldjets
      real(dp):: p(mxpart,4),pjet(mxpart,4),mu0,pt,sumptsq

      if (kcase == kW_only .or. kcase == kZ_only) then
          if (use_resummation_recoil) then
              mu0=2._dp*dot(p,3,4)+resummation_recoil**2
              mu0=sqrt(abs(mu0))
              return
          else
              mu0=2._dp*dot(p,3,4)
              mu0=sqrt(mu0)
              return
          endif
      endif

      if((kcase==kW_1jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==khjetma) .or.
     &   (kcase==kW_2jet) .or.
     &   (kcase==kZ_2jet) .or.
     &   (kcase==kggfus2) .or.
     &   (kcase==kgmgmjt) .or.
     &   (kcase==kgmgmjj)) then

c--- first work out whether this point is real radiation or not
        if (abs(p(npart+2,4)) > 1.e-8_dp) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif

c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets

        if ((kcase==kgmgmjt) .or.(kcase==kgmgmjj)) then
c--- if a photon process, call specific jet clustering routine
          call genclustphotons(p,Rcut,pjet,isub)
        else
c--- otherwise, generic routine
          call genclust2(p,Rcut,pjet,isub)
        endif

c        write(6,*) 'partons:'
c        if (p(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,p)
c      write(6,*) 'jets:'
c        if (pjet(5,4) >= 1.e-8_dp) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) >= 1.e-8_dp) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)

c--- work out sum of pt-squared
        if     (jets == 1) then
          sumptsq=pt(5,pjet)**2
        elseif (jets == 2) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2
        elseif (jets == 3) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2+pt(7,pjet)**2
        else
          sumptsq=0._dp
        endif

c--- restore old value of jets
        jets=oldjets

        mu0=2._dp*dot(p,3,4)+sumptsq
        mu0=sqrt(abs(mu0))

      else
        write(6,*) 'dynamicscale sqrt(m(34)^2+sumptj^2)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end

