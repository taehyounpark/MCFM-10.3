!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scaleset_Msqptgamm1sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt5^2), where M is the mass of the particle (34)
      include 'mxpart.f'
      include 'kprocess.f'
      include 'breit.f'
      include 'kpart.f'
      include 'npart.f'
      logical:: dummy,passedcuts_wgamma_ew,passedcuts_wgammajet_ew
      integer:: isub
      real(dp):: p(mxpart,4),mu0,pt
      real(dp):: ptgamm1,ptlep,mtranselgam
      common/ew_observables/ptgamm1,ptlep,mtranselgam
!$omp threadprivate(/ew_observables/)

      if((kcase==kWgamma) .or.
     &   (kcase==kWgajew) .or.
     &   (kcase==kWga_ew)) then

        if ((kpart == klord) .or. (kpart == kvirt)) then

          mu0=mass3**2+pt(5,p)**2
          mu0=sqrt(abs(mu0))

        else

c--- first work out whether this point is real radiation or not
          if (abs(p(npart+2,4))  >  1.e-8_dp) then
            isub=0  ! real term
          else
            isub=1  ! subtraction term
          endif

          if (kcase == kWgajew) then
            dummy=passedcuts_wgammajet_ew(isub,p)
          else
            dummy=passedcuts_wgamma_ew(isub,p)
          endif

          mu0=mass3**2+ptgamm1**2
          mu0=sqrt(abs(mu0))

        endif

      elseif((kcase==kWgaj_a) .or. (kcase==kWgajja)) then

        mu0=mass3**2+pt(5,p)**2
        mu0=sqrt(abs(mu0))

      else
        write(6,*) 'dynamicscale sqrt(M^2+ptgamm1^2)'//
     &             ' not supported for this process.'
        stop
      endif

      return
      end


