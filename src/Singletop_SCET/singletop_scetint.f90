!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_scetint
      use singletop_int
      implicit none

      private

      public :: scetint

      contains

      function scetint(r,wgt)
          use ieee_arithmetic
          use Superhisto, only : shtmpreset, shtmpcommit
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          implicit none
          include 'nf.f'
          include 'mxpart.f'
          include 'npart.f'
          include 'x1x2.f'
          include 'kprocess.f'
          include 'nproc.f'
          include 'ipsgen.f'
          include 'constants.f'
          include 'energy.f'
          include 'debug.f'
          include 'xmin.f'
          include 'maxwt.f'
          include 'taucut.f'! usescet
          include 'vegas_common.f'! ndim, mxdim
          include 'beamtype.f'

          logical :: bin
          common/bin/bin

          real(dp) :: scetint 
          real(dp), intent(in) :: r(mxdim), wgt

          ! external functions
          logical :: includedipole

          real(dp) :: p(mxpart,4), pjet(mxpart,4)
          real(dp) :: pswt, z1,z2, QB(2)

          real(dp) :: scetint_b, tmp_b

          logical :: bcontribs

          integer :: truescalevar

          integer :: j,ipdf

          p(:,:) = 0._dp
          pjet(:,:) = 0._dp
          currentPDF = 0
          b_contrib = 1

          if (.not. gen_singletop(r,p,pswt)) then
              scetint = 0._dp
              return
          endif

          if (any(ieee_is_nan(p(1:npart,:)))) then
              write(6,*) 'Discarding NaN or infinite phase space point'
              scetint = 0._dp
              return
          endif

          if (.not. passed_smallnew(p,npart,1d-9)) then
              scetint = 0._dp
              return
          endif

          ! Vegas is really good at finding those pesky zeroes
          ! hardfun_nnlo.f90.inc: w(28) = 1/(1 - 2*omx + omx**2)
          ! for some value x ~ 10^-9 the denominator is exactly zero in double precision
          ! the denominator is actually just x^2

!         if ( abs(massvec(-p(1,:)-p(6,:)) / mt**2) < 1d-8 .or. &
!              abs(massvec(-p(2,:)-p(6,:)) / mt**2) < 1d-8 ) then
!              scetint = 0._dp
!              return
!          endif

          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts
          if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp) &
              .or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
              scetint = 0._dp
              return
          endif

          ! only main b contribution here
          bcontribs = includedipole(0,p)

          if (bcontribs .eqv. .false.) then
              scetint = 0._dp
              return
          endif

          z1 = r(ndim-1)**2
          z2 = r(ndim)**2

          QB(1) = -2*p(1,4)
          QB(2) = -2*p(2,4)

!!! central value

          call singletop2_scale_setup(p)
          scetint_b = scetint_fillxmsq(p,xx,z1,z2,QB,pswt,central=.true.)

           if (ieee_is_nan(scetint_b) .or. (.not. ieee_is_finite(scetint_b))) then
               !write(6,*) 'Discarding NaN matrix element!'
               scetint = 0._dp
               return
           endif

          if (scetint_b == 0._dp) then
              scetint = 0._dp
              return
          endif

          if (doScalevar .and. bin) then
              if (use_DDIS) then
                  truescalevar = (maxscalevar-1)/2
              else
                  truescalevar = maxscalevar
              endif

              do j=1,maxscalevar
                  if (j <= truescalevar) then
                      call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j), mult_in_fac=facscalevarmult(j))
                  elseif (j == (truescalevar + 1)) then
                      call singletop2_scale_setup(p, mult_in_ren=1d0, mult_in_fac=1d0, forcemt=.true.)
                  else
                      call singletop2_scale_setup(p, mult_in_ren=scalevarmult(j-truescalevar-1), &
                          mult_in_fac=facscalevarmult(j-truescalevar-1), forcemt=.true.)
                  endif
                  ! central = .false. makes sure scetreweight is not updated
                  tmp_b = scetint_fillxmsq(p,xx,z1,z2,QB,pswt,central=.false.)

                  scalereweight(j) = tmp_b / scetint_b
              enddo

              ! restore environment
              call singletop2_scale_setup(p)
          endif

          if (maxPDFsets > 0 .and. bin) then
              do ipdf=1,maxPDFsets
                  currentPDF = ipdf
                  if (doPDFAlphas) then
                      call singletop2_scale_setup(p)
                  endif
                  tmp_b = scetint_fillxmsq(p,xx,z1,z2,QB,pswt,central=.false.)

                  pdfreweight(ipdf) = (scetint_b - tmp_b)*wgt
              enddo

              ! restore environment, up to call to singletop2_scale_setup
              currentPDF = 0
          endif

          scetint = scetint_b

          if (bin) then
              currentNd = 0
              call getptildejet(currentNd,pjet)
              includeTaucutgrid(currentNd) = .true.

              if (doScalevar) then
                  continue ! already set
              endif
              if (maxPDFsets > 0) then
                  continue ! already set
              endif

              call nplotter_new(pjet,scetint*wgt)
              call threadStorageOp(shtmpcommit)
          endif

          if (abs(scetint*wgt) > wtmax) then
              wtmax = abs(scetint*wgt)
          endif

      end function

end module
