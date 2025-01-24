!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ptint(r,wgt)
          use ieee_arithmetic
          use SCET
          use types
          use Scalevar
          use PDFerrors
          use MCFMStorage
          use m_gencuts, only : enable_reweight_user, reweight_user
          use SafetyCuts, only : passed_smallnew
          use MCFMSetupPlots, only: nplotter_new
          use MCFMSettings
          use ptveto, only: useBNR, jetptveto
      implicit none
      real(dp):: ptint
      include 'constants.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'energy.f'
      include 'npart.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'taucut.f'
      include 'first.f'
      include 'scalevar.f'
      include 'scale.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'x1x2.f'
      include 'xmin.f'
      include 'debug.f'
      include 'newpoint.f'
      include 'Lnu.f'
      include 'Lw.f'
      include 'mpicommon.f'
      real(dp):: savescale,savefacscale
      real(dp):: p(mxpart,4),pjet(mxpart,4),r(mxdim),W,xmsq,
     & val,val2,pswt,wgt,z1,z2,flux,BrnRat,pt_xmsq,scalein,facscalein
c      real(dp):: tauc,Kfac,Kqq ! for calling ptont bits
      integer j
      logical:: bin,includedipole
      real(dp) :: xjac
      common/bin/bin
      common/BrnRat/BrnRat

      if (first) then
!$omp master
        if (rank == 0) then
          if (useBNR) then
            write(6,*) 'Entering ptint: using BNR pt implementation'
          else
            write(6,*) 'Entering ptint: using BBS pt implementation'
          endif
        endif
!$omp end master
        Lnu=20._dp    ! Lnu = log (mu/nu)
        first=.false.
      endif

      ptint=0._dp

      W=sqrts**2

      currentPDF = 0

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      call gen_lops(r,p,pswt,*999)
      call dotem(npart+2,p,s)

c small safety cuts
      if (.not. passed_smallnew(p,npart)) then
        goto 999
      endif

c---  see whether this point will pass cuts - if it will not, do not
c---  bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

      z1=r(ndim-1)**2
      z2=r(ndim)**2
      xjac = four*sqrt(z1*z2)

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

      if ( (xx(1) > one)  .or. (xx(2) > one)
     & .or.(xx(1) < xmin) .or. (xx(2) < xmin)) goto 999

      newpoint = .true.

      ! CENTRAL VALUE
      if (dynamicscale) call scaleset(initscale,initfacscale,p)
c      write(6,*) 'calling us'
      xmsq = pt_xmsq(z1,z2,p,.true.)

c      write(6,*) 'calling Ye Li'
c      call CalcParam(5)
c      call LoadAllPDF(xx(1),xx(2),z1,z2,facscale)
c      if (dynamictau) then
c        tauc=getdynamictau(p)
c      else
c        tauc=taucut
c      endif
c      Kfac=Kqq(-2,1,-2,xx(1),xx(2),z1,z2,as,scale,facscale,p(1:4,:),tauc**2,1)
c      write(6,*) 'Kfac',Kfac
c      pause

      flux = fbGeV2/(two*xx(1)*xx(2)*W)
      ptint = flux*xjac*pswt*xmsq/BrnRat

      call getptildejet(0,pjet)
      val=ptint*wgt
      val2=val**2

      if (ieee_is_nan(val) .or. (.not. ieee_is_finite(val))) then
          if (debug) then
              write(6,*) 'Discarded NaN, val=',val
          endif
          goto 999
      endif

      ! SCALE VARIATION WITH CENTRAL PDF
      if (doScalevar .and. bin .and. xmsq /= 0._dp) then
          if (dynamicscale) then
              call scaleset(initscale,initfacscale,p)
          endif

          savescale = scale
          savefacscale = facscale

          do j=1,maxscalevar
              if (vetoscalevar) then
                if (scalevarmult(j) > 1._dp) then
                  scalein=max(savescale,jetptveto)*scalevarmult(j)
                else
                  scalein=min(savescale,jetptveto)*scalevarmult(j)
                endif
                if (facscalevarmult(j) > 1._dp) then
                  facscalein=max(savefacscale,jetptveto)*facscalevarmult(j)
                else
                  facscalein=min(savefacscale,jetptveto)*facscalevarmult(j)
                endif
              else
                scalein=savescale*scalevarmult(j)
                facscalein=savefacscale*facscalevarmult(j)
              endif
              call usescales(scalein,facscalein)
              scalereweight(j) = pt_xmsq(z1,z2,p,.false.)/xmsq

          enddo

          ! restore
          call usescales(savescale, savefacscale)
      endif

      ! PDF VARIATION WITH CENTRAL SCALE

      if (doPDFerrors .and. bin) then
          do j=1,maxPDFsets
              currentPDF = j
              if (doPDFAlphas) then
                  if (dynamicscale) then
                      call scaleset(initscale,initfacscale,p)
                  else
                      call usescales(initscale,savefacscale)
                  endif
                  call updateAlphas(scale)
              endif
              pdfreweight(currentPDF) = (ptint - flux*xjac*pswt*pt_xmsq(z1,z2,p,.false.)/BrnRat)*wgt
          enddo
      endif

      if (bin) then
          includeTaucutgrid(0) = .true.
          if (newStyleHistograms) then
              call nplotter_new(pjet,val)
          else
            call nplotter(pjet,val,val2,0)
          endif
      endif

      if (enable_reweight_user) then
          ptint = ptint * reweight_user(pjet)
      endif

      return

 999  continue

      ptint = 0._dp
      end


