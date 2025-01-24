!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function scet_xmsq(z1,z2,p,central)
          use types
          use, intrinsic :: iso_fortran_env, only : error_unit
          use nnlo_z1jet
          use nnlo_w1jet
          use nnlo_h1jet
          implicit none
          include 'mxpart.f'
          include 'constants.f'
          real(dp) :: scet_xmsq
          real(dp), intent(in) :: z1, z2
          real(dp), intent(in) :: p(mxpart,4)
          real(dp) :: QB(2)
          logical, intent(in) :: central
          include 'first.f'
          include 'taucut.f'
          include 'kprocess.f'
          include 'kpart.f'
          include 'x1x2.f'
          include 'energy.f'
          include 'npart.f'
          integer, save :: iorder
!$omp threadprivate(iorder)

          integer :: j
          real(dp) :: Q(4), pboost(mxpart,4), Qrest(4), Qsq, ppin(4), ppout(4)

          real(dp) :: xmsq

          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts

          QB(1)=-two*p(1,4)
          QB(2)=-two*p(2,4)
          
          if (tauboost) then
            if (ntau == 0) then
! provide beam energies in singlet c.o.m. instead if required
              QB(1)=sqrt(QB(1)*QB(2))
              QB(2)=QB(1)
! boost to Born frame (more general, for ntau > 0)
            else
              Q(:)=-p(1,:)-p(2,:)
              Qsq=Q(4)**2-Q(1)**2-Q(2)**2-Q(3)**2
              Qrest(:)=zip
              Qrest(4)=sqrt(Qsq)
              do j=1,npart+2
                ppin(:)=p(j,:)
                call boostx(ppin,Q,Qrest,ppout)
                pboost(j,:)=ppout(:)
              enddo
              QB(1)=-two*pboost(1,4)
              QB(2)=-two*pboost(2,4)
            endif
          else
            pboost(1:npart+2,:) = p(1:npart+2,:)
          endif

c--- d    etermine order of calculation on first call
          if (first) then
            first=.false.
            if (kpart==knnlo) then
              iorder=2
            elseif (kpart==ksnlo) then
              iorder=1
            else
              write(6,*) 'Error in scetint: kpart=',kpart
              stop
            endif
          endif

c--- Calculate the required matrix elements
          if     (kcase==kW_only) then
            call lumxmsq_w(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kZ_only) then
            call lumxmsq_z(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase == kZ_1jet) then
            call lumxmsq_z1jet(pboost,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase == kggfus1) then
            call lumxmsq_h1jet(pboost,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase == kW_1jet) then
            call lumxmsq_w1jet(pboost,xx,z1,z2,QB,iorder,xmsq,central)
          elseif ((kcase==kggfus0) .or. (kcase==kHigaga)) then
            call lumxmsq_h(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kHi_Zga) then
            call lumxmsq_h_Zga(p,xx,z1,z2,QB,iorder,xmsq)
          elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga)
     &       .or. (kcase==kWH__WW)) then
            call lumxmsq_wh(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga)
     &       .or. (kcase==kZH__WW)) then
            call lumxmsq_zh(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kgamgam) then
             call lumxmsq_gaga(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kZgamma) then
             call set_anomcoup(p)
             call lumxmsq_zgamma(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kWgamma) then
             call lumxmsq_wgamma(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kWWqqbr) then
#ifdef WITH_VVAMP
            call lumxmsq_vv(p,xx,z1,z2,QB,iorder,xmsq,central)
#else
            error stop "Please recompile with VVamp support to run WW production"
#endif
          elseif (kcase==kWZbbar) then
#ifdef WITH_VVAMP
            call lumxmsq_vv(p,xx,z1,z2,QB,iorder,xmsq,central)
#else
            error stop "Please recompile with VVamp support to run WZ production"

#endif
          elseif (kcase==kZZlept) then
#ifdef WITH_VVAMP
            call lumxmsq_vv(p,xx,z1,z2,QB,iorder,xmsq,central)
#else
            error stop "Please recompile with VVamp support to run ZZ production"
#endif
          else
            error stop 'Process not yet available in jettiness formalism'
          endif

          scet_xmsq = xmsq

      end function

