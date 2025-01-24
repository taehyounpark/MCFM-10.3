!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function pt_xmsq(z1,z2,p,central)
      use SCET
      use ptveto, only: useBNR
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'taucut.f'
      real(dp) :: pt_xmsq
      real(dp), intent(in) :: z1, z2
      real(dp), intent(in) :: p(mxpart,4)
      logical, intent(in) :: central
      real(dp) :: QB(2),xx(2),Q,R,ptveto
      include 'first.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'energy.f'
      include 'Rcut.f'
      integer, save :: iorder
!$omp threadprivate(iorder)
      real(dp) :: xmsq
      external hard_DY,hard_H,hard_Vgam,hard_gamgam
#ifdef WITH_VVAMP
      external hard_VV
#endif

!     dirty business
!     use taucut common block to represent ptveto
      if (dynamictau) then
        ptveto=getdynamictau(p)
      else
        ptveto=taucut
      endif
      R=Rcut

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts
      Q=sqrt(xx(1)*xx(2))*sqrts
      QB(1)=-two*p(1,4)
      QB(2)=-two*p(2,4)

c--- determine order of calculation on first call
      if (first) then
        first=.false.
        if (kpart==knnlo) then
          iorder=2
        elseif (kpart==ksnlo) then
          iorder=1
        else
          write(6,*) 'Error in ptint: kpart=',kpart
          stop
        endif
      endif

c--- Calculate the required matrix elements
      if     ((kcase==kW_only) .or. (kcase==kZ_only)) then
        if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_DY,xmsq)
        else
          call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_DY,xmsq)
        endif
!        call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_DY,xmsq)
!        write(6,*) 'W/Z BNR',xmsq
!        call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_DY,xmsq)
!        write(6,*) 'W/Z BBS',xmsq
!        pause
      elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga)
     &       .or. (kcase==kWH__WW) .or. (kcase==kZHbbar)
     &       .or. (kcase==kZHgaga) .or. (kcase==kZH__WW)) then
        if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_DY,xmsq)
        else
          call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_DY,xmsq)
        endif
      elseif (kcase==kggfus0) then
        if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_H,xmsq)
        else
          call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_H,xmsq)
        endif
!        call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_H,xmsq)
!        write(6,*) 'H BNR',xmsq
!        call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_H,xmsq)
!        write(6,*) 'H BBS',xmsq
!        pause
      elseif ((kcase==kWWqqbr)
     &       .or. (kcase==kWZbbar)
     &       .or. (kcase==kZZlept)) then
#ifdef WITH_VVAMP
         if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_VV,xmsq)
         else
           call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_VV,xmsq)
         endif
#else
        write(6,*) "error:Please recompile with VVamp support"
        write(6,*) "error:to run WW,WZ and ZZ production."
        stop 
#endif
      elseif ((kcase==kWgamma) .or. (kcase==kZgamma)) then
         if (kcase==kZgamma) then
           call set_anomcoup(p)
         endif
         if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_Vgam,xmsq)
         else
           call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_Vgam,xmsq)
         endif
      elseif (kcase==kgamgam) then
         if (useBNR) then
          call BNRptlumxmsq(p,xx,z1,z2,ptveto,Q,R,iorder,central,hard_gamgam,xmsq)
         else
           call ptlumxmsq(p,xx,z1,z2,QB,iorder,central,hard_gamgam,xmsq)
         endif
      else
        write(6,*) 'Process not available in pt (or ptveto) formalism'
        stop
      endif

c These processes are not yet implemented in the pt formalism; this code
c is left here as a reminder of how it could be upgraded in the future
c      elseif (kcase==kHi_Zga) then
c        call ptlumxmsq_h_Zga(p,xx,z1,z2,QB,iorder,xmsq)

      pt_xmsq = xmsq

      end function pt_xmsq

