!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dk1qqb_ttw_gs(p,msqc)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     May, 2012.                                                       *
c     calculate the subtraction term for radiation in the              *
c     anti-top quark decay for the process                             *
c                                                                      *
c     q(-p1)+qbar(-p2) = nu(p3)+e+(p4)+b(p5)+e(p9)+nu(p10)             *
c                        +bbar(p6)+e-(p7)+nubar(p8)+g(p11)             *
c                                                                      *
c     Top is kept strictly on-shell although all spin correlations     *
c     are retained. B-quark is taken to be either massive or massless. *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'alfacut.f'
      include 'incldip.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     & p(mxpart,4),q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb,pwsq,xr,
     & y,ymax
      integer:: j,k

      ndmax=1

      do j=-nf,nf
      do k=-nf,nf
        msqc(1,j,k)=0._dp
        incldip(1)=.true.
      enddo
      enddo

c--- special dipole for radiation in top decay
      call wtransform_generic(p,3,4,5,11,q,pbDpg,ptDpg,ptDpb)

      pwsq=2._dp*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))

c--- form of subtraction depends on whether b-quark in decay is massless or not
      if (mb < 1.e-6_dp) then
c----- massless case
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)
        z=1._dp-omz
        xr=sqrt(pwsq/mt**2)
        ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
        y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
        if ((z < 1._dp-aff) .and. (y > aff*ymax)) then
          incldip(1)=.false.
          return
        endif
        fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)
      else
c----- massive case
c-----  (no alpha-dependence at present)
        fac=gsq*cf*((mt**2+mb**2-pwsq)/(ptDpg*pbDpg)
     &             -(mt/ptDpg)**2-(mb/pbDpg)**2)
      endif

      call qqb_ttw(q,msq)

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo

      return
      end

