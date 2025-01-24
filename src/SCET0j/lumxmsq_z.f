!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine lumxmsq_z(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          use LHAPDF
      implicit none
      include 'types.f'
c----Matrix element for W production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
c For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'taucut.f'
      include 'zcouple_cms.f'
      include 'beamtype.f'
      integer:: m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2), getxmsq_z
      complex(dp):: prop,qqb,qbq
      include 'cplx.h'

      logical, intent(in) :: central
      real(dp) :: origtaucut

      fac=4._dp*abs(zesq)**2*xn

      call spinoru(4,p,za,zb)

c--   calculate propagator
      fac=aveqq*fac/s(3,4)**2
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
c---case dbar-u or ubar-d
      qqb=za(2,3)*zb(4,1)
      qbq=za(1,3)*zb(4,2)

      call softqqbis(order,soft1,soft2)
      call hardqq(s(1,2),musq,hard)

      if (order >= 0) then
      call fdist(ih1,xx(1),facscale,beama0,1)
      call fdist(ih2,xx(2),facscale,beamb0,2)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,z1,xx(1),QB(1),beama1,1)
      call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1,2)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,z1,xx(1),QB(1),beama2,1)
      call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2,2)
      endif

      xmsq = getxmsq_z(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = getxmsq_z(p,xx,order,soft1,soft2,hard,
     &                beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif


      return
      end


      function getxmsq_z(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)
          use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'taucut.f'
      include 'zcouple_cms.f'
      real(dp):: getxmsq_z
      integer:: j,k,order
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & bit,assemble,msqlo, tauc
      real(dp):: Qh, dot, powc(-5:5,-5:5)
      complex(dp):: prop,qqb,qbq

      tauc=getdynamictau(p,taucut)

c compute power corrections if required
      if ((incpowcorr) .or. (onlypowcorr)) then
        Qh=sqrt(two*dot(p,1,2))
        call powcorr_qa(order,tauc,xx(1),xx(2),Qh,beama0,beamb0,powc)
      endif

      getxmsq_z=zip
      do j=-nf,nf
      k=-j
      if (j == 0) cycle

      if (onlypowcorr) then
        bit=zip
      else
        bit=assemble(order,tauc,
     &   beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &   beama2(j,:),beamb2(k,:),soft1,soft2,hard)
      endif

      if (j > 0) then
        msqlo=fac*(
     &     abs((Q(j)*q1+zL(j)*zl1*prop)*qqb)**2
     &    +abs((Q(j)*q1+zR(j)*zr1*prop)*qqb)**2
     &    +abs((Q(j)*q1+zL(j)*zr1*prop)*qbq)**2
     &    +abs((Q(j)*q1+zR(j)*zl1*prop)*qbq)**2)
      elseif (j < 0) then
        msqlo=fac*(
     &     abs((Q(k)*q1+zL(k)*zl1*prop)*qbq)**2
     &    +abs((Q(k)*q1+zR(k)*zr1*prop)*qbq)**2
     &    +abs((Q(k)*q1+zL(k)*zr1*prop)*qqb)**2
     &    +abs((Q(k)*q1+zR(k)*zl1*prop)*qqb)**2)
      else
        msqlo=zip
      endif

      if ((incpowcorr) .or. (onlypowcorr)) then
        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
      endif

      getxmsq_z=getxmsq_z+bit*msqlo

      enddo

      return
      end

