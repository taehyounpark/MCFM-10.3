!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine lumxmsq_h(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          use LHAPDF
      implicit none
      include 'types.f'
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
c----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  b(p3)+b(p4))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'beamtype.f'
      include 'taucut.f'
      integer:: j,k,m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),
     & msq(-nf:nf,-nf:nf),getxmsq_h
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      logical, intent(in) :: central
      real(dp) :: origtaucut

      call gg_h(p,msq)

      call softggbis(order,soft1,soft2)
      call hardgg(s(1,2),musq,hard)

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

      xmsq=getxmsq_h(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,msq(0,0))

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = getxmsq_h(p,xx,order,soft1,soft2,hard,
     &                  beama0,beamb0,beama1,beamb1,beama2,beamb2,msq(0,0))
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end

      function getxmsq_h(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,msqgg)
          use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'taucut.f'
      integer:: j,order
      real(dp):: getxmsq_h
      real(dp):: p(mxpart,4),
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & msqgg,assemble,Qh,powc(-5:5,-5:5),tauc

      tauc=getdynamictau(p,taucut)

c compute power corrections if required
      if ((incpowcorr) .or. (onlypowcorr)) then
         Qh=sqrt((p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     &          -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2)
         call powcorr_gg(order,tauc,xx(1),xx(2),Qh,beama0,beamb0,powc)
      endif

      if (onlypowcorr) then
        getxmsq_h=zip
      else
        getxmsq_h=msqgg*assemble(order,tauc,
     & beama0(0),beamb0(0),beama1(0,:),beamb1(0,:),
     & beama2(0,:),beamb2(0,:),soft1,soft2,hard)
      endif

c---- power corrections
      if ((incpowcorr) .or. (onlypowcorr)) then
         do j=-nf,nf
            if(j /= 0) then
               getxmsq_h=getxmsq_h+msqgg*(powc(j,0)+powc(0,j))
            elseif(j==0) then
               getxmsq_h=getxmsq_h+msqgg*powc(0,0)
            endif
         enddo
      endif

      return
      end


