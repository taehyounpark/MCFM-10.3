!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine lumxmsq_gaga(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          use LHAPDF
      implicit none
      include 'types.f'
c====== C. Williams July 2015
c----Matrix element for Ga Ga production
c----averaged over initial colours and spins
c===== based upon similar routines for W and Z
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'facscale.f'
      include 'ewcharge.f'
      include 'taucut.f'
      include 'zcouple_cms.f'
      include 'beamtype.f'
      integer:: j,k,m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,qqb,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),getxmsq_gamgam
      real(dp)::statfac
      parameter(statfac=0.5_dp)
      real(dp):: qsum,tlrem
      real(dp):: gggaga,msqgggaga,facgg

      logical, intent(in) :: central
      real(dp) :: origtaucut

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      fac=xn*abs(zesq)**2*statfac
      gggaga=zip
      qsum=zip
      do j=1,nf
         qsum=qsum+Q(j)**2
      enddo
c======= compute Matrix elements
c      write(6,*) 'old routine'
c      call gamgamampsq(order,p,1,2,3,4,qqb,hard,tlrem)
c      write(6,*) 'hard(1),hard(2)',hard(1),hard(2)
c      write(6,*) 'new routine'
      call gamgamampsq_new(order,p,1,2,3,4,qqb,hard,tlrem)
c      write(6,*) 'hard(1),hard(2)',hard(1),hard(2)
c      pause

      qqb=fac*aveqq*qqb


      call softqqbis(order,soft1,soft2)

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
      facgg=4._dp*abs(zesq)*gsq/(16._dp*pisq)*Qsum
      gggaga=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
c      gggaga=zip
      endif

      xmsq=getxmsq_gamgam(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qsum,qqb,tlrem,gggaga)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = getxmsq_gamgam(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qsum,qqb,tlrem,gggaga)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


      function getxmsq_gamgam(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qsum,qqb,tlrem,gggaga)
          use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'taucut.f'
      real(dp):: getxmsq_gamgam
      integer:: j,k,order
      real(dp):: p(mxpart,4),fac,qqb,qbq,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),bit,assemble
      real(dp):: Qh, dot, powc(-5:5,-5:5)
      real(dp):: qsum,tlrem,bitrem
      real(dp):: gggaga,tauc

      qbq=qqb

      tauc=getdynamictau(p,taucut)

c compute power corrections if required
c note that, for this process, the power corrections computed here
c are not expected to provide any improvement due to t-channel diagrams
      if ((incpowcorr) .or. (onlypowcorr)) then
        Qh=sqrt(two*dot(p,1,2))
        call powcorr_qa(order,tauc,xx(1),xx(2),Qh,beama0,beamb0,powc)
      endif

      getxmsq_gamgam=zip
      do j=-nf,nf
         k=-j

         if (j*k > 0) cycle    ! skip, qq, aa

         if (onlypowcorr) then
           bit=zip
           bitrem=zip
         else
           bit=assemble(order,tauc,
     &          beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &          beama2(j,:),beamb2(k,:),soft1,soft2,hard)
           bitrem=fac*aveqq*tlrem*beama0(j)*beamb0(k)*ason2pi**2
         endif

         if ((incpowcorr) .or. (onlypowcorr)) then
           bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
         endif

         if ((j > 0) .and. (k < 0)) then
            bit=bit*qqb*Q(j)**4
            if(order > 1) then
c======add on two-loop A functions which go like sum over quark charges
               bit=bit+bitrem*qsum*Q(j)**2
            endif
         elseif ((j < 0) .and. (k > 0)) then
            bit=bit*qbq*Q(k)**4
            if(order > 1) then
c======add on two-loop A functions which go like sum over quark charges
               bit=bit+bitrem*qsum*Q(k)**2
            endif
         elseif ((j==0).and.(k==0).and.(order >= 2 )) then
            if (onlypowcorr) then
              bit=zip
            else
              bit=beama0(j)*beamb0(k)*gggaga
            endif
         else
            bit=zip
         endif

         getxmsq_gamgam=getxmsq_gamgam+bit

      enddo
      return
      end
