!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine lumxmsq_wh(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          use LHAPDF
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'hdecaymode.f'
      include 'beamtype.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'taucut.f'
      integer:: j,k,m,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,qqbWH,qbqWH,hdecay,prop,sH,
     & msqhtautau,msqhbb,msqhgamgam,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),getxmsq_wh
      real(dp):: qqb_wh_HtopEFT,qqbWHtoploop,qbqWHtoploop

      logical, intent(in) :: central
      real(dp) :: origtaucut

      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---calculate the two W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=xn*gwsq**3*wmass**2/prop
c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+2._dp*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+2._dp*mb**2
          hdecay=msqhbb(sH)
c--- adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
            fac=fac*GamHbb/GamHbb0
          endif
      elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
      elseif (hdecaymode == 'wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)
      else
          write(6,*) 'Unimplemented decay mode in lumxmsq_wh.f'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)

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
      qqbWHtoploop=aveqq*qqb_wh_HtopEFT(1,2,3,4,p)*hdecay
      qbqWHtoploop=aveqq*qqb_wh_HtopEFT(2,1,3,4,p)*hdecay
      if ((FixBrHbb).and.(hdecaymode=='bqba')) then
         qqbWHtoploop=qqbWHtoploop*GamHbb/GamHbb0
         qbqWHtoploop=qbqWHtoploop*GamHbb/GamHbb0
      endif
      endif

      if (toponly) then
        qqbWH=zip
        qbqWH=zip
      endif

      xmsq=getxmsq_wh(p,xx,order,soft1,soft2,hard,
     &               beama0,beamb0,beama1,beamb1,beama2,beamb2,
     &               qqbWH,qbqWH,qqbWHtoploop,qbqWHtoploop)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = getxmsq_wh(p,xx,order,soft1,soft2,hard,
     &               beama0,beamb0,beama1,beamb1,beama2,beamb2,
     &               qqbWH,qbqWH,qqbWHtoploop,qbqWHtoploop)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif

      return
      end


      function getxmsq_wh(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,
     &              qqbWH,qbqWH,qqbWHtoploop,qbqWHtoploop)
          use SCET
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ckm.f'
      include 'taucut.f'
      real(dp):: getxmsq_wh
      integer:: j,k,order
      real(dp):: p(mxpart,4),qqbWH,qbqWH,qqbWHtoploop,qbqWHtoploop,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & bit,assemble,tauc,msqlo,msqtop,
     & Qh,dot,powc(-5:5,-5:5)

      tauc=getdynamictau(p,taucut)

c compute power corrections if required
      if ((incpowcorr) .or. (onlypowcorr)) then
        Qh=sqrt(two*dot(p,1,2))
        call powcorr_qa(order,tauc,xx(1),xx(2),Qh,beama0,beamb0,powc)
      endif

      getxmsq_wh=zip
      do j=-nf,nf
      do k=-nf,nf
      if (j*k >= 0) cycle ! skip gluons, qq, aa

      if (onlypowcorr) then
        bit=zip
      else
        bit=assemble(order,tauc,
     &   beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &   beama2(j,:),beamb2(k,:),soft1,soft2,hard)
      endif

      msqtop=0d0

      if ((j > 0) .and. (k < 0)) then
        msqlo=Vsq(j,k)*qqbWH
        if(order >= 2) then
           msqtop=Vsq(j,k)*beama0(j)*beamb0(k)*qqbWHtoploop
        endif
      elseif ((j < 0) .and. (k > 0)) then
        msqlo=Vsq(j,k)*qbqWH
        if(order >= 2) then
           msqtop=Vsq(j,k)*beama0(j)*beamb0(k)*qbqWHtoploop
        endif
      else
        msqlo=zip
      endif

c---- power corrections
      if ((incpowcorr) .or. (onlypowcorr)) then
        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
      endif
      if (onlypowcorr) msqtop=zip

      getxmsq_wh=getxmsq_wh+bit*msqlo+msqtop

      enddo
      enddo

      return
      end


