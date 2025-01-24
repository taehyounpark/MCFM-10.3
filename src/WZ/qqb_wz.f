!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wz(p,msq)
      implicit none
      include 'types.f'

c----Author R.K.Ellis December 1998
c----modified by JMC to include supplementary diagrams February 1999
c----Matrix element for WZ production
c----in the notation of DKS
c----averaged over initial oolours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->mu^-(p5)+mu^+(p6)+n(p3)+e^+(p4)
c For nwz=-1
c     d(-p1)+ubar(-p2)-->mu^-(p5)+mu^+(p6)+e^-(p3)+nbar(p4)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'zcouple_cms.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'srdiags.f'
      include 'plabel.f'
      include 'nwz.f'
      include 'interference.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),q1
      complex(dp):: AWZM,AWZP,cprop,a6treea,AWZM_save(-nf:nf,-nf:nf),AWZP_save(-nf:nf,-nf:nf)
      complex(dp):: prop34,prop56,prop12
      complex(dp):: Fa123456,Fa126543,Fb123456_z,Fb123456_g
      complex(dp):: Fa123465,Fa125643,Fb123465_z,Fb123465_g
      complex(dp):: Fa213456,Fa216543,Fb213456_z,Fb213456_g
      complex(dp):: Fa213465,Fa215643,Fb213465_z,Fb213465_g
      complex(dp):: Fa346512,Fa342156,Fa652143
      complex(dp):: Fa345612,Fa342165,Fa653421
      complex(dp):: Fa346521,Fa341256,Fa651243
      complex(dp):: Fa345621,Fa341265,Fa653412
      complex(dp):: ZgL(nf),ZgR(nf),A6b_1,A6b_2,A6b_3,A6b_4
      complex(dp):: v2(2),cl1,cl2,cotw,en1,en2
      real(dp):: ave,wwflag,xfac
      complex(dp):: FAC,FACM
      integer:: j,k,ii,nmax,i3(2),i4(2),i5(2),i6(2)
      parameter(ave=0.25_dp/xn)

      q1=rq1
      FAC=-two*(zesq)**2/zxw

      if ((nwz==1) .or. (nwz == -1)) then
      FACM=nwz*FAC
      else
      write(6,*) 'nwz  /=  +1 or -1'
      stop
      endif
      if     (nwz==-1) then
        cl1=cone
        cl2=czip
        en1=zle
        en2=zln
      elseif (nwz==+1) then
        cl1=czip
        cl2=cone
        en1=zln
        en2=zle
      endif

c--- wwflag=1 for most cases, indicating presence of diagram with 2 W's
      wwflag=1._dp
c--- but for Z -> bbbar this diagram contains |V_tb|**2 which we take 0
      if (plabel(5) == 'bq') then
        wwflag=0._dp
      endif

c-- if Z -> neutrinos, we need to switch c1 and c2
      if (plabel(5) == 'nl') then
        cl1=cone-cl1
        cl2=cone-cl2
      endif

      v2(1)=zl1
      v2(2)=zr1
      cotw=sqrt((one-zxw)/zxw)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
      enddo
      enddo

c----Change the momenta to DKS notation
c   We have --- d(-p1)+ubar(-p2)-->nu(p3)+e^+(p4)+mu^-(p5)+mu^+(p6)
c   DKS have--- u( q2)+dbar( q1)-->nu(q3)+e^+(q4)+mu^-(q6)+mu^+(q5)

      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      enddo

      call spinoru(6,qdks,za,zb)

      i3(1)=3
      i4(1)=4
      i5(1)=5
      i6(1)=6

      if (interference) then
         nmax=2
         if (nwz == +1) then
           i3(2)=3
           i4(2)=5
           i5(2)=4
           i6(2)=6
         else
           i3(2)=6
           i4(2)=4
           i5(2)=5
           i6(2)=3
         endif
      else
         nmax=1
      endif

      do ii=1,nmax

c--   s returned from sprod (common block) is 2*dot product
c--   calculate propagators

c      if     (zerowidth  .eqv. .true.) then
      prop12=s(1,2)/cplx2(s(1,2)-wmass**2,wmass*wwidth)
      prop34=s(i3(ii),i4(ii))/cplx2(s(i3(ii),i4(ii))-wmass**2,wmass*wwidth)
      prop56=s(i5(ii),i6(ii))/cplx2(s(i5(ii),i6(ii))-zmass**2,zmass*zwidth)
      cprop=cone
c      elseif (zerowidth .neqv. .true.) then
c      prop12=cplx2(s(1,2)/(s(1,2)-wmass**2),zip)
c      prop34=cplx2(s(3,4)/(s(3,4)-wmass**2),zip)
c      prop56=cplx2(s(5,6)/(s(5,6)-zmass**2),zip)
c      props=(s(1,2)-wmass**2)/cplx2(s(1,2)-wmass**2,wmass*wwidth)
c      propw=(s(3,4)-wmass**2)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
c      propz=(s(5,6)-zmass**2)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
c      cprop=props*propw*propz
c      endif

c--- DEBUG to compare with Madgraph
c      prop12=s(1,2)/cplx2(s(1,2)-wmass**2,wmass*wwidth)
c      prop34=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
c      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
c      cprop=cone
c--- DEBUG to compare with Madgraph

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s(1,2)/(tevscale*1d3)**2)**2
      else
        xfac=1._dp
      endif
      xdelg1_z=xfac*delg1_z
      xdelg1_g=xfac*delg1_g
      xdelk_z=xfac*delk_z
      xdelk_g=xfac*delk_g
      xlambda_z=xfac*lambda_z
      xlambda_g=xfac*lambda_g

c---case dbar-u or ubar-d
      call A6treeb_anom_wz(1,2,i3(ii),i4(ii),i5(ii),i6(ii),za,zb,A6b_1,A6b_2,A6b_3,A6b_4)
      Fb123456_z=A6b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A6b_3*2._dp*(1._dp+xdelg1_z)
     &          +A6b_4*xlambda_z/wmass**2
      Fb123456_g=A6b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A6b_3*2._dp*(1._dp+xdelg1_g)
     &          +A6b_4*xlambda_g/wmass**2
      Fa123456=A6treea(1,2,i3(ii),i4(ii),i5(ii),i6(ii),za,zb)
      Fa126543=A6treea(1,2,i6(ii),i5(ii),i4(ii),i3(ii),za,zb)

      call A6treeb_anom_wz(1,2,i3(ii),i4(ii),i6(ii),i5(ii),za,zb,A6b_1,A6b_2,A6b_3,A6b_4)
      Fb123465_z=A6b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A6b_3*2._dp*(1._dp+xdelg1_z)
     &          +A6b_4*xlambda_z/wmass**2
      Fb123465_g=A6b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A6b_3*2._dp*(1._dp+xdelg1_g)
     &          +A6b_4*xlambda_g/wmass**2
      Fa123465=A6treea(1,2,i3(ii),i4(ii),i6(ii),i5(ii),za,zb)
      Fa125643=A6treea(1,2,i5(ii),i6(ii),i4(ii),i3(ii),za,zb)

c---case u-dbar or d-ubar
      call A6treeb_anom_wz(2,1,i3(ii),i4(ii),i5(ii),i6(ii),za,zb,A6b_1,A6b_2,A6b_3,A6b_4)
      Fb213456_z=A6b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A6b_3*2._dp*(1._dp+xdelg1_z)
     &          +A6b_4*xlambda_z/wmass**2
      Fb213456_g=A6b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A6b_3*2._dp*(1._dp+xdelg1_g)
     &          +A6b_4*xlambda_g/wmass**2
      Fa213456=A6treea(2,1,i3(ii),i4(ii),i5(ii),i6(ii),za,zb)
      Fa216543=A6treea(2,1,i6(ii),i5(ii),i4(ii),i3(ii),za,zb)

      call A6treeb_anom_wz(2,1,i3(ii),i4(ii),i6(ii),i5(ii),za,zb,A6b_1,A6b_2,A6b_3,A6b_4)
      Fb213465_z=A6b_1*(2._dp+xdelg1_z+xdelk_z+xlambda_z*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_z+xdelk_z+xlambda_z)
     &          +A6b_3*2._dp*(1._dp+xdelg1_z)
     &          +A6b_4*xlambda_z/wmass**2
      Fb213465_g=A6b_1*(2._dp+xdelg1_g+xdelk_g+xlambda_g*s(1,2)/wmass**2)
     &          +A6b_2*(2._dp+xdelg1_g+xdelk_g+xlambda_g)
     &          +A6b_3*2._dp*(1._dp+xdelg1_g)
     &          +A6b_4*xlambda_g/wmass**2
      Fa213465=A6treea(2,1,i3(ii),i4(ii),i6(ii),i5(ii),za,zb)
      Fa215643=A6treea(2,1,i5(ii),i6(ii),i4(ii),i3(ii),za,zb)

      if (srdiags) then
c---for supplementary diagrams.
      Fa346512=A6treea(i3(ii),i4(ii),i6(ii),i5(ii),1,2,za,zb)
      Fa342156=A6treea(i3(ii),i4(ii),2,1,i5(ii),i6(ii),za,zb)
      Fa652143=A6treea(i6(ii),i5(ii),2,1,i4(ii),i3(ii),za,zb)
      Fa345612=A6treea(i3(ii),i4(ii),i5(ii),i6(ii),1,2,za,zb)
      Fa342165=A6treea(i3(ii),i4(ii),2,1,i6(ii),i5(ii),za,zb)
      Fa346521=A6treea(i3(ii),i4(ii),i6(ii),i5(ii),2,1,za,zb)
      Fa341256=A6treea(i3(ii),i4(ii),1,2,i5(ii),i6(ii),za,zb)
      Fa651243=A6treea(i6(ii),i5(ii),1,2,i4(ii),i3(ii),za,zb)
      Fa345621=A6treea(i3(ii),i4(ii),i5(ii),i6(ii),2,1,za,zb)
      Fa341265=A6treea(i3(ii),i4(ii),1,2,i6(ii),i5(ii),za,zb)
      Fa653412=A6treea(i6(ii),i5(ii),i3(ii),i4(ii),1,2,za,zb)
      Fa653421=A6treea(i6(ii),i5(ii),i3(ii),i4(ii),2,1,za,zb)
      endif

c---set up left/right handed couplings for both Z and gamma*
c---note that L/R labels the LEPTON coupling v2, NOT the quarks (all L)
      do j=1,nf
        ZgL(j)=zl(j)*v2(1)*prop56+Q(j)*q1
        ZgR(j)=zl(j)*v2(2)*prop56+Q(j)*q1
      enddo

      do j=-nf,nf
      do k=-nf,nf
c--no point in wasting time if it gives zero anyway
      if (Vsq(j,k)  /=  0._dp) then
          if ((j > 0) .and. (k < 0)) then
            AWZM=(FAC*(ZgL(+j)*Fa213456+ZgL(-k)*Fa216543)
     &           +FACM*(v2(1)*cotw*prop56*Fb213456_z
     &                                +q1*Fb213456_g)*prop12)*prop34
            AWZP=(FAC*(ZgR(+j)*Fa213465+ZgR(-k)*Fa215643)
     &           +FACM*(v2(2)*cotw*prop56*Fb213465_z
     &                                +q1*Fb213465_g)*prop12)*prop34
          elseif ((j < 0) .and. (k > 0)) then
            AWZM=(FAC*(ZgL(+k)*Fa123456+ZgL(-j)*Fa126543)
     &           +FACM*(v2(1)*cotw*prop56*Fb123456_z
     &                                +q1*Fb123456_g)*prop12)*prop34
            AWZP=(FAC*(ZgR(+k)*Fa123465+ZgR(-j)*Fa125643)
     &           +FACM*(v2(2)*cotw*prop56*Fb123465_z
     &                                +q1*Fb123465_g)*prop12)*prop34
          endif
          if (srdiags) then
c---we need supplementary diagrams for gauge invariance.
c---now also assume that we have lepton decay products for W
c---so that v2(1)=le, v2(2)=re
c---1st term is diagram where Z couples to electron
c---2nd term is diagram where Z couples to neutrino
c---3rd term is diagram where gamma* couples to electron
c---4th term (l-h only) contains two W propagators
          if ((j > 0) .and. (k < 0)) then
             AWZM=AWZM+FAC*prop12*(
     &          (en1*Fa346512+en2*Fa342156)*v2(1)*prop56
     &          +q1*(-1._dp)*(cl1*Fa346512+cl2*Fa342156)
     &          +wwflag*0.5_dp/zxw*prop34*(cl1*Fa652143+cl2*Fa653412))
            AWZP=AWZP+FAC*prop12*(
     &          (en1*Fa345612+en2*Fa342165)*v2(2)*prop56
     &          +q1*(-1._dp)*(cl1*Fa345612+cl2*Fa342165))
          elseif ((j < 0) .and. (k > 0)) then
            AWZM=AWZM+FAC*prop12*(
     &          (en1*Fa346521+en2*Fa341256)*v2(1)*prop56
     &          +q1*(-1._dp)*(cl1*Fa346521+cl2*Fa341256)
     &          +wwflag*0.5_dp/zxw*prop34*(cl1*Fa651243+cl2*Fa653421))
            AWZP=AWZP+FAC*prop12*(
     &          (en1*Fa345621+en2*Fa341265)*v2(2)*prop56
     &          +q1*(-1._dp)*(cl1*Fa345621+cl2*Fa341265))
          endif
          endif

c-- Inclusion of width for W,Z a la Baur and Zeppenfeld
      AWZM=cprop*AWZM
      AWZP=cprop*AWZP

c for interference calculation, only generate (34) resonance, so instead of
c |A1+A2|^2, accummulate |A1+A2|^2+A1^2-A2^2

      if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=Vsq(j,k)*ave*vsymfact*(abs(AWZM)**2+abs(AWZP)**2)
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          AWZM_save(j,k)=AWZM
          AWZP_save(j,k)=AWZP
        else
          msq(j,k)=Vsq(j,k)*ave*vsymfact*(abs(AWZM_save(j,k)-AWZM)**2
     &                                   +abs(AWZM_save(j,k))**2-abs(AWZM)**2
     &                                   +two*abs(AWZP_save(j,k))**2)
        endif
      endif

      endif

      enddo
      enddo

      enddo   ! end ii loop

      return
      end
