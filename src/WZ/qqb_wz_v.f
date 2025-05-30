!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wz_v(p,msqv)
      implicit none
      include 'types.f'
c----Author R.K.Ellis December 1998
c----modified by JMC to include supplementary diagrams February 1999
c---  calculate the virtual matrix element squared and endpoint
c     subtraction terms for the process
c For nwz=+1
c     u(-p1)+dbar(-p2)-->mu^-(p5)+mu^+(p6)+n(p3)+e^+(p4)
c For nwz=-1
c     d(-p1)+ubar(-p2)-->mu^-(p5)+mu^+(p6)+e^-(p3)+nbar(p4)
c---  averaged(summed) over initial(final) colours and spins

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'zerowidth.f'
      include 'ewcharge.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'srdiags.f'
      include 'nwz.f'
      include 'plabel.f'
      include 'zcouple_cms.f'
      include 'blha.f'
      include 'ewinput.f'
      include 'interference.f'
      integer:: j,k,ii,nmax,i3(2),i4(2),i5(2),i6(2)
      real(dp):: p(mxpart,4),qdks(mxpart,4)
      real(dp):: msqv(-nf:nf,-nf:nf)
      real(dp):: facnlo,ave,wwflag
      real(dp):: FAC,FACM
      complex(dp):: AWZM,AWZP,BWZM,BWZP,Vpole,Vpole12,suppl
      complex(dp):: AWZM_save(-nf:nf,-nf:nf),AWZP_save(-nf:nf,-nf:nf),
     & BWZM_save(-nf:nf,-nf:nf),BWZP_save(-nf:nf,-nf:nf)
      complex(dp):: prop12,prop34,prop56
      complex(dp):: props,propw,propz,cprop,A6b_1,A6b_2,A6b_3,A6b_4
      complex(dp):: a6treea,a6loopa
      complex(dp):: Fa123456,Fa126543,Fb123456_z,Fb123456_g
      complex(dp):: Fa123465,Fa125643,Fb123465_z,Fb123465_g
      complex(dp):: Fa213456,Fa216543,Fb213456_z,Fb213456_g
      complex(dp):: Fa213465,Fa215643,Fb213465_z,Fb213465_g
      complex(dp):: La123456,La126543,Lb123456_z,Lb123456_g
      complex(dp):: La123465,La125643,Lb123465_z,Lb123465_g
      complex(dp):: La213456,La216543,Lb213456_z,Lb213456_g
      complex(dp):: La213465,La215643,Lb213465_z,Lb213465_g
      complex(dp):: Fa346512,Fa342156,Fa652143
      complex(dp):: Fa345612,Fa342165,Fa653421
      complex(dp):: Fa346521,Fa341256,Fa651243
      complex(dp):: Fa345621,Fa341265,Fa653412
c      complex(dp):: Fa561243,Fa562143
      complex(dp):: v2(2),cl1,cl2,en1,en2,cotw
      real(dp):: xfac
      complex(dp):: ZgL(1:nf),ZgR(1:nf)
      parameter(ave=0.25_dp/xn)

      scheme='dred'
      FAC=-2._dp*gwsq*esq
      FACNLO=ason2pi*cf
      if ((nwz==1) .or. (nwz == -1)) then
      FACM=nwz*FAC
      else
      write(6,*) 'nwz ne +1 or -1'
      stop
      endif
      if     (nwz==-1) then
        cl1=1._dp
        cl2=0._dp
        en1=zle
        en2=zln
      elseif (nwz==+1) then
        cl1=0._dp
        cl2=1._dp
        en1=zln
        en2=zle
      endif
      v2(1)=zl1
      v2(2)=zr1
      cotw=sqrt((one-zxw)/zxw)

c--- wwflag=1 for most cases, indicating presence of diagram with 2 W's
      wwflag=1._dp
c--- but for Z -> bbbar this diagram contains |V_tb|**2 which we take 0
      if (plabel(5) == 'bq') then
        wwflag=0._dp
      endif

c-- if Z -> neutrinos, we need to switch c1 and c2
      if (plabel(5) == 'nl') then
        cl1=1._dp-cl1
        cl2=1._dp-cl2
        wwflag=sqrt(3._dp)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

c---Change the momenta to DKS notation
c   We have --- d(-p1)+ubar(-p2)-->nu(p6)+e^+(p7)+mu^-(p4)+mu^+(p5)
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
      if     ((zerowidth  .eqv. .true.) .or. (ewscheme > 3)) then
      prop12=s(1,2)/cplx2(s(1,2)-wmass**2,wmass*wwidth)
      prop34=s(i3(ii),i4(ii))/cplx2(s(i3(ii),i4(ii))-wmass**2,wmass*wwidth)
      prop56=s(i5(ii),i6(ii))/cplx2(s(i5(ii),i6(ii))-zmass**2,zmass*zwidth)
      cprop=cone
      elseif (zerowidth .neqv. .true.) then
      prop12=cplx1(s(1,2)/(s(1,2)-wmass**2))
      prop34=cplx1(s(i3(ii),i4(ii))/(s(i3(ii),i4(ii))-wmass**2))
      prop56=cplx1(s(i5(ii),i6(ii))/(s(i5(ii),i6(ii))-zmass**2))
      props=(s(1,2)-wmass**2)/cplx2(s(1,2)-wmass**2,wmass*wwidth)
      propw=(s(i3(ii),i4(ii))-wmass**2)/cplx2(s(i3(ii),i4(ii))-wmass**2,wmass*wwidth)
      propz=(s(i5(ii),i6(ii))-zmass**2)/cplx2(s(i5(ii),i6(ii))-zmass**2,zmass*zwidth)
      cprop=props*propw*propz
      endif

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

c--- loop for a6b is simply tree*Vpole
      Vpole12=Vpole(s(1,2))

      if (useblha == 0) then
c---case dbar-u
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

      Lb123456_z=Vpole12*Fb123456_z
      Lb123456_g=Vpole12*Fb123456_g
      La123456=A6loopa(1,2,i3(ii),i4(ii),i5(ii),i6(ii),za,zb)
      La126543=A6loopa(1,2,i6(ii),i5(ii),i4(ii),i3(ii),za,zb)

      Lb123465_z=Vpole12*Fb123465_z
      Lb123465_g=Vpole12*Fb123465_g
      La123465=A6loopa(1,2,i3(ii),i4(ii),i6(ii),i5(ii),za,zb)
      La125643=A6loopa(1,2,i5(ii),i6(ii),i4(ii),i3(ii),za,zb)
      endif

c---case u-dbar
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

      Lb213456_z=Vpole12*Fb213456_z
      Lb213456_g=Vpole12*Fb213456_g
      La213456=A6loopa(2,1,i3(ii),i4(ii),i5(ii),i6(ii),za,zb)
      La216543=A6loopa(2,1,i6(ii),i5(ii),i4(ii),i3(ii),za,zb)

      Lb213465_z=Vpole12*Fb213465_z
      Lb213465_g=Vpole12*Fb213465_g
      La213465=A6loopa(2,1,i3(ii),i4(ii),i6(ii),i5(ii),za,zb)
      La215643=A6loopa(2,1,i5(ii),i6(ii),i4(ii),i3(ii),za,zb)

      if (srdiags) then
c---for supplementary diagrams.
      Fa346512=A6treea(i3(ii),i4(ii),i6(ii),i5(ii),1,2,za,zb)
      Fa342156=A6treea(i3(ii),i4(ii),2,1,i5(ii),i6(ii),za,zb)
      Fa652143=A6treea(i6(ii),i5(ii),2,1,i4(ii),i3(ii),za,zb)
      Fa345612=A6treea(i3(ii),i4(ii),i5(ii),i6(ii),1,2,za,zb)
      Fa342165=A6treea(i3(ii),i4(ii),2,1,i6(ii),i5(ii),za,zb)
c      Fa562143=A6treea(i5(ii),i6(ii),2,1,i4(ii),i3(ii),za,zb)
      Fa346521=A6treea(i3(ii),i4(ii),i6(ii),i5(ii),2,1,za,zb)
      Fa341256=A6treea(i3(ii),i4(ii),1,2,i5(ii),i6(ii),za,zb)
      Fa651243=A6treea(i6(ii),i5(ii),1,2,i4(ii),i3(ii),za,zb)
      Fa345621=A6treea(i3(ii),i4(ii),i5(ii),i6(ii),2,1,za,zb)
      Fa341265=A6treea(i3(ii),i4(ii),1,2,i6(ii),i5(ii),za,zb)
c      Fa561243=A6treea(i5(ii),i6(ii),1,2,i4(ii),i3(ii),za,zb)
      Fa653412=A6treea(i6(ii),i5(ii),i3(ii),i4(ii),1,2,za,zb)
      Fa653421=A6treea(i6(ii),i5(ii),i3(ii),i4(ii),2,1,za,zb)
c---loop diagrams just tree*Vpole since they're all triangle-type
      Vpole12=Vpole(s(1,2))
      endif

c---set up left/right handed couplings for both Z and gamma*
c---note that L/R labels the LEPTON coupling v2, NOT the quarks (all L)
      do j=1,nf
        ZgL(j)=zL(j)*v2(1)*prop56+Q(j)*q1
        ZgR(j)=zL(j)*v2(2)*prop56+Q(j)*q1
      enddo

      do j=-nf,nf
      do k=-nf,nf
      if (Vsq(j,k)  /=  0._dp) then
        if     ((j > 0) .and. (k < 0)) then
            AWZM=(FAC*(ZgL(+j)*Fa213456+ZgL(-k)*Fa216543)
     &           +FACM*(v2(1)*cotw*prop56*Fb213456_z
     &                                +q1*Fb213456_g)*prop12)*prop34
            AWZP=(FAC*(ZgR(+j)*Fa213465+ZgR(-k)*Fa215643)
     &           +FACM*(v2(2)*cotw*prop56*Fb213465_z
     &                                +q1*Fb213465_g)*prop12)*prop34
            BWZM=(FAC*(ZgL(+j)*La213456+ZgL(-k)*La216543)
     &           +FACM*(v2(1)*cotw*prop56*Lb213456_z
     &                                +q1*Lb213456_g)*prop12)*prop34
            BWZP=(FAC*(ZgR(+j)*La213465+ZgR(-k)*La215643)
     &           +FACM*(v2(2)*cotw*prop56*Lb213465_z
     &                                +q1*Lb213465_g)*prop12)*prop34
        elseif ((j < 0) .and. (k > 0)) then
            AWZM=(FAC*(ZgL(+k)*Fa123456+ZgL(-j)*Fa126543)
     &           +FACM*(v2(1)*cotw*prop56*Fb123456_z
     &                                +q1*Fb123456_g)*prop12)*prop34
            AWZP=(FAC*(ZgR(+k)*Fa123465+ZgR(-j)*Fa125643)
     &           +FACM*(v2(2)*cotw*prop56*Fb123465_z
     &                                +q1*Fb123465_g)*prop12)*prop34
            BWZM=(FAC*(ZgL(+k)*La123456+ZgL(-j)*La126543)
     &           +FACM*(v2(1)*cotw*prop56*Lb123456_z
     &                                +q1*Lb123456_g)*prop12)*prop34
            BWZP=(FAC*(ZgR(+k)*La123465+ZgR(-j)*La125643)
     &           +FACM*(v2(2)*cotw*prop56*Lb123465_z
     &                                +q1*Lb123465_g)*prop12)*prop34
        endif
        if (srdiags) then
c---we need supplementary diagrams for gauge invariance.
c---now also assume that we have lepton decay products for W
c---so that v2(1)=le, v2(2)=re
c---1st term is diagram where Z couples to electron
c---2nd term is diagram where Z couples to neutrino
c---3rd term is diagram where gamma* couples to electron
c---4th term (l-h only) contains two W propagators
        if     ((j > 0) .and. (k < 0)) then
            suppl=FAC*prop12*(
     &          (en1*Fa346512+en2*Fa342156)*v2(1)*prop56
     &          +q1*(-1._dp)*(cl1*Fa346512+cl2*Fa342156)
     &          +wwflag*0.5_dp/zxw*prop34*(cl1*Fa652143+cl2*Fa653412))
            AWZM=AWZM+suppl
            BWZM=BWZM+suppl*Vpole12
            suppl=FAC*prop12*(
     &          (en1*Fa345612+en2*Fa342165)*v2(2)*prop56
     &          +q1*(-1._dp)*(cl1*Fa345612+cl2*Fa342165))
            AWZP=AWZP+suppl
            BWZP=BWZP+suppl*Vpole12
        elseif ((j < 0) .and. (k > 0)) then
            suppl=FAC*prop12*(
     &          (en1*Fa346521+en2*Fa341256)*v2(1)*prop56
     &          +q1*(-1._dp)*(cl1*Fa346521+cl2*Fa341256)
     &          +wwflag*0.5_dp/zxw*prop34*(cl1*Fa651243+cl2*Fa653421))
            AWZM=AWZM+suppl
            BWZM=BWZM+suppl*Vpole12
            suppl=FAC*prop12*(
     &          (en1*Fa345621+en2*Fa341265)*v2(2)*prop56
     &          +q1*(-1._dp)*(cl1*Fa345621+cl2*Fa341265))
            AWZP=AWZP+suppl
            BWZP=BWZP+suppl*Vpole12
        endif
        endif

c-- Inclusion of width for W,Z a la Baur and Zeppenfeld
      AWZM=cprop*AWZM
      AWZP=cprop*AWZP
      BWZM=cprop*BWZM
      BWZP=cprop*BWZP

      if (interference .eqv. .false.) then
c--- normal case
      msqv(j,k)=FACNLO*Vsq(j,k)*ave*vsymfact*two*real(conjg(AWZM)*BWZM+conjg(AWZP)*BWZP)
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          AWZM_save(j,k)=AWZM
          AWZP_save(j,k)=AWZP
          BWZM_save(j,k)=BWZM
          BWZP_save(j,k)=BWZP
        else
          if (useblha == 1) then
            msqv(j,k)=FACNLO*Vsq(j,k)*ave*vsymfact*two*real(
     &        conjg(AWZM_save(j,k)-AWZM)*(BWZM_save(j,k)-BWZM)
     &       +conjg(AWZP)*BWZP
     &       +conjg(AWZP_save(j,k))*BWZP_save(j,k))
          else
            msqv(j,k)=FACNLO*Vsq(j,k)*ave*vsymfact*two*real(
     &        conjg(AWZM_save(j,k)-AWZM)*(BWZM_save(j,k)-BWZM)
     &       +conjg(AWZM_save(j,k))*BWZM_save(j,k)-conjg(AWZM)*BWZM
     &       +two*conjg(AWZP_save(j,k))*BWZP_save(j,k))
          endif
        endif
      endif


      endif

      enddo
      enddo

      enddo   ! end ii loop

      return
      end

