!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wzg_v(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nf.f'
      include 'masses.f'
      include 'scheme.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'plabel.f'
      include 'nwz.f'
      include 'ckm.f'
      include 'blha.f'
      integer j,k,minus,mplus,polg,polz
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),fac,fac1,facm,s127,ave
      complex(dp):: cprop,prop12,prop34,prop56,A0(2,2),A1(2,2)
      complex(dp)::
     & ubd_da_0(2,2),ubd_ua_0(2,2),ubd_b_0(2,2),ubd_sr_0(2,2),
     & ubd_da_1(2,2),ubd_ua_1(2,2),ubd_b_1(2,2),ubd_sr_1(2,2),ubd_q_1(2,2),
     & dub_da_0(2,2),dub_ua_0(2,2),dub_b_0(2,2),dub_sr_0(2,2),
     & dub_da_1(2,2),dub_ua_1(2,2),dub_b_1(2,2),dub_sr_1(2,2),dub_q_1(2,2),
     & gd_da_0(2,2) ,gd_ua_0(2,2) ,gd_b_0(2,2) ,gd_sr_0(2,2) ,
     & gd_da_1(2,2) ,gd_ua_1(2,2) ,gd_b_1(2,2) ,gd_sr_1(2,2) ,gd_q_1(2,2),
     & dg_da_0(2,2) ,dg_ua_0(2,2) ,dg_b_0(2,2) ,dg_sr_0(2,2) ,
     & dg_da_1(2,2) ,dg_ua_1(2,2) ,dg_b_1(2,2) ,dg_sr_1(2,2) ,dg_q_1(2,2),
     & ubg_da_0(2,2),ubg_ua_0(2,2),ubg_b_0(2,2),ubg_sr_0(2,2),
     & ubg_da_1(2,2),ubg_ua_1(2,2),ubg_b_1(2,2),ubg_sr_1(2,2),ubg_q_1(2,2),
     & gub_da_0(2,2),gub_ua_0(2,2),gub_b_0(2,2),gub_sr_0(2,2),
     & gub_da_1(2,2),gub_ua_1(2,2),gub_b_1(2,2),gub_sr_1(2,2),gub_q_1(2,2)
      complex(dp):: v2(2),cl1,cl2,zcotw,ZglR(nf,2),c1(2),c2(2)
      parameter(minus=1,mplus=2)

      scheme='tH-V'
c--- note that scheme = "tH-V" is also a valid choice in this case
c--- (finite terms relating the two schemes are included in, e.g. AWWjeta.f)

      msq(:,:)=0d0

c      fac=-2._dp*gwsq*esq
      fac=-2._dp*esq**2/abs(zxw)
      fac1=two*gsq*cf

      if ((nwz==1) .or. (nwz == -1)) then
        facm=nwz*fac
      else
        write(6,*) 'nwz  /=  +1 or -1 in qqb_wz_g.f'
        stop
      endif

      if     (nwz==-1) then
        cl1=1._dp
        cl2=0._dp
      elseif (nwz==+1) then
        cl1=0._dp
        cl2=1._dp
      endif
      v2(1)=zl1
      v2(2)=zr1

c-- if Z -> neutrinos, we need to switch c1 and c2
      if (plabel(5) == 'nl') then
        cl1=1._dp-cl1
        cl2=1._dp-cl2
      endif

      call spinoru(7,p,za,zb)
c--   s returned from sprodx (common block) is 2*dot product

c--   calculate propagators
      zcotw=sqrt((one-zxw)/zxw)
      s127=s(1,2)+s(1,7)+s(2,7)
      prop12=s127/cplx2(s127-wmass**2,wmass*wwidth)
      prop34=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
      cprop=cone

c--- Fill leading order and virtual amplitudes
      if (useblha == 0) then
c--- gluon-quark and gluon-antiquark
      call amp_WZg_v(p,2,7,3,4,6,5,1,gub_da_0,gub_ua_0,gub_b_0,gub_sr_0,gub_da_1,gub_ua_1,gub_b_1,gub_q_1,gub_sr_1)
      call amp_WZg_v(p,7,2,3,4,6,5,1,gd_da_0,gd_ua_0,gd_b_0,gd_sr_0,gd_da_1,gd_ua_1,gd_b_1,gd_q_1,gd_sr_1)

c--- quark-gluon and antiquark-gluon
      call amp_WZg_v(p,1,7,3,4,6,5,2,ubg_da_0,ubg_ua_0,ubg_b_0,ubg_sr_0,ubg_da_1,ubg_ua_1,ubg_b_1,ubg_q_1,ubg_sr_1)
      call amp_WZg_v(p,7,1,3,4,6,5,2,dg_da_0,dg_ua_0,dg_b_0,dg_sr_0,dg_da_1,dg_ua_1,dg_b_1,dg_q_1,dg_sr_1)

c--- quark-antiquark and antiquark-quark
      call amp_WZg_v(p,1,2,3,4,6,5,7,ubd_da_0,ubd_ua_0,ubd_b_0,ubd_sr_0,ubd_da_1,ubd_ua_1,ubd_b_1,ubd_q_1,ubd_sr_1)
      endif
      call amp_WZg_v(p,2,1,3,4,6,5,7,dub_da_0,dub_ua_0,dub_b_0,dub_sr_0,dub_da_1,dub_ua_1,dub_b_1,dub_q_1,dub_sr_1)

c---set up left/right handed couplings for both Z and gamma*
c---note that the second label corresponds to the helicity
c---of the LEPTON coupling v2, NOT the quarks (all L)
      do j=1,nf
        ZgLR(j,minus)=zL(j)*v2(1)*prop56+Q(j)*q1
        ZgLR(j,mplus)=zL(j)*v2(2)*prop56+Q(j)*q1
      enddo

      do polz=1,2
      if (nwz == 1) then
        c1(polz)=ZgLR(2,polz)
        c2(polz)=ZgLR(1,polz)
      else
        c1(polz)=ZgLR(1,polz)
        c2(polz)=ZgLR(2,polz)
      endif
      enddo

      do j=-nf,nf
      if (((j==+1).or.(j==+3).or.(j==+5)
     & .or.(j==-2).or.(j==-4)) .and. (nwz == +1))
     & go to 20
      if (((j==-1).or.(j==-3).or.(j==-5)
     & .or.(j==+2).or.(j==+4)) .and. (nwz == -1))
     & go to 20
      do k=-nf,nf
      if (((k==+1).or.(k==+3).or.(k==+5)
     & .or.(k==-2).or.(k==-4)) .and. (nwz == +1))
     & go to 19
      if (((k==-1).or.(k==-3).or.(k==-5)
     & .or.(k==+2).or.(k==+4)) .and. (nwz == -1))
     & go to 19

      if     ((j < 0) .and. (k > 0)) then

c---case ub-d
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((ZgLR(+k,polz)*ubd_da_0(polg,polz)
     &                 +ZgLR(-j,polz)*ubd_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*ubd_b_0(polg,polz)*prop12*FACM)*prop34
     &                -ubd_sr_0(polg,polz)*FACM

        A1(polg,polz)=((ZgLR(+k,polz)*ubd_da_1(polg,polz)
     &                 +ZgLR(-j,polz)*ubd_ua_1(polg,polz)
     &                 +v2(polz)*prop56*ubd_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*ubd_b_1(polg,polz)*prop12*FACM)*prop34
     &                -ubd_sr_1(polg,polz)*FACM

        enddo
        enddo
        ave=xn*aveqq*Vsq(j,k)

      elseif ((j > 0) .and. (k < 0)) then
c---case d-ub
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((ZgLR(+j,polz)*dub_da_0(polg,polz)
     &                 +ZgLR(-k,polz)*dub_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*dub_b_0(polg,polz)*prop12*FACM)*prop34
     &                -dub_sr_0(polg,polz)*FACM

        A1(polg,polz)=((ZgLR(+j,polz)*dub_da_1(polg,polz)
     &                 +ZgLR(-k,polz)*dub_ua_1(polg,polz)
     &                 +v2(polz)*prop56*dub_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*dub_b_1(polg,polz)*prop12*FACM)*prop34
     &                -dub_sr_1(polg,polz)*FACM
        enddo
        enddo
        ave=xn*aveqq*Vsq(j,k)

      elseif ((j > 0) .and. (k == 0)) then
c---case d-g
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((c1(polz)*dg_da_0(polg,polz)
     &                 +c2(polz)*dg_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*dg_b_0(polg,polz)*prop12*FACM)*prop34
     &                -dg_sr_0(polg,polz)*FACM

        A1(polg,polz)=((c1(polz)*dg_da_1(polg,polz)
     &                 +c2(polz)*dg_ua_1(polg,polz)
     &                 +v2(polz)*prop56*dg_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*dg_b_1(polg,polz)*prop12*FACM)*prop34
     &                -dg_sr_1(polg,polz)*FACM
        enddo
        enddo
        ave=xn*aveqg*Vsum(j)

      elseif ((j < 0) .and. (k == 0)) then
c---case ub-g
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((c1(polz)*ubg_da_0(polg,polz)
     &                 +c2(polz)*ubg_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*ubg_b_0(polg,polz)*prop12*FACM)*prop34
     &                -ubg_sr_0(polg,polz)*FACM

        A1(polg,polz)=((c1(polz)*ubg_da_1(polg,polz)
     &                 +c2(polz)*ubg_ua_1(polg,polz)
     &                 +v2(polz)*prop56*ubg_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*ubg_b_1(polg,polz)*prop12*FACM)*prop34
     &                -ubg_sr_1(polg,polz)*FACM
        enddo
        enddo
        ave=xn*aveqg*Vsum(j)

      elseif ((j == 0) .and. (k > 0)) then
c---case g-d
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((c1(polz)*gd_da_0(polg,polz)
     &                 +c2(polz)*gd_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*gd_b_0(polg,polz)*prop12*FACM)*prop34
     &                -gd_sr_0(polg,polz)*FACM

        A1(polg,polz)=((c1(polz)*gd_da_1(polg,polz)
     &                 +c2(polz)*gd_ua_1(polg,polz)
     &                 +v2(polz)*prop56*gd_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*gd_b_1(polg,polz)*prop12*FACM)*prop34
     &                -gd_sr_1(polg,polz)*FACM
        enddo
        enddo
        ave=xn*aveqg*Vsum(k)

      elseif ((j == 0) .and. (k < 0)) then
c---case g-ub
        do polg=1,2
        do polz=1,2
        A0(polg,polz)=((c1(polz)*gub_da_0(polg,polz)
     &                 +c2(polz)*gub_ua_0(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*gub_b_0(polg,polz)*prop12*FACM)*prop34
     &                -gub_sr_0(polg,polz)*FACM

        A1(polg,polz)=((c1(polz)*gub_da_1(polg,polz)
     &                 +c2(polz)*gub_ua_1(polg,polz)
     &                 +v2(polz)*prop56*gub_q_1(polg,polz))*FAC
     &                +(zcotw*v2(polz)*prop56+q1)*gub_b_1(polg,polz)*prop12*FACM)*prop34
     &                -gub_sr_1(polg,polz)*FACM
        enddo
        enddo
        ave=xn*aveqg*Vsum(k)

      else
          ave=0._dp
      endif

      if (ave > 0._dp) then
c--- (check of LO)
c        msq(j,k)=fac1*ave*abs(cprop)**2*real(
c     &   +A0(1,1)*conjg(A0(1,1))+A0(1,2)*conjg(A0(1,2))
c     &   +A0(2,1)*conjg(A0(2,1))+A0(2,2)*conjg(A0(2,2)),dp)
c--- (virtual)
        msq(j,k)=ason2pi*fac1*ave*abs(cprop)**2*real(
     &   +A1(1,1)*conjg(A0(1,1))+A1(1,2)*conjg(A0(1,2))
     &   +A1(2,1)*conjg(A0(2,1))+A1(2,2)*conjg(A0(2,2)),dp)
      endif

 19   continue
      enddo
 20   continue
      enddo

      return
      end

