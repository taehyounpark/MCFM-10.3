!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_vv(p,msq)
      implicit none
      include 'types.f'

c----Author R.K.Ellis December 1998
c----modified by JMC to include supplementary diagrams February 1999
c----modified by RKE (following suggestion of GZ)
c----to correct supplementary diagrams April 2011
c----Matrix element for ZZ production
c----NB: we also include virtual photons
c    in the notation of DKS
c    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->e^-(p3) + e^+(p4)   + \mu^-(p5)+ \mu^+(p6)
c    q(-p1)+qbar(-p2)-->l'(p3)  + lbar'(p4) + l(p5)    + lbar(p6)
c    with Z-leptons couplings l1 for (5,6) and l2 for(3,4)
c          and lepton charges q2 for (5,6) and q1 for (3,4)
c----No statistical factor of 1/2 included.

c--- Includes both ZZ and WW contributions, taken from qqb_zz.f & qqb_ww.f
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      include 'zcouple_cms.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'

      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),
     & ave,rescale1,rescale2

      complex(dp):: v1(2),v2(2)
      complex(dp):: qqb(2,2,2),qbq(2,2,2),q_qb,qb_q
      complex(dp):: qqb1(2,2,2),qbq1(2,2,2),qqb2(2,2,2),qbq2(2,2,2)
      complex(dp):: a6trees,prop12,prop34,prop56

      complex(dp):: FACWW
      integer:: j,k,polq,pol1,pol2
      integer:: ii,nmax
      integer,parameter::i4(2)=(/4,5/),i5(2)=(/5,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      complex(dp):: AWW(2),a6treea,A6b_1,A6b_2,A6b_3
      complex(dp):: prop34ww,prop56ww
      complex(dp):: Fa123456,Fa213456,Fb123456_z,Fb213456_z
      complex(dp):: Fa126543,Fa216543,Fb126543_z,Fb216543_z
      complex(dp):: Fb123456_g,Fb213456_g,Fb126543_g,Fb216543_g
      complex(dp):: Fa341256,Fa653421,Fa346521,Fa651243
      complex(dp):: Fa342156,Fa653412,Fa346512,Fa652143
      complex(dp):: cs_z(2,2),cs_g(2,2),cgamz(2,2),cz(2,2)
      complex(dp):: azz_pol(2,2,2),aww_pol(2,2,2)
      real(dp):: fac
      real(dp), parameter :: mp(nf)=(/-1._dp,+1._dp,-1._dp,+1._dp,-1._dp/)
      integer:: jk,tjk,minus,mplus
      parameter(minus=1,mplus=2)

      fac=-four*abs(zesq)**2
      ave=aveqq*xn

c Fix couplings by removing factors of 1/sqrt(3) (and restore at end)
      zl2=zln
      zr2=zrn
      rescale1=one
      rescale2=one

      v1(1)=zl1
      v1(2)=zr1
      v2(1)=zl2
      v2(2)=zr2

c--set msq=0 to initalize
      msq(:,:)=zip

c----Change the momenta to DKS notation
c   We have --- q(-p1)+qbar(-p2)-->l(p3)+lbar(p4) + l'(p5)+lbar'(p6)
c   DKS have--- q(q2) +qbar(q1) -->mu^-(q3)+mu^+(q4)+e^-(q6)+e^+(q5)

         do j=1,4
            qdks(1,j)=p(1,j)
            qdks(2,j)=p(2,j)
            qdks(3,j)=p(3,j)
            qdks(4,j)=p(4,j)
            qdks(5,j)=p(6,j)
            qdks(6,j)=p(5,j)
         enddo

      call spinoru(6,qdks,za,zb)

      nmax=1

      do ii=1,nmax


c--   s returned from sprod (common block) is 2*dot product
c--   calculate propagators
      prop12=s(1,2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(i5(ii),6)/cplx2(s(i5(ii),6)-zmass**2,zmass*zwidth)

c-- here the labels correspond to the polarizations of the
c-- quark, lepton 4 and lepton 6 respectively

      qbq(1,1,1)=A6trees(1,2,6,i5(ii),i4(ii),3,za,zb)
      qbq(1,1,2)=A6trees(1,2,6,i5(ii),3,i4(ii),za,zb)
      qbq(1,2,1)=A6trees(1,2,i5(ii),6,i4(ii),3,za,zb)
      qbq(1,2,2)=A6trees(1,2,i5(ii),6,3,i4(ii),za,zb)

      qqb(1,1,1)=A6trees(2,1,6,i5(ii),i4(ii),3,za,zb)
      qqb(1,1,2)=A6trees(2,1,6,i5(ii),3,i4(ii),za,zb)
      qqb(1,2,1)=A6trees(2,1,i5(ii),6,i4(ii),3,za,zb)
      qqb(1,2,2)=A6trees(2,1,i5(ii),6,3,i4(ii),za,zb)

      if (srdiags) then
c---for supplementary diagrams.
      qbq1(1,1,1)=+A6trees(3,i4(ii),1,2,i5(ii),6,za,zb)
      qbq2(1,1,1)=+A6trees(6,i5(ii),1,2,i4(ii),3,za,zb)
      qbq1(1,1,2)=-A6trees(i4(ii),3,1,2,i5(ii),6,za,zb)
      qbq2(1,1,2)=+A6trees(6,i5(ii),1,2,3,i4(ii),za,zb)
      qbq1(1,2,1)=+A6trees(3,i4(ii),1,2,6,i5(ii),za,zb)
      qbq2(1,2,1)=-A6trees(i5(ii),6,1,2,i4(ii),3,za,zb)
      qbq1(1,2,2)=-A6trees(i4(ii),3,1,2,6,i5(ii),za,zb)
      qbq2(1,2,2)=-A6trees(i5(ii),6,1,2,3,i4(ii),za,zb)

      qqb1(1,1,1)=-A6trees(3,i4(ii),2,1,i5(ii),6,za,zb)
      qqb2(1,1,1)=-A6trees(6,i5(ii),2,1,i4(ii),3,za,zb)
      qqb1(1,1,2)=+A6trees(i4(ii),3,2,1,i5(ii),6,za,zb)
      qqb2(1,1,2)=-A6trees(6,i5(ii),2,1,3,i4(ii),za,zb)
      qqb1(1,2,1)=-A6trees(3,i4(ii),2,1,6,i5(ii),za,zb)
      qqb2(1,2,1)=+A6trees(i5(ii),6,2,1,i4(ii),3,za,zb)
      qqb1(1,2,2)=+A6trees(i4(ii),3,2,1,6,i5(ii),za,zb)
      qqb2(1,2,2)=+A6trees(i5(ii),6,2,1,3,i4(ii),za,zb)
      endif

      do j=1,2
      do k=1,2
      qbq(2,j,k)=-qqb(1,j,k)
      qqb(2,j,k)=-qbq(1,j,k)
      qbq1(2,j,k)=-qqb1(1,j,k)
      qqb1(2,j,k)=-qbq1(1,j,k)
      qbq2(2,j,k)=-qqb2(1,j,k)
      qqb2(2,j,k)=-qbq2(1,j,k)
      enddo
      enddo

c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK

      FACWW=abs(zesq)**2/zxw**2 ! sqrt wrt. qqb_ww_v.f and averaging (xn*aveqq) removed
c---multiply by factor for c-sbar+u-dbar hadronic decay


c----Change the momenta to DKS notation
c   swapped possibility if we want to swap momenta
c   We have --- f(p1) + f'(p2)-->mu^-(p3)+nubar(p4)+e^+(p6)+nu(p5)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)
c----
c   or normal configuration
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)


c----all other cases
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j) ! interchanged 3<->5 wrt. qqb_ww.f
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(5,j) ! interchanged 3<->5 wrt. qqb_ww.f
      enddo

c-- s returned from sprod (common block) is 2*dot product
      call spinoru(6,qdks,za,zb)

c--   calculate propagators
      prop34ww=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      prop56ww=s(5,6)/cplx2(s(5,6)-wmass**2,wmass*wwidth)

c-- couplings with or without photon pole
      do j=1,2
      cs_z(minus,j)=+mp(j)*zl(j)*zsin2w*prop12
      cs_z(mplus,j)=-mp(j)*2._dp*Q(j)*zxw*prop12
      cs_g(minus,j)=+mp(j)*2._dp*Q(j)*zxw
      cs_g(mplus,j)=+mp(j)*2._dp*Q(j)*zxw
      if (srdiags) then
      cz(minus,j)=two*zxw*zln*zL(j)*prop12
      cz(mplus,j)=two*zxw*zln*zR(j)*prop12
      cgamz(minus,j)=two*zxw*(-Q(j)+zle*zL(j)*prop12)
      cgamz(mplus,j)=two*zxw*(-Q(j)+zle*zR(j)*prop12)
      endif
      enddo

      xdelg1_z=delg1_z
      xdelg1_g=delg1_g
      xdelk_z=delk_z
      xdelk_g=delk_g
      xlambda_z=lambda_z
      xlambda_g=lambda_g

c---case dbar-d and d-dbar

      Fa126543=A6treea(1,2,6,5,4,3,za,zb)
      Fa216543=A6treea(2,1,6,5,4,3,za,zb)
      Fa123456=A6treea(1,2,3,4,5,6,za,zb)
      Fa213456=A6treea(2,1,3,4,5,6,za,zb)

      call A6treeb_anom(1,2,3,4,5,6,za,zb,A6b_1,A6b_2,A6b_3)
      Fb123456_z=A6b_1*(2._dp+xdelg1_z+xdelk_z)
     &          +A6b_2*(2._dp*(1._dp+xdelg1_z))
     &          +A6b_3*(xlambda_z/wmass**2)
      Fb123456_g=A6b_1*(2._dp+xdelg1_g+xdelk_g)
     &          +A6b_2*(2._dp*(1._dp+xdelg1_g))
     &          +A6b_3*(xlambda_g/wmass**2)
      Fb126543_z=-Fb123456_z
      Fb126543_g=-Fb123456_g
      call A6treeb_anom(2,1,3,4,5,6,za,zb,A6b_1,A6b_2,A6b_3)
      Fb213456_z=A6b_1*(2._dp+xdelg1_z+xdelk_z)
     &          +A6b_2*(2._dp*(1._dp+xdelg1_z))
     &          +A6b_3*(xlambda_z/wmass**2)
      Fb213456_g=A6b_1*(2._dp+xdelg1_g+xdelk_g)
     &          +A6b_2*(2._dp*(1._dp+xdelg1_g))
     &          +A6b_3*(xlambda_g/wmass**2)
      Fb216543_z=-Fb213456_z
      Fb216543_g=-Fb213456_g

c      Fb123456=A6treeb(1,2,3,4,5,6,za,zb)
c      Fb126543=-Fb123456
c      Fb213456=A6treeb(2,1,3,4,5,6,za,zb)
c      Fb216543=-Fb213456

      if (srdiags) then
c---for supplementary diagrams.
      Fa341256=A6treea(3,4,1,2,5,6,za,zb)
      Fa653421=A6treea(6,5,3,4,2,1,za,zb)
      Fa346521=A6treea(3,4,6,5,2,1,za,zb)
      Fa651243=A6treea(6,5,1,2,4,3,za,zb)
      Fa342156=A6treea(3,4,2,1,5,6,za,zb)
      Fa653412=A6treea(6,5,3,4,1,2,za,zb)
      Fa346512=A6treea(3,4,6,5,1,2,za,zb)
      Fa652143=A6treea(6,5,2,1,4,3,za,zb)
      endif

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      aww_pol(:,:,:)=czip

      do j=-2,2
      k=-j
      jk=max(j,k)

      if (j==0) go to 20

      if ((j > 0).and.(k < 0)) then
      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq == 1) then
       q_qb=(prop56*v2(pol1)*zl(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zl(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (srdiags) then
          q_qb=q_qb-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(j)+q2*q(j))*qqb2(polq,pol1,pol2))

        endif
      elseif (polq == 2) then
       q_qb=(prop56*v2(pol1)*zr(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zr(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (srdiags) then
         q_qb=q_qb-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(j)+q2*q(j))*qqb2(polq,pol1,pol2))

        endif
      endif

      azz_pol(polq,pol1,pol2)=FAC*q_qb

      enddo
      enddo
      enddo

      elseif ((j < 0).and.(k > 0)) then

      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq == 1) then
       qb_q=(prop56*v2(pol1)*zl(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zl(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(k)+q2*q(k))*qbq2(polq,pol1,pol2)

        endif
      elseif (polq == 2) then
       qb_q=(prop56*v2(pol1)*zr(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zr(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(k)+q2*q(k))*qbq2(polq,pol1,pol2)

        endif
      endif

      azz_pol(polq,pol1,pol2)=FAC*qb_q

      enddo !end pol2 loop
      enddo !end pol1 loop
      enddo !end polq loop

      endif

c!! WW block
c!! WW block

c--assign values
c---Remember that base process is ub-u so this has the natural 123456 order
      if (j > 0) then
          if         (tau(jk) == +1._dp) then
            AWW(minus)=(Fa213456+cs_z(minus,2)*Fb213456_z
     &                          +cs_g(minus,2)*Fb213456_g)*prop56ww*prop34ww
            AWW(mplus)=(cs_z(mplus,2)*Fb123456_z
     &                 +cs_g(mplus,2)*Fb123456_g)*prop56ww*prop34ww
          elseif     (tau(jk) == -1._dp) then
            AWW(minus)=(Fa216543+cs_z(minus,1)*Fb216543_z
     &                          +cs_g(minus,1)*Fb216543_g)*prop56ww*prop34ww
            AWW(mplus)=(cs_z(mplus,1)*Fb126543_z
     &                 +cs_g(mplus,1)*Fb126543_g)*prop56ww*prop34ww
          endif
      elseif (j < 0) then
          if     (tau(jk) == +1._dp) then
            AWW(minus)=(Fa123456+cs_z(minus,2)*Fb123456_z
     &                          +cs_g(minus,2)*Fb123456_g)*prop56ww*prop34ww
            AWW(mplus)=(cs_z(mplus,2)*Fb213456_z
     &                 +cs_g(mplus,2)*Fb213456_g)*prop56ww*prop34ww
          elseif (tau(jk) == -1._dp) then
            AWW(minus)=(Fa126543+cs_z(minus,1)*Fb126543_z
     &                          +cs_g(minus,1)*Fb126543_g)*prop56ww*prop34ww
            AWW(mplus)=(cs_z(mplus,1)*Fb216543_z
     &                 +cs_g(mplus,1)*Fb216543_g)*prop56ww*prop34ww
          endif
      endif

      if (srdiags) then
c---we need supplementary diagrams for gauge invariance.
c---tjk is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)
      if (j > 0) then
          AWW(minus)=AWW(minus)
     &              +cgamz(minus,tjk)*(Fa342156*prop56ww+Fa653412*prop34ww)
     &                 +cz(minus,tjk)*(Fa346512*prop56ww+Fa652143*prop34ww)
          AWW(mplus)=AWW(mplus)
     &              +cgamz(mplus,tjk)*(Fa341256*prop56ww+Fa653421*prop34ww)
     &                 +cz(mplus,tjk)*(Fa346521*prop56ww+Fa651243*prop34ww)
      elseif (j < 0) then
          AWW(minus)=AWW(minus)
     &             +cgamz(minus,tjk)*(Fa341256*prop56ww+Fa653421*prop34ww)
     &                +cz(minus,tjk)*(Fa346521*prop56ww+Fa651243*prop34ww)
          AWW(mplus)=AWW(mplus)
     &             +cgamz(mplus,tjk)*(Fa342156*prop56ww+Fa653412*prop34ww)
     &                +cz(mplus,tjk)*(Fa346512*prop56ww+Fa652143*prop34ww)
      endif
      endif

      aww_pol(1,1,1)=facWW*AWW(minus)
      aww_pol(2,1,1)=facWW*AWW(mplus)

c!!
c!!

      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      msq(j,k)=msq(j,k)+ave*abs(azz_pol(polq,pol1,pol2)+aww_pol(polq,pol1,pol2))**2
      enddo
      enddo
      enddo

  20   continue
      enddo !end j-loop
      enddo   ! end ii loop

c----extend matrix element to full flavor range
      do j=-nf,nf
      k=-j
      if (j==0) go to 21
      msq(j,k)=msq(jkswitch(j),jkswitch(k))
 21   continue
      enddo

c      msq(:,:)=msq(:,:)*3d0 ! artificially inflate by factor of three

      return
      end
