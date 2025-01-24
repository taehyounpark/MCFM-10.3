!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wzg_g(p,msq)
c---  Author: R.K. Ellis
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  + e^-(p3)+ve~(p4)+mu^-(p5)+mu^+(p6)+ g(p7)+ g(p8)
      implicit none
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'pchoice.f'
      include 'srdiags.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'nflav.f'

      integer, parameter::jtype=7
      integer polg1,polg2,iww
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),ave,s3456,fac,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq
      complex(dp):: v56(2),en(4:5),qq(4:5),zcotw

      complex(dp)::prop34,prop56,prop3456,
     & AB(2,2,2),BA(2,2,2),ZgLR(2,2),
     & c1(2),c2(2),
     & d_ubAB(jtype,2,2,2),d_ubBA(jtype,2,2,2),
     & ub_dAB(jtype,2,2,2),ub_dBA(jtype,2,2,2),
     & ub_gAB(jtype,2,2,2),ub_gBA(jtype,2,2,2),
     & g_ubAB(jtype,2,2,2),g_ubBA(jtype,2,2,2),
     & d_gAB(jtype,2,2,2),d_gBA(jtype,2,2,2),
     & g_dAB(jtype,2,2,2),g_dBA(jtype,2,2,2),
     & g_gAB(jtype,2,2,2),g_gBA(jtype,2,2,2)

      msq(:,:)=0d0

c      fac=V*xn*(2d0*esq*gw**2*gsq)**2
      fac=V*xn*(2d0*esq**2/abs(zxw)*gsq)**2

      call spinoru(8,p,za,zb)
c--   s returned from sprodx (common block) is 2*dot product

c      go to 199

c--   calculate propagators
      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      prop34=cplx1(s(3,4))/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      prop56=cplx1(s(5,6))/cplx2(s(5,6)-zmass**2,zmass*zwidth)
      prop3456=cplx1(s3456)/cplx2(s3456-wmass**2,wmass*wwidth)

c---case d-ubar
      call wzggamps(1,2,3,4,5,6,7,8,d_ubAB,d_ubBA)
c---case ubar-d
      call wzggamps(2,1,3,4,5,6,7,8,ub_dAB,ub_dBA)
c---case d-g
      call wzggamps(1,7,3,4,5,6,2,8,d_gAB,d_gBA)
c---case g-d
      call wzggamps(2,7,3,4,5,6,1,8,g_dAB,g_dBA)
c---case g-ubar
      call wzggamps(7,2,3,4,5,6,1,8,g_ubAB,g_ubBA)
c---case ubar-g
      call wzggamps(7,1,3,4,5,6,2,8,ub_gAB,ub_gBA)
c---case gg
      call wzggamps(8,7,3,4,5,6,1,2,g_gAB,g_gBA)

      v56(1)=zle
      v56(2)=zre
      zcotw=zL(2)-zL(1)
      if     (nwz==-1) then
        iww=6
        en(4)=zle
        en(5)=zln
        qq(4)=qe
        qq(5)=0._dp
      elseif (nwz==+1) then
        iww=7
        en(4)=zln
        en(5)=zle
        qq(4)=0._dp
        qq(5)=qe
      endif

      do j=1,2
        ZgLR(j,:)=Q(j)*q1+zL(j)*prop56*v56(:)
      enddo
      if     (nwz == +1) then
        c1(:)=ZgLR(2,:)
        c2(:)=ZgLR(1,:)
       elseif (nwz == -1) then
        c1(:)=ZgLR(1,:)
        c2(:)=ZgLR(2,:)
      endif
c--- first populate the basic set of matrix elements
c--- (for the first generation of quarks)
      do j=-2,2
      if (((j==+1).or.(j==-2)) .and. (nwz == +1)) cycle
      if (((j==-1).or.(j==+2)) .and. (nwz == -1)) cycle
      do k=-2,2
      if (((k==+1).or.(k==-2)) .and. (nwz == +1)) cycle
      if (((k==-1).or.(k==+2)) .and. (nwz == -1)) cycle
      if ((j == k) .and. (abs(j) > 0)) cycle

      AB(:,:,:)=0
      BA(:,:,:)=0

      do polg1=1,2
      do polg2=1,2

c---case d-ubar
      if ((j > 0) .and. (k < 0)) then
c Identical particle factor of 1/2 for two gluons
      ave=half*aveqq*Vsq(j,k)
      AB(polg1,polg2,:)=prop34*(
     & +d_ubAB(1,polg1,polg2,:)*ZgLR(+j,:)
     & +d_ubAB(2,polg1,polg2,:)*ZgLR(-k,:)
     &-nwz*d_ubAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +d_ubBA(1,polg1,polg2,:)*ZgLR(+j,:)
     & +d_ubBA(2,polg1,polg2,:)*ZgLR(-k,:)
     &-nwz*d_ubBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +d_ubAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +d_ubAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +d_ubAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +d_ubBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +d_ubBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +d_ubBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif

c---case ubar-d
      elseif ((j < 0) .and. (k > 0)) then
c Identical particle factor of 1/2 for two gluons
      ave=half*aveqq*Vsq(j,k)
      AB(polg1,polg2,:)=prop34*(
     & +ub_dAB(1,polg1,polg2,:)*ZgLR(+k,:)
     & +ub_dAB(2,polg1,polg2,:)*ZgLR(-j,:)
     &-nwz*ub_dAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +ub_dBA(1,polg1,polg2,:)*ZgLR(+k,:)
     & +ub_dBA(2,polg1,polg2,:)*ZgLR(-j,:)
     &-nwz*ub_dBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +ub_dAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_dAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_dAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +ub_dBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_dBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_dBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif

c---case ubar-g
      elseif ((j < 0) .and. (k == 0)) then
      ave=aveqg*Vsum(j)
      AB(polg1,polg2,:)=prop34*(
     & +ub_gAB(1,polg1,polg2,:)*c1(:)
     & +ub_gAB(2,polg1,polg2,:)*c2(:)
     &-nwz*ub_gAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +ub_gBA(1,polg1,polg2,:)*c1(:)
     & +ub_gBA(2,polg1,polg2,:)*c2(:)
     &-nwz*ub_gBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +ub_gAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_gAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_gAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +ub_gBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_gBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_gBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif

c---case g-ubar
      elseif ((j == 0) .and. (k < 0)) then
      ave=aveqg*Vsum(k)
      AB(polg1,polg2,:)=prop34*(
     & +g_ubAB(1,polg1,polg2,:)*c1(:)
     & +g_ubAB(2,polg1,polg2,:)*c2(:)
     &-nwz*g_ubAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +g_ubBA(1,polg1,polg2,:)*c1(:)
     & +g_ubBA(2,polg1,polg2,:)*c2(:)
     &-nwz*g_ubBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +g_ubAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_ubAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_ubAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +g_ubBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_ubBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_ubBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif
c---case d-g
      elseif ((j > 0) .and. (k == 0)) then
      ave=aveqg*Vsum(j)
      AB(polg1,polg2,:)=prop34*(
     & +d_gAB(1,polg1,polg2,:)*c1(:)
     & +d_gAB(2,polg1,polg2,:)*c2(:)
     &-nwz*d_gAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +d_gBA(1,polg1,polg2,:)*c1(:)
     & +d_gBA(2,polg1,polg2,:)*c2(:)
     &-nwz*d_gBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +d_gAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +d_gAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +d_gAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +d_gBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +d_gBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +d_gBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif
c---case g-d
      elseif ((j == 0) .and. (k > 0)) then
      ave=aveqg*Vsum(k)
      AB(polg1,polg2,:)=prop34*(
     & +g_dAB(1,polg1,polg2,:)*c1(:)
     & +g_dAB(2,polg1,polg2,:)*c2(:)
     &-nwz*g_dAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +g_dBA(1,polg1,polg2,:)*c1(:)
     & +g_dBA(2,polg1,polg2,:)*c2(:)
     &-nwz*g_dBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +g_dAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_dAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_dAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +g_dBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_dBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_dBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif
c---case gg
      elseif ((j == 0) .and. (k == 0)) then
c factor two for udb + csb
      ave=two*avegg
      AB(polg1,polg2,:)=prop34*(
     & +g_gAB(1,polg1,polg2,:)*c1(:)
     & +g_gAB(2,polg1,polg2,:)*c2(:)
     &-nwz*g_gAB(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      BA(polg1,polg2,:)=prop34*(
     & +g_gBA(1,polg1,polg2,:)*c1(:)
     & +g_gBA(2,polg1,polg2,:)*c2(:)
     &-nwz*g_gBA(3,polg1,polg2,:)*prop3456*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      AB(polg1,polg2,:)=AB(polg1,polg2,:)+prop3456*(
     & +g_gAB(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_gAB(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_gAB(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      BA(polg1,polg2,:)=BA(polg1,polg2,:)+prop3456*(
     & +g_gBA(4,polg1,polg2,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_gBA(5,polg1,polg2,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_gBA(iww,polg1,polg2,:)*prop34/(2._dp*zxw))
      endif

      endif
      enddo
      enddo

      msq(j,k)=fac*ave*(
     & +cdabs(AB(1,1,1))**2+cdabs(AB(1,1,2))**2
     & +cdabs(AB(1,2,1))**2+cdabs(AB(1,2,2))**2
     & +cdabs(AB(2,2,2))**2+cdabs(AB(2,2,1))**2
     & +cdabs(AB(2,1,2))**2+cdabs(AB(2,1,1))**2

     & +cdabs(BA(1,1,1))**2+cdabs(BA(1,1,2))**2
     & +cdabs(BA(1,2,1))**2+cdabs(BA(1,2,2))**2
     & +cdabs(BA(2,2,2))**2+cdabs(BA(2,2,1))**2
     & +cdabs(BA(2,1,2))**2+cdabs(BA(2,1,1))**2
     & -(
     & +cdabs(AB(1,1,1)+BA(1,1,1))**2
     & +cdabs(AB(1,1,2)+BA(1,1,2))**2
     & +cdabs(AB(1,2,1)+BA(1,2,1))**2
     & +cdabs(AB(1,2,2)+BA(1,2,2))**2
     & +cdabs(AB(2,2,2)+BA(2,2,2))**2
     & +cdabs(AB(2,2,1)+BA(2,2,1))**2
     & +cdabs(AB(2,1,2)+BA(2,1,2))**2
     & +cdabs(AB(2,1,1)+BA(2,1,1))**2)/xn**2)

      enddo
      enddo

      msq(+3,-4)=msq(+1,-2)
      msq(-4,+3)=msq(-2,+1)
      msq(-3,+4)=msq(-1,+2)
      msq(+4,-3)=msq(+2,-1)
      do j=3,4
      msq(j,0)=msq((j-2),0)
      msq(0,j)=msq(0,(j-2))
      msq(-j,0)=msq(-(j-2),0)
      msq(0,-j)=msq(0,-(j-2))
      enddo

c      return

c---*****************4 quark diagrams********************

      call a8WZfourqsq(1,2,3,4,5,6,7,8,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq)

      if (nwz == -1) then
      msq(1,1)=ddxudsq
      msq(1,2)=half*duxuusq
      msq(1,3)=dsxussq+dsxdcsq
      msq(1,4)=dcxucsq
      msq(2,1)=half*udxuusq
      msq(2,3)=usxucsq

      msq(3,1)=msq(1,3)
      msq(3,2)=msq(1,4)
      msq(3,3)=msq(1,1)
      msq(3,4)=msq(1,2)
      msq(4,1)=msq(2,3)
      msq(4,3)=msq(2,1)

      if (nflav == 5) then
      msq(1,5)=dsxussq
      msq(3,5)=msq(1,5)
      msq(5,1)=dsxdcsq
      msq(5,3)=msq(5,1)
      endif
      elseif (nwz == +1) then
      msq(2,2)=ddxudsq
      msq(2,1)=half*duxuusq
      msq(2,4)=dsxussq+dsxdcsq
      msq(2,3)=dcxucsq
      msq(1,2)=half*udxuusq
      msq(1,4)=usxucsq

      msq(4,2)=msq(2,4)
      msq(4,1)=msq(2,3)
      msq(4,4)=msq(2,2)
      msq(4,3)=msq(2,1)
      msq(3,2)=msq(1,4)
      msq(3,4)=msq(1,2)

      if (nflav == 5) then
      msq(5,2)=msq(3,2)
      msq(5,4)=msq(3,2)
      msq(2,5)=msq(2,3)
      msq(4,5)=msq(2,3)
      endif
      endif

      call a8WZfourqsq(7,2,3,4,5,6,1,8,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq)

      if (nwz == -1) then
      msq(-1,1)=ddxdusq+dsxdc2sq
      msq(-1,3)=dsxdcsq
      msq(-2,4)=dcxucsq
      msq(-2,1)=msq(-2,1)+udxuusq+ddxudsq+dsxus2sq+usxuc2sq
      msq(-2,2)=duxuusq+dcxuc2sq
      msq(-2,3)=usxucsq+dsxussq

      msq(-4,1)=msq(-2,3)
      msq(-4,2)=msq(-2,4)
      msq(-4,3)=msq(-2,1)
      msq(-4,4)=msq(-2,2)
      msq(-3,1)=msq(-1,3)
      msq(-3,3)=msq(-1,1)

      if (nflav == 5) then
      msq(-5,1)=msq(-3,1)
      msq(-5,3)=msq(-3,1)
      msq(-2,5)=dsxussq
      msq(-4,5)=msq(-2,5)
      msq(-2,1)=msq(-2,1)+dsxus2sq
      msq(-4,3)=msq(-2,1)
      msq(-5,5)=2._dp*dsxdc2sq
      endif
      elseif (nwz == +1) then

      msq(-2,2)=ddxdusq+dsxdc2sq
      msq(-2,4)=dsxdcsq
      msq(-1,3)=dcxucsq
      msq(-1,2)=msq(-1,2)+udxuusq+ddxudsq+dsxus2sq+usxuc2sq
      msq(-1,1)=duxuusq+dcxuc2sq
      msq(-1,4)=usxucsq+dsxussq

      msq(-3,2)=msq(-1,4)
      msq(-3,1)=msq(-1,3)
      msq(-3,4)=msq(-1,2)
      msq(-3,3)=msq(-1,1)
      msq(-4,2)=msq(-2,4)
      msq(-4,4)=msq(-2,2)

      if (nflav == 5) then
      msq(-1,5)=msq(-1,3)
      msq(-3,5)=msq(-1,5)
      msq(-5,2)=usxucsq
      msq(-5,4)=msq(-5,2)
      msq(-1,2)=msq(-1,2)+usxuc2sq
      msq(-3,4)=msq(-1,2)
      msq(-5,5)=2._dp*dcxuc2sq
      endif
      endif

      call a8WZfourqsq(1,8,3,4,5,6,7,2,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq)

      if (nwz == -1) then
      msq(1,-4)=dcxucsq+dsxdcsq
      msq(1,-3)=dsxussq
      msq(1,-2)=msq(1,-2)+duxuusq+ddxdusq+dcxuc2sq+dsxdc2sq
      msq(1,-1)=ddxudsq+dsxus2sq
      msq(2,-2)=udxuusq+usxuc2sq
      msq(2,-4)=usxucsq

      msq(3,-4)=msq(1,-2)
      msq(3,-3)=msq(1,-1)
      msq(3,-2)=msq(1,-4)
      msq(3,-1)=msq(1,-3)
      msq(4,-2)=msq(2,-4)
      msq(4,-4)=msq(2,-2)

      if (nflav == 5) then
      msq(1,-5)=msq(1,-3)
      msq(3,-5)=msq(1,-3)
      msq(5,-2)=dsxdcsq
      msq(5,-4)=msq(5,-2)
      msq(1,-2)=msq(1,-2)+dsxdc2sq
      msq(3,-4)=msq(1,-2)
      msq(5,-5)=2._dp*dsxus2sq
      endif
      elseif (nwz == +1) then
      msq(2,-3)=dcxucsq+dsxdcsq
      msq(2,-4)=dsxussq
      msq(2,-1)=msq(2,-1)+duxuusq+ddxdusq+dcxuc2sq+dsxdc2sq
      msq(2,-2)=ddxudsq+dsxus2sq
      msq(1,-1)=udxuusq+usxuc2sq
      msq(1,-3)=usxucsq

      msq(4,-3)=msq(2,-1)
      msq(4,-4)=msq(2,-2)
      msq(4,-1)=msq(2,-3)
      msq(4,-2)=msq(2,-4)
      msq(3,-1)=msq(1,-3)
      msq(3,-3)=msq(1,-1)

      if (nflav == 5) then
      msq(5,-1)=msq(3,-1)
      msq(5,-3)=msq(3,-1)
      msq(4,-5)=dcxucsq
      msq(2,-5)=msq(4,-5)
      msq(2,-1)=msq(2,-1)+dcxuc2sq
      msq(4,-3)=msq(2,-1)
      msq(5,-5)=2._dp*usxuc2sq
      endif
      endif

      call a8WZfourqsq(7,8,3,4,5,6,1,2,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq)

      if (nwz == -1) then
      msq(-4,-3)=dsxdcsq
      msq(-2,-4)=dcxucsq+usxucsq
      msq(-2,-3)=dsxussq
      msq(-2,-2)=duxuusq
      msq(-2,-1)=half*ddxudsq
      msq(-1,-4)=dsxdcsq
      msq(-1,-2)=half*ddxdusq

      msq(-4,-4)=msq(-2,-2)
      msq(-4,-3)=msq(-2,-1)
      msq(-4,-2)=msq(-2,-4)
      msq(-4,-1)=msq(-2,-3)
      msq(-3,-4)=msq(-1,-2)
      msq(-3,-2)=msq(-1,-4)

      if (nflav == 5) then
      msq(-5,-4)=msq(-3,-2)
      msq(-5,-2)=msq(-3,-2)
      msq(-4,-5)=msq(-2,-3)
      msq(-2,-5)=msq(-2,-3)
      endif
      elseif (nwz == +1) then

      msq(-3,-4)=dsxdcsq
      msq(-1,-3)=dcxucsq+usxucsq
      msq(-1,-4)=dsxussq
      msq(-1,-1)=duxuusq
      msq(-1,-2)=half*ddxudsq
      msq(-2,-3)=dsxdcsq
      msq(-2,-1)=half*ddxdusq

      msq(-3,-3)=msq(-1,-1)
      msq(-3,-4)=msq(-1,-2)
      msq(-3,-1)=msq(-1,-3)
      msq(-3,-2)=msq(-1,-4)
      msq(-4,-3)=msq(-2,-1)
      msq(-4,-1)=msq(-2,-3)

      if (nflav == 5) then
      msq(-5,-1)=usxucsq
      msq(-5,-3)=msq(-5,-1)
      msq(-1,-5)=dcxucsq
      msq(-3,-5)=msq(-1,-5)
      endif
      endif

      return
      end

