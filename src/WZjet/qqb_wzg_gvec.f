!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wzg_gvec(p,n,in,msq)
      implicit none
c----Matrix element for WZ production
c----averaged over initial colours and spins
c    contracted with the vector v(mu)
c For nwz=-1
c  d(-p1)+ubar(-p2)--> W^-(e-(p3)+ve~(p4))+ Z(mu^-(p5)+mu^+(p6))+g(p7)
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'pchoice.f'
      include 'srdiags.f'
      include 'nwz.f'
      integer in,polz,iww
      integer,parameter::minus=1,mplus=2,jtype=7
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),
     & s127,fac,fac1
      complex(dp):: v56(2),en(4:5),qq(4:5),zcotw
      complex(dp):: amp(2),ZgLR(2,2),c1(2),c2(2),
     & d_ub(jtype,2),ub_d(jtype,2),d_g(jtype,2),g_d(jtype,2),
     & ub_g(jtype,2),g_ub(jtype,2),prop34,prop56,prop127

      fac1=2._dp*V*gsq*(esq**2/abs(zxw))**2
      v56(1)=zl1
      v56(2)=zr1
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

      call spinoru(7,p,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop34=cplx1(s(3,4))/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      prop56=cplx1(s(5,6))/cplx2(s(5,6)-zmass**2,zmass*zwidth)
      prop127=cplx1(s127)/cplx2(s127-wmass**2,wmass*wwidth)

      d_ub(:,:)=czip
      ub_d(:,:)=czip
      g_ub(:,:)=czip
      g_d(:,:)=czip
      d_g(:,:)=czip
      ub_g(:,:)=czip

      if     (in  ==  7) then
        call a7WZn(1,2,3,4,5,6,7,p,n,d_ub)
        call a7WZn(2,1,3,4,5,6,7,p,n,ub_d)
      elseif (in  ==  1) then
        call a7WZn(7,2,3,4,5,6,1,p,n,g_ub)
        call a7WZn(2,7,3,4,5,6,1,p,n,g_d)
      elseif (in  ==  2) then
        call a7WZn(1,7,3,4,5,6,2,p,n,d_g)
        call a7WZn(7,1,3,4,5,6,2,p,n,ub_g)
      else
        write(6,*) 'Invalid value of in=',in,' in qqb_wzg_gvec.f!'
        stop
      endif

      do j=1,2
        ZgLR(j,minus)=Q(j)*q1+zL(j)*prop56*v56(1)
        ZgLR(j,mplus)=Q(j)*q1+zL(j)*prop56*v56(2)
      enddo

      do polz=1,2
      if     (nwz == +1) then
        c1(polz)=ZgLR(2,polz)
        c2(polz)=ZgLR(1,polz)
       elseif (nwz == -1) then
        c1(polz)=ZgLR(1,polz)
        c2(polz)=ZgLR(2,polz)
      endif
      enddo

      msq(:,:)=zip
      do j=-2,2
      if (((j==+1).or.(j==-2)) .and. (nwz == +1)) cycle
      if (((j==-1).or.(j==+2)) .and. (nwz == -1)) cycle
      do k=-2,2
      if (((k==+1).or.(k==-2)) .and. (nwz == +1)) cycle
      if (((k==-1).or.(k==+2)) .and. (nwz == -1)) cycle
      if (j == k)  cycle

c--- only do g-q and g-qb for in=1
      if ((in == 1) .and. (j  /=  0)) cycle

c--- only do q-g and qb-g for in=2
      if ((in == 2) .and. (k  /=  0)) cycle

c--- only do q-qb and qb-q for in=7
      if ((in == 7) .and. ((j == 0) .or. (k == 0))) cycle

      if     ((j > 0) .and. (k < 0)) then
      fac=aveqq*fac1
      amp(:)=prop34*(
     & +d_ub(1,:)*ZgLR(j,:)
     & +d_ub(2,:)*ZgLR(-k,:)
     & -nwz*d_ub(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     & d_ub(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     &+d_ub(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     &+d_ub(iww,:)*prop34/(2._dp*zxw))
      endif

      elseif ((j < 0) .and. (k > 0)) then
      fac=aveqq*fac1
      amp(:)=prop34*(
     & +ub_d(1,:)*ZgLR(k,:)
     & +ub_d(2,:)*ZgLR(-j,:)
     & -nwz*ub_d(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     &  ub_d(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_d(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_d(iww,:)*prop34/(2._dp*zxw))
      endif

      elseif ((j > 0) .and. (k == 0)) then
      fac=aveqg*fac1
      amp(:)=prop34*(
     & +d_g(1,:)*c1(:)
     & +d_g(2,:)*c2(:)
     & -nwz*d_g(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     &  d_g(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +d_g(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +d_g(iww,:)*prop34/(2._dp*zxw))
      endif

      elseif ((j < 0) .and. (k == 0)) then
      fac=aveqg*fac1
      amp(:)=prop34*(
     & +ub_g(1,:)*c1(:)
     & +ub_g(2,:)*c2(:)
     & -nwz*ub_g(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     &  ub_g(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +ub_g(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +ub_g(iww,:)*prop34/(2._dp*zxw))
      endif

      elseif ((j == 0) .and. (k > 0)) then
      fac=aveqg*fac1
      amp(:)=prop34*(
     & +g_d(1,:)*c1(:)
     & +g_d(2,:)*c2(:)
     & -nwz*g_d(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     &  g_d(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_d(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_d(iww,:)*prop34/(2._dp*zxw))
      endif

      elseif ((j == 0) .and. (k < 0)) then
      fac=aveqg*fac1
      amp(:)=prop34*(
     & +g_ub(1,:)*c1(:)
     & +g_ub(2,:)*c2(:)
     & -nwz*g_ub(3,:)*prop127*(q1+zcotw*prop56*v56(:)))
      if (srdiags) then
      amp(:)=amp(:)+prop127*(
     &  g_ub(4,:)*(qq(4)*q1+en(4)*prop56*v56(:))
     & +g_ub(5,:)*(qq(5)*q1+en(5)*prop56*v56(:))
     & +g_ub(iww,:)*prop34/(2._dp*zxw))
      endif

      endif

      msq(j,k)=fac*(cdabs(amp(minus))**2+cdabs(amp(mplus))**2)

      enddo
      enddo

c Extend to other flavours
      do j=3,4
      msq(j,0)=msq((j-2),0)
      msq(0,j)=msq(0,(j-2))
      msq(-j,0)=msq(-(j-2),0)
      msq(0,-j)=msq(0,-(j-2))
      enddo
      msq(+3,-4)=msq(+1,-2)
      msq(-4,+3)=msq(-2,+1)
      msq(+4,-3)=msq(+2,-1)
      msq(-3,+4)=msq(-1,+2)
      return
      end

