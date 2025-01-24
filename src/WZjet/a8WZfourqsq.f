!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8WZfourqsq(j1,j2,j3,j4,j5,j6,j7,j8,
     & dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,ddxdusq,
     & dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq)

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'srdiags.f'
      include 'nwz.f'
c     Performs the square and the sum over polarization of the diagrams
c     which have the WZ radiated off the external line.
      integer j1,j2,j3,j4,j5,j6,j7,j8,h17,h28,h56,iq,iww,j,id,iu
      real(dp)::dcxucsq,dsxussq,dsxdcsq,usxucsq,duxuusq,udxuusq,ddxudsq,
     & ddxdusq,dcxuc2sq,dsxus2sq,dsxdc2sq,usxuc2sq,s3456,fac
      complex(dp):: Zq(2,2),v56(2),zcotw,en(4:5),qq(4:5)
      complex(dp):: prop34,prop56,prop3456,ZgLR(2,2),c1(2),c2(2),
     & coup17x28(2,2,2),coup28x17(2,2,2),
     & upper1(3,2,2),lower1(3,2,2),upper2(3,2,2),lower2(3,2,2),
     & updown1(2,2),updown2(2,2),downup1(2,2),downup2(2,2),
     & uppsr1(4:7,2,2),dwnsr1(4:7,2,2),
     & uppsr2(4:7,2,2),dwnsr2(4:7,2,2),
     & dcxuc(2,2,2),dcxuc2(2,2,2),dsxus(2,2,2),dsxdc(2,2,2),
     & usxuc(2,2,2),usxuc2(2,2,2),dsxus2(2,2,2),dsxdc2(2,2,2)
c      dcxuc(h17,h28,h56) etc.

      s3456=s(j3,j4)+s(j3,j5)+s(j3,j6)+s(j4,j5)+s(j4,j6)+s(j5,j6)
      prop3456=cmplx(s3456,kind=dp)
     &        /cmplx(s3456-wmass**2,wmass*wwidth,dp)
      prop34=cmplx(s(j3,j4),kind=dp)
     &      /cmplx(s(j3,j4)-wmass**2,wmass*wwidth,dp)
      prop56=cmplx(s(j5,j6),kind=dp)
     &      /cmplx(s(j5,j6)-zmass**2,zmass*zwidth,dp)

c     amplitude a(h17,h28,h34,h56)
c--- call all basic amplitudes
      call a8WZ(j1,j2,j3,j4,j5,j6,j7,j8,upper1)
      call a8WZ(j2,j1,j3,j4,j5,j6,j8,j7,lower1)

      call a8WZ(j1,j2,j3,j4,j5,j6,j8,j7,upper2)
      call a8WZ(j2,j1,j3,j4,j5,j6,j7,j8,lower2)

      call a8WZsr(j1,j2,j3,j4,j5,j6,j7,j8,uppsr1)
      call a8WZsr(j2,j1,j3,j4,j5,j6,j8,j7,dwnsr1)

      call a8WZsr(j1,j2,j3,j4,j5,j6,j8,j7,uppsr2)
      call a8WZsr(j2,j1,j3,j4,j5,j6,j7,j8,dwnsr2)

      call a8WxZ(j1,j2,j3,j4,j5,j6,j7,j8,updown1)
      call a8WxZ(j2,j1,j3,j4,j5,j6,j8,j7,downup1)

      call a8WxZ(j1,j2,j3,j4,j5,j6,j8,j7,updown2)
      call a8WxZ(j2,j1,j3,j4,j5,j6,j7,j8,downup2)

c--- initialize couplings of Z/gamma to quarks and lepton pairs.
      v56(1)=zle
      v56(2)=zre
      zcotw=zL(2)-zL(1)
      if     (nwz==-1) then
        iww=6
        en(4)=zle
        en(5)=zln
        qq(4)=qe
        qq(5)=0._dp
        id=1
        iu=2
      elseif (nwz==+1) then
        iww=7
        en(4)=zln
        en(5)=zle
        qq(4)=0._dp
        qq(5)=qe
        id=2
        iu=1
      endif

      v56(1)=zl1
      v56(2)=zr1

      do iq=1,2
      Zq(iq,1)=zL(iq)
      Zq(iq,2)=zR(iq)
      enddo

      do iq=1,2
      do h17=1,2
      do h56=1,2
      coup28x17(iq,h17,h56)=
     & (Q(iq)*q1+Zq(iq,h17)*v56(h56)*prop56)*prop34
      enddo
      enddo
      enddo

      do iq=1,2
      do h28=1,2
      do h56=1,2
      coup17x28(iq,h28,h56)=
     & (Q(iq)*q1+Zq(iq,h28)*prop56*v56(h56))*prop34
      enddo
      enddo
      enddo

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

      dcxuc(:,:,:)=czip
      dsxus(:,:,:)=czip
      dsxdc(:,:,:)=czip
      usxuc(:,:,:)=czip
      dsxdc2(:,:,:)=czip
      dcxuc2(:,:,:)=czip
      usxuc2(:,:,:)=czip
      dsxus2(:,:,:)=czip


c     dsxus

      h17=1
      do h28=1,2
      do h56=1,2
      dsxus(h17,h28,h56)=(
     & +upper1(2,h28,h56)*c2(h56)
     & +upper1(1,h28,h56)*c1(h56)
     & -nwz*upper1(3,h28,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456)*prop34
     & +updown1(h28,h56)*coup17x28(id,h28,h56)
      if (srdiags) then
      dsxus(h17,h28,h56)=dsxus(h17,h28,h56)+prop3456*(
     & +uppsr1(4,h28,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +uppsr1(5,h28,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +uppsr1(iww,h28,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

c     dsxus2 for ddxud
c     W on lower line with 7-8 swap
      h28=1
      do h17=1,2
      do h56=1,2
      dsxus2(h17,h28,h56)=(
     & -nwz*lower2(3,h17,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +lower2(2,h17,h56)*c2(h56)
     & +lower2(1,h17,h56)*c1(h56))*prop34
     & +downup2(h17,h56)*coup28x17(id,h17,h56)
      if (srdiags) then
      dsxus2(h17,h28,h56)=dsxus2(h17,h28,h56)+prop3456*(
     & +dwnsr2(4,h17,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +dwnsr2(5,h17,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +dwnsr2(iww,h17,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

c     do h17=1,2
c     do h28=1,2
c     do h56=1,2
c     write(6,*) 'dsxus',h17,h28,h56,dsxus(h17,h28,h56)
c     enddo
c     enddo
c     enddo
c     do h17=1,2
c     do h28=1,2
c     do h56=1,2
c     write(6,*) 'dsxus2',h17,h28,h56,dsxus2(h17,h28,h56)
c     enddo
c     enddo
c     enddo

c     dsxdc
      h28=1
      do h17=1,2
      do h56=1,2
      dsxdc(h17,h28,h56)=(
     & -nwz*lower1(3,h17,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +lower1(2,h17,h56)*c2(h56)
     & +lower1(1,h17,h56)*c1(h56))*prop34
     & +downup1(h17,h56)*coup28x17(id,h17,h56)
      if (srdiags) then
      dsxdc(h17,h28,h56)=dsxdc(h17,h28,h56)+prop3456*(
     & +dwnsr1(4,h17,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +dwnsr1(5,h17,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +dwnsr1(iww,h17,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

c     This is the new one
c     W on upperline with 78 swap
c     dsxdc2 for ddxdu
      h17=1
      do h28=1,2
      do h56=1,2
      dsxdc2(h17,h28,h56)=(
     & -nwz*upper2(3,h28,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +upper2(2,h28,h56)*c2(h56)
     & +upper2(1,h28,h56)*c1(h56))*prop34
     & +updown2(h28,h56)*coup17x28(id,h28,h56)

      if (srdiags) then
      dsxdc2(h17,h28,h56)=dsxdc2(h17,h28,h56)+prop3456*(
     & +uppsr2(4,h28,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +uppsr2(5,h28,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +uppsr2(iww,h28,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

c     do h17=1,2
c     do h28=1,2
c     do h56=1,2
c     write(6,*) 'dsxdc',h17,h28,h56,dsxdc(h17,h28,h56)
c     enddo
c     enddo
c     enddo!
c     do h17=1,2
c     do h28=1,2
c     do h56=1,2
c     write(6,*) 'dsxdc2',h17,h28,h56,dsxdc2(h17,h28,h56)
c     enddo
c     enddo
c     enddo

c     usxuc
      h28=1
      do h17=1,2
      do h56=1,2
      usxuc(h17,h28,h56)=(
     & -nwz*lower1(3,h17,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +lower1(2,h17,h56)*c2(h56)
     & +lower1(1,h17,h56)*c1(h56))*prop34
     & +downup1(h17,h56)*coup28x17(iu,h17,h56)

      if (srdiags) then
      usxuc(h17,h28,h56)=usxuc(h17,h28,h56)+prop3456*(
     & +dwnsr1(4,h17,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +dwnsr1(5,h17,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +dwnsr1(iww,h17,h56)*prop34/(2._dp*zxw))
      endif

      enddo
      enddo

c     usxuc2 for udxuu (7 and 8 swapped)
      h28=1
      do h17=1,2
      do h56=1,2
      usxuc2(h17,h28,h56)=(
     & -nwz*lower2(3,h17,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +lower2(2,h17,h56)*c2(h56)
     & +lower2(1,h17,h56)*c1(h56))*prop34
     & +downup2(h17,h56)*coup28x17(iu,h17,h56)

      if (srdiags) then
      usxuc2(h17,h28,h56)=usxuc2(h17,h28,h56)+prop3456*(
     & +dwnsr2(4,h17,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +dwnsr2(5,h17,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +dwnsr2(iww,h17,h56)*prop34/(2._dp*zxw))
      endif

      enddo
      enddo

c     dcxuc
      h17=1
      do h28=1,2
      do h56=1,2
      dcxuc(h17,h28,h56)=(
     & -nwz*upper1(3,h28,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +upper1(2,h28,h56)*c2(h56)
     & +upper1(1,h28,h56)*c1(h56))*prop34
     & +updown1(h28,h56)*coup17x28(iu,h28,h56)

      if (srdiags) then
      dcxuc(h17,h28,h56)=dcxuc(h17,h28,h56)+prop3456*(
     & +uppsr1(4,h28,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +uppsr1(5,h28,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +uppsr1(iww,h28,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

c     dcsuc2 for duxuu (7 and 8 swapped)
      h17=1
      do h28=1,2
      do h56=1,2
      dcxuc2(h17,h28,h56)=(
     & -nwz*upper2(3,h28,h56)*(q1+zcotw*v56(h56)*prop56)*prop3456
     & +upper2(2,h28,h56)*c2(h56)
     & +upper2(1,h28,h56)*c1(h56))*prop34
     & +updown2(h28,h56)*coup17x28(iu,h28,h56)

      if (srdiags) then
      dcxuc2(h17,h28,h56)=dcxuc2(h17,h28,h56)+prop3456*(
     & +uppsr2(4,h28,h56)*(qq(4)*q1+en(4)*prop56*v56(h56))
     & +uppsr2(5,h28,h56)*(qq(5)*q1+en(5)*prop56*v56(h56))
     & +uppsr2(iww,h28,h56)*prop34/(2._dp*zxw))
      endif
      enddo
      enddo

      dcxucsq=0._dp
      dsxussq=0._dp
      dsxdcsq=0._dp
      usxucsq=0._dp

      dcxuc2sq=0._dp
      dsxus2sq=0._dp
      dsxdc2sq=0._dp
      usxuc2sq=0._dp

      duxuusq=0._dp
      udxuusq=0._dp
      ddxudsq=0._dp
      ddxdusq=0._dp

      do h17=1,2
      do h28=1,2
      do h56=1,2
      dcxucsq=dcxucsq+abs(dcxuc(h17,h28,h56))**2
      dsxussq=dsxussq+abs(dsxus(h17,h28,h56))**2
      dsxdcsq=dsxdcsq+abs(dsxdc(h17,h28,h56))**2
      usxucsq=usxucsq+abs(usxuc(h17,h28,h56))**2

      dcxuc2sq=dcxuc2sq+abs(dcxuc2(h17,h28,h56))**2
      dsxus2sq=dsxus2sq+abs(dsxus2(h17,h28,h56))**2
      dsxdc2sq=dsxdc2sq+abs(dsxdc2(h17,h28,h56))**2
      usxuc2sq=usxuc2sq+abs(usxuc2(h17,h28,h56))**2

      duxuusq=duxuusq
     & +abs(dcxuc(h17,h28,h56))**2+abs(dcxuc2(h17,h28,h56))**2
      if (h17 == h28) duxuusq=duxuusq
     & +dble(dcxuc(h17,h28,h56)*dconjg(dcxuc2(h17,h28,h56)))*2._dp/xn

      udxuusq=udxuusq
     & +abs(usxuc(h17,h28,h56))**2+abs(usxuc2(h17,h28,h56))**2
      if (h17 == h28) udxuusq=udxuusq
     & +dble(usxuc(h17,h28,h56)*dconjg(usxuc2(h17,h28,h56)))*2._dp/xn

      ddxudsq=ddxudsq
     & +abs(dsxus(h17,h28,h56))**2+abs(dsxus2(h17,h28,h56))**2
      if (h17 == h28) ddxudsq=ddxudsq
     & +dble(dsxus(h17,h28,h56)*dconjg(dsxus2(h17,h28,h56)))*2._dp/xn

      ddxdusq=ddxdusq
     & +abs(dsxdc(h17,h28,h56))**2+abs(dsxdc2(h17,h28,h56))**2
      if (h17 == h28) ddxdusq=ddxdusq
     & +dble(dsxdc(h17,h28,h56)*dconjg(dsxdc2(h17,h28,h56)))*2._dp/xn


      enddo
      enddo
      enddo

c      fac=aveqq*V*(2._dp*gsq*esq*gwsq)**2
      fac=aveqq*V*(2._dp*gsq*esq**2/abs(zxw))**2

      dcxucsq=fac*dcxucsq
      dsxussq=fac*dsxussq
      dsxdcsq=fac*dsxdcsq
      usxucsq=fac*usxucsq

      dcxuc2sq=fac*dcxuc2sq
      dsxus2sq=fac*dsxus2sq
      dsxdc2sq=fac*dsxdc2sq
      usxuc2sq=fac*usxuc2sq

      duxuusq=fac*duxuusq
      udxuusq=fac*udxuusq
      ddxudsq=fac*ddxudsq
      ddxdusq=fac*ddxdusq

      return
      end
