!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a8ZZfourqsq(j1,j2,j3,j4,j5,j6,j7,j8,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_du2,ampsq_ds2,ampsq_ud2,ampsq_uc2)
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
c     Performs the square and the sum over polarization of the diagrams
c     which have the ZZ radiated off the external line.
      integer j1,j2,j3,j4,j5,j6,j7,j8
      real(dp)::ampsq_uc,ampsq_uu,ampsq_ds,ampsq_dd,ampsq_ud,ampsq_du,
     & ampsq_uc2,ampsq_ds2,ampsq_ud2,ampsq_du2,s3456,fac
      complex(dp)::Zq(2,2),Zl34(2),Zl56(2)
      complex(dp):: prop34,prop56,prop3456,
     & coup17(2,2,2,2),coup3456(2,2,2,2),coup5634(2,2,2,2),
     & coup17x28(2,2,2,2,2,2),coup28x17(2,2,2,2,2,2),
     & upper1(2,2,2,2),upper2(2,2,2,2),lower1(2,2,2,2),lower2(2,2,2,2),
     & updown1(2,2,2,2),updown2(2,2,2,2),
     & downup1(2,2,2,2),downup2(2,2,2,2),
     & uc(2,2,2,2),uc2(2,2,2,2),ds(2,2,2,2),ds2(2,2,2,2),
     & ud(2,2,2,2),ud2(2,2,2,2),du(2,2,2,2),du2(2,2,2,2),
     & u3456(2,2,2,2),l3456(2,2,2,2),u5634(2,2,2,2),l5634(2,2,2,2),
     & v3456(2,2,2,2),m3456(2,2,2,2),v5634(2,2,2,2),m5634(2,2,2,2)
       integer::h17,h28,h34,h56,iq1,iq2,h,i

       s3456=s(j3,j4)+s(j3,j5)+s(j3,j6)+s(j4,j5)+s(j4,j6)+s(j5,j6)
       prop3456=cmplx(s3456,kind=dp)
     &         /cmplx(s3456-zmass**2,zmass*zwidth,kind=dp)
       prop34 = cmplx(s(j3,j4),kind=dp)
     &         /cmplx(s(j3,j4)-zmass**2,zmass*zwidth,kind=dp)
       prop56 = cmplx(s(j5,j6),kind=dp)
     &         /cmplx(s(j5,j6)-zmass**2,zmass*zwidth,kind=dp)

c--- initialize couplings of Z/gamma to quarks and lepton pairs.
c      srdiags=.true.
c      write(6,*) 'srdiags',srdiags

      do h=1,2
      Zl34(h)=(2-h)*zl1+(h-1)*zr1
      Zl56(h)=(2-h)*zl2+(h-1)*zr2
      do i=1,2
      Zq(i,h)=(2-h)*zL(i)+(h-1)*zR(i)
      enddo
      enddo

      do h34=1,2
      do h56=1,2

      do iq1=1,2
      do h17=1,2
      do iq2=1,2
      do h28=1,2
      coup17x28(iq1,h17,iq2,h28,h34,h56)=
     &  (Q(iq1)*q1+Zq(iq1,h17)*Zl34(h34)*prop34)
     & *(Q(iq2)*q2+Zq(iq2,h28)*Zl56(h56)*prop56)
      enddo
      enddo
      enddo
      enddo

      do iq1=1,2
      do h17=1,2
      do iq2=1,2
      do h28=1,2
      coup28x17(iq1,h17,iq2,h28,h34,h56)=
     &  (Q(iq2)*q1+Zq(iq2,h28)*Zl34(h34)*prop34)
     & *(Q(iq1)*q2+Zq(iq1,h17)*Zl56(h56)*prop56)
      enddo
      enddo
      enddo
      enddo

      do i=1,2
      do h=1,2
      coup17(i,h,h34,h56)=
     &  (Q(i)*q1+Zq(i,h)*Zl34(h34)*prop34)
     & *(Q(i)*q2+Zq(i,h)*Zl56(h56)*prop56)
      coup3456(i,h,h34,h56)=
     &  (Q(i)*q1+Zq(i,h)*Zl34(h34)*prop3456)
     & *(q1*q2+Zl34(h34)*Zl56(h56)*prop56)
      coup5634(i,h,h34,h56)=
     &  (Q(i)*q2+Zq(i,h)*Zl56(h56)*prop3456)
     & *(q2*q1+Zl56(h56)*Zl34(h34)*prop34)
      enddo
      enddo

      enddo
      enddo

cdebug
c      srdiags=.true.
cdebug

c     amplitude a(h17,h28,h34,h56)
c--- call all basic amplitudes
      call a8ZZ(j1,j2,j3,j4,j5,j6,j7,j8,upper1)
      call a8ZZ(j2,j1,j3,j4,j5,j6,j8,j7,lower1)

      call a8ZZ(j1,j2,j3,j4,j5,j6,j8,j7,upper2)
      call a8ZZ(j2,j1,j3,j4,j5,j6,j7,j8,lower2)

      call a8ZxZ(j1,j2,j3,j4,j5,j6,j7,j8,updown1)
      call a8ZxZ(j1,j2,j5,j6,j3,j4,j7,j8,downup1)

      call a8ZxZ(j1,j2,j3,j4,j5,j6,j8,j7,updown2)
      call a8ZxZ(j1,j2,j5,j6,j3,j4,j8,j7,downup2)

      if (srdiags) then
      call a8sr4q(j1,j2,j3,j4,j5,j6,j7,j8,u3456,l3456,u5634,l5634)
      call a8sr4q(j1,j2,j3,j4,j5,j6,j8,j7,v3456,m3456,v5634,m5634)
      endif

c      write(6,*) 'sriags',srdiags
      do h17=1,2
      do h28=1,2
      do h34=1,2
      do h56=1,2
      uc(h17,h28,h34,h56)=
     & +coup17(2,h17,h34,h56)*upper1(h17,h28,h34,h56)
     & +coup17(2,h28,h34,h56)*lower1(h28,h17,h34,h56)
     & +coup17x28(2,h17,2,h28,h34,h56)*updown1(h17,h28,h34,h56)
     & +coup28x17(2,h17,2,h28,h34,h56)*downup1(h17,h28,h56,h34)
      if (srdiags) then
      uc(h17,h28,h34,h56)=uc(h17,h28,h34,h56)
     & +coup3456(2,h17,h34,h56)*u3456(h17,h28,h34,h56)
     & +coup5634(2,h17,h34,h56)*u5634(h17,h28,h34,h56)
     & +coup3456(2,h28,h34,h56)*l3456(h17,h28,h34,h56)
     & +coup5634(2,h28,h34,h56)*l5634(h17,h28,h34,h56)
      endif
      uc2(h17,h28,h34,h56)=
     & +coup17(2,h17,h34,h56)*upper2(h17,h28,h34,h56)
     & +coup17(2,h28,h34,h56)*lower2(h28,h17,h34,h56)
     & +coup17x28(2,h17,2,h28,h34,h56)*updown2(h17,h28,h34,h56)
     & +coup28x17(2,h17,2,h28,h34,h56)*downup2(h17,h28,h56,h34)
      if (srdiags) then
      uc2(h17,h28,h34,h56)=uc2(h17,h28,h34,h56)
     & +coup3456(2,h17,h34,h56)*v3456(h17,h28,h34,h56)
     & +coup5634(2,h17,h34,h56)*v5634(h17,h28,h34,h56)
     & +coup3456(2,h28,h34,h56)*m3456(h17,h28,h34,h56)
     & +coup5634(2,h28,h34,h56)*m5634(h17,h28,h34,h56)
      endif

      ds(h17,h28,h34,h56)=
     & +coup17(1,h17,h34,h56)*upper1(h17,h28,h34,h56)
     & +coup17(1,h28,h34,h56)*lower1(h28,h17,h34,h56)
     & +coup17x28(1,h17,1,h28,h34,h56)*updown1(h17,h28,h34,h56)
     & +coup28x17(1,h17,1,h28,h34,h56)*downup1(h17,h28,h56,h34)
      if (srdiags) then
      ds(h17,h28,h34,h56)=ds(h17,h28,h34,h56)
     & +coup3456(1,h17,h34,h56)*u3456(h17,h28,h34,h56)
     & +coup5634(1,h17,h34,h56)*u5634(h17,h28,h34,h56)
     & +coup3456(1,h28,h34,h56)*l3456(h17,h28,h34,h56)
     & +coup5634(1,h28,h34,h56)*l5634(h17,h28,h34,h56)
      endif
      ds2(h17,h28,h34,h56)=
     & +coup17(1,h17,h34,h56)*upper2(h17,h28,h34,h56)
     & +coup17(1,h28,h34,h56)*lower2(h28,h17,h34,h56)
     & +coup17x28(1,h17,1,h28,h34,h56)*updown2(h17,h28,h34,h56)
     & +coup28x17(1,h17,1,h28,h34,h56)*downup2(h17,h28,h56,h34)
      if (srdiags) then
      ds2(h17,h28,h34,h56)=ds2(h17,h28,h34,h56)
     & +coup3456(1,h17,h34,h56)*v3456(h17,h28,h34,h56)
     & +coup5634(1,h17,h34,h56)*v5634(h17,h28,h34,h56)
     & +coup3456(1,h28,h34,h56)*m3456(h17,h28,h34,h56)
     & +coup5634(1,h28,h34,h56)*m5634(h17,h28,h34,h56)
      endif

      ud(h17,h28,h34,h56)=
     & +coup17(2,h17,h34,h56)*upper1(h17,h28,h34,h56)
     & +coup17(1,h28,h34,h56)*lower1(h28,h17,h34,h56)
     & +coup17x28(2,h17,1,h28,h34,h56)*updown1(h17,h28,h34,h56)
     & +coup28x17(2,h17,1,h28,h34,h56)*downup1(h17,h28,h56,h34)
      if (srdiags) then
      ud(h17,h28,h34,h56)=ud(h17,h28,h34,h56)
     & +coup3456(2,h17,h34,h56)*u3456(h17,h28,h34,h56)
     & +coup5634(2,h17,h34,h56)*u5634(h17,h28,h34,h56)
     & +coup3456(1,h28,h34,h56)*l3456(h17,h28,h34,h56)
     & +coup5634(1,h28,h34,h56)*l5634(h17,h28,h34,h56)
      endif
      ud2(h17,h28,h34,h56)=
     & +coup17(2,h17,h34,h56)*upper2(h17,h28,h34,h56)
     & +coup17(1,h28,h34,h56)*lower2(h28,h17,h34,h56)
     & +coup17x28(2,h17,1,h28,h34,h56)*updown2(h17,h28,h34,h56)
     & +coup28x17(2,h17,1,h28,h34,h56)*downup2(h17,h28,h56,h34)
      if (srdiags) then
      ud2(h17,h28,h34,h56)=ud2(h17,h28,h34,h56)
     & +coup3456(2,h17,h34,h56)*v3456(h17,h28,h34,h56)
     & +coup5634(2,h17,h34,h56)*v5634(h17,h28,h34,h56)
     & +coup3456(1,h28,h34,h56)*m3456(h17,h28,h34,h56)
     & +coup5634(1,h28,h34,h56)*m5634(h17,h28,h34,h56)
      endif

      du(h17,h28,h34,h56)=
     & +coup17(1,h17,h34,h56)*upper1(h17,h28,h34,h56)
     & +coup17(2,h28,h34,h56)*lower1(h28,h17,h34,h56)
     & +coup17x28(1,h17,2,h28,h34,h56)*updown1(h17,h28,h34,h56)
     & +coup28x17(1,h17,2,h28,h34,h56)*downup1(h17,h28,h56,h34)
      if (srdiags) then
      du(h17,h28,h34,h56)=du(h17,h28,h34,h56)
     & +coup3456(1,h17,h34,h56)*u3456(h17,h28,h34,h56)
     & +coup5634(1,h17,h34,h56)*u5634(h17,h28,h34,h56)
     & +coup3456(2,h28,h34,h56)*l3456(h17,h28,h34,h56)
     & +coup5634(2,h28,h34,h56)*l5634(h17,h28,h34,h56)
      endif
      du2(h17,h28,h34,h56)=
     & +coup17(1,h17,h34,h56)*upper2(h17,h28,h34,h56)
     & +coup17(2,h28,h34,h56)*lower2(h28,h17,h34,h56)
     & +coup17x28(1,h17,2,h28,h34,h56)*updown2(h17,h28,h34,h56)
     & +coup28x17(1,h17,2,h28,h34,h56)*downup2(h17,h28,h56,h34)
      if (srdiags) then
      du2(h17,h28,h34,h56)=du2(h17,h28,h34,h56)
     & +coup3456(1,h17,h34,h56)*v3456(h17,h28,h34,h56)
     & +coup5634(1,h17,h34,h56)*v5634(h17,h28,h34,h56)
     & +coup3456(2,h28,h34,h56)*m3456(h17,h28,h34,h56)
     & +coup5634(2,h28,h34,h56)*m5634(h17,h28,h34,h56)
      endif
      enddo
      enddo
      enddo
      enddo

      ampsq_uc=0._dp
      ampsq_ud=0._dp
      ampsq_du=0._dp
      ampsq_ds=0._dp
      ampsq_uu=0._dp
      ampsq_dd=0._dp
      ampsq_uc2=0._dp
      ampsq_ds2=0._dp
      ampsq_ud2=0._dp
      ampsq_du2=0._dp

      do h17=1,2
      do h28=1,2
      do h34=1,2
      do h56=1,2
      ampsq_uc=ampsq_uc+abs(uc(h17,h28,h34,h56))**2
      ampsq_ds=ampsq_ds+abs(ds(h17,h28,h34,h56))**2
      ampsq_ud=ampsq_ud+abs(ud(h17,h28,h34,h56))**2
      ampsq_du=ampsq_du+abs(du(h17,h28,h34,h56))**2
      ampsq_uc2=ampsq_uc2+abs(uc2(h17,h28,h34,h56))**2
      ampsq_ds2=ampsq_ds2+abs(ds2(h17,h28,h34,h56))**2
      ampsq_ud2=ampsq_ud2+abs(ud2(h17,h28,h34,h56))**2
      ampsq_du2=ampsq_du2+abs(du2(h17,h28,h34,h56))**2
      ampsq_uu=ampsq_uu
     & +abs(uc(h17,h28,h34,h56))**2+abs(uc2(h17,h28,h34,h56))**2
      ampsq_dd=ampsq_dd
     & +abs(ds(h17,h28,h34,h56))**2+abs(ds2(h17,h28,h34,h56))**2
      if (h17 == h28) then
      ampsq_uu=ampsq_uu
     & +(dble(uc(h17,h28,h34,h56)*dconjg(uc2(h17,h28,h34,h56)))
     &  +dble(dconjg(uc(h17,h28,h34,h56))*uc2(h17,h28,h34,h56)))/xn
      ampsq_dd=ampsq_dd
     & +(dble(ds(h17,h28,h34,h56)*dconjg(ds2(h17,h28,h34,h56)))
     &  +dble(dconjg(ds(h17,h28,h34,h56))*ds2(h17,h28,h34,h56)))/xn
      endif
      enddo
      enddo
      enddo
      enddo
      fac=aveqq*V*(4._dp*gsq*esq**2)**2
      ampsq_uc=fac*ampsq_uc
      ampsq_du=fac*ampsq_du
      ampsq_ud=fac*ampsq_ud
      ampsq_ds=fac*ampsq_ds
      ampsq_uu=fac*ampsq_uu
      ampsq_dd=fac*ampsq_dd
      ampsq_uc2=fac*ampsq_uc2
      ampsq_ds2=fac*ampsq_ds2
      ampsq_ud2=fac*ampsq_ud2
      ampsq_du2=fac*ampsq_du2
c      write(6,*) 'ampsq_ds',ampsq_ds
c      write(6,*) 'ampsq_uc',ampsq_uc
c      write(6,*) 'ampsq_ud',ampsq_ud
c      write(6,*) 'ampsq_du',ampsq_du
c      write(6,*) 'ampsq_uu',ampsq_uu
c      write(6,*) 'ampsq_dd',ampsq_dd
c      write(6,*) 'ampsq_uc2',ampsq_uc2
c      write(6,*) 'ampsq_ds2',ampsq_ds2
c      write(6,*) 'ampsq_ud2',ampsq_ud2
c      write(6,*) 'ampsq_du2',ampsq_du2

      return
      end
