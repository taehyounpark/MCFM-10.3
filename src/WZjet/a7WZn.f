!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a7WZn(p1,p2,p3,p4,p5,p6,p7,p,n,a7n)
      implicit none
c---- Matrix element for WZ radiation from line 17
c---- including W Z interchange
c  d(-p1) +c(-p2) --> e^-(p3)+ve~^+(p4)+mu^-(p5)+mu^+(p6)+u(p7)+c(p8)

c                                         5-----<-- 6 3-----<--4
c                                                \      /
c                                               Z \    / W
c                                                  \  /
c        5-----<-- 6   3-----<--4                   \/
c            |Z            |W                        |W
c   2 ----<--|-------------|----1      2 ----<-------|-------------1
c                 0                               0
c                 0                               0
c                 0                               0
c              jtype=1                         jtype=3

c        3-----<-- 4   5-----<--6
c            |W            |Z
c   2 ----<--|-------------|----1
c                 0
c                 0
c                 0
c                jtype=2


c     with line 7 contracted with n
c     jtype=1,2 for a-type diagrams
c     jtype=3 for b-type diagrams
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'srdiags.f'
      integer i,j,p1,p2,p3,p4,p5,p6,p7
      integer,parameter:: jtype=7
      real(dp):: s3,n(4),p(mxpart,4),s34,s56,s17,s27,
     & s127,s134,s156,s234,s256
      complex(dp):: a7n(jtype,2),vecm(mxpart,mxpart),zab2,
     & aj4m,aj4p,aj5m,aj5p
c     statement functions
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      aj4m(p1,p2,p3,p4,p5,p6,p7)=
     & (za(p3,p5)*zb(p1,p4)/(za(p2,p7)*zb(p2,p7))
     & *(zab2(p2,p3,p5,p6)*vecm(p2,p2)+zab2(p7,p3,p5,p6)*vecm(p2,p7))
     & -(zb(p1,p4)*vecm(p1,p1)+zb(p7,p4)*vecm(p7,p1))
     & *za(p3,p5)*zab2(p2,p3,p5,p6)/(za(p1,p7)*zb(p1,p7)))
     & /(s(p5,p6)*s3(p3,p5,p6)*s3(p1,p2,p7))
      aj4p(p1,p2,p3,p4,p5,p6,p7)=
     & (za(p3,p6)*zb(p1,p4)/za(p2,p7)/zb(p2,p7)
     & *(zab2(p2,p3,p6,p5)*vecm(p2,p2)+zab2(p7,p3,p6,p5)*vecm(p2,p7))
     & -(zb(p1,p4)*vecm(p1,p1)+zb(p7,p4)*vecm(p7,p1))
     & *za(p3,p6)/(za(p1,p7)*zb(p1,p7))*zab2(p2,p3,p6,p5))
     & /(s(p5,p6)*s3(p3,p5,p6)*s3(p1,p2,p7))

      aj5m(p1,p2,p3,p4,p5,p6,p7)=
     & ((za(p3,p2)*vecm(p2,p2)+za(p3,p7)*vecm(p2,p7))
     & *zb(p4,p6)/(za(p2,p7)*zb(p2,p7))*zab2(p5,p4,p6,p1)
     & -(zab2(p5,p4,p6,p1)*vecm(p1,p1)+zab2(p5,p4,p6,p7)*vecm(p7,p1))
     & *za(p3,p2)*zb(p4,p6)/(za(p1,p7)*zb(p1,p7)))
     & /(s(p5,p6)*s3(p4,p5,p6)*s3(p1,p2,p7))
      aj5p(p1,p2,p3,p4,p5,p6,p7)=
     & ((za(p3,p2)*vecm(p2,p2)+za(p3,p7)*vecm(p2,p7))
     & *zb(p4,p5)/(za(p2,p7)*zb(p2,p7))*zab2(p6,p4,p5,p1)
     & -za(p3,p2)*zb(p4,p5)/(za(p1,p7)*zb(p1,p7))
     & *(zab2(p6,p4,p5,p1)*vecm(p1,p1)+zab2(p6,p4,p5,p7)*vecm(p7,p1)))
     & /(s(p5,p6)*s3(p4,p5,p6)*s3(p1,p2,p7))
c     endstatement functions

      call checkndotp(p,n,p7)
      s17=s(p1,p7)
      s27=s(p2,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s127=s3(p1,p2,p7)
      s134=s3(p1,p3,p4)
      s156=s3(p1,p5,p6)
      s234=s3(p2,p3,p4)
      s256=s3(p2,p5,p6)

      do i=1,7
      do j=i,7
      call ndveccur(i,j,n,p,vecm)
      enddo
      enddo

      a7n(:,:)=czip

      a7n(1,1)= + s17**(-1) * (  - za(p2,p3)*zb(p1,p6)*zab2(p5,p2,p3,p4
     &    )*vecm(p1,p1)*s234**(-1) + za(p2,p3)*zb(p6,p7)*zab2(p5,p2,p3,
     &    p4)*vecm(p7,p1)*s234**(-1) )
      a7n(1,1) = a7n(1,1) + s27**(-1) * (  - za(p2,p3)*zb(p1,p6)*zab2(
     &    p5,p1,p6,p4)*vecm(p2,p2)*s156**(-1) + za(p3,p7)*zb(p1,p6)*
     &    zab2(p5,p1,p6,p4)*vecm(p2,p7)*s156**(-1) )
      a7n(1,1) = a7n(1,1) + za(p1,p5)*za(p2,p3)*zb(p1,p6)*zb(p2,p4)*
     & vecm(p2,p1)*s156**(-1)*s234**(-1) + za(p1,p5)*za(p2,p3)*zb(p1,p6
     &    )*zb(p3,p4)*vecm(p3,p1)*s156**(-1)*s234**(-1) - za(p2,p3)*za(
     &    p5,p6)*zb(p1,p6)*zb(p2,p4)*vecm(p2,p6)*s156**(-1)*s234**(-1)
     &     - za(p2,p3)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*vecm(p3,p6)*
     &    s156**(-1)*s234**(-1)
      a7n(1,2)= + s17**(-1) * (  - za(p2,p3)*zb(p1,p5)*zab2(p6,p2,p3,p4
     &    )*vecm(p1,p1)*s234**(-1) + za(p2,p3)*zb(p5,p7)*zab2(p6,p2,p3,
     &    p4)*vecm(p7,p1)*s234**(-1) )
      a7n(1,2) = a7n(1,2) + s27**(-1) * (  - za(p2,p3)*zb(p1,p5)*zab2(
     &    p6,p1,p5,p4)*vecm(p2,p2)*s156**(-1) + za(p3,p7)*zb(p1,p5)*
     &    zab2(p6,p1,p5,p4)*vecm(p2,p7)*s156**(-1) )
      a7n(1,2) = a7n(1,2) + za(p1,p6)*za(p2,p3)*zb(p1,p5)*zb(p2,p4)*
     & vecm(p2,p1)*s156**(-1)*s234**(-1) + za(p1,p6)*za(p2,p3)*zb(p1,p5
     &    )*zb(p3,p4)*vecm(p3,p1)*s156**(-1)*s234**(-1) + za(p2,p3)*za(
     &    p5,p6)*zb(p1,p5)*zb(p2,p4)*vecm(p2,p5)*s156**(-1)*s234**(-1)
     &     + za(p2,p3)*za(p5,p6)*zb(p1,p5)*zb(p3,p4)*vecm(p3,p5)*
     &    s156**(-1)*s234**(-1)

      a7n(2,1)= + s17**(-1) * (  - za(p2,p5)*zb(p1,p4)*zab2(p3,p2,p5,p6
     &    )*vecm(p1,p1)*s256**(-1) + za(p2,p5)*zb(p4,p7)*zab2(p3,p2,p5,
     &    p6)*vecm(p7,p1)*s256**(-1) )
      a7n(2,1) = a7n(2,1) + s27**(-1) * (  - za(p2,p5)*zb(p1,p4)*zab2(
     &    p3,p1,p4,p6)*vecm(p2,p2)*s134**(-1) + za(p5,p7)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p6)*vecm(p2,p7)*s134**(-1) )
      a7n(2,1) = a7n(2,1) + za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p2,p6)*
     & vecm(p2,p1)*s256**(-1)*s134**(-1) + za(p1,p3)*za(p2,p5)*zb(p1,p4
     &    )*zb(p5,p6)*vecm(p5,p1)*s256**(-1)*s134**(-1) - za(p2,p5)*za(
     &    p3,p4)*zb(p1,p4)*zb(p2,p6)*vecm(p2,p4)*s256**(-1)*s134**(-1)
     &     - za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p5,p6)*vecm(p5,p4)*
     &    s256**(-1)*s134**(-1)
      a7n(2,2)= + s17**(-1) * (  - za(p2,p6)*zb(p1,p4)*zab2(p3,p2,p6,p5
     &    )*vecm(p1,p1)*s256**(-1) + za(p2,p6)*zb(p4,p7)*zab2(p3,p2,p6,
     &    p5)*vecm(p7,p1)*s256**(-1) )
      a7n(2,2) = a7n(2,2) + s27**(-1) * (  - za(p2,p6)*zb(p1,p4)*zab2(
     &    p3,p1,p4,p5)*vecm(p2,p2)*s134**(-1) + za(p6,p7)*zb(p1,p4)*
     &    zab2(p3,p1,p4,p5)*vecm(p2,p7)*s134**(-1) )
      a7n(2,2) = a7n(2,2) + za(p1,p3)*za(p2,p6)*zb(p1,p4)*zb(p2,p5)*
     & vecm(p2,p1)*s256**(-1)*s134**(-1) - za(p1,p3)*za(p2,p6)*zb(p1,p4
     &    )*zb(p5,p6)*vecm(p6,p1)*s256**(-1)*s134**(-1) - za(p2,p6)*za(
     &    p3,p4)*zb(p1,p4)*zb(p2,p5)*vecm(p2,p4)*s256**(-1)*s134**(-1)
     &     + za(p2,p6)*za(p3,p4)*zb(p1,p4)*zb(p5,p6)*vecm(p6,p4)*
     &    s256**(-1)*s134**(-1)

      a7n(3,1)= + s17**(-1) * ( za(p2,p3)*zb(p1,p4)*zab2(p5,p3,p4,p6)*
     &    vecm(p1,p1) - za(p2,p3)*zb(p4,p7)*zab2(p5,p3,p4,p6)*vecm(p7,
     &    p1) - za(p2,p5)*zb(p1,p6)*zab2(p3,p5,p6,p4)*vecm(p1,p1) + za(
     &    p2,p5)*zb(p6,p7)*zab2(p3,p5,p6,p4)*vecm(p7,p1) - za(p3,p5)*
     &    zb(p4,p6)*zab2(p2,p3,p4,p1)*vecm(p1,p1) - za(p3,p5)*zb(p4,p6)
     &    *zab2(p2,p3,p4,p7)*vecm(p7,p1) )
      a7n(3,1) = a7n(3,1) + s27**(-1) * (  - za(p2,p3)*zb(p1,p4)*zab2(
     &    p5,p3,p4,p6)*vecm(p2,p2) + za(p2,p5)*zb(p1,p6)*zab2(p3,p5,p6,
     &    p4)*vecm(p2,p2) + za(p3,p5)*zb(p4,p6)*zab2(p2,p3,p4,p1)*vecm(
     &    p2,p2) + za(p3,p5)*zb(p4,p6)*zab2(p7,p3,p4,p1)*vecm(p2,p7) +
     &    za(p3,p7)*zb(p1,p4)*zab2(p5,p3,p4,p6)*vecm(p2,p7) - za(p5,p7)
     &    *zb(p1,p6)*zab2(p3,p5,p6,p4)*vecm(p2,p7) )
      a7n(3,2)= + s17**(-1) * ( za(p2,p3)*zb(p1,p4)*zab2(p6,p3,p4,p5)*
     &    vecm(p1,p1) - za(p2,p3)*zb(p4,p7)*zab2(p6,p3,p4,p5)*vecm(p7,
     &    p1) - za(p2,p6)*zb(p1,p5)*zab2(p3,p5,p6,p4)*vecm(p1,p1) + za(
     &    p2,p6)*zb(p5,p7)*zab2(p3,p5,p6,p4)*vecm(p7,p1) - za(p3,p6)*
     &    zb(p4,p5)*zab2(p2,p3,p4,p1)*vecm(p1,p1) - za(p3,p6)*zb(p4,p5)
     &    *zab2(p2,p3,p4,p7)*vecm(p7,p1) )
      a7n(3,2) = a7n(3,2) + s27**(-1) * (  - za(p2,p3)*zb(p1,p4)*zab2(
     &    p6,p3,p4,p5)*vecm(p2,p2) + za(p2,p6)*zb(p1,p5)*zab2(p3,p5,p6,
     &    p4)*vecm(p2,p2) + za(p3,p6)*zb(p4,p5)*zab2(p2,p3,p4,p1)*vecm(
     &    p2,p2) + za(p3,p6)*zb(p4,p5)*zab2(p7,p3,p4,p1)*vecm(p2,p7) +
     &    za(p3,p7)*zb(p1,p4)*zab2(p6,p3,p4,p5)*vecm(p2,p7) - za(p6,p7)
     &    *zb(p1,p5)*zab2(p3,p5,p6,p4)*vecm(p2,p7) )

      a7n(1,:)=a7n(1,:)/(s34*s56)
      a7n(2,:)=a7n(2,:)/(s34*s56)
      a7n(3,:)=a7n(3,:)/(s34*s56*s127)

      if (srdiags .eqv. .false.) return
c     singly resonant diagrams

c     d(-p1) +u~(-p2) --> mu^-(p3)+vmu~(p4)+e^-(p5)+e^+(p6)+g(p7)

c      5--|-<---6                                     5--|---<---6
c         |Z                                             |Z
c      3--|----<-----4                      3----<--|----|----4
c            |W-               |W-                  |W-
c   2 ----<--|-------------------1        2 ----<---|-----------1
c                 0                                      0
c                 0                                      0
c                 0                                      0
c               7 0                                    7 0
c      jtype=4 (Zcoupling to electron)     jtype=5 (Z coupling to neutrino)

      a7n(4,1)=aj4m(p1,p2,p3,p4,p5,p6,p7)
      a7n(4,2)=aj4p(p1,p2,p3,p4,p5,p6,p7)
      a7n(5,1)=aj5m(p1,p2,p3,p4,p5,p6,p7)
      a7n(5,2)=aj5p(p1,p2,p3,p4,p5,p6,p7)

c  d(-p1) +u~(-p2) --> mu^-(p3)+vmu~(p4)+e^-(p5)+e^+(p6)+g(p7)

c             3----<---4
c                 |W-
c        5-----<--|---6
c            |W-
c   2 ----<--|-------------------1
c                 0
c                 0
c                 0
c               7 0
c  Overall factor of  2*i_*gs*e^2*gw^2*(T^C7*T^C8)_{i2,i1}/2
c             jtype=6

      a7n(6,1)=aj5m(p1,p2,p5,p6,p3,p4,p7)
      a7n(7,1)=aj4m(p1,p2,p5,p6,p3,p4,p7)
c     No contribution for RH 56 line (2 W's)
      a7n(6,2)=czip
      a7n(7,2)=czip

      return
      end


