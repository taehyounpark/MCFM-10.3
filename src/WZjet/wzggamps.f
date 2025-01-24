!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine wzggamps(j1,j2,j3,j4,j5,j6,j7,j8,dubAB,dubBA)
c     Amplitudes for the process
c     0 --> qbar(j1) + q(j2)+ve(j3)+e^+(p4)+mu-(p5)+mu+(p6)+g(p7)+g(p8)
c     only s34,s56,propagators included.

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'srdiags.f'
      integer,parameter::itype=7
      integer::jtype,j1,j2,j3,j4,j5,j6,j7,j8,
     & polg1,polg2,p3,p4,p5,p6
      complex(dp):: ab8m(2,2),ab8p(2,2),ba8m(2,2),ba8p(2,2),
     & nab8m(2,2),nab8p(2,2),nba8m(2,2),nba8p(2,2),
     & dubAB(itype,2,2,2),dubBA(itype,2,2,2),zab2,zab3,
     & AB4mmm,AB4mpm,AB4pmm,AB4ppm,AB5mmm,AB5mpm,AB5pmm,AB5ppm
      real(dp)::s3,s4
      integer,parameter::minus=1,pplus=2
c     statement functions
      s3(j1,j2,j3)=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s4(j1,j2,j3,j4)=s(j1,j2)+s(j1,j3)+s(j1,j4)
     &               +s(j2,j3)+s(j2,j4)+s(j3,j4)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=
     & za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)

c Statement functions for colour ordered amplitude for singly resonant
c diagrams

c d(-p1) +u~(-p2) --> mu^-(p3)+vmu~(p4)+e^-(p5)+e^+(p6)+g(p7)+g(p8)

c      5--|-<---6                                     5--|---<---6
c         |Z                                             |Z
c      3--|----<-----4                      3----<--|----|----4
c            |W-               |W-                  |W-
c   2 ----<--|-------------------1        2 ----<---|-----------1
c                 0    0                                 0   0
c                 0    0                                 0   0
c                 0    0                                 0   0
c               7 0    0 8                             7 0   0 8
c      jtype=4 (Zcoupling to electron)  jtype=5 (Z coupling to neutrino)

      AB4mmm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & (zab2(j8,j2,j7,j1)/(zb(j1,j7)*zb(j1,j8)*zb(j2,j7))
     & +za(j2,j8)/(zb(j1,j7)*zb(j7,j8))+za(j2,j7)/(zb(j1,j8)*zb(j7,j8)))
     & *za(j3,j5)*zb(j1,j4)
     & *(zb(j3,j6)*zab3(j3,j2,j7,j8,j1)+zb(j5,j6)*zab3(j5,j2,j7,j8,j1)))
     & /(s3(j2,j7,j8)*s(j5,j6)*s3(j3,j5,j6)*s4(j3,j4,j5,j6))

      AB4mpm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -za(j3,j5)*zab2(j7,j1,j8,j4)
     & *(zb(j3,j6)*zab2(j3,j2,j7,j8)+zb(j5,j6)*zab2(j5,j2,j7,j8))
     & /(za(j1,j8)*zb(j2,j7)*s(j7,j8))
     & +za(j1,j7)*zab2(j2,j3,j5,j6)*za(j3,j5)*zb(j1,j8)*zab2(j7,j1,j8,j4)
     & /(za(j1,j8)*s(j7,j8)*s3(j1,j7,j8))
     & +za(j3,j5)*zb(j1,j4)*zb(j2,j8)*za(j2,j7)
     & *(zb(j3,j6)*zab2(j3,j2,j7,j8)+zb(j5,j6)*zab2(j5,j2,j7,j8))
     & /(zb(j2,j7)*s(j7,j8)*s3(j2,j7,j8)))
     & /(s(j5,j6)*s3(j3,j5,j6)*s4(j3,j4,j5,j6))


      AB4pmm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -zab2(j2,j3,j5,j6)*za(j2,j8)*za(j3,j5)*zb(j1,j4)*zb(j1,j7)
     & /(za(j2,j7)*zb(j1,j8)*s(j7,j8))
     & +za(j2,j8)**2*za(j3,j5)*zb(j1,j4)
     & *(zb(j3,j6)*zab3(j3,j2,j7,j8,j7)+zb(j5,j6)*zab3(j5,j2,j7,j8,j7))
     & /(za(j2,j7)*s(j7,j8)*s3(j2,j7,j8))
     & +zab2(j2,j3,j5,j6)*za(j3,j5)*zb(j1,j7)**2*zab2(j8,j1,j7,j4)
     & /(zb(j1,j8)*s(j7,j8)*s3(j1,j7,j8)))
     & /(s(j5,j6)*s3(j3,j5,j6)*s4(j3,j4,j5,j6))


      AB4ppm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & +zab2(j2,j1,j8,j7)/za(j1,j8)/(za(j2,j7)*za(j2,j8))
     & -zb(j1,j8)/(za(j2,j7)*za(j7,j8))-zb(j1,j7)/(za(j2,j8)*za(j7,j8)))
     & *za(j3,j5)*zab2(j2,j3,j5,j6)*zab3(j2,j1,j7,j8,j4)
     & /(s3(j1,j7,j8)*s(j5,j6)*s3(j3,j5,j6)*s4(j3,j4,j5,j6))

      AB5mmm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -zab2(j8,j2,j7,j1)/(zb(j1,j7)*zb(j1,j8)*zb(j2,j7))
     & -za(j2,j8)/(zb(j1,j7)*zb(j7,j8))-za(j2,j7)/(zb(j1,j8)*zb(j7,j8)))
     & *zb(j4,j6)*zab2(j5,j4,j6,j1)*zab3(j3,j2,j7,j8,j1)
     & /(s3(j2,j7,j8)*s(j5,j6)*s3(j4,j5,j6)*s4(j3,j4,j5,j6))

      AB5mpm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -za(j1,j7)*za(j2,j3)*zb(j1,j8)*zb(j4,j6)
     & *(zab2(j7,j1,j8,j4)*za(j4,j5)+zab2(j7,j1,j8,j6)*za(j6,j5))
     & /(za(j1,j8)*s(j7,j8)*s3(j1,j7,j8))
     & +(zab2(j7,j1,j8,j4)*za(j4,j5)+zab2(j7,j1,j8,j6)*za(j6,j5))
     & *zb(j4,j6)*zab2(j3,j2,j7,j8)
     & /(za(j1,j8)*zb(j2,j7)*s(j7,j8))
     & -zab2(j5,j4,j6,j1)*za(j2,j7)*zb(j2,j8)*zb(j4,j6)
     & *zab2(j3,j2,j7,j8)/(zb(j2,j7)*s(j7,j8)*s3(j2,j7,j8)))
     & /(s(j5,j6)*s3(j4,j5,j6)*s4(j3,j4,j5,j6))

      AB5pmm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -za(j2,j8)**2*zab2(j5,j4,j6,j1)*zb(j4,j6)*zab3(j3,j2,j7,j8,j7)
     & /(za(j2,j7)*s(j7,j8)*s3(j2,j7,j8))
     & -za(j2,j3)*zb(j1,j7)**2*zb(j4,j6)
     & *(zab2(j8,j1,j7,j4)*za(j4,j5)+zab2(j8,j1,j7,j6)*za(j6,j5))
     & /(zb(j1,j8)*s(j7,j8)*s3(j1,j7,j8))
     & +za(j2,j3)*za(j2,j8)*zab2(j5,j4,j6,j1)*zb(j1,j7)*zb(j4,j6)
     & /(za(j2,j7)*zb(j1,j8)*s(j7,j8)))
     & /(s(j5,j6)*s3(j4,j5,j6)*s4(j3,j4,j5,j6))

      AB5ppm(j1,j2,j3,j4,j5,j6,j7,j8)=(
     & -zab2(j2,j1,j8,j7)/za(j1,j8)/(za(j2,j7)*za(j2,j8))
     & +zb(j1,j8)/(za(j2,j7)*za(j7,j8))+zb(j1,j7)/(za(j2,j8)*za(j7,j8)))
     & *za(j2,j3)*zb(j4,j6)
     & *(zab3(j2,j1,j7,j8,j4)*za(j4,j5)+zab3(j2,j1,j7,j8,j6)*za(j6,j5))
     & /(s3(j1,j7,j8)*s(j5,j6)*s3(j4,j5,j6)*s4(j3,j4,j5,j6))


c---Remember that the a8tree's are in DKS notation,
c-  so we swap (1<->2) and (5<-->6) to get to our standard notation

c     function a8treea as written in DKS has poles in s156 and s234
c     and corresponds to the amplitude <6 |5]  <1| |2] <3 |4]
c     (and dependence on gluon momenta fixed by polarizations)
c     for our basic process we swap 1<-->2,5<-->6
c     so we would have poles is s256, so after swap the amplitude
c     in our standard notation is,
c                                         5-----<--6  3-----<--4
c                                                \      /
c                                                 \    /
c                                                  \  /
c        5-----<-- 6   3-----<--4                   \/
c            |Z            |W                        |W
c   2 ----<--|-------------|----1      2 --<---------|-------------1
c                 0   0                       0   0
c              7  0   0 8                   7 0   0 8
c                 0   0                       0   0
c                 0   0                       0   0

c              jtype=2                         jtype=3

c        3-----<-- 4   5-----<--6
c            |W            |Z
c   2 ----<--|-------------|----1
c                 0   0
c                 0   0
c               7 0   0 8
c                 0   0

c                jtype=1

c     Now calculate both orderings of bosons
      do jtype=1,2
      if (jtype == 1) then
c swapped ordering
      p3=j5
      p4=j6
      p5=j3
      p6=j4
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,za,zb,ab8m)
      call a8treea(j2,j1,p4,p3,p6,p5,j7,j8,za,zb,ab8p)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,za,zb,ba8m)
      call a8treea(j2,j1,p4,p3,p6,p5,j8,j7,za,zb,ba8p)
      elseif (jtype == 2) then
c standard ordering
      p3=j3
      p4=j4
      p5=j5
      p6=j6
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,za,zb,ab8m)
      call a8treea(j2,j1,p3,p4,p5,p6,j7,j8,za,zb,ab8p)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,za,zb,ba8m)
      call a8treea(j2,j1,p3,p4,p5,p6,j8,j7,za,zb,ba8p)

      call a8treeb(j2,j1,p3,p4,p6,p5,j7,j8,za,zb,nab8m)
      call a8treeb(j2,j1,p3,p4,p5,p6,j7,j8,za,zb,nab8p)

      call a8treeb(j2,j1,p3,p4,p6,p5,j8,j7,za,zb,nba8m)
      call a8treeb(j2,j1,p3,p4,p5,p6,j8,j7,za,zb,nba8p)
      endif

      do polg1=1,2
      do polg2=1,2
      dubAB(jtype,polg1,polg2,minus)=+ab8m(polg1,polg2)
      dubAB(jtype,polg1,polg2,pplus)=+ab8p(polg1,polg2)
      dubBA(jtype,polg1,polg2,minus)=+ba8m(polg2,polg1)
      dubBA(jtype,polg1,polg2,pplus)=+ba8p(polg2,polg1)

      if (jtype == 2) then
      dubAB(3,polg1,polg2,minus)=-nab8m(polg1,polg2)
      dubAB(3,polg1,polg2,pplus)=-nab8p(polg1,polg2)
      dubBA(3,polg1,polg2,minus)=-nba8m(polg2,polg1)
      dubBA(3,polg1,polg2,pplus)=-nba8p(polg2,polg1)
      endif

      enddo !end polg2
      enddo !end polg1


      enddo ! end jtype loop

      if (srdiags .eqv. .false.) return

      dubAB(4,1,1,1)=AB4mmm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(4,2,1,1)=AB4pmm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(4,1,2,1)=AB4mpm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(4,2,2,1)=AB4ppm(j1,j2,j3,j4,j5,j6,j7,j8)

      dubAB(4,1,1,2)=AB4mmm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(4,2,1,2)=AB4pmm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(4,1,2,2)=AB4mpm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(4,2,2,2)=AB4ppm(j1,j2,j3,j4,j6,j5,j7,j8)

      dubAB(5,1,1,1)=AB5mmm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(5,2,1,1)=AB5pmm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(5,1,2,1)=AB5mpm(j1,j2,j3,j4,j5,j6,j7,j8)
      dubAB(5,2,2,1)=AB5ppm(j1,j2,j3,j4,j5,j6,j7,j8)

      dubAB(5,1,1,2)=AB5mmm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(5,2,1,2)=AB5pmm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(5,1,2,2)=AB5mpm(j1,j2,j3,j4,j6,j5,j7,j8)
      dubAB(5,2,2,2)=AB5ppm(j1,j2,j3,j4,j6,j5,j7,j8)


      dubBA(4,1,1,1)=AB4mmm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(4,2,1,1)=AB4mpm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(4,1,2,1)=AB4pmm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(4,2,2,1)=AB4ppm(j1,j2,j3,j4,j5,j6,j8,j7)

      dubBA(4,1,1,2)=AB4mmm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(4,2,1,2)=AB4mpm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(4,1,2,2)=AB4pmm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(4,2,2,2)=AB4ppm(j1,j2,j3,j4,j6,j5,j8,j7)

      dubBA(5,1,1,1)=AB5mmm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(5,2,1,1)=AB5mpm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(5,1,2,1)=AB5pmm(j1,j2,j3,j4,j5,j6,j8,j7)
      dubBA(5,2,2,1)=AB5ppm(j1,j2,j3,j4,j5,j6,j8,j7)

      dubBA(5,1,1,2)=AB5mmm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(5,2,1,2)=AB5mpm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(5,1,2,2)=AB5pmm(j1,j2,j3,j4,j6,j5,j8,j7)
      dubBA(5,2,2,2)=AB5ppm(j1,j2,j3,j4,j6,j5,j8,j7)

c d(-p1) +u~(-p2) --> mu^-(p3)+vmu~(p4)+e^-(p5)+e^+(p6)+g(p7)+g(p8)

c             3----<---4                    3----<---4
c                 |W-                           |W-
c        5-----<--|---6                     5-----<--|---6
c            |W-                                     |W-
c   2 ----<--|-------------------1        2 ----<----|-----------------1
c                 0      0                                 0      0
c                 0      0                                 0      0
c                 0      0                                 0      0
c               7 0      0 8                             7 0      0 8
c             jtype=6                                 jtype=7


      dubAB(6,1,1,1)=AB5mmm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(6,2,1,1)=AB5pmm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(6,1,2,1)=AB5mpm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(6,2,2,1)=AB5ppm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(6,:,:,2)=czip

      dubAB(7,1,1,1)=AB4mmm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(7,2,1,1)=AB4pmm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(7,1,2,1)=AB4mpm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(7,2,2,1)=AB4ppm(j1,j2,j5,j6,j3,j4,j7,j8)
      dubAB(7,:,:,2)=czip

      dubBA(6,1,1,1)=AB5mmm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(6,2,1,1)=AB5mpm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(6,1,2,1)=AB5pmm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(6,2,2,1)=AB5ppm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(6,:,:,2)=czip

      dubBA(7,1,1,1)=AB4mmm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(7,2,1,1)=AB4mpm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(7,1,2,1)=AB4pmm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(7,2,2,1)=AB4ppm(j1,j2,j5,j6,j3,j4,j8,j7)
      dubBA(7,:,:,2)=czip


      return
      end

