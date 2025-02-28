!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c***************************************************************
c   Color ordered amplitudes for:
c   0 -> q(-p1) + qb(-p4) + a(p2) + a(p3) + lb(p5) + l(p6)
c***************************************************************
      subroutine zaa_a60h(j1,j2,j3,j4,j5,j6,za,zb,a60h)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'ipsgen.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treega,a6treeQLslc
      complex(dp):: a60h(4,16)
      integer:: st1,st2,st3,st4,st5
      integer:: i,j

c---- helicity stamps
c     'q+qb-g+g+'=1
c     'q+qb-g+g-'=2
c     'q+qb-g-g+'=3
c     'q+qb-g+g+lb-l+'=4
c     'q+qb-g+g-lb-l+'=5

c---- helicity stamps
      st1=1
      st2=2
      st3=3
      st4=4
      st5=5

c-----initialize a60h
      do i=1,4
      do j=1,8
      a60h(i,j)=czip
      enddo
      enddo

c-----QQ contributions
      if (ipsgen  /=  2) then
      a60h(1,1) = a6treega(st1,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(1,2) = a6treega(st2,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(1,3) = a6treega(st3,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(1,4) = a6treega(st1,j4,j1,j2,j3,j6,j5,zb,za)
      a60h(1,5) = a6treega(st1,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(1,6) = a6treega(st2,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(1,7) = a6treega(st3,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(1,8) = a6treega(st1,j4,j1,j2,j3,j5,j6,zb,za)
      a60h(1,9) = a6treega(st1,j1,j4,j2,j3,j5,j6,zb,za)
      a60h(1,10) = a6treega(st3,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(1,11) = a6treega(st2,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(1,12) = a6treega(st1,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(1,13) = a6treega(st1,j1,j4,j2,j3,j6,j5,zb,za)
      a60h(1,14) = a6treega(st3,j4,j1,j2,j3,j5,j6,za,zb)
      a60h(1,15) = a6treega(st2,j4,j1,j2,j3,j5,j6,za,zb)
      a60h(1,16) = a6treega(st1,j4,j1,j2,j3,j5,j6,za,zb)
      endif
c-----LL contributions
      if (ipsgen  /=  3) then
      a60h(2,1) = a6treega(st1,j6,j5,j2,j3,j4,j1,za,zb)
      a60h(2,2) = a6treega(st2,j6,j5,j2,j3,j4,j1,za,zb)
      a60h(2,3) = a6treega(st3,j6,j5,j2,j3,j4,j1,za,zb)
      a60h(2,4) = a6treega(st1,j5,j6,j2,j3,j1,j4,zb,za)
      a60h(2,16) = a6treega(st1,j6,j5,j2,j3,j1,j4,za,zb)
      a60h(2,15) = a6treega(st2,j6,j5,j2,j3,j1,j4,za,zb)
      a60h(2,14) = a6treega(st3,j6,j5,j2,j3,j1,j4,za,zb)
      a60h(2,13) = a6treega(st1,j5,j6,j2,j3,j4,j1,zb,za)
      a60h(2,9) = a6treega(st1,j6,j5,j2,j3,j4,j1,zb,za)
      a60h(2,10) = a6treega(st3,j5,j6,j2,j3,j1,j4,za,zb)
      a60h(2,11) = a6treega(st2,j5,j6,j2,j3,j1,j4,za,zb)
      a60h(2,12) = a6treega(st1,j5,j6,j2,j3,j1,j4,za,zb)
      a60h(2,8) = a6treega(st1,j6,j5,j2,j3,j1,j4,zb,za)
      a60h(2,7) = a6treega(st3,j5,j6,j2,j3,j4,j1,za,zb)
      a60h(2,6) = a6treega(st2,j5,j6,j2,j3,j4,j1,za,zb)
      a60h(2,5) = a6treega(st1,j5,j6,j2,j3,j4,j1,za,zb)
      endif
c-----QL contributions
      if ((ipsgen == 2) .or. (ipsgen == 3)) then
      a60h(3,1) = -a6treeQLslc(st4,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(3,9) = -a6treeQLslc(st4,j1,j4,j2,j3,j5,j6,zb,za)
      a60h(3,2) = -a6treeQLslc(st5,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(3,10)= -a6treeQLslc(st5,j1,j4,j2,j3,j5,j6,zb,za)
      a60h(3,3) = -a6treeQLslc(st5,j4,j1,j2,j3,j6,j5,zb,za)
      a60h(3,11)= -a6treeQLslc(st5,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(3,4) = -a6treeQLslc(st4,j4,j1,j2,j3,j6,j5,zb,za)
      a60h(3,12)= -a6treeQLslc(st4,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(3,5) = -a6treeQLslc(st4,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(3,13)= -a6treeQLslc(st4,j1,j4,j2,j3,j6,j5,zb,za)
      a60h(3,6) = -a6treeQLslc(st5,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(3,14)= -a6treeQLslc(st5,j1,j4,j2,j3,j6,j5,zb,za)
      a60h(3,7) = -a6treeQLslc(st5,j4,j1,j2,j3,j5,j6,zb,za)
      a60h(3,15)= -a6treeQLslc(st5,j4,j1,j2,j3,j5,j6,za,zb)
      a60h(3,8) = -a6treeQLslc(st4,j4,j1,j2,j3,j5,j6,zb,za)
      a60h(3,16)= -a6treeQLslc(st4,j4,j1,j2,j3,j5,j6,za,zb)
      endif

      if ((ipsgen == 3) .or. (ipsgen == 4)) then
      a60h(4,1) = -a6treeQLslc(st4,j1,j4,j3,j2,j5,j6,za,zb)
      a60h(4,9) = -a6treeQLslc(st4,j1,j4,j3,j2,j5,j6,zb,za)
      a60h(4,2) = -a6treeQLslc(st5,j4,j1,j3,j2,j6,j5,zb,za)
      a60h(4,10)= -a6treeQLslc(st5,j4,j1,j3,j2,j6,j5,za,zb)
      a60h(4,3) = -a6treeQLslc(st5,j1,j4,j3,j2,j5,j6,za,zb)
      a60h(4,11)= -a6treeQLslc(st5,j1,j4,j3,j2,j5,j6,zb,za)
      a60h(4,4) = -a6treeQLslc(st4,j4,j1,j3,j2,j6,j5,zb,za)
      a60h(4,12)= -a6treeQLslc(st4,j4,j1,j3,j2,j6,j5,za,zb)
      a60h(4,5) = -a6treeQLslc(st4,j1,j4,j3,j2,j6,j5,za,zb)
      a60h(4,13)= -a6treeQLslc(st4,j1,j4,j3,j2,j6,j5,zb,za)
      a60h(4,6) = -a6treeQLslc(st5,j4,j1,j3,j2,j5,j6,zb,za)
      a60h(4,14)= -a6treeQLslc(st5,j4,j1,j3,j2,j5,j6,za,zb)
      a60h(4,7) = -a6treeQLslc(st5,j1,j4,j3,j2,j6,j5,za,zb)
      a60h(4,15)= -a6treeQLslc(st5,j1,j4,j3,j2,j6,j5,zb,za)
      a60h(4,8) = -a6treeQLslc(st4,j4,j1,j3,j2,j5,j6,zb,za)
      a60h(4,16)= -a6treeQLslc(st4,j4,j1,j3,j2,j5,j6,za,zb)
      endif
c-----done here
      return
      end


