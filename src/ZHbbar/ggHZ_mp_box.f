!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c---- CW Feb 2015
c---- amplitude for ++ gluons => Z(i3,i4) + H

      function ggHZ_mp_box(i1,i2,i3,i4,za,zb,mt2)
c---- this is helicity amplitude for
c===== g(i1)^+ + g(i2)^+ + (Z=>(l1(i3)^- + l2(i3)^+)) + H
c===== and the left-handed proejection of the Z is chosen
      implicit none
      include 'types.f'
      complex(dp)::  ggHZ_mp_box
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4
      real(dp):: mt2
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=6)

      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15,c12_Z_H
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5,c12_Z_H=6)
c--- nb I swapped the notation wrt to KC for 3m triangle, so that I
c---- can unify two routines and save calls to QCDLoop.
      complex(dp)::  ggHZ_mp_2me_b1,ggHZ_mp_2mh_b1
      external  ggHZ_mp_2me_b1,ggHZ_mp_2mh_b1
      complex(dp)::  ggHZ_mp_tri1,ggHZ_mp_tri2,ggHZ_mp_3mtri
      external  ggHZ_mp_tri1,ggHZ_mp_tri2,ggHZ_mp_3mtri

      complex(dp):: di(Nbox),ci(Ntri)
      complex(dp):: D0(Nbox),C0(Ntri)
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

      di(:)=czip
      ci(:)=czip


c----- fill integrals from QCD loop : boxes

c      D0(d25_12)=loopI4(mZsq,zip,zip,mHsq,s25,s12,mt2,mt2,mt2,mt2,musq,0)
c      D0(d15_12)=loopI4(mZsq,zip,zip,mHsq,s15,s12,mt2,mt2,mt2,mt2,musq,0)
c      D0(d15_25)=loopI4(mHsq,zip,mZsq,zip,s15,s25,mt2,mt2,mt2,mt2,musq,0)

c----- fill integrals from QCD loop : triangles
c      C0(c25_Z)=loopI3(s25,zip,mZsq,mt2,mt2,mt2,musq,0)
c      C0(cH_25)=loopI3(mHsq,zip,s25,mt2,mt2,mt2,musq,0)
c      C0(c12)=loopI3(s12,zip,zip,mt2,mt2,mt2,musq,0)
c      C0(c15_H)=loopI3(s15,zip,mHsq,mt2,mt2,mt2,musq,0)
c      C0(cZ_15)=loopI3(mZsq,zip,s15,mt2,mt2,mt2,musq,0)
c      C0(c12_Z_H)=loopI3(s12,mZsq,mHsq,mt2,mt2,mt2,musq,0)

c------ box coefficients
      di(d25_12)= ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,mt2)
      di(d15_12)= -ggHZ_mp_2mh_b1(i2,i1,i4,i3,zb,za,mt2)
      di(d15_25)= ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,mt2)


c==== triangle coefficients
      ci(c25_Z)=ggHZ_mp_tri2(i1,i2,i3,i4,za,zb)
      ci(cZ_15)=-ggHZ_mp_tri2(i2,i1,i4,i3,zb,za)
      ci(cH_25)=-ggHZ_mp_tri1(i1,i2,i3,i4,za,zb)
      ci(c15_H)=ggHZ_mp_tri1(i2,i1,i4,i3,zb,za)
      ci(c12_Z_H)=ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb)

      ci(c12)=(-(za(i1,i3)*(za(i1,i4)*za(i2,i3)
     & +za(i1,i3)*za(i2,i4))*zb(i2,i1)
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))*zb(i4,i3))
     & +za(i1,i2)**2*za(i3,i4)*zb(i2,i1)
     & *(-(za(i2,i3)*zb(i2,i1))+za(i3,i4)*zb(i4,i1))*zb(i4,i3)+za(i1,i2)
     & *(2*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)
     & +2*za(i2,i4)**2*za(i3,i4)*zb(i4,i1)
     & *zb(i4,i2)**2-za(i2,i3)**2*zb(i2,i1)
     & *(2*za(i2,i4)*zb(i3,i2)*zb(i4,i2)
     & +za(i1,i4)*zb(i2,i1)*zb(i4,i3))
     & +za(i1,i3)*zb(i2,i1)*(2*za(i1,i4)*za(i2,i4)*zb(i4,i1)**2
     & +2*za(i2,i4)**2*zb(i4,i1)*zb(i4,i2)
     & -za(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i4,i3)
     & +za(i2,i4)*(-(za(i2,i3)*zb(i2,i1))
     & +2*za(i3,i4)*zb(i4,i1))*zb(i4,i3))
     & +za(i2,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)
     & *(-2*za(i2,i4)*zb(i4,i2)+za(i3,i4)*zb(i4,i3))
     & -2*za(i2,i4)*zb(i4,i2)*(za(i2,i4)*zb(i2,i1)*zb(i4,i2)+za(i3,i4)
     & *(-(zb(i3,i1)*zb(i4,i2))+2*zb(i2,i1)*zb(i4,i3))))))
     & /(2*za(i2,i4)*za(i3,i4)
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))**2*zb(i4,i3))


c===== make amplitude
      ggHZ_mp_box=czip
      do i=1,Nbox
           ggHZ_mp_box=ggHZ_mp_box+D0(i)*di(i)
      enddo
      do i=1,Ntri
           ggHZ_mp_box=ggHZ_mp_box+C0(i)*ci(i)
      enddo

      return
      end


      function ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(dp):: ggHZ_mp_2mh_b1
c---- coefficient of 2mh box
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4
      real(dp):: mt2,t
      complex(dp):: zab2
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      ggHZ_mp_2mh_b1=(-2*za(i1,i3)*(2*mt2*zb(i4,i1)
     & *zab2(i2,i3,i4,i1)
     & *zab2(i1,i3,i4,i2)
     & +t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)
     & *(t(i1,i3,i4)*zb(i4,i1)
     & +zab2(i2,i3,i4,i1)*zb(i4,i2)))
     & +zb(i4,i2)*(-(t(i1,i3,i4)*za(i1,i2)*za(i3,i4)
     & *zab2(i2,i3,i4,i1)*zb(i4,i2))
     & +2*za(i2,i3)*(t(i1,i3,i4)**2*za(i1,i2)*zb(i2,i1)
     & +2*mt2*zab2(i2,i3,i4,i1)
     & *zab2(i1,i3,i4,i2)))
     & +t(i1,i3,i4)*za(i1,i3)**2*zb(i2,i1)*zab2(i2,i3,i4,i1)*zb(i4,i3))
     & /(2*s(i3,i4)*zab2(i2,i3,i4,i1)**2)
      return
      end

      function ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,mt2)
      implicit none
      include 'types.f'
      complex(dp):: ggHZ_mp_2me_b1
c---- coefficient of 2mh box
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4
      real(dp):: mt2
      complex(dp)::zab2,denom,z(9)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      denom=2*za(i1,i2)**2*za(i3,i4)*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i3)
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))**3
      z(1) = za(i1,i2)**2*za(i1,i3)*zb(i2,i1)*(za(i2,i3)*za(i2,i4)*zb(i3
     &   ,i1)*zb(i4,i1)*((-za(i3,i4)*((-zb(i3,i1)*zb(i4,i2)*(zb(i2,i1)*z
     &   b(i4,i3)*(8*mt2+7*s(i3,i4))+2*zb(i3,i2)*zb(i4,i1)*(4*mt2-s(i3,i
     &   4))))+zb(i2,i1)*zb(i4,i3)*(zb(i3,i2)*zb(i4,i1)*(4*mt2+9*s(i3,i4
     &   ))+zb(i2,i1)*za(i3,i4)*zb(i4,i3)**2)+6*zb(i3,i1)**2*za(i3,i4)*z
     &   b(i4,i2)**2*zb(i4,i3)))+s(i2,i4)*(zb(i2,i1)*zb(i3,i2)*zb(i4,i1)
     &   *(4*mt2+s(i3,i4))+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*(2*zb(i3,i1)*zb
     &   (i4,i2)+zb(i3,i2)*zb(i4,i1)))+4*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*(
     &   zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3)-2*zb(i3,i1)*zb(i4,i2)*m
     &   t2))+za(i2,i3)**2*zb(i3,i1)**2*(s(i2,i4)*(zb(i2,i1)*(zb(i3,i2)*
     &   zb(i4,i1)*(8*mt2+s(i3,i4))-2*zb(i3,i1)*za(i3,i4)*zb(i4,i2)*zb(i
     &   4,i3))+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*(zb(i3,i1)*zb(i4,i2)+2*zb(
     &   i3,i2)*zb(i4,i1)))+za(i3,i4)*(zb(i3,i1)*zb(i4,i2)*(zb(i2,i1)*zb
     &   (i4,i3)*(4*mt2+s(i3,i4))+zb(i3,i2)*zb(i4,i1)*(4*mt2-3*s(i3,i4))
     &   )+2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)*(s(i3,i4)-4*mt2)-zb
     &   (i3,i1)**2*za(i3,i4)*zb(i4,i2)**2*zb(i4,i3))-2*za(i1,i4)*zb(i2,
     &   i1)*zb(i4,i1)*(zb(i3,i1)*zb(i4,i2)*(2*mt2+s(i3,i4))-2*zb(i3,i2)
     &   *za(i3,i4)*zb(i4,i1)*zb(i4,i3)))+za(i2,i4)**2*zb(i4,i1)**2*(za(
     &   i3,i4)*(zb(i2,i1)*zb(i4,i3)*(zb(i3,i2)*zb(i4,i1)*(4*mt2+s(i3,i4
     &   )-s(i2,i4))-zb(i2,i1)*za(i3,i4)*zb(i4,i3)**2)+zb(i3,i1)*zb(i4,i
     &   2)*(zb(i3,i2)*zb(i4,i1)*(4*mt2+s(i3,i4))+2*zb(i2,i1)*zb(i4,i3)*
     &   (2*mt2-3*s(i3,i4)+s(i2,i4)))+zb(i3,i1)**2*(s(i2,i4)-5*s(i3,i4))
     &   *zb(i4,i2)**2)+2*za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)*zb(i4,
     &   i2)*(s(i3,i4)-2*mt2))+za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)*zb(i4
     &   ,i2)*(zb(i2,i1)*(4*mt2-s(i3,i4))+zb(i3,i1)*za(i3,i4)*zb(i4,i2))
     &   )
      z(2) = za(i1,i2)**2*za(i1,i4)*zb(i2,i1)*(za(i2,i3)*za(i2,i4)*zb(i4
     &   ,i1)**2*(s(i2,i4)*(zb(i2,i1)*(zb(i3,i1)*zb(i4,i2)*(4*mt2-7*s(i3
     &   ,i4))+2*zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3))+zb(i3,i1)*za(i
     &   3,i4)*zb(i4,i2)*(2*zb(i3,i1)*zb(i4,i2)+zb(i3,i2)*zb(i4,i1)))-2*
     &   za(i3,i4)*(zb(i3,i1)*zb(i4,i2)*zb(i4,i3)*(zb(i2,i1)*(2*mt2-5*s(
     &   i3,i4))+2*zb(i3,i2)*za(i3,i4)*zb(i4,i1))+2*zb(i3,i1)**2*zb(i4,i
     &   2)**2*(s(i3,i4)-2*mt2)+3*zb(i2,i1)*zb(i3,i2)*za(i3,i4)*zb(i4,i1
     &   )*zb(i4,i3)**2))+za(i2,i3)**2*zb(i3,i1)*zb(i4,i1)*(2*za(i3,i4)*
     &   ((-2*zb(i3,i1)*zb(i4,i2)*zb(i4,i3)*(2*zb(i2,i1)*mt2+zb(i3,i2)*z
     &   a(i3,i4)*zb(i4,i1)))+2*zb(i3,i1)**2*zb(i4,i2)**2*mt2+3*zb(i2,i1
     &   )*zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3)**2)+za(i2,i4)*zb(i3,i
     &   1)*zb(i4,i2)**2*(zb(i2,i1)*(8*mt2-7*s(i3,i4))+za(i3,i4)*(zb(i3,
     &   i1)*zb(i4,i2)+2*zb(i3,i2)*zb(i4,i1))))+za(i2,i3)**3*zb(i3,i1)**
     &   2*zb(i4,i2)*(zb(i2,i1)*(zb(i3,i1)*zb(i4,i2)*(4*mt2-s(i3,i4))-2*
     &   zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3))+zb(i3,i1)*zb(i3,i2)*za
     &   (i3,i4)*zb(i4,i1)*zb(i4,i2))+za(i2,i4)**2*za(i3,i4)*zb(i4,i1)**
     &   3*zb(i4,i2)*(zb(i3,i1)*zb(i4,i2)*(4*mt2-4*s(i3,i4)+s(i2,i4))-zb
     &   (i2,i1)*zb(i4,i3)*((-4*mt2)+2*s(i3,i4)+s(i2,i4))))
      z(3) = -za(i1,i2)**2*za(i1,i3)**2*zb(i2,i1)**2*zb(i3,i1)*(za(i2,i4
     &   )*zb(i4,i1)+za(i2,i3)*zb(i3,i1))*(za(i2,i4)*zb(i4,i1)*(zb(i3,i2
     &   )*zb(i4,i1)*(4*mt2+s(i3,i4))-3*zb(i3,i1)*za(i3,i4)*zb(i4,i2)*zb
     &   (i4,i3))+za(i2,i3)*zb(i3,i1)*(zb(i3,i2)*zb(i4,i1)*(4*mt2-3*s(i3
     &   ,i4))+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*zb(i4,i3)))
      z(4) = za(i1,i2)*za(i1,i3)*za(i1,i4)*(za(i2,i3)*za(i2,i4)**2*zb(i4
     &   ,i1)**2*(zb(i2,i1)*(zb(i3,i1)**2*zb(i4,i2)**2*((-2*(mt2+9*s(i3,
     &   i4)))+s(i2,i4)+2*s(i1,i4))+zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,
     &   i2)*(2*mt2+23*s(i3,i4)-6*s(i2,i4)-2*s(i1,i4))+zb(i3,i2)**2*zb(i
     &   4,i1)**2*((-4*mt2)-2*s(i3,i4)+s(i2,i4)))+2*zb(i2,i1)**2*zb(i4,i
     &   3)*(zb(i3,i1)*zb(i4,i2)*mt2+zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4
     &   ,i3))+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*(6*zb(i3,i1)**2*zb(i4,i2)**
     &   2+7*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)-zb(i3,i2)**2*zb(i4,
     &   i1)**2))+za(i2,i3)**2*za(i2,i4)*zb(i3,i1)*zb(i4,i1)*(zb(i2,i1)*
     &   ((-zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*(8*mt2+27*s(i3,i4)))
     &   +zb(i3,i1)**2*zb(i4,i2)**2*(8*mt2+s(i2,i4)+2*s(i1,i4))+zb(i3,i2
     &   )**2*zb(i4,i1)**2*(4*mt2+18*s(i3,i4)-s(i2,i4)-2*s(i1,i4)))+zb(i
     &   2,i1)**2*zb(i4,i3)*(zb(i3,i1)*zb(i4,i2)*(4*mt2-s(i3,i4))+2*zb(i
     &   3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3))+zb(i3,i1)*za(i3,i4)*zb(i4,
     &   i2)*(zb(i3,i1)**2*zb(i4,i2)**2+9*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     &   zb(i4,i2)+2*zb(i3,i2)**2*zb(i4,i1)**2))-za(i2,i4)**3*zb(i4,i1)*
     &   *3*zb(i4,i2)*(zb(i2,i1)*(zb(i3,i1)*zb(i4,i2)*(8*mt2-6*s(i3,i4)+
     &   s(i2,i4))+zb(i3,i2)*zb(i4,i1)*((-4*mt2)+s(i3,i4)+s(i2,i4)))-zb(
     &   i2,i1)**2*za(i3,i4)*zb(i4,i3)**2+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*
     &   (zb(i3,i2)*zb(i4,i1)-5*zb(i3,i1)*zb(i4,i2)))+za(i2,i3)**3*zb(i3
     &   ,i1)**2*(zb(i2,i1)*(zb(i3,i1)**2*zb(i4,i2)**2*(2*mt2-s(i2,i4))+
     &   zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*((-6*mt2)-3*s(i3,i4)+6*
     &   s(i2,i4)+2*s(i1,i4))-zb(i3,i2)**2*zb(i4,i1)**2*((-8*mt2)+4*s(i3
     &   ,i4)+s(i2,i4)+2*s(i1,i4)))+2*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i2)*z
     &   b(i4,i3)*mt2+zb(i3,i1)*zb(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i2)*
     &   (zb(i3,i1)*zb(i4,i2)+3*zb(i3,i2)*zb(i4,i1)))+zb(i2,i1)*za(i2,i3
     &   )**4*zb(i3,i1)**3*zb(i3,i2)*zb(i4,i2)*(zb(i3,i1)*zb(i4,i2)+zb(i
     &   3,i2)*zb(i4,i1)))
      z(5) = za(i1,i2)*za(i1,i4)**2*zb(i4,i1)*(za(i2,i3)*za(i2,i4)**2*zb
     &   (i4,i1)**2*zb(i4,i2)*(zb(i2,i1)*(zb(i3,i1)*zb(i4,i2)*(8*mt2-10*
     &   s(i3,i4)+s(i2,i4))+zb(i3,i2)*zb(i4,i1)*((-4*mt2)+4*s(i3,i4)+s(i
     &   2,i4)))+2*zb(i3,i1)*za(i3,i4)*zb(i4,i2)*(zb(i3,i1)*zb(i4,i2)+2*
     &   zb(i3,i2)*zb(i4,i1)))+za(i2,i3)**3*zb(i3,i1)*(2*zb(i3,i1)*zb(i3
     &   ,i2)**2*za(i3,i4)*zb(i4,i1)**2*zb(i4,i2)-zb(i2,i1)*((-zb(i3,i1)
     &   *zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*(8*mt2+7*s(i2,i4)))+zb(i3,i1)**2
     &   *zb(i4,i2)**2*(4*mt2+s(i2,i4))+6*zb(i3,i2)**2*za(i3,i4)*zb(i4,i
     &   1)**2*zb(i4,i3)))+za(i2,i3)**2*za(i2,i4)*zb(i4,i1)*(zb(i2,i1)*(
     &   zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*(4*mt2-20*s(i3,i4)+7*s(
     &   i2,i4))+2*zb(i3,i1)**2*zb(i4,i2)**2*(2*mt2+s(i3,i4))+zb(i3,i2)*
     &   *2*(6*s(i3,i4)-s(i2,i4))*zb(i4,i1)**2)+2*zb(i3,i1)*zb(i3,i2)*za
     &   (i3,i4)*zb(i4,i1)*zb(i4,i2)*(2*zb(i3,i1)*zb(i4,i2)+zb(i3,i2)*zb
     &   (i4,i1)))+2*za(i2,i4)**3*zb(i3,i1)*za(i3,i4)*zb(i4,i1)**3*zb(i4
     &   ,i2)**3+zb(i2,i1)*za(i2,i3)**4*zb(i3,i1)**2*zb(i3,i2)*zb(i4,i2)
     &   *(zb(i3,i1)*zb(i4,i2)+zb(i3,i2)*zb(i4,i1)))
      z(6) = -za(i1,i2)*za(i1,i3)**2*(za(i2,i3)**2*za(i2,i4)*zb(i3,i1)**
     &   2*(zb(i2,i1)*(2*zb(i3,i2)**2*zb(i4,i1)**2*(5*mt2-4*s(i3,i4)+s(i
     &   1,i4))-zb(i3,i1)**2*zb(i4,i2)**2*((-4*mt2)+s(i2,i4)+3*s(i1,i4))
     &   +zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*((-10*mt2)+5*s(i3,i4)+
     &   7*s(i2,i4)+s(i1,i4)))+zb(i2,i1)**2*zb(i4,i3)*(2*zb(i3,i2)*zb(i4
     &   ,i1)*(mt2-s(i3,i4))+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*zb(i4,i3))-zb
     &   (i3,i1)**2*za(i3,i4)*zb(i4,i2)**2*(zb(i3,i1)*zb(i4,i2)+5*zb(i3,
     &   i2)*zb(i4,i1)))+za(i2,i3)*za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)*(zb(
     &   i2,i1)*(2*zb(i3,i1)**2*zb(i4,i2)**2*(4*mt2+2*s(i3,i4)-s(i1,i4))
     &   -zb(i3,i2)**2*zb(i4,i1)**2*((-8*mt2)-4*s(i3,i4)+s(i2,i4)+s(i1,i
     &   4))+zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*((-8*mt2)-17*s(i3,i
     &   4)+7*s(i2,i4)+3*s(i1,i4)))-2*zb(i2,i1)**2*zb(i3,i2)*zb(i4,i1)*z
     &   b(i4,i3)*(s(i3,i4)-2*mt2)+zb(i3,i1)*za(i3,i4)*zb(i4,i2)*((-4*zb
     &   (i3,i1)**2*zb(i4,i2)**2)-3*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,
     &   i2)+zb(i3,i2)**2*zb(i4,i1)**2))+za(i2,i4)**3*zb(i4,i1)**2*(zb(i
     &   2,i1)*(zb(i3,i1)**2*zb(i4,i2)**2*(4*mt2-8*s(i3,i4)+s(i2,i4)+s(i
     &   1,i4))+zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)*(2*mt2+3*s(i3,i4
     &   )+s(i2,i4)-s(i1,i4))-2*zb(i3,i2)**2*zb(i4,i1)**2*mt2)+zb(i2,i1)
     &   **2*zb(i4,i3)*(2*zb(i3,i2)*zb(i4,i1)*mt2-zb(i3,i1)*za(i3,i4)*zb
     &   (i4,i2)*zb(i4,i3))+zb(i3,i1)**2*za(i3,i4)*zb(i4,i2)**2*(zb(i3,i
     &   2)*zb(i4,i1)-3*zb(i3,i1)*zb(i4,i2)))+za(i2,i3)**3*zb(i3,i1)**3*
     &   zb(i3,i2)*(zb(i4,i2)*(zb(i2,i1)*za(i2,i4)*(zb(i3,i1)*zb(i4,i2)+
     &   zb(i3,i2)*zb(i4,i1))-2*zb(i3,i1)*zb(i3,i2)*za(i3,i4)*zb(i4,i1))
     &   +3*za(i1,i4)*zb(i2,i1)*zb(i4,i1)*(zb(i3,i2)*zb(i4,i1)-zb(i3,i1)
     &   *zb(i4,i2))))
      z(7) = za(i1,i2)*za(i1,i3)**3*zb(i2,i1)*(s(i2,i4)+s(i2,i3))*zb(i3,
     &   i1)**2*(za(i2,i3)**2*zb(i3,i1)**2-za(i2,i4)**2*zb(i4,i1)**2)*(z
     &   b(i3,i1)*zb(i4,i2)-zb(i3,i2)*zb(i4,i1))
      z(8) = -zab2(i1,i3,i4,i1)*(s(i2,i4)+s(i2,i3))*(za(i1,i3)*s(i1,i4)*
     &   ((-za(i2,i3)*za(i2,i4)**2*zb(i4,i1)*((-zb(i3,i2)*zb(i4,i1)*(zb(
     &   i2,i1)*zb(i4,i3)+8*zb(i3,i1)*zb(i4,i2)))+zb(i3,i1)*zb(i4,i2)*(z
     &   b(i2,i1)*zb(i4,i3)+5*zb(i3,i1)*zb(i4,i2))+zb(i3,i2)**2*zb(i4,i1
     &   )**2))+za(i2,i3)**2*za(i2,i4)*zb(i3,i1)*(zb(i3,i2)*zb(i4,i1)*(z
     &   b(i2,i1)*zb(i4,i3)-10*zb(i3,i1)*zb(i4,i2))+zb(i3,i1)*zb(i3,i2)*
     &   zb(i4,i1)*zb(i4,i2)+7*zb(i3,i2)**2*zb(i4,i1)**2)+2*za(i2,i4)**3
     &   *zb(i3,i1)*zb(i4,i1)**2*zb(i4,i2)**2-2*za(i2,i3)**3*zb(i3,i1)**
     &   3*zb(i3,i2)*zb(i4,i2))+za(i1,i3)**2*za(i2,i4)*zb(i3,i1)*(za(i2,
     &   i3)*za(i2,i4)*zb(i4,i1)*(zb(i3,i2)*zb(i4,i1)*(zb(i2,i1)*zb(i4,i
     &   3)+6*zb(i3,i1)*zb(i4,i2))+zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i
     &   2)-3*zb(i3,i2)**2*zb(i4,i1)**2)+za(i2,i3)**2*zb(i3,i1)*(zb(i2,i
     &   1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)+zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*
     &   zb(i4,i2)+zb(i3,i2)**2*zb(i4,i1)**2)-2*za(i2,i4)**2*zb(i4,i1)**
     &   2*zb(i4,i2)*(zb(i3,i2)*zb(i4,i1)-2*zb(i3,i1)*zb(i4,i2)))-2*za(i
     &   1,i4)**2*za(i2,i3)*zb(i4,i1)**2*((-za(i2,i3)*za(i2,i4)*(zb(i3,i
     &   1)**2*zb(i4,i2)**2-4*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)+zb
     &   (i3,i2)**2*zb(i4,i1)**2))+za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)*zb(i
     &   4,i2)**2+za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)**2*zb(i4,i1)))
      z(9) = -za(i1,i2)**3*zb(i2,i1)**2*za(i3,i4)**2*zb(i4,i3)*(za(i2,i4
     &   )*zb(i4,i1)**2*(s(i2,i4)*(zb(i2,i1)*zb(i4,i3)+zb(i3,i1)*zb(i4,i
     &   2))-2*zb(i4,i3)*(za(i3,i4)*(zb(i2,i1)*zb(i4,i3)+zb(i3,i1)*zb(i4
     &   ,i2))-za(i1,i3)*zb(i2,i1)*zb(i3,i1)))+2*za(i2,i3)*zb(i3,i1)*zb(
     &   i4,i1)*(zb(i4,i3)*(za(i3,i4)*(zb(i2,i1)*zb(i4,i3)-zb(i3,i1)*zb(
     &   i4,i2))+za(i1,i3)*zb(i2,i1)*zb(i3,i1))+za(i2,i4)*zb(i3,i1)*zb(i
     &   4,i2)**2)+s(i2,i3)*za(i2,i3)*zb(i3,i1)**2*zb(i4,i1)*zb(i4,i2))
      ggHZ_mp_2me_b1=
     & (z(1)+z(2)+z(3)+z(4)+z(5)+z(6)+z(7)+z(8)+z(9))/denom
      return
      end

      function ggHZ_mp_tri1(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp)::ggHZ_mp_tri1
c---- coefficient of 2mh box
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4

      ggHZ_mp_tri1=(-(za(i1,i2)**3*za(i3,i4)*zb(i2,i1)**2*zb(i4,i1)
     & *zb(i4,i2))+za(i1,i2)**2*zb(i2,i1)*zb(i4,i1)
     & *(2*za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)+za(i1,i3)
     & *(2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)
     & +(3*za(i2,i4)*zb(i2,i1)-za(i3,i4)*zb(i3,i1))*zb(i4,i2))
     & -zb(i4,i2)*(3*za(i1,i4)*za(i2,i3)*zb(i2,i1)+2*za(i3,i4)
     & *(za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2))))
     & -za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)
     & *zb(i4,i2)*(za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2))*zb(i4,i3)
     & +za(i1,i3)**2*zb(i3,i1)*(za(i2,i4)**2*zb(i4,i1)*zb(i4,i2)**2
     & -2*za(i2,i3)**2*zb(i2,i1)*zb(i3,i2)*zb(i4,i3)
     & +za(i2,i3)*za(i2,i4)*zb(i4,i2)
     & *(zb(i3,i2)*zb(i4,i1)-3*zb(i2,i1)*zb(i4,i3)))+za(i1,i3)
     & *(za(i2,i4)*za(i3,i4)*zb(i4,i1)
     & *zb(i4,i2)*(za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2))*zb(i4,i3)
     & -za(i1,i4)*za(i2,i3)*(za(i2,i4)*zb(i4,i1)*zb(i4,i2)
     & *(zb(i3,i2)*zb(i4,i1)+3*zb(i2,i1)*zb(i4,i3))
     & +za(i2,i3)*(zb(i3,i2)**2*zb(i4,i1)**2
     & +2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)
     & -zb(i2,i1)**2*zb(i4,i3)**2)))
     & -za(i1,i2)*(za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)
     & *(-3*za(i2,i4)*zb(i4,i1)*zb(i4,i2)
     & +za(i2,i3)*(-2*zb(i3,i2)*zb(i4,i1)
     & +2*zb(i2,i1)*zb(i4,i3)))+zb(i4,i2)*(za(i3,i4)
     & *(2*za(i2,i3)*zb(i2,i1)+za(i3,i4)*zb(i4,i1))
     & *(za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2))*zb(i4,i3)
     & +2*za(i1,i4)*za(i2,i3)*zb(i2,i1)*(2*za(i2,i4)*zb(i4,i1)*zb(i4,i2)
     & +za(i2,i3)*(2*zb(i3,i2)*zb(i4,i1)
     & -zb(i2,i1)*zb(i4,i3))))+za(i1,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)
     & *(-2*za(i2,i4)*zb(i4,i1)*zb(i4,i2)+za(i2,i3)
     & *(-(zb(i3,i2)*zb(i4,i1))+3*zb(i2,i1)*zb(i4,i3)))
     & +zb(i4,i2)*(za(i2,i4)
     & *(-4*za(i2,i4)*zb(i2,i1)+za(i3,i4)*zb(i3,i1))*zb(i4,i1)
     & *zb(i4,i2)+za(i2,i3)
     & *(za(i3,i4)*zb(i3,i1)*(zb(i3,i2)*zb(i4,i1)
     & -zb(i2,i1)*zb(i4,i3))+2*za(i2,i4)*zb(i2,i1)
     & *(-2*zb(i3,i2)*zb(i4,i1)+zb(i2,i1)*zb(i4,i3)))))))/
     & (2*s(i3,i4)*za(i1,i2)*zb(i2,i1)
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))**2)
      return
      end

      function ggHZ_mp_tri2(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: ggHZ_mp_tri2
c---- coefficient of 2mh box
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4

      ggHZ_mp_tri2=(zb(i4,i1)*(2*za(i1,i3)**3*zb(i2,i1)*zb(i3,i1)*
     & (za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))+
     & 2*za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)*
     & zb(i4,i2)*(za(i1,i2)*zb(i2,i1)+za(i2,i3)*zb(i3,i2)+za(i2,i4)
     & *zb(i4,i2))+2*za(i1,i3)*za(i3,i4)*(za(i1,i2)*zb(i2,i1)*
     & (-(za(i1,i4)*zb(i4,i1)**2)+za(i2,i3)*zb(i3,i1)*zb(i4,i2))+
     & za(i2,i3)*(zb(i3,i1)*zb(i4,i2)*
     & (za(i2,i3)*zb(i3,i2)+za(i2,i4)*zb(i4,i2))+za(i1,i4)*zb(i4,i1)*
     & (-(zb(i3,i2)*zb(i4,i1))+zb(i3,i1)*zb(i4,i2))))+za(i1,i3)**2*
     & (2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))
     & +za(i3,i4)*(-2*za(i1,i2)*zb(i2,i1)*zb(i3,i1)*
     & zb(i4,i1)+za(i2,i4)*zb(i4,i1)*(zb(i3,i2)*zb(i4,i1)
     & -zb(i3,i1)*zb(i4,i2)+zb(i2,i1)*zb(i4,i3))+za(i2,i3)*zb(i3,i1)*
     & (-(zb(i3,i2)*zb(i4,i1))+zb(i3,i1)*zb(i4,i2)
     & +zb(i2,i1)*zb(i4,i3))))))
     & /(2*za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     & (za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))**2*zb(i4,i3))
      return
      end

      function ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: ggHZ_mp_3mtri
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4
      real(dp):: ssum,root,lamp,den,den2,den3,bit
      complex(dp):: zab2,denom
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

c==== part 1
      ssum=0.5_dp*(s(i1,i3)+s(i2,i3)+s(i1,i4)+s(i2,i4))
      root=sqrt(ssum**2-s(i1,i2)*s(i3,i4))
      lamp=ssum+root
      den=(1._dp-(s(i1,i2)*s(i3,i4))/lamp**2)
      den2=s(i1,i2)*(1._dp-(s(i1,i3)+s(i1,i4))/lamp)
      den3=(s(i1,i3)+s(i1,i4)-(s(i1,i2)*s(i3,i4))/lamp)
      bit=-2*s(i3,i4)*(s(i1,i2)*den2*den3)**2*(lamp**2-s(i1,i2)*s(i3,i4))
      denom=cmplx(bit,kind=dp)*zab2(i2,i3,i4,i1)**4

      ggHZ_mp_3mtri=(-2*lamp**2*den**8*((lamp*den2*(s(i1,i2)*((zab2(i2,i3,i4,i1)*zab2(i1,i3,i4,i2))/den**2
     & -(den3*(s(i2,i3)+s(i2,i4)-(s(i1,i2)*s(i3,i4))/lamp))/den**2)
     & *((s(i1,i2)**2*zab2(i2,i3,i4,i1)**2*(-(za(i1,i2)*zb(i4,i2))+(s(i1,i2)*za(i1,i3)*zb(i4,i3))/lamp)*den3
     & **2*(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp))/(lamp**2*den**6)
     & +(zab2(i2,i3,i4,i1)**2*(-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp)*den2
     & **2*(-(za(i1,i3)*zb(i4,i3))+(za(i1,i2)*zb(i4,i2)*s(i3,i4))/lamp))/den**6)-(s(i1,i2)*zab2(i2,i3,i4,i1)**2*den3
     & *(za(i1,i2)*zb(i4,i2)*((zab2(i2,i3,i4,i1)*zab2(i1,i3,i4,i2)*(-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp
     & )*(s(i1,i2)-(s(i1,i2)*(s(i1,i3)+s(i1,i4)))/lamp))/den**4
     & +(s(i1,i2)**2*zab2(i2,i3,i4,i1)*zab2(i1,i3,i4,i2)*den3
     & *(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp))/(lamp**2*den**4))
     & +za(i1,i3)*zb(i2,i1)*(-((s(i1,i2)*zab2(i2,i3,i4,i1)
     & *(-(za(i1,i2)*zb(i4,i2))+(s(i1,i2)*za(i1,i3)*zb(i4,i3))/lamp)*den3**2)/(lamp*den**4))
     & +(zab2(i2,i3,i4,i1)*den2**2*(-(za(i1,i3)*zb(i4,i3))+(za(i1,i2)*zb(i4,i2)*s(i3,i4))/lamp))/den**4)))/den**3))/den
     & +s(i3,i4)*(-((s(i1,i2)**4*zab2(i2,i3,i4,i1)**3*zab2(i1,i3,i4,i2)*(s(i1,i2)+lamp)
     & *(-(za(i1,i2)*zb(i4,i2))+(s(i1,i2)*za(i1,i3)*zb(i4,i3))/lamp)*den3**3
     & *(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp))/(lamp**4*den**9))
     & +(s(i1,i2)*zab2(i2,i3,i4,i1)**3*(-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp)
     & *(s(i1,i2)-(s(i1,i2)*(s(i1,i3)+s(i1,i4)))/lamp)**2*den3
     & *((za(i1,i2)*s(i1,i2)*zb(i4,i2)*zab2(i1,i3,i4,i2))/den-(s(i1,i2)*zab2(i1,i3,i4,i2)*(s(i1,i2)+lamp)
     & *(-(za(i1,i3)*zb(i4,i3))+(za(i1,i2)*zb(i4,i2)*s(i3,i4))/lamp))/(lamp*den
     & **2)))/(lamp*den**7)+(zab2(i2,i3,i4,i1)**2*den2**3*(-(za(i1,i3)*zb(i4,i3))+(za(i1,i2)*zb(i4,i2)*s(i3,i4))/lamp
     & )*((s(i1,i2)*zab2(i2,i3,i4,i1)*zab2(i1,i3,i4,i2)
     & *(-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp))/den**3+(den3*(((s(i1,i2)+lamp)
     & *(-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp)*(s(i1,i2)-(s(i1,i2)*(s(i2,i3)+s(i2,i4)))/lamp
     & ))/den**2-s(i1,i2)*(-((s(i1,i2)*za(i1,i3)*zb(i2,i1)*zab2(i2,i3,i4,i1))/(lamp*den
     & ))+((-(za(i2,i3)*zb(i2,i1))-(s(i1,i2)*za(i3,i4)*zb(i4,i1))/lamp)*(s(i2,i3)+s(i2,i4)-(s(i1,i2)*s(i3,i4))/lamp
     & ))/den**2)))/den))/den**6+(s(i1,i2)**2*zab2(i2,i3,i4,i1)**2*den2*den3
     & **2*((za(i1,i2)*s(i1,i2)**2*zab2(i2,i3,i4,i1)*zb(i4,i2)*zab2(i1,i3,i4,i2)
     & *(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp
     & ))/(lamp*den**3)+((-(za(i1,i2)*zb(i4,i2))+(s(i1,i2)*za(i1,i3)*zb(i4,i3))/lamp
     & )*((s(i1,i2)*zab2(i2,i3,i4,i1)*zab2(i1,i3,i4,i2)*(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp
     & ))/den**3+(den3*(((s(i1,i2)+lamp)*(s(i1,i2)-(s(i1,i2)*(s(i2,i3)+s(i2,i4)))/lamp
     & )*(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp
     & ))/den**2-s(i1,i2)*((za(i1,i3)*zb(i2,i1)*zab2(i2,i3,i4,i1))/den+((s(i2,i3)+s(i2,i4)-(s(i1,i2)*s(i3,i4))/lamp
     & )*(za(i3,i4)*zb(i4,i1)+(za(i2,i3)*zb(i2,i1)*s(i3,i4))/lamp))/den**2)))/den))/den))/(lamp**2*den**5))))
      ggHZ_mp_3mtri=ggHZ_mp_3mtri/denom
      return
      end
