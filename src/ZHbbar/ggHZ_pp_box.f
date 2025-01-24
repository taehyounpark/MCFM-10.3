!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c---- CW Feb 2015
c---- amplitude for ++ gluons => Z(i3,i4) + H

      function ggHZ_pp_box(i1,i2,i3,i4,za,zb,mt2)
c---- this is helicity amplitude for
c===== g(i1)^+ + g(i2)^+ + (Z=>(l1(i3)^- + l2(i3)^+)) + H
c===== and the left-handed proejection of the Z is chosen
      implicit none
      include 'types.f'
      complex(dp)::  ggHZ_pp_box
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4
      real(dp):: mt2,t
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=6)
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5)
      complex(dp):: di(Nbox),ci(Ntri),z(3),denom

      complex(dp):: D0(Nbox),C0(Ntri)
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)

      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      di(:)=czip
      ci(:)=czip

c====== integrals are now filled elsewhere and stored in common block
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


c------- coefficients of integrals

      di(d25_12)=
     & (zb(i2,i1)*(2*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(-2*mt2*za(i2,i3) + za(i1,i2)*za(i3,i4)*zb(i4,i1))
     & *zb(i4,i2)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & - za(i1,i3)*(za(i2,i3)*zb(i3,i1)
     & +za(i2,i4)*zb(i4,i1))
     & *(t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i2)
     & +2*(2*mt2 - za(i1,i2)*zb(i2,i1))*zb(i4,i1)
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))
     & -t(i1,i3,i4)*za(i1,i2)
     & *(za(i2,i3)*zb(i2,i1)*zb(i4,i1)
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & +za(i3,i4)*((za(i2,i3)*zb(i3,i1)
     & +za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2
     & -zb(i4,i1)**2
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/
     & (2.*s(i3,i4)*za(i1,i2)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))

      di(d15_12)= (zb(i2,i1)*(2*zb(i4,i1)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & *(-2*mt2*za(i1,i3) - za(i1,i2)*za(i3,i4)*zb(i4,i2))
     & -za(i2,i3)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & *(t(i2,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i1)
     & +2*(2*mt2 - za(i1,i2)*zb(i2,i1))
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *zb(i4,i2))
     & +t(i2,i3,i4)*za(i1,i2)
     & *(-(za(i1,i3)*zb(i2,i1)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *zb(i4,i2))
     & +za(i3,i4)*(-((za(i2,i3)*zb(i3,i1)
     & +za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2)
     & +zb(i4,i1)**2
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/
     & (2.*s(i3,i4)*za(i1,i2)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))

      z(1) = -4*za(i1,i2)*zb(i2,i1)*(za(i1,i3)*zb(i4,i2)*(za(i2,i3)*(za(
     &   i2,i4)*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))+za(i1,i4)*zb
     &   (i3,i1)*zb(i4,i1))+za(i2,i4)*zb(i4,i1)*(za(i2,i4)*zb(i4,i2)+za(
     &   i1,i4)*zb(i4,i1))+za(i2,i3)**2*zb(i3,i1)*zb(i3,i2))-za(i1,i2)*z
     &   a(i3,i4)*(za(i2,i4)*zb(i4,i1)+za(i2,i3)*zb(i3,i1))*zb(i4,i2)**2
     &   +za(i1,i3)**2*zb(i3,i2)*zb(i4,i1)*(za(i2,i4)*zb(i4,i1)+za(i2,i3
     &   )*zb(i3,i1)))*mt2
      z(2) = -2*za(i1,i3)*(za(i1,i3)*(za(i2,i3)*za(i2,i4)*(za(i1,i4)*zb(
     &   i4,i1)*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))**2-za(i3,i4)
     &   *zb(i4,i3)*(zb(i2,i1)**2*zb(i4,i3)**2+zb(i2,i1)*zb(i3,i2)*zb(i4
     &   ,i1)*zb(i4,i3)-zb(i3,i2)**2*zb(i4,i1)**2)+za(i2,i4)*zb(i3,i2)*z
     &   b(i4,i1)*zb(i4,i2)*(zb(i2,i1)*zb(i4,i3)+zb(i3,i2)*zb(i4,i1)))+z
     &   a(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)*(za(i1,i4)*(zb(i2,i1)*zb(i4,
     &   i3)+2*zb(i3,i2)*zb(i4,i1))+zb(i3,i2)*za(i3,i4)*zb(i4,i3))+za(i2
     &   ,i3)**2*zb(i3,i1)*zb(i3,i2)*(za(i2,i4)*zb(i4,i2)+za(i1,i4)*zb(i
     &   4,i1))*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))+za(i2,i3)**3
     &   *zb(i3,i1)*zb(i3,i2)**2*(zb(i2,i1)*zb(i4,i3)+zb(i3,i2)*zb(i4,i1
     &   )))+za(i1,i4)*(za(i2,i3)**2*(za(i2,i4)*zb(i4,i2)*(zb(i2,i1)**2*
     &   zb(i4,i3)**2+5*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)+5*zb(i3,
     &   i2)**2*zb(i4,i1)**2)+za(i3,i4)*zb(i4,i3)*(zb(i2,i1)**2*zb(i4,i3
     &   )**2+zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)-zb(i3,i2)**2*zb(i4
     &   ,i1)**2)+za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2*(zb(i2,i1)*zb(i4,i3)+
     &   zb(i3,i2)*zb(i4,i1)))+za(i2,i3)**3*zb(i3,i2)*(zb(i2,i1)**2*zb(i
     &   4,i3)**2+3*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)+2*zb(i3,i2)*
     &   *2*zb(i4,i1)**2)+za(i2,i3)*za(i2,i4)*zb(i4,i1)*zb(i4,i2)*(2*za(
     &   i2,i4)*zb(i4,i2)*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))+za
     &   (i1,i4)*zb(i4,i1)*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))-z
     &   b(i3,i2)*za(i3,i4)*zb(i4,i1)*zb(i4,i3))+za(i2,i4)**2*zb(i4,i1)*
     &   *2*zb(i4,i2)**2*(za(i2,i4)*zb(i4,i2)+za(i1,i4)*zb(i4,i1)))+za(i
     &   1,i3)**2*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*(za(i2,i4)*zb(i4,i1)+za(
     &   i2,i3)*zb(i3,i1))*(za(i2,i4)*zb(i4,i2)+za(i2,i3)*zb(i3,i2)))

      z(3) = -za(i1,i2)*(za(i1,i4)*zb(i4,i2)*((-za(i2,i3)**2*(za(i3,i4)*
     &   (zb(i2,i1)*zb(i4,i3)+zb(i3,i2)*zb(i4,i1))**2-za(i1,i4)*zb(i2,i1
     &   )*zb(i3,i1)*zb(i4,i1)*zb(i4,i2)))+za(i2,i3)*zb(i4,i1)*zb(i4,i2)
     &   *(za(i1,i4)*(zb(i2,i1)*za(i2,i4)-zb(i3,i1)*za(i3,i4))*zb(i4,i1)
     &   -2*za(i2,i4)*za(i3,i4)*(zb(i2,i1)*zb(i4,i3)+zb(i3,i2)*zb(i4,i1)
     &   ))-za(i2,i4)*za(i3,i4)*zb(i4,i1)**2*zb(i4,i2)*(za(i2,i4)*zb(i4,
     &   i2)+za(i1,i4)*zb(i4,i1)))+za(i1,i3)*(za(i2,i3)*(2*za(i3,i4)*(za
     &   (i3,i4)*zb(i4,i3)*(zb(i2,i1)**2*zb(i4,i3)**2+zb(i2,i1)*zb(i3,i2
     &   )*zb(i4,i1)*zb(i4,i3)-zb(i3,i2)**2*zb(i4,i1)**2)+za(i2,i4)*zb(i
     &   3,i2)**2*zb(i4,i1)**2*zb(i4,i2))+za(i1,i4)*zb(i4,i1)*zb(i4,i2)*
     &   (zb(i2,i1)*za(i2,i4)*(zb(i2,i1)*zb(i4,i3)+zb(i3,i1)*zb(i4,i2)+3
     &   *zb(i3,i2)*zb(i4,i1))-2*zb(i3,i1)**2*za(i3,i4)*zb(i4,i2)))+za(i
     &   2,i3)**2*(za(i3,i4)*(zb(i3,i2)**3*zb(i4,i1)**2-zb(i2,i1)**2*zb(
     &   i3,i2)*zb(i4,i3)**2)+za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i2)*(z
     &   b(i2,i1)*zb(i4,i3)+3*zb(i3,i2)*zb(i4,i1)))+za(i2,i4)*zb(i4,i1)*
     &   *2*zb(i4,i2)*(zb(i3,i2)*za(i3,i4)*(za(i2,i4)*zb(i4,i2)-2*za(i3,
     &   i4)*zb(i4,i3))+za(i1,i4)*(zb(i2,i1)*za(i2,i4)-2*zb(i3,i1)*za(i3
     &   ,i4))*zb(i4,i2)))+za(i1,i3)**2*zb(i3,i2)*(za(i2,i4)*zb(i4,i1)+z
     &   a(i2,i3)*zb(i3,i1))*(zb(i4,i1)*(zb(i2,i1)*za(i2,i4)*zb(i4,i2)-z
     &   a(i3,i4)*(2*zb(i2,i1)*zb(i4,i3)+zb(i3,i2)*zb(i4,i1)))+zb(i2,i1)
     &     *za(i2,i3)*(zb(i2,i1)*zb(i4,i3)+2*zb(i3,i2)*zb(i4,i1))))

      denom=2*s(i3,i4)*za(i1,i2)**2
     & *(za(i2,i3)*zb(i3,i1)+za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2)+za(i1,i4)*zb(i4,i2))
      di(d15_25)=(z(1)+z(2)+z(3))/denom


c---- triangle coefficients
      ci(c25_Z)=(za(i1,i3)*(za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2)*zb(i4,i1)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & +zb(i4,i2)*(za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)
     & +za(i1,i2)*zb(i2,i1)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & +za(i2,i4)*zb(i4,i1)
     & *(za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2))
     & +za(i2,i3)*(za(i1,i4)*zb(i3,i1)*zb(i4,i1)
     & +za(i2,i4)
     & *(2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))))))
     & /(za(i1,i2)**2*za(i3,i4)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      ci(cZ_15)= (za(i2,i3)*(za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))
     & *(za(i2,i3)*zb(i3,i1)*zb(i4,i2)
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & +zb(i4,i1)*(za(i1,i3)**2*zb(i3,i1)*zb(i3,i2)
     & +za(i1,i2)*zb(i2,i1)
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))
     & +za(i1,i4)*zb(i4,i2)
     & *(za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2))
     & +za(i1,i3)*(za(i2,i4)*zb(i3,i2)*zb(i4,i2)
     & +za(i1,i4)
     & *(2*zb(i3,i1)*zb(i4,i2) - zb(i2,i1)*zb(i4,i3))))))
     & /(za(i1,i2)**2*za(i3,i4)
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      ci(cH_25)=-((za(i2,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i1)
     & +(za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))
     & *(za(i1,i3)*zb(i3,i1)*zb(i4,i1)
     & +za(i1,i4)*zb(i4,i1)**2
     & +(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & *zb(i4,i2))
     & +za(i1,i2)*zb(i2,i1)
     & *(za(i1,i3)*zb(i3,i1)*zb(i4,i1)
     & +za(i1,i4)*zb(i4,i1)**2
     & +2*za(i2,i3)*zb(i3,i1)*zb(i4,i2)
     & +2*za(i2,i4)*zb(i4,i1)*zb(i4,i2)
     & -za(i2,i3)*zb(i2,i1)*zb(i4,i3))))/
     & (s(i3,i4)*za(i1,i2)**2
     & *(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))))

      ci(c15_H)= -((za(i1,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i2)
     & +(za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))
     & *(za(i2,i3)*zb(i3,i2)*zb(i4,i2)
     & +za(i2,i4)*zb(i4,i2)**2
     & +zb(i4,i1)*(za(i1,i3)*zb(i3,i2)
     & +za(i1,i4)*zb(i4,i2)))
     & +za(i1,i2)*zb(i2,i1)
     & *(2*za(i1,i3)*zb(i3,i2)*zb(i4,i1)
     & +za(i2,i3)*zb(i3,i2)*zb(i4,i2)
     & +2*za(i1,i4)*zb(i4,i1)*zb(i4,i2)
     & +za(i2,i4)*zb(i4,i2)**2
     & +za(i1,i3)*zb(i2,i1)*zb(i4,i3))))/
     & (s(i3,i4)*za(i1,i2)**2
     & *(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))

      ci(c12)=  -((zb(i2,i1)*((zb(i4,i1)
     & *(-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1)))/
     & (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))
     & -(zb(i4,i2)*(za(i1,i3)*zb(i2,i1)
     & +za(i3,i4)*zb(i4,i2)))/
     & (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))/
     & s(i3,i4))

      ggHZ_pp_box=czip
c===== make amplitude
        do i=1,Nbox
           ggHZ_pp_box=ggHZ_pp_box+D0(i)*di(i)
        enddo
        do i=1,Ntri
           ggHZ_pp_box=ggHZ_pp_box+C0(i)*ci(i)
        enddo

      return
      end
