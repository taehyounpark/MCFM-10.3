!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine GLYI2(z,I2)
c     Elaborated using the form results of 1403.6451v2
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nfl.f'
      include 'zeta.f'
      include 'transitionlabels.f'
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     -2=finite,-1=delta=-1,0=L0(z),1=L1(z)
      real(dp):: I2(0:6,dmin:dmax,0:4)
c      I2, index one transition flavor
c      I2, index2,distribution type
c      I2, index3,power of Lb
      real(dp):: tp0(0:6,dmin:dmax),mp0(0:6,dmin:dmax)
      real(dp):: z,mz,Li2,Li3,lnz,omz,opz,lnomz,lnopz,
     & Li2z,Li2mz,Li2omz,Li3z,Li3mz,Li3omz,Li3zonopz

c     initialize to zero
      I2(:,:,:)=0
      omz=1._dp-z
      opz=1._dp+z
      lnomz=log(omz); lnopz=log(opz); lnz=log(z)
      Li2z=Li2(z); Li2mz=Li2(-z); Li2omz=Li2(omz)
      Li3z=Li3(z); Li3mz=Li3(-z); Li3omz=Li3(omz)
      Li3zonopz=Li3(z/opz)

      call tildep0(+z,tp0)
      mz=-z
      call tildep0(mz,mp0)

cgg
      I2(gg,rglr,2)=I2(gg,rglr,2) + CA**2 * ( 24.D0 - 88.D0/3.D0*z**(-1)
     &     + 16*z**(-1)*lnomz - 32*lnomz - 8*lnz*omz**(-1) - 16*lnz -
     &    24*z + 16*z*lnomz - 16*z*lnz + 88.D0/3.D0*z**2 - 16*z**2*
     &    lnomz )
      I2(gg,rglr,2) = I2(gg,rglr,2) + CF*TF * ( 16.D0/3.D0*z**(-1)*nfl
     &     + 4*nfl + 8*nfl*lnz - 4*z*nfl + 8*z*nfl*lnz - 16.D0/3.D0*
     &    z**2*nfl )
      I2(gg,rglr,1)=I2(gg,rglr,1) + CA*TF * ( 104.D0/9.D0*z**(-1)*nfl
     &     - 8*nfl + 16.D0/3.D0*nfl*lnz + 8*z*nfl + 16.D0/3.D0*z*nfl*
     &    lnz - 104.D0/9.D0*z**2*nfl )
      I2(gg,rglr,1) = I2(gg,rglr,1) + CA**2 * (  - 54.D0 + 268.D0/9.D0*
     &    z**(-1) + 100.D0/3.D0*lnz + 16*lnz*lnomz*omz**(-1) - 4*lnz**2
     &    *omz**(-1) - 16*lnz**2 + 54*z - 44.D0/3.D0*z*lnz - 16*z*
     &    lnz**2 - 268.D0/9.D0*z**2 + 176.D0/3.D0*z**2*lnz )
      I2(gg,rglr,1) = I2(gg,rglr,1) + CF*TF * (  - 16.D0/3.D0*z**(-1)*
     &    nfl + 56*nfl + 24*nfl*lnz + 8*nfl*lnz**2 - 40*z*nfl + 24*z*
     &    nfl*lnz + 8*z*nfl*lnz**2 - 32.D0/3.D0*z**2*nfl )
      I2(gg,rglr,0)=I2(gg,rglr,0) + CA*TF * ( 260.D0/27.D0*z**(-1)*nfl
     &     - 8*nfl + 52.D0/9.D0*nfl*lnz + 4.D0/3.D0*nfl*lnz**2 + 8*z*
     &    nfl - 4.D0/3.D0*z*nfl*lnomz + 40.D0/9.D0*z*nfl*lnz + 4.D0/3.D0
     &    *z*nfl*lnz**2 - 332.D0/27.D0*z**2*nfl )
      I2(gg,rglr,0) = I2(gg,rglr,0) + CA**2 * ( 232.D0/3.D0 - 784.D0/9.D
     &    0*z**(-1) - 88.D0/3.D0*z**(-1)*Li2z - 88.D0/3.D0*z**(-1)*lnz*
     &    lnomz + 88.D0/3.D0*z**(-1)*zeta2 - 24*Li3z*omz**(-1) + 32*
     &    Li2z - 701.D0/9.D0*lnz + 16*lnz*Li2z*omz**(-1) + 32*lnz*lnomz
     &     + 4*lnz*lnomz**2*omz**(-1) + 25.D0/3.D0*lnz**2 + 4*lnz**2*
     &    lnomz*omz**(-1) - 2.D0/3.D0*lnz**3*omz**(-1) - 8.D0/3.D0*
     &    lnz**3 + 24*zeta3*omz**(-1) - 32*zeta2 - 248.D0/3.D0*z - 32*z
     &    *Li2z + 2.D0/3.D0*z*lnomz - 149.D0/9.D0*z*lnz - 32*z*lnz*
     &    lnomz - 11.D0/3.D0*z*lnz**2 - 8.D0/3.D0*z*lnz**3 + 32*z*zeta2
     &     + 844.D0/9.D0*z**2 + 88.D0/3.D0*z**2*Li2z - 536.D0/9.D0*z**2
     &    *lnz + 88.D0/3.D0*z**2*lnz*lnomz + 44.D0/3.D0*z**2*lnz**2 -
     &    88.D0/3.D0*z**2*zeta2 )
      I2(gg,rglr,0) = I2(gg,rglr,0) + CF*TF * (  - 8.D0/3.D0*z**(-1)*
     &    nfl + 64*nfl + 24*nfl*lnz + 6*nfl*lnz**2 + 4.D0/3.D0*nfl*
     &    lnz**3 - 64*z*nfl + 24*z*nfl*lnz + 2*z*nfl*lnz**2 + 4.D0/3.D0
     &    *z*nfl*lnz**3 + 8.D0/3.D0*z**2*nfl )

      I2(gg,rglr,3)=I2(gg,rglr,3) + tp0(gg,rglr)*CA**2 * (  - 4.D0 )
      I2(gg,rglr,2)=I2(gg,rglr,2) + tp0(gg,rglr)*CA*TF * ( 8.D0/3.D0*
     &    nfl )
      I2(gg,rglr,2) = I2(gg,rglr,2) + tp0(gg,rglr)*CA**2 * (  - 22.D0/3.
     &    D0 - 8*lnz )
      I2(gg,rglr,1)=I2(gg,rglr,1) + tp0(gg,rglr)*CA*TF * ( 80.D0/9.D0*
     &    nfl )
      I2(gg,rglr,1) = I2(gg,rglr,1) + tp0(gg,rglr)*CA**2 * (  - 268.D0/
     &    9.D0 + 16*lnz*lnomz - 4*lnz**2 + 12*zeta2 )
      I2(gg,rglr,0)=I2(gg,rglr,0) + tp0(gg,rglr)*CA*TF * ( 224.D0/27.D0
     &    *nfl )
      I2(gg,rglr,0) = I2(gg,rglr,0) + tp0(gg,rglr)*CA**2 * (  - 808.D0/
     &    27.D0 - 24*Li3z + 16*lnz*Li2z + 4*lnz*lnomz**2 + 4*lnz**2*
     &    lnomz - 2.D0/3.D0*lnz**3 + 52*zeta3 )

      I2(gg,rglr,1)=I2(gg,rglr,1) + mp0(gg,rglr)*CA**2 * ( 16*Li2mz +
     &    16*lnz*lnopz - 4*lnz**2 + 8*zeta2 )
      I2(gg,rglr,0)=I2(gg,rglr,0) + mp0(gg,rglr)*CA**2 * (  - 16*
     &    Li3zonopz + 8*Li3mz + 16*Li3z + 8.D0/3.D0*lnopz**3 - 8*lnz*
     &    Li2mz - 8*lnz*Li2z - 8*lnz*lnopz**2 + 4*lnz**2*lnopz - 2.D0/3.
     &    D0*lnz**3 + 4*zeta3 - 8*zeta2*lnopz )

      I2(gg,delt,4)=I2(gg,delt,4) + CA**2 * ( 1.D0/2.D0 )
      I2(gg,delt,3)=I2(gg,delt,3) + CA*TF * (  - 4.D0/9.D0*nfl )
      I2(gg,delt,3) = I2(gg,delt,3) + CA**2 * ( 11.D0/9.D0 )
      I2(gg,delt,2)=I2(gg,delt,2) + CA*TF * (  - 20.D0/9.D0*nfl )
      I2(gg,delt,2) = I2(gg,delt,2) + CA**2 * ( 67.D0/9.D0 - 11*zeta2 )
      I2(gg,delt,1)=I2(gg,delt,1) + CA*TF * (  - 112.D0/27.D0*nfl + 8.D0
     &    /3.D0*zeta2*nfl )
      I2(gg,delt,1) = I2(gg,delt,1) + CA**2 * ( 404.D0/27.D0 - 14*zeta3
     &     - 22.D0/3.D0*zeta2 )
      I2(gg,delt,0)=I2(gg,delt,0) + CA*TF * (  - 328.D0/81.D0*nfl + 28.D
     &    0/9.D0*zeta3*nfl + 10.D0/3.D0*zeta2*nfl )
      I2(gg,delt,0) = I2(gg,delt,0) + CA**2 * ( 1214.D0/81.D0 + 25.D0/4.
     &    D0*zeta4 - 77.D0/9.D0*zeta3 - 67.D0/6.D0*zeta2 )

      I2(gg,plus,3)=I2(gg,plus,3) + CA**2 * (  - 4 )
      I2(gg,plus,2)=I2(gg,plus,2) + CA*TF * ( 8.D0/3.D0*nfl )
      I2(gg,plus,2) = I2(gg,plus,2) + CA**2 * (  - 22.D0/3.D0 )
      I2(gg,plus,1)=I2(gg,plus,1) + CA*TF * ( 80.D0/9.D0*nfl )
      I2(gg,plus,1) = I2(gg,plus,1) + CA**2 * (  - 268.D0/9.D0 + 12*
     &    zeta2 )
      I2(gg,plus,0)=I2(gg,plus,0) + CA*TF * ( 224.D0/27.D0*nfl )
      I2(gg,plus,0) = I2(gg,plus,0) + CA**2 * (  - 808.D0/27.D0 + 28*
     &    zeta3 )

      I2(gg,lpls,2)=I2(gg,lpls,2) + CA**2 * ( 16.D0 )

cqq
      I2(qq,rglr,2)=I2(qq,rglr,2) + CF*TF*( 2.D0 + 8.D0/3.D0*z**(-1) + 4
     &    *lnz - 2*z + 4*z*lnz - 8.D0/3.D0*z**2 )
      I2(qq,rglr,2) = I2(qq,rglr,2) + CF**2 * ( - 2.D0 - 8*lnomz + 2*lnz
     &     + 2*z - 8*z*lnomz + 2*z*lnz )
      I2(qq,rglr,1)=I2(qq,rglr,1) + CF*TF *(16.D0 - 104.D0/9.D0*z**(-1)
     &     - 4*lnz + 4*lnz**2 + 8.D0/3.D0*nfl - 24*z - 12*z*lnz + 4*z*
     &    lnz**2 - 8.D0/3.D0*z*nfl + 176.D0/9.D0*z**2 - 32.D0/3.D0*z**2
     &    *lnz )
      I2(qq,rglr,1) = I2(qq,rglr,1) + CF*CA * (  - 58.D0/3.D0 - 4*lnz
     &     + 58.D0/3.D0*z - 4*z*lnz )
      I2(qq,rglr,1) = I2(qq,rglr,1) + CF**2 * ( 28.D0 - 8*lnomz + 10*lnz
     &     + 2*lnz**2 - 28*z + 8*z*lnomz + 10*z*lnz + 2*z*lnz**2 )
      I2(qq,rglr,0)=I2(qq,rglr,0) + CF*TF * (  - 70.D0/3.D0 + 344.D0/27.
     &    D0*z**(-1) + 16.D0/3.D0*z**(-1)*Li2z + 16.D0/3.D0*z**(-1)*lnz
     &    *lnomz - 16.D0/3.D0*z**(-1)*zeta2 - 8*Li2z + 28.D0/3.D0*lnz
     &     - 8*lnz*lnomz - lnz**2 + 2.D0/3.D0*lnz**3 - 4.D0/3.D0*nfl +
     &    8*zeta2 + 62.D0/3.D0*z + 8*z*Li2z - 40.D0/3.D0*z*lnz + 8*z*
     &    lnz*lnomz - z*lnz**2 + 2.D0/3.D0*z*lnz**3 + 4.D0/3.D0*z*nfl
     &     - 8*z*zeta2 - 272.D0/27.D0*z**2 - 16.D0/3.D0*z**2*Li2z + 128.
     &    D0/9.D0*z**2*lnz - 16.D0/3.D0*z**2*lnz*lnomz - 8.D0/3.D0*z**2
     &    *lnz**2 + 16.D0/3.D0*z**2*zeta2 )
      I2(qq,rglr,0) = I2(qq,rglr,0) + CF*CA * ( 44.D0/3.D0 + 4*Li2z + 2
     &    *lnz + 4*lnz*lnomz - 6*zeta2 - 44.D0/3.D0*z - 4*z*Li2z + 2*z*
     &    lnomz + 10*z*lnz - 4*z*lnz*lnomz - 2*z*lnz**2 + 6*z*zeta2 )
      I2(qq,rglr,0) = I2(qq,rglr,0) + CF**2 *(-22.D0 - 8*Li2z + 2*lnz
     &     - 12*lnz*lnomz + 3.D0/2.D0*lnz**2 + 1.D0/3.D0*lnz**3 + 6*
     &    zeta2 + 22*z + 8*z*Li2z - 2*z*lnomz - 24*z*lnz + 12*z*lnz*
     &    lnomz + 7.D0/2.D0*z*lnz**2 + 1.D0/3.D0*z*lnz**3 - 6*z*zeta2 )

      I2(qq,rglr,3)=I2(qq,rglr,3) + tp0(qq,rglr)*CF**2 * (  - 2 )
      I2(qq,rglr,2)=I2(qq,rglr,2) + CF**2 * (  - 8*lnz*omz**(-1) )
      I2(qq,rglr,2) = I2(qq,rglr,2) + tp0(qq,rglr)*CF*TF * ( 4.D0/3.D0*
     &    nfl )
      I2(qq,rglr,2) = I2(qq,rglr,2) + tp0(qq,rglr)*CF*CA * (  - 11.D0/3.
     &    D0 )
      I2(qq,rglr,2) = I2(qq,rglr,2) + tp0(qq,rglr)*CF**2 * (  - 4*lnz )
      I2(qq,rglr,1)=I2(qq,rglr,1) + CF*TF * ( 16.D0/3.D0*nfl*lnz*
     &    omz**(-1) )
      I2(qq,rglr,1) = I2(qq,rglr,1) + CF*CA * (  - 44.D0/3.D0*lnz*
     &    omz**(-1) - 4*lnz**2*omz**(-1) )
      I2(qq,rglr,1) = I2(qq,rglr,1) + CF**2 * ( 12*lnz*omz**(-1) + 16*
     &    lnz*lnomz*omz**(-1) )
      I2(qq,rglr,1) = I2(qq,rglr,1) + tp0(qq,rglr)*CF*TF * ( 40.D0/9.D0
     &    *nfl + 8.D0/3.D0*nfl*lnz )
      I2(qq,rglr,1) = I2(qq,rglr,1) + tp0(qq,rglr)*CF*CA * (  - 134.D0/
     &    9.D0 - 22.D0/3.D0*lnz - 2*lnz**2 + 4*zeta2 )
      I2(qq,rglr,1) = I2(qq,rglr,1) + tp0(qq,rglr)*CF**2 * ( 6*lnz + 8*
     &    lnz*lnomz + 2*zeta2 )
      I2(qq,rglr,0)=I2(qq,rglr,0) + CF*TF * ( 40.D0/9.D0*nfl*lnz*
     &    omz**(-1) + 4.D0/3.D0*nfl*lnz**2*omz**(-1) )
      I2(qq,rglr,0) = I2(qq,rglr,0) + CF*CA * (  - 8*Li3omz*omz**(-1)
     &     + 16*Li3z*omz**(-1) + 8*lnomz*Li2omz*omz**(-1) - 152.D0/9.D0
     &    *lnz*omz**(-1) - 8*lnz*Li2z*omz**(-1) - 11.D0/3.D0*lnz**2*
     &    omz**(-1) - 2.D0/3.D0*lnz**3*omz**(-1) - 16*zeta3*omz**(-1) )
      I2(qq,rglr,0) = I2(qq,rglr,0) + CF**2 * ( 8*Li3omz*omz**(-1) - 40
     &    *Li3z*omz**(-1) - 8*lnomz*Li2omz*omz**(-1) + 16*lnz*omz**(-1)
     &     + 24*lnz*Li2z*omz**(-1) + 4*lnz*lnomz**2*omz**(-1) + 3*
     &    lnz**2*omz**(-1) + 4*lnz**2*lnomz*omz**(-1) + 40*zeta3*
     &    omz**(-1) )
      I2(qq,rglr,0) = I2(qq,rglr,0) + tp0(qq,rglr)*CF*TF * ( 112.D0/27.D
     &    0*nfl + 20.D0/9.D0*nfl*lnz + 2.D0/3.D0*nfl*lnz**2 )
      I2(qq,rglr,0) = I2(qq,rglr,0) + tp0(qq,rglr)*CF*CA * (  - 404.D0/
     &    27.D0 - 4*Li3omz + 8*Li3z + 4*lnomz*Li2omz - 76.D0/9.D0*lnz
     &     - 4*lnz*Li2z - 11.D0/6.D0*lnz**2 - 1.D0/3.D0*lnz**3 + 6*
     &    zeta3 )
      I2(qq,rglr,0) = I2(qq,rglr,0) + tp0(qq,rglr)*CF**2 * ( 4*Li3omz
     &     - 20*Li3z - 4*lnomz*Li2omz + 8*lnz + 12*lnz*Li2z + 2*lnz*
     &    lnomz**2 + 3.D0/2.D0*lnz**2 + 2*lnz**2*lnomz + 20*zeta3 )

      I2(qq,delt,4)=I2(qq,delt,4) + CF**2 * ( 1.D0/2.D0 )
      I2(qq,delt,3)=I2(qq,delt,3) + CF*TF * (  - 4.D0/9.D0*nfl )
      I2(qq,delt,3) = I2(qq,delt,3) + CF*CA * ( 11.D0/9.D0 )
      I2(qq,delt,2)=I2(qq,delt,2) + CF*TF * (  - 20.D0/9.D0*nfl )
      I2(qq,delt,2) = I2(qq,delt,2) + CF*CA * ( 67.D0/9.D0 - 2*zeta2 )
      I2(qq,delt,2) = I2(qq,delt,2) + CF**2 * (  - 9*zeta2 )
      I2(qq,delt,1)=I2(qq,delt,1) + CF*TF * (  - 112.D0/27.D0*nfl + 8.D0
     &    /3.D0*zeta2*nfl )
      I2(qq,delt,1) = I2(qq,delt,1) + CF*CA * ( 404.D0/27.D0 - 14*zeta3
     &     - 22.D0/3.D0*zeta2 )
      I2(qq,delt,0)=I2(qq,delt,0) + CF*TF * (  - 328.D0/81.D0*nfl + 28.D
     &    0/9.D0*zeta3*nfl + 10.D0/3.D0*zeta2*nfl )
      I2(qq,delt,0) = I2(qq,delt,0) + CF*CA * ( 1214.D0/81.D0 + 5*zeta4
     &     - 77.D0/9.D0*zeta3 - 67.D0/6.D0*zeta2 )
      I2(qq,delt,0) = I2(qq,delt,0) + CF**2 * ( 5.D0/4.D0*zeta4 )

      I2(qq,plus,3)=I2(qq,plus,3) + CF**2 * (  - 4 )
      I2(qq,plus,2)=I2(qq,plus,2) + CF*TF * ( 8.D0/3.D0*nfl )
      I2(qq,plus,2) = I2(qq,plus,2) + CF*CA * (  - 22.D0/3.D0 )
      I2(qq,plus,1)=I2(qq,plus,1) + CF*TF * ( 80.D0/9.D0*nfl )
      I2(qq,plus,1) = I2(qq,plus,1) + CF*CA * (  - 268.D0/9.D0 + 8*
     &    zeta2 )
      I2(qq,plus,1) = I2(qq,plus,1) + CF**2 * ( 4*zeta2 )
      I2(qq,plus,0)=I2(qq,plus,0) + CF*TF * ( 224.D0/27.D0*nfl )
      I2(qq,plus,0) = I2(qq,plus,0) + CF*CA * (  - 808.D0/27.D0 + 28*
     &    zeta3 )

      I2(qq,lpls,2)=I2(qq,lpls,2) + CF**2 * ( 16 )

cqg
      I2(qg,rglr,2)=I2(qg,rglr,2) + CA*TF * ( 37.D0/3.D0 + 8.D0/3.D0*
     &    z**(-1) + 4*lnz - 14.D0/3.D0*z + 16*z*lnz )
      I2(qg,rglr,2) = I2(qg,rglr,2) + CF*TF * (1.D0 + 2*lnz + 4*z - 4*z*
     &    lnz )
      I2(qg,rglr,1)=I2(qg,rglr,1) + CA*TF * (  - 88.D0/9.D0 - 104.D0/9.D
     &    0*z**(-1) + 76.D0/3.D0*lnz + 4*lnz**2 - 76.D0/9.D0*z - 176.D0/
     &    3.D0*z*lnz + 8*z*lnz**2 + 16*z*zeta2 )
      I2(qg,rglr,1) = I2(qg,rglr,1) + CF*TF *(8.D0-8*lnomz + 6*lnz + 2
     &    *lnz**2 + 10*z - 4*z*lnz**2 )
      I2(qg,rglr,0)=I2(qg,rglr,0) + CA*TF * (  - 332.D0/27.D0 + 344.D0/
     &    27.D0*z**(-1) + 16.D0/3.D0*z**(-1)*Li2z + 16.D0/3.D0*z**(-1)*
     &    lnz*lnomz - 16.D0/3.D0*z**(-1)*zeta2 - 4*Li2mz + 20.D0/3.D0*
     &    Li2z + 4*lnomz - 2*lnomz**2 - 52.D0/9.D0*lnz - 4*lnz*lnopz +
     &    20.D0/3.D0*lnz*lnomz + 19.D0/3.D0*lnz**2 + 2.D0/3.D0*lnz**3
     &     - 20.D0/3.D0*zeta2 + 178.D0/27.D0*z + 32*z*Li3z + 8.D0/3.D0*
     &    z*Li2z - 2*z*lnomz + 152.D0/9.D0*z*lnz - 16*z*lnz*Li2z + 8.D0/
     &    3.D0*z*lnz*lnomz - 32.D0/3.D0*z*lnz**2 + 4.D0/3.D0*z*lnz**3
     &     + 8*z*zeta3 + 16.D0/3.D0*z*zeta2 )
      I2(qg,rglr,0) = I2(qg,rglr,0) + CF*TF * (23.D0 - 4*lnomz + 2*
     &    lnomz**2 + 12*lnz - 4*lnz*lnomz + 5.D0/2.D0*lnz**2 + 1.D0/3.D0
     &    *lnz**3 - 6*zeta2 + 3*z + 2*z*lnomz + 7*z*lnz + 2*z*lnz**2 -
     &    2.D0/3.D0*z*lnz**3 )

      I2(qg,rglr,3)=I2(qg,rglr,3) + tp0(qg,rglr)*CF*TF * (  - 2 )
      I2(qg,rglr,2)=I2(qg,rglr,2) + tp0(qg,rglr)*CA*TF * (  - 31.D0/3.D0
     &     + 4*lnomz )
      I2(qg,rglr,2) = I2(qg,rglr,2) + tp0(qg,rglr)*CF*TF * (-8.D0 + 4*
     &    lnomz - 4*lnz )
      I2(qg,rglr,1)=I2(qg,rglr,1) + tp0(qg,rglr)*CA*TF * ( 232.D0/9.D0
     &     + 4*lnomz**2 - 88.D0/3.D0*lnz )
      I2(qg,rglr,1) = I2(qg,rglr,1) + tp0(qg,rglr)*CF*TF *(-28.D0 + 8*
     &    lnomz - 4*lnomz**2 - 8*lnz + 8*lnz*lnomz - 4*lnz**2 + 10*
     &    zeta2 )
      I2(qg,rglr,0)=I2(qg,rglr,0) + tp0(qg,rglr)*CA*TF * (  - 298.D0/27.
     &    D0 + 4*Li3omz - 44.D0/3.D0*Li2z - 4*lnomz - 4*lnomz*Li2omz +
     &    2*lnomz**2 + 2.D0/3.D0*lnomz**3 + 136.D0/9.D0*lnz - 44.D0/3.D0
     &    *lnz*lnomz - 22.D0/3.D0*lnz**2 - 4*zeta3 + 44.D0/3.D0*zeta2 )
      I2(qg,rglr,0) = I2(qg,rglr,0) + tp0(qg,rglr)*CF*TF *(-36.D0 - 4*
     &    Li3omz - 4*Li3z + 4*lnomz + 4*lnomz*Li2omz - 2*lnomz**2 - 2.D0
     &    /3.D0*lnomz**3 - 4*lnz + 4*lnz*Li2z + 4*lnz*lnomz + 2*lnz*
     &    lnomz**2 - 2*lnz**2 + 2*lnz**2*lnomz - 2.D0/3.D0*lnz**3 + 32*
     &    zeta3 + 6*zeta2 )

      I2(qg,plus,3)=I2(qg,plus,3) + tp0(qg,plus)*CF*TF * (  - 2.D0 )
      I2(qg,plus,2)=I2(qg,plus,2) + tp0(qg,plus)*CA*TF * (  - 31.D0/3.D0
     &     + 4*lnomz )
      I2(qg,plus,2) = I2(qg,plus,2) + tp0(qg,plus)*CF*TF * (-8.D0 + 4*
     &    lnomz - 4*lnz )
      I2(qg,plus,1)=I2(qg,plus,1) + tp0(qg,plus)*CA*TF * ( 232.D0/9.D0
     &     + 4*lnomz**2 - 88.D0/3.D0*lnz )
      I2(qg,plus,1) = I2(qg,plus,1) + tp0(qg,plus)*CF*TF *(-28.D0 + 8*
     &    lnomz - 4*lnomz**2 - 8*lnz + 8*lnz*lnomz - 4*lnz**2 + 10*
     &    zeta2 )
      I2(qg,plus,0)=I2(qg,plus,0) + tp0(qg,plus)*CA*TF * (  - 298.D0/27.
     &    D0 + 4*Li3omz - 44.D0/3.D0*Li2z - 4*lnomz - 4*lnomz*Li2omz +
     &    2*lnomz**2 + 2.D0/3.D0*lnomz**3 + 136.D0/9.D0*lnz - 44.D0/3.D0
     &    *lnz*lnomz - 22.D0/3.D0*lnz**2 - 4*zeta3 + 44.D0/3.D0*zeta2 )
      I2(qg,plus,0) = I2(qg,plus,0) + tp0(qg,plus)*CF*TF *(-36.D0 - 4*
     &    Li3omz - 4*Li3z + 4*lnomz + 4*lnomz*Li2omz - 2*lnomz**2 - 2.D0
     &    /3.D0*lnomz**3 - 4*lnz + 4*lnz*Li2z + 4*lnz*lnomz + 2*lnz*
     &    lnomz**2 - 2*lnz**2 + 2*lnz**2*lnomz - 2.D0/3.D0*lnz**3 + 32*
     &    zeta3 + 6*zeta2 )

      I2(qg,rglr,1)=I2(qg,rglr,1) + mp0(qg,rglr)*CA*TF * ( 8*Li2mz + 8*
     &    lnz*lnopz )
      I2(qg,rglr,0)=I2(qg,rglr,0) + mp0(qg,rglr)*CA*TF * (  - 8*
     &    Li3zonopz + 4*Li3mz + 4*Li2mz + 4.D0/3.D0*lnopz**3 - 4*lnz*
     &    Li2mz + 4*lnz*lnopz - 4*lnz*lnopz**2 + 2*lnz**2*lnopz - 4*
     &    zeta2*lnopz )



cgq
      I2(gq,rglr,2)=I2(gq,rglr,2) + CF*CA * (  - 14.D0/3.D0 - 16*lnz +
     &    43.D0/3.D0*z - 4*z*lnz + 8.D0/3.D0*z**2 )
      I2(gq,rglr,2) = I2(gq,rglr,2) + CF**2*(4.D0 + 4*lnz - z - 2*z*lnz
     &     )
      I2(gq,rglr,1)=I2(gq,rglr,1) + CF*CA * (  - 112.D0/9.D0 + 48*lnz
     &     - 8*lnz**2 - 16*zeta2 + 2.D0/9.D0*z - 8*z*lnomz + 20*z*lnz
     &     - 4*z*lnz**2 - 176.D0/9.D0*z**2 + 32.D0/3.D0*z**2*lnz )
      I2(gq,rglr,1) = I2(gq,rglr,1) + CF**2*(6.D0 - 8*lnz + 4*lnz**2 +
     &    12*z - 10*z*lnz - 2*z*lnz**2 )
      I2(gq,rglr,0)=I2(gq,rglr,0) + CF*TF * (  - 40.D0/9.D0*z*nfl - 8.D0
     &    /3.D0*z*nfl*lnomz )
      I2(gq,rglr,0) = I2(gq,rglr,0) + CF*CA * ( 4.D0/27.D0 - 32*Li3z +
     &    8.D0/3.D0*Li2z - 166.D0/3.D0*lnz + 16*lnz*Li2z + 8.D0/3.D0*
     &    lnz*lnomz + 12*lnz**2 - 4.D0/3.D0*lnz**3 - 8*zeta3 - 8.D0/3.D0
     &    *zeta2 + 508.D0/27.D0*z + 4*z*Li2mz + 20.D0/3.D0*z*Li2z + 22.D
     &    0/3.D0*z*lnomz + 2*z*lnomz**2 + 4.D0/3.D0*z*lnz + 4*z*lnz*
     &    lnopz + 8.D0/3.D0*z*lnz*lnomz + 3*z*lnz**2 - 2.D0/3.D0*z*
     &    lnz**3 - 26.D0/3.D0*z*zeta2 + 608.D0/27.D0*z**2 + 16.D0/3.D0*
     &    z**2*Li2z - 176.D0/9.D0*z**2*lnz + 16.D0/3.D0*z**2*lnz*lnomz
     &     + 8.D0/3.D0*z**2*lnz**2 - 16.D0/3.D0*z**2*zeta2 )
      I2(gq,rglr,0) = I2(gq,rglr,0) + CF**2 * (10.D0 - 15*lnz - 2*lnz**2
     &     + 2.D0/3.D0*lnz**3 - z - 6*z*lnomz - 2*z*lnomz**2 + 5*z*lnz
     &     - 3.D0/2.D0*z*lnz**2 - 1.D0/3.D0*z*lnz**3 )

      I2(gq,rglr,3)=I2(gq,rglr,3) + tp0(gq,rglr)*CF*CA * (  - 2.D0 )
      I2(gq,rglr,2)=I2(gq,rglr,2) + tp0(gq,rglr)*CF*TF * ( 8.D0/3.D0*
     &    nfl )
      I2(gq,rglr,2) = I2(gq,rglr,2) + tp0(gq,rglr)*CF*CA * (  - 53.D0/3.
     &    D0 + 4*lnomz - 4*lnz )
      I2(gq,rglr,2) = I2(gq,rglr,2) + tp0(gq,rglr)*CF**2 * ( 4*lnomz )
      I2(gq,rglr,1)=I2(gq,rglr,1) + tp0(gq,rglr)*CF*TF * ( 80.D0/9.D0*
     &    nfl + 16.D0/3.D0*nfl*lnomz )
      I2(gq,rglr,1) = I2(gq,rglr,1) + tp0(gq,rglr)*CF*CA *(-2.D0 - 44.D
     &    0/3.D0*lnomz - 4*lnomz**2 + 8*lnz*lnomz + 2*zeta2 )
      I2(gq,rglr,1) = I2(gq,rglr,1) + tp0(gq,rglr)*CF**2 * ( 12*lnomz
     &     + 4*lnomz**2 )
      I2(gq,rglr,0)=I2(gq,rglr,0) + tp0(gq,rglr)*CF*TF * ( 224.D0/27.D0
     &    *nfl + 40.D0/9.D0*nfl*lnomz + 4.D0/3.D0*nfl*lnomz**2 )
      I2(gq,rglr,0) = I2(gq,rglr,0) + tp0(gq,rglr)*CF*CA * (  - 1580.D0/
     &    27.D0 - 20*Li3z - 44.D0/3.D0*Li2z - 152.D0/9.D0*lnomz - 11.D0/
     &    3.D0*lnomz**2 - 2.D0/3.D0*lnomz**3 + 12*lnz*Li2z - 44.D0/3.D0
     &    *lnz*lnomz + 2*lnz*lnomz**2 + 2*lnz**2*lnomz + 24*zeta3 + 44.D
     &    0/3.D0*zeta2 )
      I2(gq,rglr,0) = I2(gq,rglr,0) + tp0(gq,rglr)*CF**2 * ( 16*lnomz
     &     + 3*lnomz**2 + 2.D0/3.D0*lnomz**3 )

      I2(gq,plus,3)=I2(gq,plus,3) + tp0(gq,plus)*CF*TF * ( 8.D0/3.D0*
     &    nfl )
      I2(gq,plus,3) = I2(gq,plus,3) + tp0(gq,plus)*CF*CA * (  - 53.D0/3.
     &    D0 + 4*lnomz - 4*lnz )
      I2(gq,plus,3) = I2(gq,plus,3) + tp0(gq,plus)*CF**2 * ( 4*lnomz )
      I2(gq,plus,2)=I2(gq,plus,2) + tp0(gq,plus)*CF*TF * ( 8.D0/3.D0*
     &    nfl )
      I2(gq,plus,2) = I2(gq,plus,2) + tp0(gq,plus)*CF*CA * (  - 53.D0/3.
     &    D0 + 4*lnomz - 4*lnz )
      I2(gq,plus,2) = I2(gq,plus,2) + tp0(gq,plus)*CF**2 * ( 4*lnomz )
      I2(gq,plus,1)=I2(gq,plus,1) + tp0(gq,plus)*CF*TF * ( 80.D0/9.D0*
     &    nfl + 16.D0/3.D0*nfl*lnomz )
      I2(gq,plus,1) = I2(gq,plus,1) + tp0(gq,plus)*CF*CA *(-2.D0 - 44.D
     &    0/3.D0*lnomz - 4*lnomz**2 + 8*lnz*lnomz + 2*zeta2 )
      I2(gq,plus,1) = I2(gq,plus,1) + tp0(gq,plus)*CF**2 * ( 12*lnomz
     &     + 4*lnomz**2 )
      I2(gq,plus,0)=I2(gq,plus,0) + tp0(gq,plus)*CF*TF * ( 224.D0/27.D0
     &    *nfl + 40.D0/9.D0*nfl*lnomz + 4.D0/3.D0*nfl*lnomz**2 )
      I2(gq,plus,0) = I2(gq,plus,0) + tp0(gq,plus)*CF*CA * (  - 1580.D0/
     &    27.D0 - 20*Li3z - 44.D0/3.D0*Li2z - 152.D0/9.D0*lnomz - 11.D0/
     &    3.D0*lnomz**2 - 2.D0/3.D0*lnomz**3 + 12*lnz*Li2z - 44.D0/3.D0
     &    *lnz*lnomz + 2*lnz*lnomz**2 + 2*lnz**2*lnomz + 24*zeta3 + 44.D
     &    0/3.D0*zeta2 )
      I2(gq,plus,0) = I2(gq,plus,0) + tp0(gq,plus)*CF**2 * ( 16*lnomz
     &     + 3*lnomz**2 + 2.D0/3.D0*lnomz**3 )

      I2(gq,rglr,1)=I2(gq,rglr,1) + mp0(gq,rglr)*CF*CA * ( 8*Li2mz + 8*
     &    lnz*lnopz )
      I2(gq,rglr,0)=I2(gq,rglr,0) + mp0(gq,rglr)*CF*CA * (  - 8*
     &    Li3zonopz + 4*Li3mz + 4.D0/3.D0*lnopz**3 - 4*lnz*Li2mz - 4*
     &    lnz*lnopz**2 + 2*lnz**2*lnopz - 4*zeta2*lnopz )


cqbq
      I2(qbq,rglr,2)=I2(qbq,rglr,2) + CF*TF * (2.D0 +8.D0/3.D0*z**(-1)
     &     + 4*lnz - 2*z + 4*z*lnz - 8.D0/3.D0*z**2 )
      I2(qbq,rglr,1)=I2(qbq,rglr,1) + CF*TF * ( 16.D0 - 104.D0/9.D0*
     &    z**(-1) - 4*lnz + 4*lnz**2 - 24*z - 12*z*lnz + 4*z*lnz**2 +
     &    176.D0/9.D0*z**2 - 32.D0/3.D0*z**2*lnz )
      I2(qbq,rglr,1) = I2(qbq,rglr,1) + CF*CA*(8.D0 + 4*lnz - 8*z + 4*z
     &    *lnz )
      I2(qbq,rglr,1) = I2(qbq,rglr,1) + CF**2 *(-16.D0 - 8*lnz + 16*z
     &     - 8*z*lnz )
      I2(qbq,rglr,0)=I2(qbq,rglr,0) + CF*TF * (  - 70.D0/3.D0 + 344.D0/
     &    27.D0*z**(-1) + 16.D0/3.D0*z**(-1)*Li2z + 16.D0/3.D0*z**(-1)*
     &    lnz*lnomz - 16.D0/3.D0*z**(-1)*zeta2 - 8*Li2z + 28.D0/3.D0*
     &    lnz - 8*lnz*lnomz - lnz**2 + 2.D0/3.D0*lnz**3 + 8*zeta2 + 62.D
     &    0/3.D0*z + 8*z*Li2z - 40.D0/3.D0*z*lnz + 8*z*lnz*lnomz - z*
     &    lnz**2 + 2.D0/3.D0*z*lnz**3 - 8*z*zeta2 - 272.D0/27.D0*z**2
     &     - 16.D0/3.D0*z**2*Li2z + 128.D0/9.D0*z**2*lnz - 16.D0/3.D0*
     &    z**2*lnz*lnomz - 8.D0/3.D0*z**2*lnz**2 + 16.D0/3.D0*z**2*
     &    zeta2 )
      I2(qbq,rglr,0) = I2(qbq,rglr,0) + CF*CA * (-15.D0 + 4*Li2mz - 4*
     &    Li2z - 3*lnz + 4*lnz*lnopz - 4*lnz*lnomz + 6*zeta2 + 15*z + 4
     &    *z*Li2mz + 4*z*Li2z - 11*z*lnz + 4*z*lnz*lnopz + 4*z*lnz*
     &    lnomz - 2*z*zeta2 )
      I2(qbq,rglr,0) = I2(qbq,rglr,0) + CF**2*(30.D0 - 8*Li2mz + 8*Li2z
     &     + 6*lnz - 8*lnz*lnopz + 8*lnz*lnomz - 12*zeta2 - 30*z - 8*z*
     &    Li2mz - 8*z*Li2z + 22*z*lnz - 8*z*lnz*lnopz - 8*z*lnz*lnomz
     &     + 4*z*zeta2 )

      I2(qbq,rglr,1)=I2(qbq,rglr,1) + mp0(qq,rglr)*CF*CA * (  - 8*Li2mz
     &     - 8*lnz*lnopz + 2*lnz**2 - 4*zeta2 )
      I2(qbq,rglr,1) = I2(qbq,rglr,1) + mp0(qq,rglr)*CF**2 * ( 16*Li2mz
     &     + 16*lnz*lnopz - 4*lnz**2 + 8*zeta2 )
      I2(qbq,rglr,0)=I2(qbq,rglr,0) + mp0(qq,rglr)*CF*CA * ( 8*
     &    Li3zonopz - 4*Li3mz - 8*Li3z - 4.D0/3.D0*lnopz**3 + 4*lnz*
     &    Li2mz + 4*lnz*Li2z + 4*lnz*lnopz**2 - 2*lnz**2*lnopz + 1.D0/3.
     &    D0*lnz**3 - 2*zeta3 + 4*zeta2*lnopz )
      I2(qbq,rglr,0) = I2(qbq,rglr,0) + mp0(qq,rglr)*CF**2 * (  - 16*
     &    Li3zonopz + 8*Li3mz + 16*Li3z + 8.D0/3.D0*lnopz**3 - 8*lnz*
     &    Li2mz - 8*lnz*Li2z - 8*lnz*lnopz**2 + 4*lnz**2*lnopz - 2.D0/3.
     &    D0*lnz**3 + 4*zeta3 - 8*zeta2*lnopz )

cqpq
      I2(qpq,rglr,2)=I2(qpq,rglr,2) + CF*TF * ( 2.D0 + 8.D0/3.D0*z**(-1)
     &     + 4*lnz - 2*z + 4*z*lnz - 8.D0/3.D0*z**2 )
      I2(qpq,rglr,1)=I2(qpq,rglr,1) + CF*TF * ( 16.D0 - 104.D0/9.D0*
     &    z**(-1) - 4*lnz + 4*lnz**2 - 24*z - 12*z*lnz + 4*z*lnz**2 +
     &    176.D0/9.D0*z**2 - 32.D0/3.D0*z**2*lnz )
      I2(qpq,rglr,0)=I2(qpq,rglr,0) + CF*TF * (  - 70.D0/3.D0 + 344.D0/
     &    27.D0*z**(-1) + 16.D0/3.D0*z**(-1)*Li2z + 16.D0/3.D0*z**(-1)*
     &    lnz*lnomz - 16.D0/3.D0*z**(-1)*zeta2 - 8*Li2z + 28.D0/3.D0*
     &    lnz - 8*lnz*lnomz - lnz**2 + 2.D0/3.D0*lnz**3 + 8*zeta2 + 62.D
     &    0/3.D0*z + 8*z*Li2z - 40.D0/3.D0*z*lnz + 8*z*lnz*lnomz - z*
     &    lnz**2 + 2.D0/3.D0*z*lnz**3 - 8*z*zeta2 - 272.D0/27.D0*z**2
     &     - 16.D0/3.D0*z**2*Li2z + 128.D0/9.D0*z**2*lnz - 16.D0/3.D0*
     &    z**2*lnz*lnomz - 8.D0/3.D0*z**2*lnz**2 + 16.D0/3.D0*z**2*
     &    zeta2 )


cqbpq
      I2(qbpq,rglr,2)=I2(qbpq,rglr,2)+CF*TF*(2.D0+ 8.D0/3.D0*z**(-1)
     &     + 4*lnz - 2*z + 4*z*lnz - 8.D0/3.D0*z**2 )
      I2(qbpq,rglr,1)=I2(qbpq,rglr,1) + CF*TF * ( 16.D0 - 104.D0/9.D0*
     &    z**(-1) - 4*lnz + 4*lnz**2 - 24*z - 12*z*lnz + 4*z*lnz**2 +
     &    176.D0/9.D0*z**2 - 32.D0/3.D0*z**2*lnz )
      I2(qbpq,rglr,0)=I2(qbpq,rglr,0) + CF*TF * (  - 70.D0/3.D0 + 344.D0
     &    /27.D0*z**(-1) + 16.D0/3.D0*z**(-1)*Li2z + 16.D0/3.D0*z**(-1)
     &    *lnz*lnomz - 16.D0/3.D0*z**(-1)*zeta2 - 8*Li2z + 28.D0/3.D0*
     &    lnz - 8*lnz*lnomz - lnz**2 + 2.D0/3.D0*lnz**3 + 8*zeta2 + 62.D
     &    0/3.D0*z + 8*z*Li2z - 40.D0/3.D0*z*lnz + 8*z*lnz*lnomz - z*
     &    lnz**2 + 2.D0/3.D0*z*lnz**3 - 8*z*zeta2 - 272.D0/27.D0*z**2
     &     - 16.D0/3.D0*z**2*Li2z + 128.D0/9.D0*z**2*lnz - 16.D0/3.D0*
     &    z**2*lnz*lnomz - 8.D0/3.D0*z**2*lnz**2 + 16.D0/3.D0*z**2*
     &    zeta2 )

      return
      end
