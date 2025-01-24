!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine DeltaI2R(x,LQ,R,DeltaI2)
          use ptveto, only: xlim
!     implementation of DeltaI piece of Eq.(3.2)
!     of https://arxiv.org/abs/2207.07037v1
!     written by writeDeltaI2R.frm
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      include 'transitionlabels.f'
!     dmin=-1,dmax=2,delt=-1,plus=0,lpls=1,rglr=2
      include 'distributions.f'
      real(dp):: x,LQ,R,DeltaI2(0:6,dmin:dmax)
!      DeltaI2, index one transition flavor
!      DeltaI2, index2,distribution type
      include 'DelI2bit.f'
      real(dp)::BoundaryConditionGGCACA,
     & BoundaryConditionGGCATF,BoundaryConditionGGCFTF,
     & BoundaryConditionQQCACF,
     & BoundaryConditionQQCFCF,BoundaryConditionQQCFTF,
     & BoundaryConditionGQCFCA,BoundaryConditionGQCFCF,
     & BoundaryConditionQGCATF,BoundaryConditionQGCFTF,
     & BoundaryConditionQQBARS,BoundaryConditionQQBARNS

      real(dp)::omx,opx,opxon2,lnx,lnopxon2,ddilog,Li2omx,Li2mx,
     & lnR,LnRonpi,R10,R8,R6,R4,R2
      real(dp),parameter::ln2=0.6931471805599453_dp,
     & ln512=6.238324625039508_dp,ln4096=8.317766166719343_dp,
     & lnpi=1.1447298858494_dp
      integer,parameter:: GGCACA=1,GGCATF=2,GGCFTF=3,
     &                    QQCACF=4,QQCFCF=5,QQCFTF=6,
     & GQCFCA=7,GQCFCF=8,QGCATF=9,QGCFTF=10,QQBARS=11,QQBARNS=12
      integer::nf

      omx=one-x
      opx=one+x
      opxon2=0.5d0*opx
      lnx=log(x)
      lnopxon2=log(opxon2)
      Li2omx=ddilog(omx)
      Li2mx=ddilog(-x)
      R2=R**2
      R4=R2*R2
      R6=R4*R2
      R8=R6*R2
      R10=R8*R2
      lnR=log(R)
      lnRonpi=lnR-lnpi
      nf=nfl

c     initialize arrays to zero
      DeltaI2(:,:)=0
      DelI2bit(:,:)=0

      if (x > xlim) then
      include 'GGCACA.f'
      include 'GGCATF.f'
      include 'GGCFTF.f'
      include 'QQCACF.f'
      include 'QQCFCF.f'
      include 'QQCFTF.f'
      include 'GQCFCA.f'
      include 'GQCFCF.f'
      include 'QGCATF.f'
      include 'QGCFTF.f'
      include 'QQBARS.f'
      include 'QQBARNS.f'
      else
      DelI2bit(GGCACA,rglr) =  + CA**2*omx**(-7)*R6 * (  - 160.D0/3.D0*
     &    Li2mx + 296.D0/9.D0*lnx + 160.D0/3.D0*Li2omx - 40.D0/9.D0*
     &    pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-6)*
     & R6 * (  - 184.D0/9.D0 + 880.D0/3.D0*Li2mx - 1628.D0/9.D0*lnx + 
     &    80.D0/3.D0*lnopxon2 + 160.D0/3.D0*ln2 - 880.D0/3.D0*Li2omx + 
     &    220.D0/9.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-5)*
     & R4 * (  - 24*Li2mx + 14*lnx + 24*Li2omx - 2*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-5)*
     & R6 * ( 1040.D0/9.D0 - 2140.D0/3.D0*Li2mx + 3977.D0/9.D0*lnx - 
     &    376.D0/3.D0*lnopxon2 - 800.D0/3.D0*ln2 + 2140.D0/3.D0*Li2omx
     &     - 535.D0/9.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-4)*
     & R4 * (  - 10 + 96*Li2mx - 56*lnx + 12*lnopxon2 + 24*ln2 - 96*
     &    Li2omx + 8*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-4)*
     & R6 * (  - 7577.D0/27.D0 + 1008*Li2mx - 1892.D0/3.D0*lnx + 2306.D0
     &    /9.D0*lnopxon2 + 5260.D0/9.D0*ln2 - 1008*Li2omx + 84*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-3)*
     & R2 * (  - 16*Li2mx + 8*lnx + 16*Li2omx - 4.D0/3.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-3)*
     & R4 * ( 41 - 482.D0/3.D0*Li2mx + 289.D0/3.D0*lnx - 38*lnopxon2 - 
     &    84*ln2 + 482.D0/3.D0*Li2omx - 241.D0/18.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-3)*
     & R6 * ( 20767.D0/54.D0 - 2704.D0/3.D0*Li2mx + 1736.D0/3.D0*lnx - 
     &    2693.D0/9.D0*lnopxon2 - 6622.D0/9.D0*ln2 + 2704.D0/3.D0*
     &    Li2omx - 676.D0/9.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-2)*
     & R2 * (  - 8 + 40*Li2mx - 24*lnx + 8*lnopxon2 + 16*ln2 - 40*
     &    Li2omx + 10.D0/3.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-2)*
     & R4 * (  - 389.D0/6.D0 + 436.D0/3.D0*Li2mx - 93*lnx + 145.D0/3.D0
     &    *lnopxon2 + 362.D0/3.D0*ln2 - 436.D0/3.D0*Li2omx + 109.D0/9.D0
     &    *pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-2)*
     & R6 * (  - 17161.D0/54.D0 + 1528.D0/3.D0*Li2mx - 1070.D0/3.D0*lnx
     &     + 1928.D0/9.D0*lnopxon2 + 5152.D0/9.D0*ln2 - 1528.D0/3.D0*
     &    Li2omx + 382.D0/9.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-1)*
     & R2 * ( 16 - 24*Li2mx + 40*lnx - 16*lnopxon2 - 32*ln2 + 24*Li2omx
     &     - 2*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-1)*
     & R4 * ( 601.D0/12.D0 - 58*Li2mx + 125.D0/3.D0*lnx - 104.D0/3.D0*
     &    lnopxon2 - 91*ln2 + 58*Li2omx - 29.D0/6.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*omx**(-1)*
     & R6 * ( 853.D0/6.D0 - 164*Li2mx + 1240.D0/9.D0*lnx - 847.D0/9.D0*
     &    lnopxon2 - 2395.D0/9.D0*ln2 + 164*Li2omx - 41.D0/3.D0*pisq )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2 * (  - 26.D0
     &    /9.D0 - 524.D0/9.D0*x**(-1)*ln2 + 176.D0/3.D0*x**(-1)*ln2**2
     &     + 16.D0/3.D0*x**(-1)*pisq*ln2 + 1072.D0/9.D0*ln2 - 352.D0/3.D
     &    0*ln2**2 - 32.D0/3.D0*pisq*ln2 - 524.D0/9.D0*x*ln2 + 176.D0/3.
     &    D0*x*ln2**2 + 16.D0/3.D0*x*pisq*ln2 + 524.D0/9.D0*x**2*ln2 - 
     &    176.D0/3.D0*x**2*ln2**2 - 16.D0/3.D0*x**2*pisq*ln2 )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*R2 * (  - 
     &    94139.D0/2700.D0 + 289.D0/12.D0*x**(-1) - 22*x**(-1)*ln2 - 2*
     &    x**(-1)*pisq - 8*opx**(-1)*Li2mx + 8*opx**(-1)*Li2omx - 16*
     &    Li2mx - 24*lnx + 16*lnopxon2 + 2092.D0/45.D0*ln2 + 16*Li2omx
     &     - 2.D0/3.D0*pisq*opx**(-1) + 8.D0/3.D0*pisq + 54383.D0/1800.D
     &    0*x + 16*x*Li2mx - 8*x*lnx - 16*x*lnopxon2 - 662.D0/15.D0*x*
     &    ln2 - 16*x*Li2omx - 2.D0/3.D0*x*pisq + 989.D0/150.D0*x**2 - 
     &    24*x**2*Li2mx - 24*x**2*lnx + 24*x**2*lnopxon2 + 76.D0/5.D0*
     &    x**2*ln2 + 24*x**2*Li2omx )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*R4 * ( 
     &    2312371.D0/259200.D0 + 42959.D0/51840.D0*x**(-1) - 11.D0/108.D
     &    0*x**(-1)*ln2 - 1.D0/108.D0*x**(-1)*pisq - 49.D0/8.D0*
     &    opx**(-2) - 95.D0/8.D0*opx**(-1) + 2.D0/3.D0*Li2mx - 3*lnx + 
     &    56.D0/3.D0*lnopxon2 + 3566.D0/135.D0*ln2 - 2.D0/3.D0*Li2omx
     &     + 2.D0/27.D0*pisq - 109304063.D0/12700800.D0*x + 26.D0/3.D0*
     &    x*Li2mx + 53.D0/3.D0*x*lnx - 18*x*lnopxon2 - 96871.D0/3780.D0
     &    *x*ln2 - 26.D0/3.D0*x*Li2omx + 77.D0/108.D0*x*pisq - 86103083.
     &    D0/6350400.D0*x**2 + 8*x**2*Li2mx - 25.D0/3.D0*x**2*lnx + 25.D
     &    0/3.D0*x**2*lnopxon2 + 82213.D0/3780.D0*x**2*ln2 - 8*x**2*
     &    Li2omx + 73.D0/108.D0*x**2*pisq + 117731.D0/117600.D0*x**3 - 
     &    38.D0/35.D0*x**3*ln2 )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*R6 * (  - 
     &    8826163921.D0/435456000.D0 + 83.D0/4608.D0*x**(-1) - 61.D0/48.
     &    D0*opx**(-4) + 37.D0/48.D0*opx**(-3) - 223.D0/192.D0*
     &    opx**(-2) - 7.D0/6.D0*opx**(-1) + 64.D0/3.D0*Li2mx - 23*lnx
     &     + 65.D0/3.D0*lnopxon2 + 156641.D0/2700.D0*ln2 - 64.D0/3.D0*
     &    Li2omx + 16.D0/9.D0*pisq + 63601373.D0/12192768.D0*x + 4.D0/3.
     &    D0*x*Li2mx + 83.D0/9.D0*x*lnx - 35.D0/9.D0*x*lnopxon2 - 902.D0
     &    /189.D0*x*ln2 - 4.D0/3.D0*x*Li2omx + 1.D0/9.D0*x*pisq - 
     &    1553655541.D0/152409600.D0*x**2 + 20.D0/3.D0*x**2*Li2mx + 20.D
     &    0/3.D0*x**2*lnx - 20.D0/3.D0*x**2*lnopxon2 + 3121.D0/945.D0*
     &    x**2*ln2 - 20.D0/3.D0*x**2*Li2omx + 5.D0/9.D0*x**2*pisq - 
     &    186397847.D0/76204800.D0*x**3 + 2524.D0/945.D0*x**3*ln2 + 
     &    539159.D0/3810240.D0*x**4 - 25.D0/189.D0*x**4*ln2 )
      DelI2bit(GGCACA,rglr) = DelI2bit(GGCACA,rglr) + CA**2*lnRonpi
     &  * (  - 1048.D0/9.D0 + 524.D0/9.D0*x**(-1) - 176.D0/3.D0*x**(-1)
     &    *ln2 - 16.D0/3.D0*x**(-1)*pisq + 352.D0/3.D0*ln2 + 32.D0/3.D0
     &    *pisq + 524.D0/9.D0*x - 176.D0/3.D0*x*ln2 - 16.D0/3.D0*x*pisq
     &     - 524.D0/9.D0*x**2 + 176.D0/3.D0*x**2*ln2 + 16.D0/3.D0*x**2*
     &    pisq )

      DelI2bit(GGCATF,rglr) =  + CA*TF*nf * ( 52.D0/9.D0 + 184.D0/9.D0*
     &    x**(-1)*ln2 - 64.D0/3.D0*x**(-1)*ln2**2 - 416.D0/9.D0*ln2 + 
     &    128.D0/3.D0*ln2**2 + 184.D0/9.D0*x*ln2 - 64.D0/3.D0*x*ln2**2
     &     - 184.D0/9.D0*x**2*ln2 + 64.D0/3.D0*x**2*ln2**2 )
      DelI2bit(GGCATF,rglr) = DelI2bit(GGCATF,rglr) + CA*TF*nf*R2 * ( 
     &     - 6811.D0/1350.D0 + 55.D0/18.D0*x**(-1) - 8.D0/3.D0*x**(-1)*
     &    ln2 + 196.D0/45.D0*ln2 + 3517.D0/900.D0*x - 56.D0/15.D0*x*ln2
     &     - 742.D0/225.D0*x**2 + 44.D0/15.D0*x**2*ln2 )
      DelI2bit(GGCATF,rglr) = DelI2bit(GGCATF,rglr) + CA*TF*nf*R4 * ( 
     &    12113.D0/129600.D0 - 179.D0/1620.D0*x**(-1) + 1.D0/27.D0*
     &    x**(-1)*ln2 - 4.D0/135.D0*ln2 + 210359.D0/1587600.D0*x - 169.D
     &    0/945.D0*x*ln2 - 214691.D0/6350400.D0*x**2 + 64.D0/945.D0*
     &    x**2*ln2 + 226.D0/2205.D0*x**3 - 2.D0/21.D0*x**3*ln2 )
      DelI2bit(GGCATF,rglr) = DelI2bit(GGCATF,rglr) + CA*TF*nf*R6 * ( 
     &     - 3559313.D0/1524096000.D0 + 5.D0/6912.D0*x**(-1) + 73.D0/
     &    9450.D0*ln2 + 2091127.D0/152409600.D0*x - 22.D0/945.D0*x*ln2
     &     - 5000369.D0/76204800.D0*x**2 + 58.D0/945.D0*x**2*ln2 + 
     &    2156093.D0/38102400.D0*x**3 - 47.D0/945.D0*x**3*ln2 - 249667.D
     &    0/9525600.D0*x**4 + 22.D0/945.D0*x**4*ln2 )
      DelI2bit(GGCATF,rglr) = DelI2bit(GGCATF,rglr) + CA*TF*nf*R8 * ( 
     &     - 194032051.D0/1609445376000.D0 + 384581.D0/4180377600.D0*
     &    x**(-1) + 1.D0/194400.D0*x**(-1)*ln2 + 1.D0/136080.D0*ln2 + 
     &    46395565501.D0/17703899136000.D0*x - 50587.D0/14968800.D0*x*
     &    ln2 - 41813973583.D0/4425974784000.D0*x**2 + 161401.D0/
     &    14968800.D0*x**2*ln2 + 1473964069.D0/73766246400.D0*x**3 - 
     &    263.D0/13860.D0*x**3*ln2 - 2403887.D0/140507136.D0*x**4 + 1.D0
     &    /66.D0*x**4*ln2 + 77879803.D0/12294374400.D0*x**5 - 19.D0/
     &    3465.D0*x**5*ln2 )
      DelI2bit(GGCATF,rglr) = DelI2bit(GGCATF,rglr) + CA*TF*nf*lnRonpi
     &  * ( 368.D0/9.D0 - 184.D0/9.D0*x**(-1) + 64.D0/3.D0*x**(-1)*ln2
     &     - 128.D0/3.D0*ln2 - 184.D0/9.D0*x + 64.D0/3.D0*x*ln2 + 184.D0
     &    /9.D0*x**2 - 64.D0/3.D0*x**2*ln2 )

      DelI2bit(GGCFTF,rglr) =  + CF*TF*nf*R2 * ( 184.D0/9.D0 - 92.D0/9.D
     &    0*x**(-1) + 32.D0/3.D0*x**(-1)*ln2 - 64.D0/3.D0*ln2 - 46.D0/3.
     &    D0*x + 16*x*ln2 + 46.D0/9.D0*x**2 - 16.D0/3.D0*x**2*ln2 )
      DelI2bit(GGCFTF,rglr) = DelI2bit(GGCFTF,rglr) + CF*TF*nf*R4 * ( 
     &    23.D0/50.D0 + 1.D0/4.D0*x**(-1) - 4.D0/5.D0*ln2 - 3217.D0/
     &    1800.D0*x + 28.D0/15.D0*x*ln2 + 1939.D0/1800.D0*x**2 - 16.D0/
     &    15.D0*x**2*ln2 )
      DelI2bit(GGCFTF,rglr) = DelI2bit(GGCFTF,rglr) + CF*TF*nf*R6 * ( 
     &     - 437.D0/105840.D0 - 2.D0/63.D0*ln2 - 14059.D0/66150.D0*x + 
     &    82.D0/315.D0*x*ln2 + 56029.D0/105840.D0*x**2 - 32.D0/63.D0*
     &    x**2*ln2 - 97607.D0/264600.D0*x**3 + 104.D0/315.D0*x**3*ln2
     &     + 14863.D0/264600.D0*x**4 - 16.D0/315.D0*x**4*ln2 )
      DelI2bit(GGCFTF,rglr) = DelI2bit(GGCFTF,rglr) + CF*TF*nf*R8 * ( 
     &    9511.D0/6350400.D0 - 1.D0/3072.D0*x**(-1) - 1.D0/630.D0*ln2
     &     - 3359477.D0/203212800.D0*x + 43.D0/1890.D0*x*ln2 + 5749721.D
     &    0/67737600.D0*x**2 - 86.D0/945.D0*x**2*ln2 - 1870111.D0/
     &    12700800.D0*x**3 + 44.D0/315.D0*x**3*ln2 + 811177.D0/8467200.D
     &    0*x**4 - 3.D0/35.D0*x**4*ln2 - 91859.D0/5080320.D0*x**5 + 1.D0
     &    /63.D0*x**5*ln2 )

      DelI2bit(QQCACF,rglr) =  + CA*CF*omx**(-9)*R8 * ( 140*Li2mx - 533.
     &    D0/6.D0*lnx - 140*Li2omx + 35.D0/3.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-8)*
     & R8 * ( 307.D0/6.D0 - 980*Li2mx + 3731.D0/6.D0*lnx - 70*lnopxon2
     &     - 140*ln2 + 980*Li2omx - 245.D0/3.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-7)*
     & R6 * ( 160.D0/3.D0*Li2mx - 296.D0/9.D0*lnx - 160.D0/3.D0*Li2omx
     &     + 40.D0/9.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-7)*
     & R8 * (  - 4411.D0/12.D0 + 3060*Li2mx - 11647.D0/6.D0*lnx + 435*
     &    lnopxon2 + 910*ln2 - 3060*Li2omx + 255*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-6)*
     & R6 * ( 184.D0/9.D0 - 880.D0/3.D0*Li2mx + 1628.D0/9.D0*lnx - 80.D0
     &    /3.D0*lnopxon2 - 160.D0/3.D0*ln2 + 880.D0/3.D0*Li2omx - 220.D0
     &    /9.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-6)*
     & R8 * ( 20945.D0/18.D0 - 16840.D0/3.D0*Li2mx + 64063.D0/18.D0*lnx
     &     - 3565.D0/3.D0*lnopxon2 - 7850.D0/3.D0*ln2 + 16840.D0/3.D0*
     &    Li2omx - 4210.D0/9.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-5)*
     & R4 * ( 24*Li2mx - 14*lnx - 24*Li2omx + 2*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-5)*
     & R6 * (  - 1040.D0/9.D0 + 2092.D0/3.D0*Li2mx - 3869.D0/9.D0*lnx
     &     + 376.D0/3.D0*lnopxon2 + 800.D0/3.D0*ln2 - 2092.D0/3.D0*
     &    Li2omx + 523.D0/9.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-5)*
     & R8 * (  - 76805.D0/36.D0 + 1205999.D0/180.D0*Li2mx - 9168679.D0/
     &    2160.D0*lnx + 11239.D0/6.D0*lnopxon2 + 4375*ln2 - 1205999.D0/
     &    180.D0*Li2omx + 1205999.D0/2160.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-4)*
     & R4 * ( 10 - 96*Li2mx + 56*lnx - 12*lnopxon2 - 24*ln2 + 96*Li2omx
     &     - 8*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-4)*
     & R6 * ( 7469.D0/27.D0 - 936*Li2mx + 1730.D0/3.D0*lnx - 2234.D0/9.D
     &    0*lnopxon2 - 5116.D0/9.D0*ln2 + 936*Li2omx - 78*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-4)*
     & R8 * ( 10806163.D0/4320.D0 - 326179.D0/60.D0*Li2mx + 7433033.D0/
     &    2160.D0*lnx - 675899.D0/360.D0*lnopxon2 - 845339.D0/180.D0*
     &    ln2 + 326179.D0/60.D0*Li2omx - 326179.D0/720.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-3)*
     & R2 * ( 16*Li2mx - 8*lnx - 16*Li2omx + 4.D0/3.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-3)*
     & R4 * (  - 41 + 470.D0/3.D0*Li2mx - 274.D0/3.D0*lnx + 38*lnopxon2
     &     + 84*ln2 - 470.D0/3.D0*Li2omx + 235.D0/18.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-3)*
     & R6 * (  - 19687.D0/54.D0 + 2326.D0/3.D0*Li2mx - 8597.D0/18.D0*
     &    lnx + 2417.D0/9.D0*lnopxon2 + 6046.D0/9.D0*ln2 - 2326.D0/3.D0
     &    *Li2omx + 1163.D0/18.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-3)*
     & R8 * (  - 16843579.D0/8640.D0 + 1638799.D0/540.D0*Li2mx - 
     &    8295319.D0/4320.D0*lnx + 178631.D0/144.D0*lnopxon2 + 242287.D0
     &    /72.D0*ln2 - 1638799.D0/540.D0*Li2omx + 1638799.D0/6480.D0*
     &    pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-2)*
     & R2 * ( 8 - 40*Li2mx + 24*lnx - 8*lnopxon2 - 16*ln2 + 40*Li2omx
     &     - 10.D0/3.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-2)*
     & R4 * ( 395.D0/6.D0 - 400.D0/3.D0*Li2mx + 78*lnx - 139.D0/3.D0*
     &    lnopxon2 - 350.D0/3.D0*ln2 + 400.D0/3.D0*Li2omx - 100.D0/9.D0
     &    *pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-2)*
     & R6 * ( 7733.D0/27.D0 - 1214.D0/3.D0*Li2mx + 4501.D0/18.D0*lnx - 
     &    1541.D0/9.D0*lnopxon2 - 4294.D0/9.D0*ln2 + 1214.D0/3.D0*
     &    Li2omx - 607.D0/18.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-2)*
     & R8 * ( 26257819.D0/25920.D0 - 620363.D0/540.D0*Li2mx + 262039.D0/
     &    360.D0*lnx - 77735.D0/144.D0*lnopxon2 - 116303.D0/72.D0*ln2
     &     + 620363.D0/540.D0*Li2omx - 620363.D0/6480.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-1)*
     & R2 * (  - 16 + 32*Li2mx - 28*lnx + 16*lnopxon2 + 32*ln2 - 32*
     &    Li2omx + 8.D0/3.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-1)*
     & R4 * (  - 619.D0/12.D0 + 61*Li2mx - 229.D0/6.D0*lnx + 83.D0/3.D0
     &    *lnopxon2 + 81*ln2 - 61*Li2omx + 61.D0/12.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-1)*
     & R6 * (  - 537.D0/4.D0 + 386.D0/3.D0*Li2mx - 1457.D0/18.D0*lnx + 
     &    1157.D0/18.D0*lnopxon2 + 1840.D0/9.D0*ln2 - 386.D0/3.D0*
     &    Li2omx + 193.D0/18.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*omx**(-1)*
     & R8 * (  - 17621993.D0/51840.D0 + 67399.D0/240.D0*Li2mx - 258071.D
     &    0/1440.D0*lnx + 42911.D0/288.D0*lnopxon2 + 8059.D0/16.D0*ln2
     &     - 67399.D0/240.D0*Li2omx + 67399.D0/2880.D0*pisq )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF * (  - 26.D0
     &    /9.D0 + 286.D0/9.D0*ln2 - 88.D0/3.D0*ln2**2 - 8.D0/3.D0*pisq*
     &    ln2 + 262.D0/9.D0*x*ln2 - 88.D0/3.D0*x*ln2**2 - 8.D0/3.D0*x*
     &    pisq*ln2 )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*R2 * ( 3913.
     &    D0/675.D0 - 8*Li2mx + 12*lnx - 8*lnopxon2 - 671.D0/45.D0*ln2
     &     + 8*Li2omx - pisq - 5011.D0/1800.D0*x + 4*x*Li2mx + 4*x*lnx
     &     - 4*x*lnopxon2 - 41.D0/15.D0*x*ln2 - 4*x*Li2omx - 847.D0/300.
     &    D0*x**2 + 16.D0/5.D0*x**2*ln2 )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*R4 * ( 
     &    453509339.D0/25401600.D0 - 1.D0/8.D0*opx**(-2) - 3.D0/8.D0*
     &    opx**(-1) - 37.D0/3.D0*Li2mx + 19.D0/2.D0*lnx - 22.D0/3.D0*
     &    lnopxon2 - 185807.D0/7560.D0*ln2 + 37.D0/3.D0*Li2omx - 221.D0/
     &    216.D0*pisq + 57755723.D0/25401600.D0*x + 13.D0/6.D0*x*lnx - 
     &    13.D0/6.D0*x*lnopxon2 - 33419.D0/7560.D0*x*ln2 + 1.D0/216.D0*
     &    x*pisq - 42743.D0/50400.D0*x**2 + 14.D0/15.D0*x**2*ln2 + 3301.
     &    D0/11760.D0*x**3 - 2.D0/7.D0*x**3*ln2 )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*R6 * ( 
     &    31457362511.D0/1016064000.D0 - 1.D0/48.D0*opx**(-4) + 1.D0/
     &    144.D0*opx**(-3) - 11.D0/576.D0*opx**(-2) - 1.D0/18.D0*
     &    opx**(-1) - 62.D0/3.D0*Li2mx + 247.D0/18.D0*lnx - 217.D0/18.D0
     &    *lnopxon2 - 831293.D0/18900.D0*ln2 + 62.D0/3.D0*Li2omx - 31.D0
     &    /18.D0*pisq + 360275579.D0/101606400.D0*x - 4.D0/3.D0*x*Li2mx
     &     + 4.D0/3.D0*x*lnx - 4.D0/3.D0*x*lnopxon2 - 3229.D0/630.D0*x*
     &    ln2 + 4.D0/3.D0*x*Li2omx - 1.D0/9.D0*x*pisq - 51752123.D0/
     &    152409600.D0*x**2 + 323.D0/945.D0*x**2*ln2 + 1245577.D0/
     &    7620480.D0*x**3 - 31.D0/189.D0*x**3*ln2 - 168493.D0/3810240.D0
     &    *x**4 + 8.D0/189.D0*x**4*ln2 )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*R8 * ( 
     &    8674382583266579.D0/141631193088000.D0 - 5.D0/768.D0*
     &    opx**(-6) + 3.D0/256.D0*opx**(-5) - 13.D0/1536.D0*opx**(-4)
     &     + 11.D0/3072.D0*opx**(-3) - 1117.D0/276480.D0*opx**(-2) - 
     &    701.D0/138240.D0*opx**(-1) - 80279.D0/2160.D0*Li2mx + 26129.D0
     &    /1080.D0*lnx - 31679.D0/1440.D0*lnopxon2 - 10235195009.D0/
     &    119750400.D0*ln2 + 80279.D0/2160.D0*Li2omx - 4816739.D0/
     &    1555200.D0*pisq + 857223100927811.D0/141631193088000.D0*x - 7.
     &    D0/3.D0*x*Li2mx + 5.D0/3.D0*x*lnx - 5.D0/3.D0*x*lnopxon2 - 
     &    961756217.D0/119750400.D0*x*ln2 + 7.D0/3.D0*x*Li2omx - 302399.
     &    D0/1555200.D0*x*pisq + 3219497383.D0/98354995200.D0*x**2 - 
     &    953.D0/20790.D0*x**2*ln2 + 2684880109.D0/36883123200.D0*x**3
     &     - 1469.D0/20790.D0*x**3*ln2 - 891289519.D0/24588748800.D0*
     &    x**4 + 479.D0/13860.D0*x**4*ln2 + 1929101.D0/234178560.D0*
     &    x**5 - 1.D0/132.D0*x**5*ln2 )
      DelI2bit(QQCACF,rglr) = DelI2bit(QQCACF,rglr) + CA*CF*lnRonpi
     &  * (  - 262.D0/9.D0 + 88.D0/3.D0*ln2 + 8.D0/3.D0*pisq - 262.D0/9.
     &    D0*x + 88.D0/3.D0*x*ln2 + 8.D0/3.D0*x*pisq )

      DelI2bit(QQCFCF,rglr) =  + CF**2*omx**(-9)*R8 * (  - 280*Li2mx + 
     &    533.D0/3.D0*lnx + 280*Li2omx - 70.D0/3.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-8)*
     & R8 * (  - 307.D0/3.D0 + 1960*Li2mx - 3731.D0/3.D0*lnx + 140*
     &    lnopxon2 + 280*ln2 - 1960*Li2omx + 490.D0/3.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-7)*
     & R6 * (  - 320.D0/3.D0*Li2mx + 592.D0/9.D0*lnx + 320.D0/3.D0*
     &    Li2omx - 80.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-7)*
     & R8 * ( 4411.D0/6.D0 - 6120*Li2mx + 11647.D0/3.D0*lnx - 870*
     &    lnopxon2 - 1820*ln2 + 6120*Li2omx - 510*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-6)*
     & R6 * (  - 368.D0/9.D0 + 1760.D0/3.D0*Li2mx - 3256.D0/9.D0*lnx + 
     &    160.D0/3.D0*lnopxon2 + 320.D0/3.D0*ln2 - 1760.D0/3.D0*Li2omx
     &     + 440.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-6)*
     & R8 * (  - 20945.D0/9.D0 + 33680.D0/3.D0*Li2mx - 64063.D0/9.D0*
     &    lnx + 7130.D0/3.D0*lnopxon2 + 15700.D0/3.D0*ln2 - 33680.D0/3.D
     &    0*Li2omx + 8420.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-5)*
     & R4 * (  - 48*Li2mx + 28*lnx + 48*Li2omx - 4*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-5)*
     & R6 * ( 2080.D0/9.D0 - 4184.D0/3.D0*Li2mx + 7738.D0/9.D0*lnx - 
     &    752.D0/3.D0*lnopxon2 - 1600.D0/3.D0*ln2 + 4184.D0/3.D0*Li2omx
     &     - 1046.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-5)*
     & R8 * ( 76805.D0/18.D0 - 1205999.D0/90.D0*Li2mx + 9168679.D0/1080.
     &    D0*lnx - 11239.D0/3.D0*lnopxon2 - 8750*ln2 + 1205999.D0/90.D0
     &    *Li2omx - 1205999.D0/1080.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-4)*
     & R4 * (  - 20 + 192*Li2mx - 112*lnx + 24*lnopxon2 + 48*ln2 - 192*
     &    Li2omx + 16*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-4)*
     & R6 * (  - 14938.D0/27.D0 + 1872*Li2mx - 3460.D0/3.D0*lnx + 4468.D
     &    0/9.D0*lnopxon2 + 10232.D0/9.D0*ln2 - 1872*Li2omx + 156*pisq
     &     )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-4)*
     & R8 * (  - 10806163.D0/2160.D0 + 326179.D0/30.D0*Li2mx - 7433033.D
     &    0/1080.D0*lnx + 675899.D0/180.D0*lnopxon2 + 845339.D0/90.D0*
     &    ln2 - 326179.D0/30.D0*Li2omx + 326179.D0/360.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-3)*
     & R2 * (  - 32*Li2mx + 16*lnx + 32*Li2omx - 8.D0/3.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-3)*
     & R4 * ( 82 - 940.D0/3.D0*Li2mx + 548.D0/3.D0*lnx - 76*lnopxon2 - 
     &    168*ln2 + 940.D0/3.D0*Li2omx - 235.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-3)*
     & R6 * ( 19687.D0/27.D0 - 4652.D0/3.D0*Li2mx + 8597.D0/9.D0*lnx - 
     &    4834.D0/9.D0*lnopxon2 - 12092.D0/9.D0*ln2 + 4652.D0/3.D0*
     &    Li2omx - 1163.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-3)*
     & R8 * ( 16843579.D0/4320.D0 - 1638799.D0/270.D0*Li2mx + 8295319.D0
     &    /2160.D0*lnx - 178631.D0/72.D0*lnopxon2 - 242287.D0/36.D0*ln2
     &     + 1638799.D0/270.D0*Li2omx - 1638799.D0/3240.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-2)*
     & R2 * (  - 16 + 80*Li2mx - 48*lnx + 16*lnopxon2 + 32*ln2 - 80*
     &    Li2omx + 20.D0/3.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-2)*
     & R4 * (  - 395.D0/3.D0 + 800.D0/3.D0*Li2mx - 156*lnx + 278.D0/3.D0
     &    *lnopxon2 + 700.D0/3.D0*ln2 - 800.D0/3.D0*Li2omx + 200.D0/9.D0
     &    *pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-2)*
     & R6 * (  - 15466.D0/27.D0 + 2428.D0/3.D0*Li2mx - 4501.D0/9.D0*lnx
     &     + 3082.D0/9.D0*lnopxon2 + 8588.D0/9.D0*ln2 - 2428.D0/3.D0*
     &    Li2omx + 607.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-2)*
     & R8 * (  - 26257819.D0/12960.D0 + 620363.D0/270.D0*Li2mx - 262039.
     &    D0/180.D0*lnx + 77735.D0/72.D0*lnopxon2 + 116303.D0/36.D0*ln2
     &     - 620363.D0/270.D0*Li2omx + 620363.D0/3240.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-1)*
     & R2 * ( 32 - 64*Li2mx + 56*lnx - 32*lnopxon2 - 64*ln2 + 64*Li2omx
     &     - 16.D0/3.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-1)*
     & R4 * ( 619.D0/6.D0 - 122*Li2mx + 229.D0/3.D0*lnx - 166.D0/3.D0*
     &    lnopxon2 - 162*ln2 + 122*Li2omx - 61.D0/6.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-1)*
     & R6 * ( 537.D0/2.D0 - 772.D0/3.D0*Li2mx + 1457.D0/9.D0*lnx - 1157.
     &    D0/9.D0*lnopxon2 - 3680.D0/9.D0*ln2 + 772.D0/3.D0*Li2omx - 
     &    193.D0/9.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*omx**(-1)*
     & R8 * ( 17621993.D0/25920.D0 - 67399.D0/120.D0*Li2mx + 258071.D0/
     &    720.D0*lnx - 42911.D0/144.D0*lnopxon2 - 8059.D0/8.D0*ln2 + 
     &    67399.D0/120.D0*Li2omx - 67399.D0/1440.D0*pisq )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*R2 * (  - 
     &    26.D0/9.D0 + 16*Li2mx - 24*lnx + 16*lnopxon2 + 56.D0/3.D0*ln2
     &     - 16*Li2omx + 8.D0/3.D0*pisq + 52.D0/9.D0*x - 8*x*Li2mx - 8*
     &    x*lnx + 8*x*lnopxon2 + 8.D0/3.D0*x*ln2 + 8*x*Li2omx + 2.D0/3.D
     &    0*x*pisq + 11.D0/3.D0*x**2 - 4*x**2*ln2 )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*R4 * (  - 
     &    258629.D0/7200.D0 + 1.D0/4.D0*opx**(-2) + 3.D0/4.D0*opx**(-1)
     &     + 74.D0/3.D0*Li2mx - 19*lnx + 44.D0/3.D0*lnopxon2 + 1463.D0/
     &    30.D0*ln2 - 74.D0/3.D0*Li2omx + 37.D0/18.D0*pisq - 55901.D0/
     &    10800.D0*x - 13.D0/3.D0*x*lnx + 13.D0/3.D0*x*lnopxon2 + 809.D0
     &    /90.D0*x*ln2 + 2173.D0/3600.D0*x**2 - 11.D0/15.D0*x**2*ln2 - 
     &    239.D0/720.D0*x**3 + 1.D0/3.D0*x**3*ln2 )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*R6 * (  - 
     &    8198009.D0/132300.D0 + 1.D0/24.D0*opx**(-4) - 1.D0/72.D0*
     &    opx**(-3) + 11.D0/288.D0*opx**(-2) + 1.D0/9.D0*opx**(-1) + 
     &    124.D0/3.D0*Li2mx - 247.D0/9.D0*lnx + 217.D0/9.D0*lnopxon2 + 
     &    18479.D0/210.D0*ln2 - 124.D0/3.D0*Li2omx + 31.D0/9.D0*pisq - 
     &    3774257.D0/529200.D0*x + 8.D0/3.D0*x*Li2mx - 8.D0/3.D0*x*lnx
     &     + 8.D0/3.D0*x*lnopxon2 + 6479.D0/630.D0*x*ln2 - 8.D0/3.D0*x*
     &    Li2omx + 2.D0/9.D0*x*pisq + 2628887.D0/4233600.D0*x**2 - 373.D
     &    0/630.D0*x**2*ln2 - 32917.D0/264600.D0*x**3 + 83.D0/630.D0*
     &    x**3*ln2 + 1751.D0/37800.D0*x**4 - 2.D0/45.D0*x**4*ln2 )
      DelI2bit(QQCFCF,rglr) = DelI2bit(QQCFCF,rglr) + CF**2*R8 * (  - 
     &    995669660791.D0/8128512000.D0 + 5.D0/384.D0*opx**(-6) - 3.D0/
     &    128.D0*opx**(-5) + 13.D0/768.D0*opx**(-4) - 11.D0/1536.D0*
     &    opx**(-3) + 1117.D0/138240.D0*opx**(-2) + 701.D0/69120.D0*
     &    opx**(-1) + 80279.D0/1080.D0*Li2mx - 26129.D0/540.D0*lnx + 
     &    31679.D0/720.D0*lnopxon2 + 2871817.D0/16800.D0*ln2 - 80279.D0/
     &    1080.D0*Li2omx + 80279.D0/12960.D0*pisq - 73772512021.D0/
     &    6096384000.D0*x + 14.D0/3.D0*x*Li2mx - 10.D0/3.D0*x*lnx + 10.D
     &    0/3.D0*x*lnopxon2 + 7284787.D0/453600.D0*x*ln2 - 14.D0/3.D0*x
     &    *Li2omx + 7.D0/18.D0*x*pisq - 3032819.D0/45158400.D0*x**2 + 
     &    733.D0/7560.D0*x**2*ln2 - 44843851.D0/406425600.D0*x**3 + 773.
     &    D0/7560.D0*x**3*ln2 + 8171.D0/290304.D0*x**4 - 1.D0/36.D0*
     &    x**4*ln2 - 86623.D0/11289600.D0*x**5 + 1.D0/140.D0*x**5*ln2 )

      DelI2bit(QQCFTF,rglr) =  + CF*TF * ( 52.D0/9.D0 - 140.D0/9.D0*ln2
     &     + 32.D0/3.D0*ln2**2 - 92.D0/9.D0*x*ln2 + 32.D0/3.D0*x*ln2**2
     &     )
      DelI2bit(QQCFTF,rglr) = DelI2bit(QQCFTF,rglr) + CF*TF*R2 * (  - 
     &    2077.D0/1350.D0 + 52.D0/45.D0*ln2 + 361.D0/900.D0*x - 8.D0/15.
     &    D0*x*ln2 - 253.D0/150.D0*x**2 + 8.D0/5.D0*x**2*ln2 )
      DelI2bit(QQCFTF,rglr) = DelI2bit(QQCFTF,rglr) + CF*TF*R4 * ( 
     &    294401.D0/6350400.D0 - 53.D0/1890.D0*ln2 + 1253831.D0/6350400.
     &    D0*x - 353.D0/1890.D0*x*ln2 - 3307.D0/12600.D0*x**2 + 4.D0/15.
     &    D0*x**2*ln2 + 14863.D0/58800.D0*x**3 - 8.D0/35.D0*x**3*ln2 )
      DelI2bit(QQCFTF,rglr) = DelI2bit(QQCFTF,rglr) + CF*TF*R6 * (  - 
     &    2669531.D0/508032000.D0 + 53.D0/9450.D0*ln2 + 430939.D0/
     &    50803200.D0*x - 4.D0/315.D0*x*ln2 - 90457.D0/1555200.D0*x**2
     &     + 8.D0/135.D0*x**2*ln2 + 1565863.D0/19051200.D0*x**3 - 74.D0/
     &    945.D0*x**3*ln2 - 91859.D0/1905120.D0*x**4 + 8.D0/189.D0*x**4
     &    *ln2 )
      DelI2bit(QQCFTF,rglr) = DelI2bit(QQCFTF,rglr) + CF*TF*R8 * ( 
     &    12601498991.D0/35407798272000.D0 - 4001.D0/29937600.D0*ln2 + 
     &    55378445891.D0/35407798272000.D0*x - 67241.D0/29937600.D0*x*
     &    ln2 - 22092977.D0/3073593600.D0*x**2 + 179.D0/20790.D0*x**2*
     &    ln2 + 3036176107.D0/147532492800.D0*x**3 - 214.D0/10395.D0*
     &    x**3*ln2 - 35752141.D0/1536796800.D0*x**4 + 74.D0/3465.D0*
     &    x**4*ln2 + 2062751.D0/204906240.D0*x**5 - 2.D0/231.D0*x**5*
     &    ln2 )
      DelI2bit(QQCFTF,rglr) = DelI2bit(QQCFTF,rglr) + CF*TF*lnRonpi
     &  * ( 92.D0/9.D0 - 32.D0/3.D0*ln2 + 92.D0/9.D0*x - 32.D0/3.D0*x*
     &    ln2 )

      DelI2bit(GQCFCA,rglr) =  + CA*CF*omx**(-8)*R8 * (  - 70*Li2mx + 
     &    533.D0/12.D0*lnx + 70*Li2omx - 35.D0/6.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-7)*
     & R8 * (  - 307.D0/12.D0 + 470*Li2mx - 3553.D0/12.D0*lnx + 35*
     &    lnopxon2 + 70*ln2 - 470*Li2omx + 235.D0/6.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-6)*
     & R6 * (  - 80.D0/3.D0*Li2mx + 148.D0/9.D0*lnx + 80.D0/3.D0*Li2omx
     &     - 20.D0/9.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-6)*
     & R8 * ( 1429.D0/8.D0 - 1435*Li2mx + 7191.D0/8.D0*lnx - 415.D0/2.D0
     &    *lnopxon2 - 435*ln2 + 1435*Li2omx - 1435.D0/12.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-5)*
     & R6 * (  - 92.D0/9.D0 + 424.D0/3.D0*Li2mx - 770.D0/9.D0*lnx + 40.D
     &    0/3.D0*lnopxon2 + 80.D0/3.D0*ln2 - 424.D0/3.D0*Li2omx + 106.D0
     &    /9.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-5)*
     & R8 * (  - 40039.D0/72.D0 + 7916.D0/3.D0*Li2mx - 59197.D0/36.D0*
     &    lnx + 1661.D0/3.D0*lnopxon2 + 3670.D0/3.D0*ln2 - 7916.D0/3.D0
     &    *Li2omx + 1979.D0/9.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-4)*
     & R4 * (  - 12*Li2mx + 7*lnx + 12*Li2omx - pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-4)*
     & R6 * ( 172.D0/3.D0 - 998.D0/3.D0*Li2mx + 3581.D0/18.D0*lnx - 60*
     &    lnopxon2 - 128*ln2 + 998.D0/3.D0*Li2omx - 499.D0/18.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-4)*
     & R8 * ( 49205.D0/48.D0 - 1169219.D0/360.D0*Li2mx + 8701189.D0/
     &    4320.D0*lnx - 880*lnopxon2 - 6181.D0/3.D0*ln2 + 1169219.D0/
     &    360.D0*Li2omx - 1169219.D0/4320.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-3)*
     & R4 * (  - 5 + 48*Li2mx - 26*lnx + 6*lnopxon2 + 12*ln2 - 48*
     &    Li2omx + 4*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-3)*
     & R6 * (  - 7409.D0/54.D0 + 1384.D0/3.D0*Li2mx - 817.D0/3.D0*lnx
     &     + 1057.D0/9.D0*lnopxon2 + 2438.D0/9.D0*ln2 - 1384.D0/3.D0*
     &    Li2omx + 346.D0/9.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-3)*
     & R8 * (  - 10747843.D0/8640.D0 + 1007519.D0/360.D0*Li2mx - 498319.
     &    D0/288.D0*lnx + 663599.D0/720.D0*lnopxon2 + 829199.D0/360.D0*
     &    ln2 - 1007519.D0/360.D0*Li2omx + 1007519.D0/4320.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-2)*
     & R2 * (  - 8*Li2mx + 4*lnx + 8*Li2omx - 2.D0/3.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-2)*
     & R4 * ( 45.D0/2.D0 - 241.D0/3.D0*Li2mx + 253.D0/6.D0*lnx - 19*
     &    lnopxon2 - 42*ln2 + 241.D0/3.D0*Li2omx - 241.D0/36.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-2)*
     & R6 * ( 20363.D0/108.D0 - 1240.D0/3.D0*Li2mx + 2123.D0/9.D0*lnx
     &     - 2401.D0/18.D0*lnopxon2 - 3019.D0/9.D0*ln2 + 1240.D0/3.D0*
     &    Li2omx - 310.D0/9.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-2)*
     & R8 * ( 1985723.D0/1920.D0 - 3685193.D0/2160.D0*Li2mx + 304043.D0/
     &    288.D0*lnx - 949199.D0/1440.D0*lnopxon2 - 142231.D0/80.D0*ln2
     &     + 3685193.D0/2160.D0*Li2omx - 3685193.D0/25920.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-1)*
     & R2 * (  - 4 + 28*Li2mx - 4*lnx + 4*lnopxon2 + 8*ln2 - 28*Li2omx
     &     + 7.D0/3.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-1)*
     & R4 * (  - 449.D0/12.D0 + 80*Li2mx - 349.D0/6.D0*lnx + 133.D0/6.D0
     &    *lnopxon2 + 181.D0/3.D0*ln2 - 80*Li2omx + 20.D0/3.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-1)*
     & R6 * (  - 18103.D0/108.D0 + 240*Li2mx - 1297.D0/9.D0*lnx + 841.D0
     &    /9.D0*lnopxon2 + 2372.D0/9.D0*ln2 - 240*Li2omx + 20*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*omx**(-1)*
     & R8 * (  - 7647787.D0/12960.D0 + 86533.D0/120.D0*Li2mx - 43573.D0/
     &    96.D0*lnx + 231599.D0/720.D0*lnopxon2 + 341419.D0/360.D0*ln2
     &     - 86533.D0/120.D0*Li2omx + 86533.D0/1440.D0*pisq )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*R2 * ( 47.D0
     &    /18.D0 - 11.D0/3.D0*x**(-1) + 4*x**(-1)*ln2 - 4.D0/3.D0*
     &    x**(-1)*pisq - 44*Li2mx + 4*lnopxon2 + 2.D0/3.D0*ln2 + 44*
     &    Li2omx - 7.D0/3.D0*pisq + 451.D0/18.D0*x + 20*x*Li2mx + 4*x*
     &    lnx - 20*x*lnopxon2 - 134.D0/3.D0*x*ln2 - 20*x*Li2omx + x*
     &    pisq - 8*x**2 - 16*x**2*lnx + 16*x**2*lnopxon2 + 24*x**2*ln2
     &     )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*R4 * ( 
     &    1976299.D0/43200.D0 + 55.D0/288.D0*x**(-1) + 2.D0/3.D0*
     &    x**(-1)*ln2 - 37.D0/8.D0*opx**(-2) - 197.D0/16.D0*opx**(-1)
     &     - 109.D0/3.D0*Li2mx + 35*lnx - 17.D0/6.D0*lnopxon2 - 3109.D0/
     &    90.D0*ln2 + 109.D0/3.D0*Li2omx - 109.D0/36.D0*pisq + 264713.D0
     &    /43200.D0*x + 8*x*Li2mx + 193.D0/6.D0*x*lnx - 145.D0/6.D0*x*
     &    lnopxon2 - 4103.D0/90.D0*x*ln2 - 8*x*Li2omx + 2.D0/3.D0*x*
     &    pisq - 37379.D0/1440.D0*x**2 + 16*x**2*Li2mx + 77.D0/3.D0*
     &    x**2*ln2 - 16*x**2*Li2omx + 4.D0/3.D0*x**2*pisq - 4541.D0/
     &    3600.D0*x**3 + 22.D0/15.D0*x**3*ln2 )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*R6 * ( 
     &    40830623.D0/564480.D0 + 7.D0/384.D0*x**(-1) - 53.D0/48.D0*
     &    opx**(-4) + 133.D0/288.D0*opx**(-3) - 419.D0/576.D0*opx**(-2)
     &     - 1667.D0/1152.D0*opx**(-1) - 70*Li2mx + 913.D0/18.D0*lnx - 
     &    553.D0/18.D0*lnopxon2 - 6173.D0/63.D0*ln2 + 70*Li2omx - 35.D0/
     &    6.D0*pisq + 37979987.D0/1693440.D0*x - 28.D0/3.D0*x*Li2mx + 
     &    24*x*lnx - 44.D0/3.D0*x*lnopxon2 - 4499.D0/126.D0*x*ln2 + 28.D
     &    0/3.D0*x*Li2omx - 7.D0/9.D0*x*pisq - 1513369.D0/120960.D0*
     &    x**2 + 8*x**2*Li2mx + 28.D0/3.D0*x**2*lnx - 28.D0/3.D0*x**2*
     &    lnopxon2 + 41.D0/18.D0*x**2*ln2 - 8*x**2*Li2omx + 2.D0/3.D0*
     &    x**2*pisq - 1863271.D0/529200.D0*x**3 + 1241.D0/315.D0*x**3*
     &    ln2 + 8987.D0/29400.D0*x**4 - 32.D0/105.D0*x**4*ln2 )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*R8 * ( 
     &    1236113042791.D0/6967296000.D0 - 11441.D0/33177600.D0*x**(-1)
     &     - 115.D0/256.D0*opx**(-6) + 1271.D0/1536.D0*opx**(-5) - 875.D
     &    0/1536.D0*opx**(-4) + 1081.D0/9216.D0*opx**(-3) + 847.D0/6912.
     &    D0*opx**(-2) + 185911.D0/552960.D0*opx**(-1) - 366119.D0/2160.
     &    D0*Li2mx + 120359.D0/1080.D0*lnx - 122879.D0/1440.D0*lnopxon2
     &     - 17615669.D0/64800.D0*ln2 + 366119.D0/2160.D0*Li2omx - 
     &    366119.D0/25920.D0*pisq + 2022645417353.D0/48771072000.D0*x
     &     - 239.D0/9.D0*x*Li2mx + 203.D0/9.D0*x*lnx - 160.D0/9.D0*x*
     &    lnopxon2 - 26741467.D0/453600.D0*x*ln2 + 239.D0/9.D0*x*Li2omx
     &     - 239.D0/108.D0*x*pisq + 146255561.D0/162570240.D0*x**2 + 2.D
     &    0/9.D0*x**2*Li2mx + 47.D0/9.D0*x**2*lnx - 47.D0/9.D0*x**2*
     &    lnopxon2 - 10117.D0/1512.D0*x**2*ln2 - 2.D0/9.D0*x**2*Li2omx
     &     + 1.D0/54.D0*x**2*pisq - 241772707.D0/135475200.D0*x**3 + 
     &    3877.D0/1890.D0*x**3*ln2 + 114859393.D0/203212800.D0*x**4 - 
     &    359.D0/630.D0*x**4*ln2 )
      DelI2bit(GQCFCA,rglr) = DelI2bit(GQCFCA,rglr) + CA*CF*R8 * (  - 
     &    346247.D0/4838400.D0*x**5 + 1.D0/15.D0*x**5*ln2 )

      DelI2bit(GQCFCF,rglr) =  + CF**2 * (  - 16.D0/3.D0*x**(-1)*lnpi*
     &    ln512 - 48*x**(-1)*ln2 + 48*x**(-1)*ln2*lnpi + 48*x**(-1)*
     &    ln2**2 + 16.D0/3.D0*x**(-1)*pisq*ln2 + 16.D0/3.D0*lnpi*ln512
     &     + 48*ln2 - 48*ln2*lnpi - 48*ln2**2 - 16.D0/3.D0*pisq*ln2 - 8.
     &    D0/3.D0*x*lnpi*ln512 - 24*x*ln2 + 24*x*ln2*lnpi + 24*x*ln2**2
     &     + 8.D0/3.D0*x*pisq*ln2 )
      DelI2bit(GQCFCF,rglr) = DelI2bit(GQCFCF,rglr) + CF**2*R2 * (  - 
     &    305.D0/18.D0 + 145.D0/6.D0*x**(-1) - 22*x**(-1)*ln2 - 2.D0/3.D
     &    0*x**(-1)*pisq + 46.D0/3.D0*ln2 + 2.D0/3.D0*pisq - 245.D0/36.D
     &    0*x + 23.D0/3.D0*x*ln2 - 1.D0/3.D0*x*pisq + 8*x**2 - 8*x**2*
     &    ln2 )
      DelI2bit(GQCFCF,rglr) = DelI2bit(GQCFCF,rglr) + CF**2*R4 * (  - 
     &    18221.D0/86400.D0 + 815.D0/1152.D0*x**(-1) - 3.D0/4.D0*
     &    x**(-1)*ln2 - 1.D0/108.D0*x**(-1)*pisq + 11.D0/180.D0*ln2 + 1.
     &    D0/108.D0*pisq - 314483.D0/172800.D0*x + 773.D0/360.D0*x*ln2
     &     - 1.D0/216.D0*x*pisq + 5291.D0/1800.D0*x**2 - 44.D0/15.D0*
     &    x**2*ln2 - 8341.D0/7200.D0*x**3 + 16.D0/15.D0*x**3*ln2 )
      DelI2bit(GQCFCF,rglr) = DelI2bit(GQCFCF,rglr) + CF**2*R6 * ( 8371.
     &    D0/2822400.D0 + 1.D0/6912.D0*x**(-1) - 8.D0/315.D0*ln2 - 
     &    411161.D0/2419200.D0*x + 7.D0/30.D0*x*ln2 + 391243.D0/705600.D
     &    0*x**2 - 62.D0/105.D0*x**2*ln2 - 653143.D0/1058400.D0*x**3 + 
     &    184.D0/315.D0*x**3*ln2 + 58841.D0/264600.D0*x**4 - 62.D0/315.D
     &    0*x**4*ln2 )
      DelI2bit(GQCFCF,rglr) = DelI2bit(GQCFCF,rglr) + CF**2*R8 * (  - 
     &    13824469.D0/195084288000.D0 + 4987.D0/265420800.D0*x**(-1) - 
     &    1.D0/86400.D0*x**(-1)*ln2 - 1.D0/777600.D0*x**(-1)*pisq - 
     &    2803.D0/1814400.D0*ln2 + 1.D0/777600.D0*pisq - 4676939899.D0/
     &    390168576000.D0*x + 3131.D0/145152.D0*x*ln2 - 1.D0/1555200.D0
     &    *x*pisq + 15586147.D0/203212800.D0*x**2 - 173.D0/1890.D0*x**2
     &    *ln2 - 133793623.D0/812851200.D0*x**3 + 313.D0/1890.D0*x**3*
     &    ln2 + 598765.D0/4064256.D0*x**4 - 17.D0/126.D0*x**4*ln2 - 
     &    3174257.D0/67737600.D0*x**5 + 17.D0/420.D0*x**5*ln2 )
      DelI2bit(GQCFCF,rglr) = DelI2bit(GQCFCF,rglr) + CF**2*lnRonpi
     &  * (  - 48 + 48*x**(-1) - 16.D0/3.D0*x**(-1)*ln512 - 16.D0/3.D0*
     &    x**(-1)*pisq + 16.D0/3.D0*ln512 + 16.D0/3.D0*pisq + 24*x - 8.D
     &    0/3.D0*x*ln512 - 8.D0/3.D0*x*pisq )

      DelI2bit(QGCATF,rglr) =  + CA*TF*omx**(-8)*R8 * (  - 70*Li2mx + 
     &    533.D0/12.D0*lnx + 70*Li2omx - 35.D0/6.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-7)*
     & R8 * (  - 307.D0/12.D0 + 440*Li2mx - 844.D0/3.D0*lnx + 35*
     &    lnopxon2 + 70*ln2 - 440*Li2omx + 110.D0/3.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-6)*
     & R6 * (  - 80.D0/3.D0*Li2mx + 148.D0/9.D0*lnx + 80.D0/3.D0*Li2omx
     &     - 20.D0/9.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-6)*
     & R8 * ( 1307.D0/8.D0 - 1255*Li2mx + 6483.D0/8.D0*lnx - 385.D0/2.D0
     &    *lnopxon2 - 405*ln2 + 1255*Li2omx - 1255.D0/12.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-5)*
     & R6 * (  - 92.D0/9.D0 + 376.D0/3.D0*Li2mx - 710.D0/9.D0*lnx + 40.D
     &    0/3.D0*lnopxon2 + 80.D0/3.D0*ln2 - 376.D0/3.D0*Li2omx + 94.D0/
     &    9.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-5)*
     & R8 * (  - 8365.D0/18.D0 + 6419.D0/3.D0*Li2mx - 100859.D0/72.D0*
     &    lnx + 2863.D0/6.D0*lnopxon2 + 3175.D0/3.D0*ln2 - 6419.D0/3.D0
     &    *Li2omx + 6419.D0/36.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-4)*
     & R4 * (  - 12*Li2mx + 7*lnx + 12*Li2omx - pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-4)*
     & R6 * ( 48 - 782.D0/3.D0*Li2mx + 3041.D0/18.D0*lnx - 52*lnopxon2
     &     - 112*ln2 + 782.D0/3.D0*Li2omx - 391.D0/18.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-4)*
     & R8 * ( 37169.D0/48.D0 - 865739.D0/360.D0*Li2mx + 6948889.D0/4320.
     &    D0*lnx - 1401.D0/2.D0*lnopxon2 - 4924.D0/3.D0*ln2 + 865739.D0/
     &    360.D0*Li2omx - 865739.D0/4320.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-3)*
     & R4 * (  - 5 + 36*Li2mx - 23*lnx + 6*lnopxon2 + 12*ln2 - 36*
     &    Li2omx + 3*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-3)*
     & R6 * (  - 5177.D0/54.D0 + 310*Li2mx - 3805.D0/18.D0*lnx + 805.D0/
     &    9.D0*lnopxon2 + 1862.D0/9.D0*ln2 - 310*Li2omx + 155.D0/6.D0*
     &    pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-3)*
     & R8 * (  - 7159903.D0/8640.D0 + 82882.D0/45.D0*Li2mx - 459431.D0/
     &    360.D0*lnx + 481079.D0/720.D0*lnopxon2 + 596639.D0/360.D0*ln2
     &     - 82882.D0/45.D0*Li2omx + 41441.D0/270.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-2)*
     & R2 * (  - 8*Li2mx + 4*lnx + 8*Li2omx - 2.D0/3.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-2)*
     & R4 * ( 27.D0/2.D0 - 133.D0/3.D0*Li2mx + 199.D0/6.D0*lnx - 13*
     &    lnopxon2 - 30*ln2 + 133.D0/3.D0*Li2omx - 133.D0/36.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-2)*
     & R6 * ( 11153.D0/108.D0 - 662.D0/3.D0*Li2mx + 347.D0/2.D0*lnx - 
     &    1567.D0/18.D0*lnopxon2 - 1933.D0/9.D0*ln2 + 662.D0/3.D0*
     &    Li2omx - 331.D0/18.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-2)*
     & R8 * ( 669827.D0/1152.D0 - 2055437.D0/2160.D0*Li2mx + 1013009.D0/
     &    1440.D0*lnx - 612473.D0/1440.D0*lnopxon2 - 89217.D0/80.D0*ln2
     &     + 2055437.D0/2160.D0*Li2omx - 2055437.D0/25920.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-1)*
     & R2 * (  - 4 + 4*Li2mx - 16*lnx + 4*lnopxon2 + 8*ln2 - 4*Li2omx
     &     + 1.D0/3.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-1)*
     & R4 * (  - 143.D0/12.D0 + 21*Li2mx - 23.D0/3.D0*lnx + 79.D0/6.D0*
     &    lnopxon2 + 91.D0/3.D0*ln2 - 21*Li2omx + 7.D0/4.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-1)*
     & R6 * (  - 2783.D0/54.D0 + 86*Li2mx - 1505.D0/18.D0*lnx + 899.D0/
     &    18.D0*lnopxon2 + 1139.D0/9.D0*ln2 - 86*Li2omx + 43.D0/6.D0*
     &    pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*omx**(-1)*
     & R8 * (  - 2641049.D0/10368.D0 + 224587.D0/720.D0*Li2mx - 562813.D
     &    0/2160.D0*lnx + 255829.D0/1440.D0*lnopxon2 + 349229.D0/720.D0
     &    *ln2 - 224587.D0/720.D0*Li2omx + 224587.D0/8640.D0*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*R2 * ( 59.D0
     &    /6.D0 + 4*Li2mx + 12*lnx - 4*lnopxon2 - 14*ln2 - 4*Li2omx - 1.
     &    D0/3.D0*pisq - 92.D0/9.D0*x - 16*x*Li2mx + 8*x*lnx + 16*x*
     &    lnopxon2 + 80.D0/3.D0*x*ln2 + 16*x*Li2omx - 65.D0/18.D0*x**2
     &     + 24*x**2*Li2mx + 24*x**2*lnx - 24*x**2*lnopxon2 - 62.D0/3.D0
     &    *x**2*ln2 - 24*x**2*Li2omx + 2.D0/3.D0*x**2*pisq )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*R4 * (  - 
     &    193829.D0/14400.D0 + 37.D0/8.D0*opx**(-2) + 167.D0/16.D0*
     &    opx**(-1) - 2.D0/3.D0*Li2mx - 19.D0/2.D0*lnx - 37.D0/6.D0*
     &    lnopxon2 - 51.D0/5.D0*ln2 + 2.D0/3.D0*Li2omx - 1.D0/18.D0*
     &    pisq + 22387.D0/1800.D0*x - 26.D0/3.D0*x*Li2mx - 53.D0/3.D0*x
     &    *lnx + 18*x*lnopxon2 + 109.D0/5.D0*x*ln2 + 26.D0/3.D0*x*
     &    Li2omx - 13.D0/18.D0*x*pisq + 59917.D0/4800.D0*x**2 - 8*x**2*
     &    Li2mx + 25.D0/3.D0*x**2*lnx - 25.D0/3.D0*x**2*lnopxon2 - 106.D
     &    0/5.D0*x**2*ln2 + 8*x**2*Li2omx - 2.D0/3.D0*x**2*pisq - 8659.D
     &    0/7200.D0*x**3 + 19.D0/15.D0*x**3*ln2 )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*R6 * ( 
     &    7012483.D0/1693440.D0 + 53.D0/48.D0*opx**(-4) - 199.D0/288.D0
     &    *opx**(-3) + 45.D0/64.D0*opx**(-2) + 1133.D0/1152.D0*
     &    opx**(-1) - 40.D0/3.D0*Li2mx + 15*lnx - 41.D0/3.D0*lnopxon2
     &     - 2090.D0/63.D0*ln2 + 40.D0/3.D0*Li2omx - 10.D0/9.D0*pisq - 
     &    3977959.D0/846720.D0*x - 4.D0/3.D0*x*Li2mx - 83.D0/9.D0*x*lnx
     &     + 35.D0/9.D0*x*lnopxon2 + 59.D0/14.D0*x*ln2 + 4.D0/3.D0*x*
     &    Li2omx - 1.D0/9.D0*x*pisq + 16365011.D0/1693440.D0*x**2 - 20.D
     &    0/3.D0*x**2*Li2mx - 20.D0/3.D0*x**2*lnx + 20.D0/3.D0*x**2*
     &    lnopxon2 - 178.D0/63.D0*x**2*ln2 + 20.D0/3.D0*x**2*Li2omx - 5.
     &    D0/9.D0*x**2*pisq + 10721689.D0/4233600.D0*x**3 - 1721.D0/630.
     &    D0*x**3*ln2 - 38833.D0/529200.D0*x**4 + 23.D0/315.D0*x**4*ln2
     &     )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*R8 * ( 
     &    25570422979.D0/464486400.D0 + 115.D0/256.D0*opx**(-6) - 1373.D
     &    0/1536.D0*opx**(-5) + 499.D0/768.D0*opx**(-4) - 1817.D0/9216.D
     &    0*opx**(-3) - 931.D0/27648.D0*opx**(-2) - 87673.D0/184320.D0*
     &    opx**(-1) - 56153.D0/1080.D0*Li2mx + 112541.D0/2160.D0*lnx - 
     &    3173.D0/80.D0*lnopxon2 - 117883.D0/1080.D0*ln2 + 56153.D0/
     &    1080.D0*Li2omx - 56153.D0/12960.D0*pisq + 196103581.D0/
     &    77414400.D0*x - 3719.D0/1080.D0*x*Li2mx + 2099.D0/540.D0*x*
     &    lnx - 1093.D0/240.D0*x*lnopxon2 - 243.D0/20.D0*x*ln2 + 3719.D0
     &    /1080.D0*x*Li2omx - 3719.D0/12960.D0*x*pisq + 6487237597.D0/
     &    3251404800.D0*x**2 - 13.D0/9.D0*x**2*Li2mx - 26.D0/9.D0*x**2*
     &    lnx + 26.D0/9.D0*x**2*lnopxon2 + 944.D0/945.D0*x**2*ln2 + 13.D
     &    0/9.D0*x**2*Li2omx - 13.D0/108.D0*x**2*pisq + 1362503729.D0/
     &    812851200.D0*x**3 - 13961.D0/7560.D0*x**3*ln2 - 92317843.D0/
     &    203212800.D0*x**4 + 563.D0/1260.D0*x**4*ln2 + 8029621.D0/
     &    203212800.D0*x**5 )
      DelI2bit(QGCATF,rglr) = DelI2bit(QGCATF,rglr) + CA*TF*R8 * (  - 
     &    23.D0/630.D0*x**5*ln2 )

      DelI2bit(QGCFTF,rglr) =  + CF*TF * (  - 8.D0/3.D0*lnpi*ln512 - 24
     &    *ln2 + 24*ln2*lnpi + 24*ln2**2 + 8.D0/3.D0*pisq*ln2 + 16.D0/3.
     &    D0*x*lnpi*ln512 + 48*x*ln2 - 48*x*ln2*lnpi - 48*x*ln2**2 - 16.
     &    D0/3.D0*x*pisq*ln2 - 16.D0/3.D0*x**2*lnpi*ln512 - 48*x**2*ln2
     &     + 48*x**2*ln2*lnpi + 48*x**2*ln2**2 + 16.D0/3.D0*x**2*pisq*
     &    ln2 )
      DelI2bit(QGCFTF,rglr) = DelI2bit(QGCFTF,rglr) + CF*TF*R2 * (  - 
     &    17.D0/12.D0 + 3*ln2 - 1.D0/3.D0*pisq - 41.D0/18.D0*x - 2.D0/3.
     &    D0*x*ln2 + 2.D0/3.D0*x*pisq - 35.D0/9.D0*x**2 + 5.D0/9.D0*
     &    x**2*ln4096 - 2.D0/3.D0*x**2*pisq )
      DelI2bit(QGCFTF,rglr) = DelI2bit(QGCFTF,rglr) + CF*TF*R4 * (  - 
     &    7073.D0/57600.D0 + 43.D0/120.D0*ln2 - 1.D0/216.D0*pisq + 107.D
     &    0/384.D0*x - 7.D0/12.D0*x*ln2 + 1.D0/108.D0*x*pisq - 1693.D0/
     &    3200.D0*x**2 + 19.D0/20.D0*x**2*ln2 - 1.D0/108.D0*x**2*pisq
     &     + 541.D0/3600.D0*x**3 - 2.D0/15.D0*x**3*ln2 )
      DelI2bit(QGCFTF,rglr) = DelI2bit(QGCFTF,rglr) + CF*TF*R6 * ( 5983.
     &    D0/1881600.D0 + 4.D0/315.D0*ln2 + 74717.D0/1693440.D0*x - 5.D0
     &    /63.D0*x*ln2 - 90833.D0/604800.D0*x**2 + 17.D0/90.D0*x**2*ln2
     &     + 232759.D0/1058400.D0*x**3 - 67.D0/315.D0*x**3*ln2 - 21989.D
     &    0/264600.D0*x**4 + 23.D0/315.D0*x**4*ln2 )
      DelI2bit(QGCFTF,rglr) = DelI2bit(QGCFTF,rglr) + CF*TF*R8 * ( 
     &    806377.D0/8670412800.D0 + 953.D0/1209600.D0*ln2 - 1.D0/
     &    1555200.D0*pisq + 84496229.D0/13005619200.D0*x - 6553.D0/
     &    604800.D0*x*ln2 + 1.D0/777600.D0*x*pisq - 2297374657.D0/
     &    65028096000.D0*x**2 + 26021.D0/604800.D0*x**2*ln2 - 1.D0/
     &    777600.D0*x**2*pisq + 165091.D0/2150400.D0*x**3 - 43.D0/540.D0
     &    *x**3*ln2 - 7458643.D0/101606400.D0*x**4 + 43.D0/630.D0*x**4*
     &    ln2 + 1217729.D0/50803200.D0*x**5 - 13.D0/630.D0*x**5*ln2 )
      DelI2bit(QGCFTF,rglr) = DelI2bit(QGCFTF,rglr) + CF*TF*lnRonpi
     &  * ( 24 - 8.D0/3.D0*ln512 - 8.D0/3.D0*pisq - 48*x + 16.D0/3.D0*x
     &    *ln512 + 16.D0/3.D0*x*pisq + 48*x**2 - 16.D0/3.D0*x**2*ln512
     &     - 16.D0/3.D0*x**2*pisq )

      DelI2bit(QQBARS,rglr) =  + CF*TF*omx**(-7)*R8 * (  - 80*Li2mx + 
     &    148.D0/3.D0*lnx + 80*Li2omx - 20.D0/3.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-6)*
     & R8 * (  - 92.D0/3.D0 + 480*Li2mx - 296*lnx + 40*lnopxon2 + 80*
     &    ln2 - 480*Li2omx + 40*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-5)*
     & R6 * (  - 32*Li2mx + 56.D0/3.D0*lnx + 32*Li2omx - 8.D0/3.D0*pisq
     &     )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-5)*
     & R8 * ( 566.D0/3.D0 - 3856.D0/3.D0*Li2mx + 7136.D0/9.D0*lnx - 208
     &    *lnopxon2 - 440*ln2 + 3856.D0/3.D0*Li2omx - 964.D0/9.D0*pisq
     &     )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-4)*
     & R6 * (  - 40.D0/3.D0 + 144*Li2mx - 84*lnx + 16*lnopxon2 + 32*ln2
     &     - 144*Li2omx + 12*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-4)*
     & R8 * (  - 4577.D0/9.D0 + 6068.D0/3.D0*Li2mx - 11239.D0/9.D0*lnx
     &     + 476*lnopxon2 + 1072*ln2 - 6068.D0/3.D0*Li2omx + 1517.D0/9.D
     &    0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-3)*
     & R4 * (  - 16*Li2mx + 8*lnx + 16*Li2omx - 4.D0/3.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-3)*
     & R6 * ( 184.D0/3.D0 - 280*Li2mx + 164*lnx - 176.D0/3.D0*lnopxon2
     &     - 128*ln2 + 280*Li2omx - 70.D0/3.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-3)*
     & R8 * ( 790 - 18590.D0/9.D0*Li2mx + 3838.D0/3.D0*lnx - 1888.D0/3.D
     &    0*lnopxon2 - 1520*ln2 + 18590.D0/9.D0*Li2omx - 9295.D0/54.D0*
     &    pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-2)*
     & R4 * (  - 8 + 48*Li2mx - 24*lnx + 8*lnopxon2 + 16*ln2 - 48*
     &    Li2omx + 4*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-2)*
     & R6 * (  - 352.D0/3.D0 + 920.D0/3.D0*Li2mx - 176*lnx + 272.D0/3.D0
     &    *lnopxon2 + 656.D0/3.D0*ln2 - 920.D0/3.D0*Li2omx + 230.D0/9.D0
     &    *pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-2)*
     & R8 * (  - 167839.D0/216.D0 + 4258.D0/3.D0*Li2mx - 7966.D0/9.D0*
     &    lnx + 1582.D0/3.D0*lnopxon2 + 4136.D0/3.D0*ln2 - 4258.D0/3.D0
     &    *Li2omx + 2129.D0/18.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-1)*
     & R2 * (  - 16*Li2mx + 16*Li2omx - 4.D0/3.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-1)*
     & R4 * ( 24 - 56*Li2mx + 56*lnx - 16*lnopxon2 - 40*ln2 + 56*Li2omx
     &     - 14.D0/3.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-1)*
     & R6 * ( 380.D0/3.D0 - 604.D0/3.D0*Li2mx + 380.D0/3.D0*lnx - 232.D0
     &    /3.D0*lnopxon2 - 620.D0/3.D0*ln2 + 604.D0/3.D0*Li2omx - 151.D0
     &    /9.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*omx**(-1)*
     & R8 * ( 107741.D0/216.D0 - 5875.D0/9.D0*Li2mx + 3745.D0/9.D0*lnx
     &     - 863.D0/3.D0*lnopxon2 - 7399.D0/9.D0*ln2 + 5875.D0/9.D0*
     &    Li2omx - 5875.D0/108.D0*pisq )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*R2 * (  - 
     &    95.D0/9.D0 + 16*Li2mx + 32.D0/3.D0*ln2 - 16*Li2omx + 4.D0/3.D0
     &    *pisq - 49.D0/9.D0*x - 16*x*Li2mx + 16*x*lnopxon2 + 64.D0/3.D0
     &    *x*ln2 + 16*x*Li2omx - 4.D0/3.D0*x*pisq + 16*x**2*lnx - 16*
     &    x**2*lnopxon2 - 16*x**2*ln2 )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*R4 * (  - 
     &    111511.D0/3600.D0 + 4*opx**(-2) + 9*opx**(-1) + 24*Li2mx - 40
     &    *lnx + 8*lnopxon2 + 392.D0/15.D0*ln2 - 24*Li2omx + 2*pisq + 
     &    821.D0/300.D0*x - 8*x*Li2mx - 32*x*lnx + 24*x*lnopxon2 + 182.D
     &    0/5.D0*x*ln2 + 8*x*Li2omx - 2.D0/3.D0*x*pisq + 41959.D0/1800.D
     &    0*x**2 - 16*x**2*Li2mx - 346.D0/15.D0*x**2*ln2 + 16*x**2*
     &    Li2omx - 4.D0/3.D0*x**2*pisq + 7141.D0/3600.D0*x**3 - 32.D0/
     &    15.D0*x**3*ln2 )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*R6 * (  - 
     &    12550123.D0/211680.D0 + opx**(-4) - 7.D0/12.D0*opx**(-3) + 7.D
     &    0/12.D0*opx**(-2) + 5.D0/6.D0*opx**(-1) + 188.D0/3.D0*Li2mx
     &     - 148.D0/3.D0*lnx + 88.D0/3.D0*lnopxon2 + 5302.D0/63.D0*ln2
     &     - 188.D0/3.D0*Li2omx + 47.D0/9.D0*pisq - 41503291.D0/2116800.
     &    D0*x + 28.D0/3.D0*x*Li2mx - 24*x*lnx + 44.D0/3.D0*x*lnopxon2
     &     + 10259.D0/315.D0*x*ln2 - 28.D0/3.D0*x*Li2omx + 7.D0/9.D0*x*
     &    pisq + 2647843.D0/235200.D0*x**2 - 8*x**2*Li2mx - 28.D0/3.D0*
     &    x**2*lnx + 28.D0/3.D0*x**2*lnopxon2 - 37.D0/35.D0*x**2*ln2 + 
     &    8*x**2*Li2omx - 2.D0/3.D0*x**2*pisq + 2084269.D0/529200.D0*
     &    x**3 - 1364.D0/315.D0*x**3*ln2 - 13424.D0/33075.D0*x**4 + 124.
     &    D0/315.D0*x**4*ln2 )
      DelI2bit(QQBARS,rglr) = DelI2bit(QQBARS,rglr) + CF*TF*R8 * (  - 
     &    9331105411.D0/58060800.D0 + 5.D0/12.D0*opx**(-6) - 13.D0/16.D0
     &    *opx**(-5) + 55.D0/96.D0*opx**(-4) - 83.D0/576.D0*opx**(-3)
     &     - 127.D0/1152.D0*opx**(-2) - 491.D0/1152.D0*opx**(-1) + 485.D
     &    0/3.D0*Li2mx - 970.D0/9.D0*lnx + 245.D0/3.D0*lnopxon2 + 33946.
     &    D0/135.D0*ln2 - 485.D0/3.D0*Li2omx + 485.D0/36.D0*pisq - 
     &    8188909919.D0/203212800.D0*x + 239.D0/9.D0*x*Li2mx - 203.D0/9.
     &    D0*x*lnx + 160.D0/9.D0*x*lnopxon2 + 72539.D0/1260.D0*x*ln2 - 
     &    239.D0/9.D0*x*Li2omx + 239.D0/108.D0*x*pisq - 27275893.D0/
     &    20321280.D0*x**2 - 2.D0/9.D0*x**2*Li2mx - 47.D0/9.D0*x**2*lnx
     &     + 47.D0/9.D0*x**2*lnopxon2 + 5399.D0/756.D0*x**2*ln2 + 2.D0/
     &    9.D0*x**2*Li2omx - 1.D0/54.D0*x**2*pisq + 10904881.D0/5419008.
     &    D0*x**3 - 428.D0/189.D0*x**3*ln2 - 8157791.D0/12700800.D0*
     &    x**4 + 403.D0/630.D0*x**4*ln2 + 2982737.D0/33868800.D0*x**5
     &     - 17.D0/210.D0*x**5*ln2 )

      DelI2bit(QQBARNS,rglr) =  + CF**2*omx**(-5)*R8 * ( 16*Li2mx - 32.D
     &    0/3.D0*lnx - 16*Li2omx + 4.D0/3.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-4)
     & *R8 * ( 16.D0/3.D0 - 80*Li2mx + 160.D0/3.D0*lnx - 8*lnopxon2 - 
     &    16*ln2 + 80*Li2omx - 20.D0/3.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-3)
     & *R6 * ( 16.D0/3.D0*Li2mx - 4*lnx - 16.D0/3.D0*Li2omx + 4.D0/9.D0
     &    *pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-3)
     & *R8 * (  - 28 + 496.D0/3.D0*Li2mx - 994.D0/9.D0*lnx + 34*
     &    lnopxon2 + 72*ln2 - 496.D0/3.D0*Li2omx + 124.D0/9.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-2)
     & *R6 * ( 4.D0/3.D0 - 56.D0/3.D0*Li2mx + 14*lnx - 8.D0/3.D0*
     &    lnopxon2 - 16.D0/3.D0*ln2 + 56.D0/3.D0*Li2omx - 14.D0/9.D0*
     &    pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-2)
     & *R8 * ( 532.D0/9.D0 - 542.D0/3.D0*Li2mx + 2183.D0/18.D0*lnx - 
     &    172.D0/3.D0*lnopxon2 - 392.D0/3.D0*ln2 + 542.D0/3.D0*Li2omx
     &     - 271.D0/18.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-1)
     & *R4 * (  - 4*lnx )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-1)
     & *R6 * (  - 16.D0/3.D0 + 24*Li2mx - 56.D0/3.D0*lnx + 8*lnopxon2
     &     + 16*ln2 - 24*Li2omx + 2*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*omx**(-1)
     & *R8 * (  - 2285.D0/36.D0 + 110*Li2mx - 671.D0/9.D0*lnx + 48*
     &    lnopxon2 + 362.D0/3.D0*ln2 - 110*Li2omx + 55.D0/6.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*R2 * ( 8
     &     - 16*opx**(-1)*Li2mx + 16*opx**(-1)*Li2omx + 16*Li2mx - 8*
     &    ln2 - 16*Li2omx - 4.D0/3.D0*pisq*opx**(-1) + 4.D0/3.D0*pisq
     &     - 8*x + 8*x*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*R4 * ( 
     &     - 65.D0/36.D0 - 2*opx**(-2) + 3*opx**(-1) + 4*lnx + 2.D0/3.D0
     &    *ln2 - 23.D0/12.D0*x + 4*x*lnx - 4*x*lnopxon2 - 2*x*ln2 - 2*
     &    x**2 + 2*x**2*ln2 + 13.D0/18.D0*x**3 - 2.D0/3.D0*x**3*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*R6 * ( 
     &    158293.D0/43200.D0 - 1.D0/6.D0*opx**(-4) + 5.D0/12.D0*
     &    opx**(-3) - 2.D0/3.D0*opx**(-2) + 3.D0/4.D0*opx**(-1) - 32.D0/
     &    3.D0*Li2mx + 26.D0/3.D0*lnx - 16.D0/3.D0*lnopxon2 - 479.D0/45.
     &    D0*ln2 + 32.D0/3.D0*Li2omx - 8.D0/9.D0*pisq + 40141.D0/14400.D
     &    0*x - 8.D0/3.D0*x*Li2mx + 8.D0/3.D0*x*lnx - 8.D0/3.D0*x*
     &    lnopxon2 - 83.D0/15.D0*x*ln2 + 8.D0/3.D0*x*Li2omx - 2.D0/9.D0
     &    *x*pisq - 13.D0/3600.D0*x**2 + 1.D0/15.D0*x**2*ln2 + 364.D0/
     &    675.D0*x**3 - 23.D0/45.D0*x**3*ln2 - 541.D0/2700.D0*x**4 + 8.D
     &    0/45.D0*x**4*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CF**2*R8 * ( 
     &    136552219.D0/5080320.D0 - 1.D0/24.D0*opx**(-6) + 7.D0/48.D0*
     &    opx**(-5) - 73.D0/288.D0*opx**(-4) + 173.D0/576.D0*opx**(-3)
     &     - 29.D0/96.D0*opx**(-2) + 173.D0/576.D0*opx**(-1) - 92.D0/3.D
     &    0*Li2mx + 379.D0/18.D0*lnx - 50.D0/3.D0*lnopxon2 - 34775.D0/
     &    756.D0*ln2 + 92.D0/3.D0*Li2omx - 23.D0/9.D0*pisq + 7517387.D0/
     &    940800.D0*x - 14.D0/3.D0*x*Li2mx + 10.D0/3.D0*x*lnx - 10.D0/3.
     &    D0*x*lnopxon2 - 4769.D0/420.D0*x*ln2 + 14.D0/3.D0*x*Li2omx - 
     &    7.D0/18.D0*x*pisq + 380029.D0/627200.D0*x**2 - 781.D0/1260.D0
     &    *x**2*ln2 + 4610261.D0/50803200.D0*x**3 - 377.D0/3780.D0*x**3
     &    *ln2 - 3319.D0/22050.D0*x**4 + 29.D0/210.D0*x**4*ln2 + 1941.D0
     &    /39200.D0*x**5 - 3.D0/70.D0*x**5*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-5)
     & *R8 * (  - 8*Li2mx + 16.D0/3.D0*lnx + 8*Li2omx - 2.D0/3.D0*pisq
     &     )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-4)
     & *R8 * (  - 8.D0/3.D0 + 40*Li2mx - 80.D0/3.D0*lnx + 4*lnopxon2 + 
     &    8*ln2 - 40*Li2omx + 10.D0/3.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-3)
     & *R6 * (  - 8.D0/3.D0*Li2mx + 2*lnx + 8.D0/3.D0*Li2omx - 2.D0/9.D0
     &    *pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-3)
     & *R8 * ( 14 - 248.D0/3.D0*Li2mx + 497.D0/9.D0*lnx - 17*lnopxon2
     &     - 36*ln2 + 248.D0/3.D0*Li2omx - 62.D0/9.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-2)
     & *R6 * (  - 2.D0/3.D0 + 28.D0/3.D0*Li2mx - 7*lnx + 4.D0/3.D0*
     &    lnopxon2 + 8.D0/3.D0*ln2 - 28.D0/3.D0*Li2omx + 7.D0/9.D0*pisq
     &     )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-2)
     & *R8 * (  - 266.D0/9.D0 + 271.D0/3.D0*Li2mx - 2183.D0/36.D0*lnx
     &     + 86.D0/3.D0*lnopxon2 + 196.D0/3.D0*ln2 - 271.D0/3.D0*Li2omx
     &     + 271.D0/36.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-1)
     & *R4 * ( 2*lnx )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-1)
     & *R6 * ( 8.D0/3.D0 - 12*Li2mx + 28.D0/3.D0*lnx - 4*lnopxon2 - 8*
     &    ln2 + 12*Li2omx - pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*omx**(-1)
     & *R8 * ( 2285.D0/72.D0 - 55*Li2mx + 671.D0/18.D0*lnx - 24*
     &    lnopxon2 - 181.D0/3.D0*ln2 + 55*Li2omx - 55.D0/12.D0*pisq )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*R2 * ( 
     &     - 4 + 8*opx**(-1)*Li2mx - 8*opx**(-1)*Li2omx - 8*Li2mx + 4*
     &    ln2 + 8*Li2omx + 2.D0/3.D0*pisq*opx**(-1) - 2.D0/3.D0*pisq + 
     &    4*x - 4*x*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*R4 * ( 65.
     &    D0/72.D0 + opx**(-2) - 3.D0/2.D0*opx**(-1) - 2*lnx - 1.D0/3.D0
     &    *ln2 + 23.D0/24.D0*x - 2*x*lnx + 2*x*lnopxon2 + x*ln2 + x**2
     &     - x**2*ln2 - 13.D0/36.D0*x**3 + 1.D0/3.D0*x**3*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*R6 * ( 
     &     - 158293.D0/86400.D0 + 1.D0/12.D0*opx**(-4) - 5.D0/24.D0*
     &    opx**(-3) + 1.D0/3.D0*opx**(-2) - 3.D0/8.D0*opx**(-1) + 16.D0/
     &    3.D0*Li2mx - 13.D0/3.D0*lnx + 8.D0/3.D0*lnopxon2 + 479.D0/90.D
     &    0*ln2 - 16.D0/3.D0*Li2omx + 4.D0/9.D0*pisq - 40141.D0/28800.D0
     &    *x + 4.D0/3.D0*x*Li2mx - 4.D0/3.D0*x*lnx + 4.D0/3.D0*x*
     &    lnopxon2 + 83.D0/30.D0*x*ln2 - 4.D0/3.D0*x*Li2omx + 1.D0/9.D0
     &    *x*pisq + 13.D0/7200.D0*x**2 - 1.D0/30.D0*x**2*ln2 - 182.D0/
     &    675.D0*x**3 + 23.D0/90.D0*x**3*ln2 + 541.D0/5400.D0*x**4 - 4.D
     &    0/45.D0*x**4*ln2 )
      DelI2bit(QQBARNS,rglr) = DelI2bit(QQBARNS,rglr) + CA*CF*R8 * ( 
     &     - 136552219.D0/10160640.D0 + 1.D0/48.D0*opx**(-6) - 7.D0/96.D
     &    0*opx**(-5) + 73.D0/576.D0*opx**(-4) - 173.D0/1152.D0*
     &    opx**(-3) + 29.D0/192.D0*opx**(-2) - 173.D0/1152.D0*opx**(-1)
     &     + 46.D0/3.D0*Li2mx - 379.D0/36.D0*lnx + 25.D0/3.D0*lnopxon2
     &     + 34775.D0/1512.D0*ln2 - 46.D0/3.D0*Li2omx + 23.D0/18.D0*
     &    pisq - 7517387.D0/1881600.D0*x + 7.D0/3.D0*x*Li2mx - 5.D0/3.D0
     &    *x*lnx + 5.D0/3.D0*x*lnopxon2 + 4769.D0/840.D0*x*ln2 - 7.D0/3.
     &    D0*x*Li2omx + 7.D0/36.D0*x*pisq - 380029.D0/1254400.D0*x**2
     &     + 781.D0/2520.D0*x**2*ln2 - 4610261.D0/101606400.D0*x**3 + 
     &    377.D0/7560.D0*x**3*ln2 + 3319.D0/44100.D0*x**4 - 29.D0/420.D0
     &    *x**4*ln2 - 1941.D0/78400.D0*x**5 + 3.D0/140.D0*x**5*ln2 )

      endif


! Add boundary terms for both x < xmin and x > xmin
      DelI2bit(GGCACA,rglr)=DelI2bit(GGCACA,rglr)-8*CA**2*pi**3
     & *(x*BoundaryConditionGGCACA(x)-BoundaryConditionGGCACA(one))/omx
      DelI2bit(GGCATF,rglr)=DelI2bit(GGCATF,rglr)-8*CA*TF*nf*pi**3
     & *(x*BoundaryConditionGGCATF(x)-BoundaryConditionGGCATF(one))/omx
!   BoundaryConditionsGGCFTF(one)=0
      DelI2bit(GGCFTF,rglr)=DelI2bit(GGCFTF,rglr)-8*CF*TF*nf*pi**3
     & *x*BoundaryConditionGGCFTF(x)/omx
      DelI2bit(QQCACF,rglr)=DelI2bit(QQCACF,rglr)-8*CA*CF*pi**3
     & *(BoundaryConditionQQCACF(x)-BoundaryConditionQQCACF(one))/omx
!   BoundaryConditionsQQCFCF(one)=0
      DelI2bit(QQCFCF,rglr)=DelI2bit(QQCFCF,rglr)-8*CF**2*pi**3
     & *x*BoundaryConditionQQCFCF(x)/omx
      DelI2bit(QQCFTF,rglr)=DelI2bit(QQCFTF,rglr)-8*CF*TF*pi**3
     & *(BoundaryConditionQQCFTF(x)-BoundaryConditionQQCFTF(one))/omx
!   BoundaryConditionsGQCFCF(one)=0
      DelI2bit(GQCFCA,rglr)=DelI2bit(GQCFCA,rglr)-8*CA*CF*pi**3
     & *x*BoundaryConditionGQCFCA(x)/omx
!   BoundaryConditionsGQCFCA(one)=0
      DelI2bit(GQCFCF,rglr)=DelI2bit(GQCFCF,rglr)-8*CF**2*pi**3
     & *x*BoundaryConditionGQCFCF(x)/omx
!   BoundaryConditionsQGCATF(one)=0
      DelI2bit(QGCATF,rglr)=DelI2bit(QGCATF,rglr)-8*CA*TF*pi**3
     & *x*BoundaryConditionQGCATF(x)/omx
!   BoundaryConditionsQGCFTF(one)=0
      DelI2bit(QGCFTF,rglr)=DelI2bit(QGCFTF,rglr)-8*CF*TF*pi**3
     & *x*BoundaryConditionQGCFTF(x)/omx
!   BoundaryConditionsQQBARS(one)=0
      DelI2bit(QQBARS,rglr)=DelI2bit(QQBARS,rglr)-8*CF*TF*pi**3
     & *x*BoundaryConditionQQBARS(x)/omx
!   BoundaryConditionsQQBARNS(one)=0
      DelI2bit(QQBARNS,rglr)=DelI2bit(QQBARNS,rglr)
     & +16*CF*(CF-CA/2)*pi**3*x*BoundaryConditionQQBARNS(x)/omx


      DelI2bit(GGCACA,plus) =  + CA**2 * ( 26.D0/9.D0 - 548.D0/9.D0*ln2
     &     + 176.D0/3.D0*ln2**2 + 16.D0/3.D0*pisq*ln2 - 8*
     &    BoundaryConditionGGCACA(one)*pi**3 )
      DelI2bit(GGCACA,plus) = DelI2bit(GGCACA,plus) + CA**2*R2 * ( 1429.
     &    D0/5400.D0 + 104.D0/45.D0*ln2 - 2*pisq )
      DelI2bit(GGCACA,plus) = DelI2bit(GGCACA,plus) + CA**2*R4 * ( 
     &    16018321.D0/12700800.D0 + 587.D0/3780.D0*ln2 - 1.D0/108.D0*
     &    pisq )
      DelI2bit(GGCACA,plus) = DelI2bit(GGCACA,plus) + CA**2*R6 * ( 
     &    74801417.D0/3048192000.D0 - 23.D0/2100.D0*ln2 )
      DelI2bit(GGCACA,plus) = DelI2bit(GGCACA,plus) + CA**2*lnRonpi
     &  * ( 524.D0/9.D0 - 176.D0/3.D0*ln2 - 16.D0/3.D0*pisq )

      DelI2bit(GGCATF,plus) =  + CA*TF*nf * (  - 52.D0/9.D0 + 232.D0/9.D
     &    0*ln2 - 64.D0/3.D0*ln2**2 - 8*BoundaryConditionGGCATF(one)*
     &    pi**3 )
      DelI2bit(GGCATF,plus) = DelI2bit(GGCATF,plus) + CA*TF*nf*R2 * ( 
     &    3071.D0/2700.D0 - 28.D0/45.D0*ln2 )
      DelI2bit(GGCATF,plus) = DelI2bit(GGCATF,plus) + CA*TF*nf*R4 * ( 
     &     - 168401.D0/3175200.D0 + 53.D0/945.D0*ln2 )
      DelI2bit(GGCATF,plus) = DelI2bit(GGCATF,plus) + CA*TF*nf*R6 * ( 
     &    7001023.D0/1524096000.D0 - 11.D0/3150.D0*ln2 )
      DelI2bit(GGCATF,plus) = DelI2bit(GGCATF,plus) + CA*TF*nf*R8 * ( 
     &     - 5664846191.D0/17703899136000.D0 + 4001.D0/14968800.D0*ln2
     &     )
      DelI2bit(GGCATF,plus) = DelI2bit(GGCATF,plus) + CA*TF*nf*lnRonpi
     &  * (  - 184.D0/9.D0 + 64.D0/3.D0*ln2 )

      DelI2bit(GGCFTF,plus) =  0

      DelI2bit(QQCACF,plus) =  + CA*CF * ( 26.D0/9.D0 - 548.D0/9.D0*ln2
     &     + 176.D0/3.D0*ln2**2 + 16.D0/3.D0*pisq*ln2 - 8*
     &    BoundaryConditionQQCACF(one)*pi**3 )
      DelI2bit(QQCACF,plus) = DelI2bit(QQCACF,plus) + CA*CF*R2 * ( 1429.
     &    D0/5400.D0 + 104.D0/45.D0*ln2 + 2.D0/3.D0*pisq )
      DelI2bit(QQCACF,plus) = DelI2bit(QQCACF,plus) + CA*CF*R4 * (  - 
     &    9383279.D0/12700800.D0 + 587.D0/3780.D0*ln2 - 1.D0/108.D0*
     &    pisq )
      DelI2bit(QQCACF,plus) = DelI2bit(QQCACF,plus) + CA*CF*R6 * ( 
     &    74801417.D0/3048192000.D0 - 23.D0/2100.D0*ln2 )
      DelI2bit(QQCACF,plus) = DelI2bit(QQCACF,plus) + CA*CF*R8 * (  - 
     &    50937246539.D0/70815596544000.D0 + 28529.D0/59875200.D0*ln2
     &     - 1.D0/777600.D0*pisq )
      DelI2bit(QQCACF,plus) = DelI2bit(QQCACF,plus) + CA*CF*lnRonpi
     &  * ( 524.D0/9.D0 - 176.D0/3.D0*ln2 - 16.D0/3.D0*pisq )

      DelI2bit(QQCFCF,plus) =  + CF**2*R2 * (  - 8.D0/3.D0*pisq )
      DelI2bit(QQCFCF,plus) = DelI2bit(QQCFCF,plus) + CF**2*R4 * ( 2 )

      DelI2bit(QQCFTF,plus) =  + CF*TF * (  - 52.D0/9.D0 + 232.D0/9.D0*
     &    ln2 - 64.D0/3.D0*ln2**2 - 8*BoundaryConditionQQCFTF(one)*
     &    pi**3 )
      DelI2bit(QQCFTF,plus) = DelI2bit(QQCFTF,plus) + CF*TF*R2 * ( 3071.
     &    D0/2700.D0 - 28.D0/45.D0*ln2 )
      DelI2bit(QQCFTF,plus) = DelI2bit(QQCFTF,plus) + CF*TF*R4 * (  - 
     &    168401.D0/3175200.D0 + 53.D0/945.D0*ln2 )
      DelI2bit(QQCFTF,plus) = DelI2bit(QQCFTF,plus) + CF*TF*R6 * ( 
     &    7001023.D0/1524096000.D0 - 11.D0/3150.D0*ln2 )
      DelI2bit(QQCFTF,plus) = DelI2bit(QQCFTF,plus) + CF*TF*R8 * (  - 
     &    5664846191.D0/17703899136000.D0 + 4001.D0/14968800.D0*ln2 )
      DelI2bit(QQCFTF,plus) = DelI2bit(QQCFTF,plus) + CF*TF*lnRonpi
     &  * (  - 184.D0/9.D0 + 64.D0/3.D0*ln2 )


      DelI2bit(GGCACA,delt) =  + CA**2 * ( 3220.D0/27.D0*LQ - 524.D0/9.D
     &    0*LQ*lnR - 560.D0/9.D0*LQ*ln2 + 176.D0/3.D0*LQ*ln2*lnR - 176.D
     &    0/3.D0*LQ*ln2**2 - 44.D0/9.D0*pisq*LQ + 16.D0/3.D0*pisq*LQ*
     &    lnR - 16*zeta3*LQ )
      DelI2bit(GGCACA,delt) = DelI2bit(GGCACA,delt) + CA**2*R2 * (  - 
     &    1429.D0/5400.D0*LQ - 104.D0/45.D0*LQ*ln2 + 2*pisq*LQ )
      DelI2bit(GGCACA,delt) = DelI2bit(GGCACA,delt) + CA**2*R4 * (  - 
     &    16018321.D0/12700800.D0*LQ - 587.D0/3780.D0*LQ*ln2 + 1.D0/108.
     &    D0*pisq*LQ )
      DelI2bit(GGCACA,delt) = DelI2bit(GGCACA,delt) + CA**2*R6 * (  - 
     &    74801417.D0/3048192000.D0*LQ + 23.D0/2100.D0*LQ*ln2 )
      DelI2bit(GGCACA,delt) = DelI2bit(GGCACA,delt) + CA**2*R8 * ( 
     &    50937246539.D0/70815596544000.D0*LQ - 28529.D0/59875200.D0*LQ
     &    *ln2 + 1.D0/777600.D0*pisq*LQ )

      DelI2bit(GGCATF,delt) =  + CA*TF*nf * (  - 1256.D0/27.D0*LQ + 184.
     &    D0/9.D0*LQ*lnR + 256.D0/9.D0*LQ*ln2 - 64.D0/3.D0*LQ*ln2*lnR
     &     + 64.D0/3.D0*LQ*ln2**2 + 16.D0/9.D0*pisq*LQ )
      DelI2bit(GGCATF,delt) = DelI2bit(GGCATF,delt) + CA*TF*nf*R2 * ( 
     &     - 3071.D0/2700.D0*LQ + 28.D0/45.D0*LQ*ln2 )
      DelI2bit(GGCATF,delt) = DelI2bit(GGCATF,delt) + CA*TF*nf*R4 * ( 
     &    168401.D0/3175200.D0*LQ - 53.D0/945.D0*LQ*ln2 )
      DelI2bit(GGCATF,delt) = DelI2bit(GGCATF,delt) + CA*TF*nf*R6 * ( 
     &     - 7001023.D0/1524096000.D0*LQ + 11.D0/3150.D0*LQ*ln2 )
      DelI2bit(GGCATF,delt) = DelI2bit(GGCATF,delt) + CA*TF*nf*R8 * ( 
     &    5664846191.D0/17703899136000.D0*LQ - 4001.D0/14968800.D0*LQ*
     &    ln2 )

      DelI2bit(GGCFTF,delt) =  0

      DelI2bit(QQCACF,delt) =  + CA*CF * ( 3220.D0/27.D0*LQ - 524.D0/9.D
     &    0*LQ*lnR - 560.D0/9.D0*LQ*ln2 + 176.D0/3.D0*LQ*ln2*lnR - 176.D
     &    0/3.D0*LQ*ln2**2 - 44.D0/9.D0*pisq*LQ + 16.D0/3.D0*pisq*LQ*
     &    lnR - 16*zeta3*LQ )
      DelI2bit(QQCACF,delt) = DelI2bit(QQCACF,delt) + CA*CF*R2 * (  - 
     &    1429.D0/5400.D0*LQ - 104.D0/45.D0*LQ*ln2 - 2.D0/3.D0*pisq*LQ
     &     )
      DelI2bit(QQCACF,delt) = DelI2bit(QQCACF,delt) + CA*CF*R4 * ( 
     &    9383279.D0/12700800.D0*LQ - 587.D0/3780.D0*LQ*ln2 + 1.D0/108.D
     &    0*pisq*LQ )
      DelI2bit(QQCACF,delt) = DelI2bit(QQCACF,delt) + CA*CF*R6 * (  - 
     &    74801417.D0/3048192000.D0*LQ + 23.D0/2100.D0*LQ*ln2 )
      DelI2bit(QQCACF,delt) = DelI2bit(QQCACF,delt) + CA*CF*R8 * ( 
     &    50937246539.D0/70815596544000.D0*LQ - 28529.D0/59875200.D0*LQ
     &    *ln2 + 1.D0/777600.D0*pisq*LQ )

      DelI2bit(QQCFCF,delt) =  + CF**2*R2 * ( 8.D0/3.D0*pisq*LQ )
      DelI2bit(QQCFCF,delt) = DelI2bit(QQCFCF,delt) + CF**2*R4 * (  - 2
     &    *LQ )

      DelI2bit(QQCFTF,delt) =  + CF*TF * (  - 1256.D0/27.D0*LQ + 184.D0/
     &    9.D0*LQ*lnR + 256.D0/9.D0*LQ*ln2 - 64.D0/3.D0*LQ*ln2*lnR + 64.
     &    D0/3.D0*LQ*ln2**2 + 16.D0/9.D0*pisq*LQ )
      DelI2bit(QQCFTF,delt) = DelI2bit(QQCFTF,delt) + CF*TF*R2 * (  - 
     &    3071.D0/2700.D0*LQ + 28.D0/45.D0*LQ*ln2 )
      DelI2bit(QQCFTF,delt) = DelI2bit(QQCFTF,delt) + CF*TF*R4 * ( 
     &    168401.D0/3175200.D0*LQ - 53.D0/945.D0*LQ*ln2 )
      DelI2bit(QQCFTF,delt) = DelI2bit(QQCFTF,delt) + CF*TF*R6 * (  - 
     &    7001023.D0/1524096000.D0*LQ + 11.D0/3150.D0*LQ*ln2 )
      DelI2bit(QQCFTF,delt) = DelI2bit(QQCFTF,delt) + CF*TF*R8 * ( 
     &    5664846191.D0/17703899136000.D0*LQ - 4001.D0/14968800.D0*LQ*
     &    ln2 )


c Implementation of Eq. (3.36) of arXiV:2207.07037
      DeltaI2(qq,:)=DelI2bit(QQCFCF,:)+DelI2bit(QQCACF,:)
     &             +nf*DelI2bit(QQCFTF,:)+DelI2bit(QQBARS,:)
      DeltaI2(qbq,:)=DelI2bit(QQBARS,:)+DelI2bit(QQBARNS,:)
      DeltaI2(qbpq,:)=DelI2bit(QQBARS,:)
      DeltaI2(qg,:)=DelI2bit(QGCATF,:)+DelI2bit(QGCFTF,:)
      DeltaI2(gq,:)=DelI2bit(GQCFCA,:)+DelI2bit(GQCFCF,:)
      DeltaI2(gg,:)=DelI2bit(GGCACA,:)+DelI2bit(GGCATF,:)
     &             +DelI2bit(GGCFTF,:)
      DeltaI2(qpq,:)=DeltaI2(qbpq,:)

      return
      end subroutine DeltaI2R
