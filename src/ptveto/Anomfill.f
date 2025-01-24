!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine Anomfill(nfl,
     & beta0,beta1,beta2,beta3,beta4,
     & GammaA0,GammaA1,GammaA2,GammaA3,GammaA4,
     & gammagq0,gammagq1,gammagq2,gammagq3,
     & gamma0S,gamma1S,gamma2S,gamma3S)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      integer::nfl,k
      real(dp):: beta0,beta1,beta2,beta3,beta4
      real(dp):: CR(0:1),GammaA0(0:1),GammaA1(0:1),GammaA2(0:1),GammaA3(0:1),GammaA4(0:1)
      real(dp):: gammagq0(0:1),gammagq1(0:1),gammagq2(0:1),gammagq3(0:1)
      real(dp):: gamma0S,gamma1S,gamma2S,gamma3S
      real(dp):: d4AA,d4AF,d4FF,NA,NF,d4RA,d4RF,NR
      real(dp), parameter :: dAANA = 135._dp/8._dp
      real(dp), parameter :: dRANA = 15._dp/16._dp
      real(dp), parameter :: dRRNA = 5._dp/96._dp

!     all of these factors should be put in a separate routine and initialized once and for all
!     Taken from 9701390, Eq11
      beta0=11._dp-2._dp*nfl/3._dp
      beta1=102._dp-38*nfl/3._dp
      beta2=2857/2._dp-5033/18._dp*nfl+325/54._dp*nfl**2
      beta3=149753._dp/6._dp+3564*zeta3
     & -(1078361/162._dp+6508/27._dp*zeta3)*nfl
     & +(50065/162._dp+6472/81._dp*zeta3)*nfl**2
     & +1093/729._dp*nfl**3
      beta4 = (8157455._dp/16._dp+621885._dp/2._dp*zeta3-88209._dp/2._dp*zeta4-288090._dp*zeta5 
     & +nfl*(-336460813._dp/1944._dp-4811164._dp/81._dp*zeta3 
     & +33935._dp/6._dp*zeta4+1358995._dp/27._dp*zeta5) 
     & +nfl**2*(25960913._dp/1944._dp+698531._dp/81._dp*zeta3-10526._dp/9._dp*zeta4-381760._dp/81._dp*zeta5) 
     & +nfl**3*(-630559._dp/5832._dp-48722._dp/243._dp*zeta3+1618._dp/27._dp*zeta4+460._dp/9._dp*zeta5) 
     & +nfl**4*(1205._dp/2916._dp-152._dp/81._dp*zeta3))

      CR(0)=CA
      CR(1)=CF
! Taken from ff4l.m, supplemental material for 2002.04617,
! or equivalently 1911.10174, Eq 6.5, noting that d4 factors here
! are just sums d(abcd)d(abcd), i.e. not normalized by NA
      NA=8._dp    ! N^2-1
      NF=3._dp    ! N
      d4AA=135d0        ! d4AA=NA*N^2*(N^2+36)/24
      d4AF=15/2d0       ! d4AF=NA*N*(N^2+6)/48;
      d4FF=5/12d0       ! d4FF=NA*(N^4-6*N^2+18)/96/N^2;
      do k=0,1
      GammaA0(k)=4*CR(k)
      GammaA1(k)=16*CR(k)*(67/12d0-pisq/4d0-nfl*5/18d0)
      GammaA2(k)=64*CR(k)
     & *((735/32d0+33/8d0*zeta3-67/24d0*pisq+11/80d0*pisq**2)
     & +nfl*(-319/144d0-13/12d0*zeta3+5/36d0*pisq)
     & -nfl**2/108d0)
!      GammaA3(k)=256*CR(k)
!     & *(42139/384d0+1309/16d0*zeta3-1353/32d0*zeta5-27/16d0*zeta3**2
!     & -5525/288d0*pisq-33/16d0*pisq*zeta3+1353/640d0*pisq**2
!     & -313/3360d0*pisq**3
!     & +nfl*(-344345/20736d0+1219/144d0*zeta5-2405/108d0*zeta3
!     & +50/27d0*pisq+13/24d0*pisq*zeta3-77/1440d0*pisq**2)
!     & +nfl**2*(+17875/62208d0+65/108d0*zeta3-19/1296d0*pisq
!     & -13/4320d0*pisq**2)
!     & +nfl**3*(-1/648d0+zeta3/108d0))
!      write(6,*) 'old GammaA3',k,GammaA3(k)
! This version is correct only for k = 0
!      GammaA3(k)=256*CR(k)
!     & *(42139/384d0 - 1415/72d0*pisq + 1353/640d0*pisq**2 - 781/6720d0*pisq**3
!     &  -344345/20736d0*nfl + 1645/864d0*nfl*pisq - 77/1440d0*nfl*pisq**2
!     & + 17875/62208d0*nfl**2 - 19/1296d0*nfl**2*pisq - 13/4320d0*nfl**2*pisq**2
!     & - 1/648d0*nfl**3 - 33/2d0*zeta5 + 143/18d0*zeta5*nfl + 331/4d0*zeta3
!     & - 33/16d0*zeta3*pisq - 9665/432d0*zeta3*nfl + 13/24d0*zeta3*nfl*pisq
!     & + 65/108d0*zeta3*nfl**2 + 1d0/108d0*zeta3*nfl**3 - 81/8d0*zeta3**2)
!      write(6,*) 'new GammaA3',k,GammaA3(k)
      if (k == 0) then
        d4RA=d4AA
        d4RF=d4AF
        NR=NA
      else
        d4RA=d4AF
        d4RF=d4FF
        NR=NF
      endif
      GammaA3(k)=CR(k)
     & *( 84278/3._dp - 44200/9._dp*pisq + 2706/5._dp*pisq**2 - 2504/105._dp*pisq**3 - 
     &   344345/81._dp*nfl + 12800/27._dp*nfl*pisq - 616/45._dp*nfl*pisq**2 + 17875/243._dp*nfl**2
     &    - 304/81._dp*nfl**2*pisq - 104/135._dp*nfl**2*pisq**2 - 32/81._dp*nfl**3 - 10824*zeta5
     &    + 19504/9._dp*zeta5*nfl + 20944*zeta3 - 528*zeta3*pisq - 153920/27._dp*zeta3*nfl
     &    + 416/3._dp*zeta3*nfl*pisq + 4160/27._dp*zeta3*nfl**2
     &    + 64/27._dp*zeta3*nfl**3 - 432*zeta3**2 )
     & + d4RA/NR* (  - 64/3._dp*pisq - 992/945._dp*pisq**3 + 3520/3._dp*zeta5
     &               + 128/3._dp*zeta3 - 384*zeta3**2 )
     & + d4RF/NR * ( 128/3._dp*nfl*pisq - 1280/3._dp*zeta5*nfl - 256/3._dp*zeta3*nfl )

      enddo
! Numerical values taken from arXiv:1812.11818, Eq. (13)
      do k=0,1
!        GammaA0(k)=0.42441_dp*fourpi/CF*CR(k)
        if (nfl == 3) then
!          GammaA1(k)=0.42441_dp*0.7266_dp*fourpi**2/CF*CR(k)
!          GammaA2(k)=0.42441_dp*0.7341_dp*fourpi**3/CF*CR(k)
!          GammaA3(k)=0.42441_dp*0.665_dp*fourpi**4/CF*CR(k)
          GammaA4(k)=0.42441_dp*1.3_dp*fourpi**5/CF*CR(k)
        elseif (nfl == 4) then
!          GammaA1(k)=0.42441_dp*0.6382_dp*fourpi**2/CF*CR(k)
!          GammaA2(k)=0.42441_dp*0.5100_dp*fourpi**3/CF*CR(k)
!          GammaA3(k)=0.42441_dp*0.317_dp*fourpi**4/CF*CR(k)
          GammaA4(k)=0.42441_dp*0.8_dp*fourpi**5/CF*CR(k)
        elseif (nfl == 5) then
!          GammaA1(k)=0.42441_dp*0.5497_dp*fourpi**2/CF*CR(k)
!          GammaA2(k)=0.42441_dp*0.2840_dp*fourpi**3/CF*CR(k)
!          GammaA3(k)=0.42441_dp*0.013_dp*fourpi**4/CF*CR(k)
          GammaA4(k)=0.42441_dp*0.5_dp*fourpi**5/CF*CR(k)
        else
          write(6,*) 'Gamma4 not defined for nfl = ',nfl
          call exit(1)
      endif
      enddo
!    Taken from 0903.1126, Eqs (A4) and (A5)
      gammagq0(0)=-beta0;
      gammagq1(0)=CA**2*(-692/27._dp+11*pisq/18._dp+2*zeta3)
     & +CA*TF*nfl*(256/27._dp-2*pisq/9._dp)+4*CF*TF*nfl
      gammagq2(0)= CA**3*( - 97186/729._dp + 6109*pi**2/486._dp - 319*pi**4/270._dp 
     &  + 122/3._dp*zeta3 - 20*pi**2/9._dp*zeta3 - 16*zeta5 )
     & + CA**2*TF*nfl*( 30715/729._dp - 1198*pi**2/243._dp + 82*pi**4/135._dp + 712/27._dp*zeta3 ) 
     & + CA*CF*TF*nfl*( 2434/27._dp - 2*pi**2/3._dp - 8*pi**4/45._dp - 304/9._dp*zeta3 ) 
     &   - 2*CF**2*TF*nfl + CA*TF**2*nfl**2*( - 538/729._dp + 40*pi**2/81._dp - 224/27._dp*zeta3 ) 
     &   - 44/9._dp*CF*TF**2*nfl**2
      gammagq0(1)=-3*CF
      gammagq1(1)=CF**2*(-3/2._dp+2*pisq-24*zeta3)
     & +CF*CA*(-961/54._dp-11*pisq/6._dp+26*zeta3)
     & +CF*TF*nfl*(130/27._dp+2*pisq/3._dp)
      gammagq2(1)=CF**3*( -29/2._dp - 3*pi**2 - 8*pi**4/5._dp - 68*zeta3 + 16*pi**2/3._dp*zeta3 + 240*zeta5 )
     & + CF**2*CA*( - 151/4._dp + 205*pi**2/9._dp
     &    + 247*pi**4/135._dp - 844/3._dp*zeta3
     &    - 8*pi**2/3._dp*zeta3 - 120*zeta5 ) 
     & + CF*CA**2*( - 139345/2916._dp - 7163*pi**2/486._dp
     &    - 83*pi**4/90._dp + 3526/9._dp*zeta3
     &    - 44*pi**2/9._dp*zeta3 - 136*zeta5 ) 
     & + CF**2*TF*nfl*( 2953/27._dp - 26*pi**2/9._dp - 28*pi**4/27._dp + 512/9._dp*zeta3 )
     & + CF*CA*TF*nfl*( - 17318/729._dp + 2594*pi**2/243._dp + 22*pi**4/45._dp  - 1928/27._dp*zeta3 ) 
     & + CF*TF**2*nfl**2*( 9668/729._dp - 40*pi**2/27._dp - 32/27._dp*zeta3 )
! Numerical value taken from 2102.09725, Eq. (23)
! Note that we need an overall factor of (-1/2), c.f. comment in 2002.04617 before Eq. (8)
      gammagq3(0) = -0.5_dp*(
     &   + dAANA * ( - 2451.040712450_dp)
     &   +CA**4 * ( 1557.4287417889_dp)
     &   +nfl*dRANA * ( - 41.2080190194_dp)
     &   +nfl*CA**3 * ( - 1033.98729659_dp)
     &   +nfl*CA**2*CF * ( - 57.9377499658_dp)
     &   +nfl*CA*CF**2 * ( - 100.315097910_dp)
     &   +nfl*CF**3  * ( 46._dp)
     &   +nfl**2*dRRNA * ( 253.857645167_dp)
     &   +nfl**2*CA**2  * ( 70.7744401902_dp)
     &   +nfl**2*CA*CF  * ( 73.9372035966_dp)
     &   +nfl**2*CF**2 * ( - 21.9767440643_dp)
     &   +nfl**3*CA  * ( 0.405507202650_dp)
     &   +nfl**3*CF  * ( 1.26748971193_dp)
     &   )
! Numerical value taken from 2102.09725, Eq. (22)
! Note that we need an overall factor of (-1/2), c.f. comment in 2002.04617 before Eq. (8)
      gammagq3(1) = -0.5_dp*(
     &   + dRANA * ( - 2195.670535473_dp)
     &   + CA**3*CF * ( - 13.809312037_dp)
     &   + CA**2*CF**2 * ( 2438.569338812_dp)
     &   + CA*CF**3 * (- 1373.764650948_dp)
     &   +CF**4  * ( 392.899478384_dp)
     &   +nfl*dRRNA * ( - 425.019550390_dp)
     &   +nfl*CA**2*CF * ( - 274.147360589_dp)
     &   +nfl*CA*CF**2 * ( - 912.844845636_dp)
     &   +nfl*CF**3  * ( 151.933788877_dp)
     &   +nfl**2*CA*CF  * ( 109.081415293_dp)
     &   +nfl**2*CF**2 * ( - 12.5342425083_dp)
     &   +nfl**3*CF  * ( 4.88682798281_dp)
     &   )

!    * Relationship in 1410.1892 Eq. (I.5)
!    0 = 2*gammag-(gammat+gammaS+beta/as);
!    where
! 1410.1892 Eq. (I.2)
! gammag=gamma0(g)*ason4pi+gamma1(g)*ason4pi^2
! +gamma2(g)*ason4pi^3+gamma3(g)*ason4pi^4;

! 1410.1892 Eq. (8.33)
! beta=-2*as*(beta0*ason4pi+beta1*ason4pi^2+beta2*ason4pi^3+beta3*ason4pi^4);

! 1205.3806 Eq. (A1)
! id,gammat=-2*(beta1*ason4pi^2+2*beta2*ason4pi^3+3*beta3*ason4pi^4);
!
      gamma0S=2*gammagq0(0)+2*beta0
      gamma1S=2*gammagq1(0)+4*beta1
      gamma2S=2*gammagq2(0)+6*beta2
      gamma3S=2*gammagq3(0)+8*beta3

!      gamma0S=0._dp
!      gamma1S=(CA**2*(-160._dp/27._dp+11._dp*pisq/9._dp+4._dp*zeta3)
!     &    +CA*TF*nfl*(-208._dp/27._dp-4._dp*pisq/9._dp)-8._dp*CF*TF*nfl)
!      gamma2S=
!     &    CA**3*( 37045._dp/729._dp + 6109._dp*pi**2/243._dp
!     &   - 319._dp*pi**4/135 ._dp + ( 244._dp/3._dp - 40._dp*pi**2/9._dp)*zeta3
!     &   - 32._dp*zeta5 )
!     &  + CA**2*TF*nfl*( -167800._dp/729._dp
!     &   - 2396._dp*pi**2/243._dp + 164._dp*pi**4/135._dp 
!     &   + 1424._dp/27._dp*zeta3 ) 
!     &  + CA*CF*TF*nfl*( 1178._dp/27._dp
!     &   - 4._dp*pi**2/3._dp - 16._dp*pi**4/45._dp 
!     &   - 608._dp/9._dp*zeta3 )
!     &  + 8*CF**2*TF*nfl
!     &  + CA*TF**2*nfl**2*( 24520._dp/729._dp
!     &   + 80*pi**2/81._dp - 448._dp/27._dp*zeta3 ) 
!     &  + 176._dp/9._dp*CF*TF**2*nfl**2
!      write(6,*) gamma0S,2*gammagq0(0)+2*beta0
!      write(6,*) gamma1S,2*gammagq1(0)+4*beta1
!      write(6,*) gamma2S,2*gammagq2(0)+6*beta2
!      pause

      return
      end subroutine Anomfill
