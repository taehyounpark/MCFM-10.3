!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function fBNR(R,CB)
!     f(R) from Eqn 17/A1 of 1307.0025
!     CB=CA for Higgs, CB=CF for Drell-Yan
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::fBNR,ccA(-1:10),ccf(-1:10),R,lnR,ln2,sumCFCA,sumCFTF,CB
      integer::j

      ln2=log(2d0)
      ccA(:)=0d0
      ccf(:)=0d0

!1307.0025 Eq.A16/Eq.43
      ccA(-1)=131d0/72d0-pisq/6-11d0/6d0*ln2
      ccA(0)=-805d0/216d0+11d0*pisq/72+35d0/18d0*ln2+11d0/6d0*ln2**2
     & +zeta3/2
!1307.0025 Eq.A21
      ccA(2)=1429d0/172800d0+pisq/48d0+13d0/180d0*ln2
      ccA(4)=-9383279d0/406425600d0-pisq/3456d0+587d0/120960d0*ln2
      ccA(6)=74801417d0/97542144000d0-23d0/67200d0*ln2
      ccA(8)=-50937246539d0/2266099089408d3-pisq/24883200d0
     & +28529d0/1916006400d0*ln2
      ccA(10)=348989849431d0/243708656615424d3-3509d0/3962649600d0*ln2

!1307.0025 Eq.A16/Eq.43
      ccf(-1)=-23d0/36d0+2d0/3d0*ln2
      ccf(0)=157d0/108d0-pisq/18d0-8d0/9d0*ln2-2d0/3d0*ln2**2
!1307.0025 Eq.A20
      ccf(2)=3071d0/86400d0-7d0/360d0*ln2
      ccf(4)=-168401d0/101606400d0+53d0/30240d0*ln2
      ccf(6)=7001023d0/48771072000d0-11d0/100800d0*ln2
      ccf(8)=-5664846191d0/566524772352000d0+4001d0/479001600d0*ln2
      ccf(10)=68089272001d0/83774850711552000d0
     & -13817d0/21794572800d0*ln2

      sumCFCA=ccA(-1)*log(R)
      sumCFTF=ccf(-1)*log(R)
      do j=0,10,2
         sumCFCA=sumCFCA+ccA(j)*R**j
         sumCFTF=sumCFTF+ccf(j)*R**j
      enddo

!1307.0025 Eq.A1
      fBNR=CA*sumCFCA+TF*nfl*sumCFTF+CB*(-pisq/12*R**2+R**4/16)

      return
      end

      function fBNR4(R,CB)
!     f(R) from Eqn 17/A1 of 1307.0025
!     CB=CA for Higgs, CB=CF for Drell-Yan
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp),intent(in)::R,CB
      real(dp)::fBNR4,ccA(-1:10),ccf(-1:10),sumCFCA,sumCFTF,
     & lnR,R2,R4,R6,R8,R10
      real(dp),parameter::ln2=0.6931471805599453_dp
      integer::j

      ccA(:)=0d0
      ccf(:)=0d0
      lnR=log(R); R2=R**2; R4=R2*R2; R6=R4*R2; R8=R6*R2; R10=R8*R2
!1307.0025 Eq.A16/Eq.43
      ccA(-1)=131d0/72d0-pisq/6-11d0/6d0*ln2
      ccA(0)=-805d0/216d0+11d0*pisq/72+35d0/18d0*ln2+11d0/6d0*ln2**2
     & +zeta3/2
!1307.0025 Eq.A21
      ccA(2)=1429d0/172800d0+pisq/48d0+13d0/180d0*ln2
      ccA(4)=-9383279d0/406425600d0-pisq/3456d0+587d0/120960d0*ln2

!1307.0025 Eq.A16/Eq.43
      ccf(-1)=-23d0/36d0+2d0/3d0*ln2
      ccf(0)=157d0/108d0-pisq/18d0-8d0/9d0*ln2-2d0/3d0*ln2**2
!1307.0025 Eq.A20
      ccf(2)=3071d0/86400d0-7d0/360d0*ln2
      ccf(4)=-168401d0/101606400d0+53d0/30240d0*ln2

      sumCFCA=ccA(-1)*lnR+ccA(0)
      sumCFTF=ccf(-1)*lnR+ccf(0)
      sumCFCA=sumCFCA+ccA(2)*R2+ccA(4)*R4
      sumCFTF=sumCFTF+ccf(2)*R2+ccf(4)*R4

!1307.0025 Eq.A1
      fBNR4=CA*sumCFCA+TF*nfl*sumCFTF+CB*(-pisq/12._dp*R2+R4/16._dp)

      return
      end

      function fBNR2(R,CB)
!     f(R) from Eqn 17/A1 of 1307.0025
!     CB=CA for Higgs, CB=CF for Drell-Yan
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp),intent(in)::R,CB
      real(dp)::fBNR2,ccA(-1:10),ccf(-1:10),sumCFCA,sumCFTF,
     & lnR,R2,R4,R6,R8,R10
      real(dp),parameter::ln2=0.6931471805599453_dp
      integer::j

      ccA(:)=0d0
      ccf(:)=0d0
      lnR=log(R); R2=R**2; R4=R2*R2; R6=R4*R2; R8=R6*R2; R10=R8*R2
!1307.0025 Eq.A16/Eq.43
      ccA(-1)=131d0/72d0-pisq/6-11d0/6d0*ln2
      ccA(0)=-805d0/216d0+11d0*pisq/72+35d0/18d0*ln2+11d0/6d0*ln2**2
     & +zeta3/2
!1307.0025 Eq.A21
      ccA(2)=1429d0/172800d0+pisq/48d0+13d0/180d0*ln2

!1307.0025 Eq.A16/Eq.43
      ccf(-1)=-23d0/36d0+2d0/3d0*ln2
      ccf(0)=157d0/108d0-pisq/18d0-8d0/9d0*ln2-2d0/3d0*ln2**2
!1307.0025 Eq.A20
      ccf(2)=3071d0/86400d0-7d0/360d0*ln2

      sumCFCA=ccA(-1)*lnR+ccA(0)
      sumCFTF=ccf(-1)*lnR+ccf(0)
      sumCFCA=sumCFCA+ccA(2)*R2
      sumCFTF=sumCFTF+ccf(2)*R2

!1307.0025 Eq.A1
      fBNR2=CA*sumCFCA+TF*nfl*sumCFTF+CB*(-pisq/12._dp*R2)

      return
      end


      function fBNR0(R,CB)
!     f(R) from Eqn 17/A1 of 1307.0025
!     CB=CA for Higgs, CB=CF for Drell-Yan
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp),intent(in)::R,CB
      real(dp)::fBNR0,ccA(-1:10),ccf(-1:10),sumCFCA,sumCFTF,
     & lnR,R2,R4,R6,R8,R10
      real(dp),parameter::ln2=0.6931471805599453_dp
      integer::j

      ccA(:)=0d0
      ccf(:)=0d0
      lnR=log(R); R2=R**2; R4=R2*R2; R6=R4*R2; R8=R6*R2; R10=R8*R2
!1307.0025 Eq.A16/Eq.43
      ccA(-1)=131d0/72d0-pisq/6-11d0/6d0*ln2
      ccA(0)=-805d0/216d0+11d0*pisq/72+35d0/18d0*ln2+11d0/6d0*ln2**2
     & +zeta3/2

!1307.0025 Eq.A16/Eq.43
      ccf(-1)=-23d0/36d0+2d0/3d0*ln2
      ccf(0)=157d0/108d0-pisq/18d0-8d0/9d0*ln2-2d0/3d0*ln2**2

      sumCFCA=ccA(-1)*lnR+ccA(0)
      sumCFTF=ccf(-1)*lnR+ccf(0)

!1307.0025 Eq.A1
      fBNR0=CA*sumCFCA+TF*nfl*sumCFTF
      return
      end
