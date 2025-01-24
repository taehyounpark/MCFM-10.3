!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine I2bar(x,lnmuonpt,R,I2)
      implicit none
!     using results from https://arxiv.org/abs/2207.07037
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      include 'transitionlabels.f'
      include 'distributions.f'
!     dmin=-1,dmax=2,delt=-1,plus=0,lpls=1,rglr=2
      real(dp)::I2(0:6,dmin:dmax),I2P(0:6,dmin:dmax),I2R(0:6,dmin:dmax)
!      I2, index one transition flavor
!      I2, index2,distribution type

      real(dp):: x,lnmuonpt,R,Li2,Li3,omx,opx,lnx,lnomx,lnopx,xonopx,
     & Li2x,Li2mx,Li2omx,Li3x,Li3mx,Li3omx,Li3xonopx,lnopxon2,lnR,R2,R4,R6,R8,R10
      real(dp),parameter::picube=pi*pisq,lnpi=1.1447298858494_dp,
     & xmin=0.9_dp,
     & Li4xhalf=0.5174790616738991_dp,ln2=0.6931471805599453_dp,
     & ln64=6*ln2,ln512=9*ln2,ln4096=12*ln2
      real(dp)::
     & BoundaryConditionGGCACA,BCGGCACAx,BCGGCACA1,
     & BoundaryConditionGGCATF,BCGGCATFx,BCGGCATF1,
     & BoundaryConditionGGCFTF,BCGGCFTFx,BCGGCFTF1,

     & BoundaryConditionQQCACF,BCQQCACFx,BCQQCACF1,
     & BoundaryConditionQQCFCF,BCQQCFCFx,BCQQCFCF1,
     & BoundaryConditionQQCFTF,BCQQCFTFx,BCQQCFTF1,
     & 
     & BoundaryConditionGQCFCA,BCGQCFCAx,BCGQCFCA1,
     & BoundaryConditionGQCFCF,BCGQCFCFx,BCGQCFCF1,
     & 
     & BoundaryConditionQGCATF,BCQGCATFx,BCQGCATF1,
     & BoundaryConditionQGCFTF,BCQGCFTFx,BCQGCFTF1,
     & BoundaryConditionQQBARS,BCQQBARSx,BCQQBARS1,
     & BoundaryConditionQQBARNS,BCQQBARNSx,BCQQBARNS1

      omx=one-x
      opx=one+x
      lnx=log(x); lnomx=log(omx); lnopx=log(opx);lnopxon2=log(opx)-ln2
      xonopx=x/opx
      Li2x=Li2(x)
      Li2mx=Li2(-x)
      Li2omx=Li2(omx)
      Li3x=Li3(x)
      Li3mx=Li3(-x)
      Li3omx=Li3(omx)
      Li3xonopx=Li3(xonopx)
      
!     initialize to zero
      I2P(:,:)=0
      I2R(:,:)=0

!     I2P
      include 'I2Pgg.f'
      include 'I2Pqq.f'
      include 'I2Pgq.f'
      include 'I2Pqg.f'
      include 'I2Pqpq.f'
      include 'I2Pqbq.f'
      include 'I2Pqbpq.f'
      
! I2R
!     initialize Boundary Conditions
      lnR=log(R); R2=R**2; R4=R2*R2; R6=R4*R2; R8=R6*R2; R10=R8*R2
      
      BCGGCACAx=BoundaryConditionGGCACA(x)
      BCGGCACA1=BoundaryConditionGGCACA(one)

      BCGGCATFx=BoundaryConditionGGCATF(x)
      BCGGCATF1=BoundaryConditionGGCATF(one)

      BCGGCFTFx=BoundaryConditionGGCFTF(x)
      BCGGCFTF1=BoundaryConditionGGCFTF(one)

      BCQQCACFx=BoundaryConditionQQCACF(x)
      BCQQCACF1=BoundaryConditionQQCACF(one)

      BCQQCFCFx=BoundaryConditionQQCFCF(x)
      BCQQCFCF1=BoundaryConditionQQCFCF(one)

      BCQQCFTFx=BoundaryConditionQQCFTF(x)
      BCQQCFTF1=BoundaryConditionQQCFTF(one)

      BCGQCFCAx=BoundaryConditionGQCFCA(x)
      BCGQCFCA1=BoundaryConditionGQCFCA(one)

      BCGQCFCFx=BoundaryConditionGQCFCF(x)
      BCGQCFCF1=BoundaryConditionGQCFCF(one)

      BCQGCATFx=BoundaryConditionQGCATF(x)
      BCQGCATF1=BoundaryConditionQGCATF(one)

      BCQGCFTFx=BoundaryConditionQGCFTF(x)
      BCQGCFTF1=BoundaryConditionQGCFTF(one)

      BCQQBARSx=BoundaryConditionQQBARS(x)
      BCQQBARS1=BoundaryConditionQQBARS(one)

      BCQQBARNSx=BoundaryConditionQQBARNS(x)
      BCQQBARNS1=BoundaryConditionQQBARNS(one)

      include 'I2Rgg.f'
      include 'I2Rqq.f'
      include 'I2Rgq.f'
      include 'I2Rqg.f'
      include 'I2Rqpq.f'
      include 'I2Rqbq.f'
      include 'I2Rqbpq.f'
     
      I2(:,:)=I2p(:,:)+I2R(:,:)
      return
      end
