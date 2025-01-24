!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjet_qlooptri_sr(p1,p2,p3,p4,p5,p6,p7,ampsr)
c--- This is the amplitude for
c---   0 -> q(p1) q~(p2) e-(p3) e+(p4) mu-(p5) mu+(p) g(p7)
c--- for all helicity combinations:
c---  amp(helq,hel3,hel6,helg)
c---
c--- where the Z current contains all singly-resonant contributions
c--- such as the diagram below (but also interchanging Z's)
c---

c        l5
c           |      l3
c           | Z   /              q1
c           |~~~~~              /
c           |     \         oooo
c           |      a4     /|    \
c           ^            / |     q~2
c           |     Z     /  |
c           |~~~~~~~~~~~   |
c           |           \  |
c           |            \ |
c           |             \|
c        l6                 ooooo g7

c--- It is adapted from Eq. (A.14) of arXiv:0710.1832 (CEZ)
c--- and is multiplied by (s127-Mz^2)/s127
c---
c--- Supplemented by Higgs contributions from arXiv:1409.1897, Sec. 3
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'scale.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'scalarselect.f'
      include 'yukawas.f'
      integer p1,p2,p3,p4,p5,p6,p7,hel12,j1,j2,hel34,j3,j4,hel56,j5,j6
      real(dp):: s127,mfsq
      complex(dp):: cF1,Cmt,B12mt,B1mt,Cmb,B12mb,B1mb,
     & F1mt,F1mb,MZsq,ampsr(2,2,2,2),v1,v2,FH,
     & amp_qlooptri_sr_mmmp,amp_qlooptri_sr_mmpp

      s127=s(p1,p2)+s(p1,p7)+s(p2,p7)

c--- relevant scalar integrals
      mfsq=mt**2
      Cmt=loopI3(s(p1,p2),zip,s127,mfsq,mfsq,mfsq,musq,0)
      B12mt=loopI2(s127,mfsq,mfsq,musq,0)
      B1mt=loopI2(s(p1,p2),mfsq,mfsq,musq,0)
      mfsq=mb**2
      Cmb=loopI3(s(p1,p2),zip,s127,mfsq,mfsq,mfsq,musq,0)
      B12mb=loopI2(s127,mfsq,mfsq,musq,0)
      B1mb=loopI2(s(p1,p2),mfsq,mfsq,musq,0)

c--- Eq. (A.7)
      F1mt=1d0/4d0/(s(p1,p7)+s(p2,p7))*(2d0+4d0*mt**2*Cmt
     & +(2d0+2d0*s(p1,p2)/(s(p1,p7)+s(p2,p7)))*(B12mt-B1mt))
      F1mb=1d0/4d0/(s(p1,p7)+s(p2,p7))*(2d0+4d0*mb**2*Cmb
     & +(2d0+2d0*s(p1,p2)/(s(p1,p7)+s(p2,p7)))*(B12mb-B1mb))
c Note that these expressions are equivalent, but we need the
c scalar integrals above for the Higgs contributions anyway
c      F1mt=F1anom(s(p1,p2),s127,mt**2,musq)
c      F1mb=F1anom(s(p1,p2),s127,mb**2,musq)

c--- combination of functions appearing in amplitude
      cF1=F1mt-F1mb

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

      ampsr(:,:,:,:)=czip
      ampsr(1,1,1,2)=+amp_qlooptri_sr_mmmp(p1,p2,p3,p4,p5,p6,p7,za,zb,zl1,zl2)
      ampsr(1,1,2,2)=+amp_qlooptri_sr_mmpp(p1,p2,p3,p4,p5,p6,p7,za,zb,zl1,zr2)
      ampsr(1,2,2,2)=-amp_qlooptri_sr_mmmp(p1,p2,p4,p3,p6,p5,p7,za,zb,zr1,zr2)
      ampsr(1,2,1,2)=-amp_qlooptri_sr_mmpp(p1,p2,p4,p3,p6,p5,p7,za,zb,zr1,zl2)

      ampsr(2,1,1,2)=+amp_qlooptri_sr_mmmp(p2,p1,p3,p4,p5,p6,p7,za,zb,zl1,zl2)
      ampsr(2,1,2,2)=+amp_qlooptri_sr_mmpp(p2,p1,p3,p4,p5,p6,p7,za,zb,zl1,zr2)
      ampsr(2,2,2,2)=-amp_qlooptri_sr_mmmp(p2,p1,p4,p3,p6,p5,p7,za,zb,zr1,zr2)
      ampsr(2,2,1,2)=-amp_qlooptri_sr_mmpp(p2,p1,p4,p3,p6,p5,p7,za,zb,zr1,zl2)

      ampsr(1,1,1,1)=-amp_qlooptri_sr_mmmp(p2,p1,p4,p3,p6,p5,p7,zb,za,zl1,zl2)
      ampsr(1,1,2,1)=-amp_qlooptri_sr_mmpp(p2,p1,p4,p3,p6,p5,p7,zb,za,zl1,zr2)
      ampsr(1,2,2,1)=+amp_qlooptri_sr_mmmp(p2,p1,p3,p4,p5,p6,p7,zb,za,zr1,zr2)
      ampsr(1,2,1,1)=+amp_qlooptri_sr_mmpp(p2,p1,p3,p4,p5,p6,p7,zb,za,zr1,zl2)

      ampsr(2,1,1,1)=-amp_qlooptri_sr_mmmp(p1,p2,p4,p3,p6,p5,p7,zb,za,zl1,zl2)
      ampsr(2,1,2,1)=-amp_qlooptri_sr_mmpp(p1,p2,p4,p3,p6,p5,p7,zb,za,zl1,zr2)
      ampsr(2,2,2,1)=+amp_qlooptri_sr_mmmp(p1,p2,p3,p4,p5,p6,p7,zb,za,zr1,zr2)
      ampsr(2,2,1,1)=+amp_qlooptri_sr_mmpp(p1,p2,p3,p4,p5,p6,p7,zb,za,zr1,zl2)


      ampsr(:,:,:,:)=ampsr(:,:,:,:)*cF1/sqrt(zxw*(cone-zxw))/(s127-MZsq)

c--- additional contribution from Higgs diagrams, c.f. 1409.1897, Sec. 3
c Note: we do not account for other finite quark-mass effects in qloop
      FH=2._dp*mt*mt_yuk*(two-(s127-s(p1,p2)-4._dp*mt**2)*Cmt
     &                   +two*s(p1,p2)/(s127-s(p1,p2))*(B12mt-B1mt))
     &  +2._dp*mb*mb_yuk*(two-(s127-s(p1,p2)-4._dp*mb**2)*Cmb
     &                   +two*s(p1,p2)/(s127-s(p1,p2))*(B12mb-B1mb))
      do hel12=1,2
      if (hel12 == 1) then
        j1=p1
        j2=p2
      else
        j1=p2
        j2=p1
      endif
      do hel34=1,2
      if (hel34 == 1) then
        j3=p3
        j4=p4
        v1=zl1
      else
        j3=p4
        j4=p3
        v1=zr1
      endif
      do hel56=1,2
      if (hel56 == 1) then
        j5=p5
        j6=p6
        v2=zl2
      else
        j5=p6
        j6=p5
        v2=zr2
      endif
      ampsr(hel12,hel34,hel56,1)=ampsr(hel12,hel34,hel56,1)
     & +FH/cmplx(s127-hmass**2,hmass*hwidth,dp)/s(j1,j2)
     &  *(zb(j1,j2)*za(j1,p7)**2/(s127-s(j1,j2)))
     &  *za(j3,j5)*zb(j6,j4)/(s(j3,j4)*s(j5,j6))
     & *v1*s(j3,j4)/(s(p3,p4)-MZsq)
     & *v2*s(j5,j6)/(s(p5,p6)-MZsq)
     & /(4._dp)/(zxw*(cone-zxw))

      ampsr(hel12,hel34,hel56,2)=ampsr(hel12,hel34,hel56,2)
     & +FH/cmplx(s127-hmass**2,hmass*hwidth,dp)/s(j1,j2)
     &  *(za(j1,j2)*zb(j2,p7)**2/(s127-s(j1,j2)))
     &  *za(j3,j5)*zb(j6,j4)/(s(j3,j4)*s(j5,j6))
     & *v1*s(j3,j4)/(s(p3,p4)-MZsq)
     & *v2*s(j5,j6)/(s(p5,p6)-MZsq)
     & /(4._dp)/(zxw*(cone-zxw))

      enddo
      enddo
      enddo

      return
      end




      function amp_qlooptri_sr_mmmp(p1,p2,p3,p4,p5,p6,p7,za,zb,v1,v2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'nf.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zcouple.f'
      complex(dp):: amp_qlooptri_sr_mmmp,v1,v2,zab2,PZ,MZsq
      real(dp):: s56,s3
      integer p1,p2,p3,p4,p5,p6,p7

c--- statement functions
      zab2(p1,p3,p4,p2)=za(p1,p3)*zb(p3,p2)+za(p1,p4)*zb(p4,p2)
      PZ(q1,q2,v1,v2,s56)=q1*q2+v1*v2*s56/(s56-MZsq)
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

c essential amplitudes
      amp_qlooptri_sr_mmmp = zb(p2,p7)*(
     & + PZ(q1,q2,v1,v2,s(p5,p6))*v1/s(p5,p6) * (
     &       za(p1,p3)*zb(p4,p6)*zab2(p5,p4,p6,p7)/s3(p4,p5,p6)
     &     + za(p3,p5)*zb(p4,p7)*zab2(p1,p3,p5,p6)/s3(p3,p5,p6) )
     & - PZ(q1,q2,v1,v2,s(p3,p4))*v2/s(p3,p4) * (
     &       za(p1,p5)*zb(p4,p6)*zab2(p3,p6,p4,p7)/s3(p3,p4,p6)
     &     + za(p3,p5)*zb(p6,p7)*zab2(p1,p5,p3,p4)/s3(p3,p4,p5) ) )

      return
      end


      function amp_qlooptri_sr_mmpp(p1,p2,p3,p4,p5,p6,p7,za,zb,v1,v2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'nf.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zcouple.f'
      complex(dp):: amp_qlooptri_sr_mmpp,v1,v2,zab2,PZ,MZsq
      real(dp):: s56,s3
      integer p1,p2,p3,p4,p5,p6,p7

c--- statement functions
      zab2(p1,p3,p4,p2)=za(p1,p3)*zb(p3,p2)+za(p1,p4)*zb(p4,p2)
      PZ(q1,q2,v1,v2,s56)=q1*q2+v1*v2*s56/(s56-MZsq)
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)

      MZsq=cmplx(zmass**2,-zmass*zwidth,kind=dp)

c essential amplitudes
      amp_qlooptri_sr_mmpp = zb(p2,p7)*(
     & + PZ(q1,q2,v1,v2,s(p5,p6))*v1/s(p5,p6) * (
     &    za(p1,p3)*zb(p4,p5)*zab2(p6,p4,p5,p7)/s3(p4,p5,p6)
     &  + za(p3,p6)*zb(p4,p7)*zab2(p1,p3,p6,p5)/s3(p3,p5,p6) )
     & + PZ(q1,q2,v1,v2,s(p3,p4))*v2/s(p3,p4) * (
     &    za(p1,p6)*zb(p4,p5)*zab2(p3,p5,p4,p7)/s3(p3,p4,p5)
     &  + za(p3,p6)*zb(p5,p7)*zab2(p1,p6,p3,p4)/s3(p3,p4,p6) ) )

      return
      end
