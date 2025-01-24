!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine fillcoupfac(ftype1,ftype2,s12,s56,s78)
c Fills array coupfac(6,2,2,2,2) containing combination of couplings that
c appear in all contributions:
c coupfac(ABCF, qtype, h12, h56, h78) with ABCF = 1(A), 2(B), 3(C), 4(F)

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'nf.f'
      include 'nflav.f'
      include 'zcouple_cms.f'
      include 'masses.f'
      include 'coupfac_vv.f'
      include 'ABCF.f'
      include 'iDen.f'
      include 'first.f'
      integer:: qtype,ftype1,ftype2
      integer:: nq,h12,h56,h78
      real(dp):: s12,s56,s78,activeflavour
      complex(dp):: prop12,prop56,prop78,sinw,cosw,sin2w,cotw,
     & cz(2,2),cgamz(2,2)
      complex(dp),save:: VQ(2,2),vf(2,2) !VQ(qtype,hel)
      real(dp),parameter::eQ(2)=[-1._dp/3._dp,2._dp/3._dp]
      real(dp),parameter::I3(2)=[-0.5_dp,0.5_dp]
      real(dp),parameter::I3f(2)=[-0.5_dp,0.5_dp]
      real(dp),parameter::ef(2)=[-1._dp,0._dp]

      sinw=sqrt(zxw)
      cosw=sqrt(1._dp-zxw)
      sin2w=2*sinw*cosw
      cotw=cosw/sinw
      if (first) then  ! fill Z-coupling matrix
         do h12=1,2; do qtype=1,2
         VQ(qtype,h12)=((2-h12)*2*I3(qtype)-2*zxw*eQ(qtype))/(sin2w)
         vf(qtype,h12)=((2-h12)*2*I3f(qtype)-2*zxw*ef(qtype))/(sin2w)
         enddo; enddo
         first=.false.
      endif


      coupfac(:,:,:,:,:)=czip     ! Initialize to zero


      if  (kcase==kWWqqbr) then                ! WW
        iDen(i12)=cone/cmplx(s12-zmass**2,zwidth*zmass,dp)
        iDen(i56)=cone/cmplx(s56-wmass**2,wwidth*wmass,dp)
        iDen(i78)=cone/cmplx(s78-wmass**2,wwidth*wmass,dp)
        activeflavour=real(nflav,kind=dp)
      do qtype=1,2
      coupfac(A,qtype,1,1,1)=iDen(i56)*iDen(i78)
      coupfac(B,qtype,1,1,1)=coupfac(A,qtype,1,1,1)
      do h12=1,2
      coupfac(C,qtype,h12,1,1)=0.5_dp*activeflavour*iDen(i56)*iDen(i78)
      coupfac(F,qtype,h12,1,1)=(-sin2w*VQ(qtype,h12)*iDen(i12)
     & -2*zxw*eQ(qtype)/s12)*iDen(i56)*iDen(i78)
      enddo
      enddo

c     singly resonant diagrams
c     Note that coupfac(5,....) and coupfac(6,....) are normalized differently,
c     so define prop12 etc.
      prop12=s12*iDen(i12)
      prop56=s56*iDen(i56)
      prop78=s78*iDen(i78)

      do qtype=1,2
      do h12=1,2
      cz(qtype,h12)=2._dp*zxw*zln*VQ(qtype,h12)*prop12
      cgamz(qtype,h12)=2._dp*zxw*(zle*VQ(qtype,h12)*prop12-eQ(qtype))
      coupfac(5,qtype,h12,1,1)=cgamz(qtype,h12)*prop56
      coupfac(6,qtype,h12,1,1)=cgamz(qtype,h12)*prop78
      coupfac(7,qtype,h12,1,1)=cz(qtype,h12)*prop56
      coupfac(8,qtype,h12,1,1)=cz(qtype,h12)*prop78
      enddo
      enddo

c      pause
c Upgrade with virtual photon contributions
      elseif (kcase == kZZlept) then
        iDen(i56)=cone/cmplx(s56-zmass**2,zwidth*zmass,dp)
        iDen(i78)=cone/cmplx(s78-zmass**2,zwidth*zmass,dp)
        iDen(i12)=cone/cmplx(s12-zmass**2,zwidth*zmass,dp)
c coupfac(ABCF, qtype, h12, h56, 78) with ABCF = 1(A), 2(B), 3(C), 4(F)
      do h56=1,2; do h78=1,2
      coupfac(C,1,1,h56,h78)=czip
      do qtype=1,2 ! for C-type sum over qtypes in loop
      nq=4-qtype   ! number of active flavours for C-type
      do h12=1,2   ! for C-type sum over handedness in loop

      coupfac(A,qtype,h12,h56,h78)=
     &   (VQ(qtype,h12)*vf(ftype1,h56)*iDen(i56)+eQ(qtype)*ef(ftype1)/s56)
     &  *(VQ(qtype,h12)*vf(ftype2,h78)*iDen(i78)+eQ(qtype)*ef(ftype2)/s78)

      coupfac(C,1,1,h56,h78)=coupfac(C,1,1,h56,h78)+0.5*nq*(
     & (VQ(qtype,h12)*vf(ftype1,h56)*iden(i56)+eQ(qtype)*ef(ftype1)/s56)
     &*(VQ(qtype,h12)*vf(ftype2,h78)*iden(i78)+eQ(qtype)*ef(ftype2)/s78))
      enddo
      enddo       ! not really sure yet of the normalization for C-type
      coupfac(C,2,2,h56,h78)=coupfac(C,1,1,h56,h78)
      coupfac(C,1,2,h56,h78)=coupfac(C,1,1,h56,h78)
      coupfac(C,2,1,h56,h78)=coupfac(C,1,1,h56,h78)
      enddo; enddo
      coupfac(B,:,:,:,:)=coupfac(A,:,:,:,:)


c     singly resonant diagrams
c     Note that coupfac(5,....) and coupfac(6,....) are normalized differently
      prop12=s12*iDen(i12)
      prop56=s56*iDen(i56)
      prop78=s78*iDen(i78)
      do qtype=1,2; do h12=1,2; do h56=1,2; do h78=1,2
       coupfac(5,qtype,h12,h56,h78)=
     &   (prop12*VQ(qtype,h12)*vf(ftype1,h56)+eQ(qtype)*ef(ftype1))
     &  *(prop78*vf(ftype1,h56)*vf(ftype2,h78)+ef(ftype1)*ef(ftype2))
      coupfac(6,qtype,h12,h56,h78)=
     &   (prop12*VQ(qtype,h12)*vf(ftype2,h78)+eQ(qtype)*ef(ftype2))
     &  *(prop56*vf(ftype2,h78)*vf(ftype1,h56)+ef(ftype2)*ef(ftype1))
      enddo; enddo; enddo; enddo
      coupfac(B,:,:,:,:)=coupfac(A,:,:,:,:)

      elseif (kcase == kWZbbar) then

       iDen(i56)=cone/cmplx(s56-wmass**2,wwidth*wmass,dp)
       iDen(i78)=cone/cmplx(s78-zmass**2,zwidth*zmass,dp)
       iDen(i12)=cone/cmplx(s12-wmass**2,wwidth*wmass,dp)
       do qtype=1,2;do h78=1,2
       coupfac(A,qtype,1,1,h78)=
     & (VQ(qtype,1)*vf(ftype2,h78)*iDen(i78)+eQ(qtype)*ef(ftype2)/s78)
     & *iDen(i56)
       coupfac(F,qtype,1,1,h78)=
     & (cotw*vf(ftype2,h78)*iDen(i78)+ef(ftype2)/s78)
     & *iDen(i56)*iDen(i12)
       enddo
       enddo
       coupfac(B,:,:,:,:)=coupfac(A,:,:,:,:)
      endif

      return
      end
