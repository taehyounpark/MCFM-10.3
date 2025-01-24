!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module vv_msq
      use qT_amp
      implicit none

c Converts from amplitudes in qT-scheme (provided by qtvvamp)
c to amplitudes in N-nettiness scheme using results in schemer.frm,
c c.f. Heinrich et al.,  arXiv:1710.06294, suitably extended to
c include all renormalization group logs

      contains

      function msq_tree(qqb,qtype,ftype1,ftype2,i1,i2,i5,i6,i7,i8,qqbAj,reset)
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'newpoint.f'
      integer:: i1,i2,i5,i6,i7,i8
      integer:: qtype,ftype1,ftype2
      real(dp):: s12,s56,s78,msq_tree
      complex(dp)::resqT00,amp0(2,2,2),qqbAj(0:2,4,10)
      logical::qqb,reset
      complex(dp),save::amp0_save(2,2,2,2,2)
      logical, save:: savedamps(2,2)
!$omp threadprivate(amp0_save)
!$omp threadprivate(savedamps)

      s12=s(i1,i2)
      s56=s(i5,i6)
      s78=s(i7,i8)
      call fillcoupfac(ftype1,ftype2,s12,s56,s78)

c This cache system is based on the fact that the routine is called twice,
c first with (i1,i2) = (2,1) [qqb=F) and again with (i1,i2)=(1,2) (qqb=T)
c and also for two values of qtype (=1,2)
      if ((i1 >2) .or. (i2 > 2)) then
        write(6,*)
     & 'Cacheing system in vv_msq will not work with i1,i2 = ',i1,i2
        stop
      endif

      if (newpoint .and. reset) then
        savedamps(:,:)=.false.
        reset=.false.
      endif

      if (savedamps(qtype,i1) .eqv. .false.) then
        call MXXX(0,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp0_save(qtype,i1,:,:,:))
        savedamps(qtype,i1)=.true.
      endif
      amp0(:,:,:)=amp0_save(qtype,i1,:,:,:)

      resqT00=qTsum(amp0,amp0)
      msq_tree=real(resqT00,kind=dp)

      end function msq_tree

      function msq_onelp(qqb,qtype,ftype1,ftype2,i1,i2,i5,i6,i7,i8,qqbAj,reset)
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'newpoint.f'
      integer:: i1,i2,i5,i6,i7,i8
      integer:: qtype,ftype1,ftype2
      real(dp):: s12,s56,s78,msq_onelp
      complex(dp):: clog,lnrat
      complex(dp):: DeltaI1NQ,resqT10,amp0(2,2,2),amp1(2,2,2),qqbAj(0:2,4,10)
      logical::qqb,reset
      complex(dp),save::amp0_save(2,2,2,2,2),amp1_save(2,2,2,2,2)
      logical, save:: savedamps(2,2)
!$omp threadprivate(amp0_save,amp1_save)
!$omp threadprivate(savedamps)

      s12=s(i1,i2)
      s56=s(i5,i6)
      s78=s(i7,i8)
      call fillcoupfac(ftype1,ftype2,s12,s56,s78)

! Previous implementation for musq > 0, now extended to also include musq < 0
!      rlog=log(musq/s12)
!      DeltaI1NQ =
!     & + CF * (
!     &    - 3._dp/2*rlog
!     &    - 1._dp/2*rlog**2
!     &    - impi*rlog
!     &    + 1._dp/12*pi**2
!     &    )

      clog=lnrat(musq,-s12)

      DeltaI1NQ =
     & + CF * (
     &    - 3._dp/2*clog
     &    - 1._dp/2*clog**2
     &    + 3._dp/2._dp*impi
     &    - 5._dp/12*pi**2
     &    )

c This cache system is based on the fact that the routine is called twice,
c first with (i1,i2) = (2,1) [qqb=F) and again with (i1,i2)=(1,2) (qqb=T)
c and also for two values of qtype (=1,2)
      if ((i1 >2) .or. (i2 > 2)) then
        write(6,*)
     & 'Cacheing system in vv_msq will not work with i1,i2 = ',i1,i2
        stop
      endif

      if (newpoint .and. reset) then
        savedamps(:,:)=.false.
        reset=.false.
      endif

      if (savedamps(qtype,i1) .eqv. .false.) then
        call MXXX(0,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp0_save(qtype,i1,:,:,:))
        call MXXX(1,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp1_save(qtype,i1,:,:,:))
        savedamps(qtype,i1)=.true.
      endif
      amp0(:,:,:)=amp0_save(qtype,i1,:,:,:)
      amp1(:,:,:)=amp1_save(qtype,i1,:,:,:)

c Scheme conversion
      amp1=amp1+DeltaI1NQ*amp0

      resqT10=qTsum(amp1,amp0)
      msq_onelp=ason2pi*two*real(resqT10,kind=dp)

c to make it MCFM equivalent
c      msq_onelp=msq_onelp-ason2pi*CF*(1d0-zeta2)*qT00(qtype,i1,i2,i5,i6,i7,i8)

      end function msq_onelp

      function msq_twolp(qqb,qtype,ftype1,ftype2,i1,i2,i5,i6,i7,i8,qqbAj,reset)
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'zeta.f'
      include 'nflav.f'
      include 'newpoint.f'
      integer:: i1,i2,i5,i6,i7,i8
      integer:: qtype,ftype1,ftype2
      real(dp):: s12,s56,s78,be0,msq_twolp
      complex(dp):: clog,lnrat
      complex(dp):: DeltaI1NQ,DeltaI2NQ_0,DeltaI2NQ_1,
     & resqT11,resqT20,amp0(2,2,2),amp1(2,2,2),amp2(2,2,2),qqbAj(0:2,4,10)
      logical::qqb,reset
      complex(dp),save::amp0_save(2,2,2,2,2),amp1_save(2,2,2,2,2),amp2_save(2,2,2,2,2)
      logical, save:: savedamps(2,2)
!$omp threadprivate(amp0_save,amp1_save,amp2_save)
!$omp threadprivate(savedamps)

      s12=s(i1,i2)
      s56=s(i5,i6)
      s78=s(i7,i8)
      call fillcoupfac(ftype1,ftype2,s12,s56,s78)

! Previous implementation for musq > 0, now extended to also include musq < 0
!      rlog=log(musq/s12)
!      DeltaI1NQ =
!     & + CF * (
!     &    - 3._dp/2*rlog
!     &    - 1._dp/2*rlog**2
!     &    - impi*rlog
!     &    + 1._dp/12*pi**2
!     &    )

      clog=lnrat(musq,-s12)

      DeltaI1NQ =
     & + CF * (
     &    - 3._dp/2*clog
     &    - 1._dp/2*clog**2
     &    + 3._dp/2._dp*impi
     &    - 5._dp/12*pi**2
     &    )

      DeltaI2NQ_1 = DeltaI1NQ

!      DeltaI2NQ_0 =
!     & + CF*nflav * (
!     &    + 41._dp/81
!     &    + 65._dp/108*rlog
!     &    + 19._dp/36*rlog**2
!     &    + 1._dp/18*rlog**3
!     &    + 5._dp/9*pi*rlog*im
!     &    + 1._dp/6*pi*rlog**2*im
!     &    - 5._dp/72*pi**2
!     &    + 1._dp/18*pi**2*rlog
!     &    - 1._dp/36*pi**3*im
!     &    - 7._dp/18*zeta3
!     &    )
!
!      DeltaI2NQ_0 = DeltaI2NQ_0
!     & + CF*CA * (
!     &    - 607._dp/162
!     &    - 961._dp/216*rlog
!     &    - 233._dp/72*rlog**2
!     &    - 11._dp/36*rlog**3
!     &    - 67._dp/18*pi*rlog*im
!     &    - 11._dp/12*pi*rlog**2*im
!     &    + 67._dp/144*pi**2
!     &    - 11._dp/36*pi**2*rlog
!     &    + 1._dp/12*pi**2*rlog**2
!     &    + 11._dp/72*pi**3*im
!     &    + 1._dp/6*pi**3*rlog*im
!     &    - 1._dp/72*pi**4
!     &    + 77._dp/36*zeta3
!     &    + 13._dp/2*zeta3*rlog
!     &    )
!
!      DeltaI2NQ_0 = DeltaI2NQ_0
!     & + CF**2 * (
!     &    - 3._dp/8*rlog
!     &    + 9._dp/8*rlog**2
!     &    + 3._dp/4*rlog**3
!     &    + 1._dp/8*rlog**4
!     &    + 3._dp/2*pi*rlog**2*im
!     &    + 1._dp/2*pi*rlog**3*im
!     &    + 3._dp/8*pi**2*rlog
!     &    - 13._dp/24*pi**2*rlog**2
!     &    - 1._dp/12*pi**3*rlog*im
!     &    + 1._dp/288*pi**4
!     &    - 6*zeta3*rlog
!     &    )

      DeltaI2NQ_0 =
     & + CF*nflav * (
     &    + 41._dp/81._dp
     &    + 65._dp/108._dp*clog
     &    + 19._dp/36._dp*clog**2
     &    + 1._dp/18._dp*clog**3
     &    - 65._dp/108._dp*pi*im
     &    - 1._dp/2._dp*pi*clog*im
     &    - 1._dp/24._dp*pi**2
     &    + 2._dp/9._dp*pi**2*clog
     &    - 7._dp/36._dp*pi**3*im
     &    - 7._dp/18._dp*zeta3
     &    )

     & + CF*CA * (
     &    - 607._dp/162._dp
     &    - 961._dp/216._dp*clog
     &    - 233._dp/72._dp*clog**2
     &    - 11._dp/36._dp*clog**3
     &    + 961._dp/216._dp*pi*im
     &    + 11._dp/4._dp*pi*clog*im
     &    - 1._dp/48._dp*pi**2
     &    - 11._dp/9._dp*pi**2*clog
     &    + 1._dp/12._dp*pi**2*clog**2
     &    + 77._dp/72._dp*pi**3*im
     &    + 5._dp/72._dp*pi**4
     &    + 77._dp/36._dp*zeta3
     &    + 13._dp/2._dp*zeta3*clog
     &    - 13._dp/2._dp*zeta3*pi*im
     &    )

     & + CF**2 * (
     &    - 3._dp/8._dp*clog
     &    + 9._dp/8._dp*clog**2
     &    + 3._dp/4._dp*clog**3
     &    + 1._dp/8._dp*clog**4
     &    + 3._dp/8._dp*pi*im
     &    - 9._dp/4._dp*pi*clog*im
     &    - 3._dp/4._dp*pi*clog**2*im
     &    - 9._dp/8._dp*pi**2
     &    + 9._dp/8._dp*pi**2*clog
     &    + 5._dp/24*pi**2*clog**2
     &    - 9._dp/8._dp*pi**3*im
     &    + 25._dp/288._dp*pi**4
     &    - 6*zeta3*clog
     &    + 6*zeta3*pi*im
     &    )

c This cache system is based on the fact that the routine is called twice,
c first with (i1,i2) = (2,1) [qqb=F) and again with (i1,i2)=(1,2) (qqb=T)
c and also for two values of qtype (=1,2)
      if ((i1 >2) .or. (i2 > 2)) then
        write(6,*)
     & 'Cacheing system in vv_msq will not work with i1,i2 = ',i1,i2
        stop
      endif

      if (newpoint .and. reset) then
        savedamps(:,:)=.false.
        reset=.false.
      endif

      if (savedamps(qtype,i1) .eqv. .false.) then
        call MXXX(0,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp0_save(qtype,i1,:,:,:))
        call MXXX(1,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp1_save(qtype,i1,:,:,:))
        call MXXX(2,qqb,qType,i1,i2,i5,i6,i7,i8,qqbAj,amp2_save(qtype,i1,:,:,:))
        savedamps(qtype,i1)=.true.
      endif
      amp0(:,:,:)=amp0_save(qtype,i1,:,:,:)
      amp1(:,:,:)=amp1_save(qtype,i1,:,:,:)
      amp2(:,:,:)=amp2_save(qtype,i1,:,:,:)

c First add missing be0 term (since qtvvamp results quoted at musq = s)
      be0=(11._dp*xn-2._dp*nflav)/6._dp
      amp2=amp2+be0*real(clog,kind=dp)*amp1

c Scheme conversion
      amp2=amp2+DeltaI2NQ_1*amp1+DeltaI2NQ_0*amp0
      amp1=amp1+DeltaI1NQ*amp0

      resqT11=qTsum(amp1,amp1)
      resqT20=qTsum(amp2,amp0)
      msq_twolp=real(resqT11+two*resqT20,kind=dp)*ason2pi**2

      end function msq_twolp

      end module vv_msq

