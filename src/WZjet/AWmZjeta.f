!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWmZjeta(p,j1,j2,j3,j4,j5,j6,j7,A7dalo,A7ualo,A7dav,A7uav)
c--- compute amplitude A^a_7 virtual, including both
c--- leading and subleading pieces, for both positive
c--- and negative gluon helicities, and for both lepton helicities for Z decay
c--- labelling is A7xxx(gluehel, Zlephel)
c--- Note: includes both the usual strong coupling renormalization
c--- and also the finite renormalization in the dred scheme
c--- leading order amplitudes are also returned
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'zprods_com.f'
      include 'WWjetlabels.f'
      include 'scheme.f'
      include 'b0.f'
      include 'verbose.f'
      integer j1,j2,j3,j4,j5,j6,k5,k6,j7,iz
      real(dp):: p(mxpart,4)
      complex(dp):: A7dalo(2,2),A7dav(2,2),A7ualo(2,2),A7uav(2,2),
     & A7a,A7a_lc,A7a_slc,zlc,zslc

      call spinoru(7,p,za,zb)

      call WWjetcomputescalars(j1,j2,j3,j4,j5,j6,j7,scints)

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3d0,0d0)           ! Nc
      zslc=cplx2(-1d0/3d0,0d0)     ! -1/Nc

c flag for printing out KCheck comparison
      verbose=.false.

      do iz=1,2

      if (iz == 1) then
        k5=j5
        k6=j6
      else
        k5=j6
        k6=j5
      endif

c--- positive helicity gluon
c--- down-type amplitudes
      call AWWjeta_rat(j1,j2,j3,j4,k5,k6,j7,za,zb,coeff,A7a)
      call AWmZjeta_lc(.false.,.false.,p,j1,j2,j3,j4,k5,k6,j7,za,zb,.false.,scints,coeff,A7a,A7a_lc)
      call AWmZjeta_slc(.false.,.false.,p,j1,j2,j3,j4,k5,k6,j7,za,zb,.false.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'plus A7a_lc ',A7a_lc
c      write(6,*) 'plus A7a_slc',A7a_slc
      A7dalo(2,iz)=A7a
      A7dav(2,iz)=zlc*A7a_lc+zslc*A7a_slc

c--- up-type amplitudes
      call AWWjeta_rat(j1,j2,k6,k5,j4,j3,j7,za,zb,coeff,A7a)
      call AWmZjeta_lc(.false.,.true.,p,j1,j2,k6,k5,j4,j3,j7,za,zb,.false.,scints,coeff,A7a,A7a_lc)
      call AWmZjeta_slc(.false.,.true.,p,j1,j2,k6,k5,j4,j3,j7,za,zb,.false.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'plus A7a_lc ',A7a_lc
c      write(6,*) 'plus A7a_slc',A7a_slc
      A7ualo(2,iz)=A7a
      A7uav(2,iz)=zlc*A7a_lc+zslc*A7a_slc

c--- negative helicity gluon
c--- down-type amplitudes
      call AWWjeta_rat(j2,j1,k5,k6,j3,j4,j7,zb,za,coeff,A7a)
      call AWmZjeta_lc(.true.,.true.,p,j2,j1,k5,k6,j3,j4,j7,zb,za,.true.,scints,coeff,A7a,A7a_lc)
      call AWmZjeta_slc(.true.,.true.,p,j2,j1,k5,k6,j3,j4,j7,zb,za,.true.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'minus A7a_lc ',A7a_lc
c      write(6,*) 'minus A7a_slc',A7a_slc
      A7dalo(1,iz)=A7a
      A7dav(1,iz)=zlc*A7a_lc+zslc*A7a_slc

c--- up-type amplitudes
      call AWWjeta_rat(j2,j1,j4,j3,k6,k5,j7,zb,za,coeff,A7a)
      call AWmZjeta_lc(.true.,.false.,p,j2,j1,j4,j3,k6,k5,j7,zb,za,.true.,scints,coeff,A7a,A7a_lc)
      call AWmZjeta_slc(.true.,.false.,p,j2,j1,j4,j3,k6,k5,j7,zb,za,.true.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'minus A7a_lc ',A7a_lc
c      write(6,*) 'minus A7a_slc',A7a_slc
      A7ualo(1,iz)=A7a
      A7uav(1,iz)=zlc*A7a_lc+zslc*A7a_slc

      enddo

c--- strong coupling renormalization in dred scheme
      A7dav(:,:)=A7dav(:,:)-zlc*(b0/xn*epinv-1d0/6d0)*A7dalo(:,:)
      A7uav(:,:)=A7uav(:,:)-zlc*(b0/xn*epinv-1d0/6d0)*A7ualo(:,:)
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7dav(:,:)=A7dav(:,:)-(Cf+xn/6d0)*A7dalo(:,:)
        A7uav(:,:)=A7uav(:,:)-(Cf+xn/6d0)*A7ualo(:,:)
      endif

      return
      end

