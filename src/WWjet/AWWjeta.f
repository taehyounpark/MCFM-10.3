!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjeta(p,j1,j2,j3,j4,j5,j6,j7,Alop,Alom,A7vp,A7vm)
c--- compute amplitude A^a_7 virtual, including both
c--- leading and subleading pieces, for both positive
c--- and negative gluon helicities
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
      integer j1,j2,j3,j4,j5,j6,j7
      real(dp):: p(mxpart,4)
      complex(dp):: Alom,Alop,A7vm,A7vp,A7a_lc,A7a_slc,zlc,zslc

      call spinoru(7,p,za,zb)

      call WWjetcomputescalars(j1,j2,j3,j4,j5,j6,j7,scints)

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3d0,0d0)           ! Nc
      zslc=cplx2(-1d0/3d0,0d0)     ! -1/Nc

c flag for printing out KCheck comparison
      verbose=.false.

c--- positive helicity gluon
      call AWWjeta_rat(j1,j2,j3,j4,j5,j6,j7,za,zb,coeff,Alop)
      call AWWjeta_lc(.false.,.false.,p,j1,j2,j3,j4,j5,j6,j7,za,zb,.false.,scints,coeff,Alop,A7a_lc)
      call AWWjeta_slc(.false.,.false.,p,j1,j2,j3,j4,j5,j6,j7,za,zb,.false.,scints,coeff,Alop,A7a_slc)

c      write(6,*) 'plus A7a_lc ',A7a_lc
c      write(6,*) 'plus A7a_slc',A7a_slc

      A7vp=zlc*A7a_lc+zslc*A7a_slc

c--- negative helicity gluon
      call AWWjeta_rat(j2,j1,j5,j6,j3,j4,j7,zb,za,coeff,Alom)
      call AWWjeta_lc(.true.,.true.,p,j2,j1,j5,j6,j3,j4,j7,zb,za,.true.,scints,coeff,Alom,A7a_lc)
      call AWWjeta_slc(.true.,.true.,p,j2,j1,j5,j6,j3,j4,j7,zb,za,.true.,scints,coeff,Alom,A7a_slc)

c      write(6,*) 'minus A7a_lc ',A7a_lc
c      write(6,*) 'minus A7a_slc',A7a_slc

      A7vm=zlc*A7a_lc+zslc*A7a_slc

c--- strong coupling renormalization in dred scheme
      A7vp=A7vp-zlc*(b0/xn*epinv-1d0/6d0)*Alop
      A7vm=A7vm-zlc*(b0/xn*epinv-1d0/6d0)*Alom
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7vp=A7vp-(Cf+xn/6d0)*Alop
        A7vm=A7vm-(Cf+xn/6d0)*Alom
      endif

      return
      end

