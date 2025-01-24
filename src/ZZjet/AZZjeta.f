!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjeta(p,j1,j2,j3,j4,j5,j6,j7,A7alo,A7asr,A7av,A7avsr)
c--- compute amplitude A^a_7 and A^q_7 virtual, including both
c--- leading and subleading pieces, for both positive
c--- and negative gluon helicities, and for all quark and lepton helicities
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
      include 'srdiags.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,k1,k2,k3,k4,k5,k6,i12,i34,i56,Qid
      real(dp):: p(mxpart,4)
      complex(dp):: A7alo(2,2,2,2),A7av(2,2,2,2),
     & A7asr(2,2,2,2,2),A7avsr(2,2,2,2,2),A7a,A7a_lc,A7a_slc,zlc,zslc
      logical:: swap12
      character, parameter:: ah(2) = (/ 'm', 'p' /)
      character(len=4):: help,helm

      call spinoru(7,p,za,zb)

      call WWjetcomputescalars(j1,j2,j3,j4,j5,j6,j7,scints)

c--- color factors for leading and subleading amplitudes
      zlc=cplx2(3d0,0d0)           ! Nc
      zslc=cplx2(-1d0/3d0,0d0)     ! -1/Nc

c flag for printing out KCheck comparison
      verbose=.false.

      do i12=1,2
      do i34=1,2
      do i56=1,2

      if (i12 == 1) then
        k1=j1
        k2=j2
        swap12=.false.
      else
        k1=j2
        k2=j1
        swap12=.true.
      endif
      if (i34 == 1) then
        k3=j3
        k4=j4
      else
        k3=j4
        k4=j3
      endif
      if (i56 == 1) then
        k5=j5
        k6=j6
      else
        k5=j6
        k6=j5
      endif

      help=ah(i12)//ah(i34)//ah(i56)//'p'
      helm=ah(i12)//ah(i34)//ah(i56)//'m'

c--- positive helicity gluon
c--- down-type amplitudes
      call AZZjeta_rat(k1,k2,k3,k4,k5,k6,j7,za,zb,coeff,A7a)
      call AZZjeta_lc(help,swap12,.false.,p,k1,k2,k3,k4,k5,k6,j7,za,zb,.false.,scints,coeff,A7a,A7a_lc)
      call AZZjeta_slc(help,swap12,.false.,p,k1,k2,k3,k4,k5,k6,j7,za,zb,.false.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'plus A7a_lc ',A7a_lc
c      write(6,*) 'plus A7a_slc',A7a_slc
c      write(6,*) 'plus A7q',A7q
      A7alo(i12,i34,i56,2)=A7a
      A7av(i12,i34,i56,2)=zlc*A7a_lc+zslc*A7a_slc

      if (srdiags) then
        do Qid=1,2
          call AZZjetsr_amps(Qid,i12,i34,i56,k1,k2,k3,k4,k5,k6,j7,za,zb,A7a,A7a_lc,A7a_slc)
          A7asr(Qid,i12,i34,i56,2)=A7a
          A7avsr(Qid,i12,i34,i56,2)=zlc*A7a_lc+zslc*A7a_slc
        enddo
      else
        A7asr(:,i12,i34,i56,2)=czip
        A7avsr(:,i12,i34,i56,2)=czip
      endif

c--- negative helicity gluon
c--- down-type amplitudes
      call AZZjeta_rat(k2,k1,k5,k6,k3,k4,j7,zb,za,coeff,A7a)
      call AZZjeta_lc(helm,.not.(swap12),.true.,p,k2,k1,k5,k6,k3,k4,j7,zb,za,.true.,scints,coeff,A7a,A7a_lc)
      call AZZjeta_slc(helm,.not.(swap12),.true.,p,k2,k1,k5,k6,k3,k4,j7,zb,za,.true.,scints,coeff,A7a,A7a_slc)
c      write(6,*) 'minus A7a_lc ',A7a_lc
c      write(6,*) 'minus A7a_slc',A7a_slc
c      write(6,*) 'minus A7q',A7q
      A7alo(i12,i34,i56,1)=A7a
      A7av(i12,i34,i56,1)=zlc*A7a_lc+zslc*A7a_slc

      if (srdiags) then
        do Qid=1,2
          call AZZjetsr_amps(Qid,i12,i34,i56,k2,k1,k4,k3,k6,k5,j7,zb,za,A7a,A7a_lc,A7a_slc)
          A7asr(Qid,i12,i34,i56,1)=-A7a
          A7avsr(Qid,i12,i34,i56,1)=-(zlc*A7a_lc+zslc*A7a_slc)
        enddo
      else
        A7asr(:,i12,i34,i56,1)=czip
        A7avsr(:,i12,i34,i56,1)=czip
      endif

      enddo
      enddo
      enddo

c--- strong coupling renormalization in dred scheme
      A7av(:,:,:,:)=A7av(:,:,:,:)-zlc*(b0/xn*epinv-1d0/6d0)*A7alo(:,:,:,:)
      A7avsr(:,:,:,:,:)=A7avsr(:,:,:,:,:)-zlc*(b0/xn*epinv-1d0/6d0)*A7asr(:,:,:,:,:)
      if (scheme  ==  'tH-V') then
c--- known translation rules
        A7av(:,:,:,:)=A7av(:,:,:,:)-(Cf+xn/6d0)*A7alo(:,:,:,:)
        A7avsr(:,:,:,:,:)=A7avsr(:,:,:,:,:)-(Cf+xn/6d0)*A7asr(:,:,:,:,:)
      endif

      return
      end

