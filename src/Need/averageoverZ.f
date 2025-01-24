!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
      subroutine averageoverZ(p,msq,order)
      use nnlo_z1jet
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'zcouple_cms.f'
      integer:: kcount,jcount,order
      real(dp):: msq(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),
     & p(mxpart,4),pout(mxpart,4),statfac
      pout(:,:)=0d0
      msq(:,:)=0d0
C---sum over 6 values;
      do jcount=1,6
      call pgen(p,jcount,3,4,pout)
      call qqb_z1jet_vbis(pout,msq1,order)
!      call writeout(pout)
!      write(6,*) 'msq1 bit',msq1(2,-2)
      msq(:,:)=msq(:,:)+msq1(:,:)
      enddo
      msq(:,:)=msq(:,:)*(zwidth**2/(4d0*abs(zesq)*(zle**2+zre**2)))
      return
      end
