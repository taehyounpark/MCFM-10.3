!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c      use Hgggg_integralfill_generic
      implicit none

      include 'Inc/hgggglabels.f'
      include 'Inc/IntResults.f'
      complex(dp) :: hgggg_assemble_res,Dcoeff(dmax),Ccoeff(cmax),Bcoeff(bmax),Rat
      integer k

      hgggg_assemble_res=czip
      do k=1,dmax
        hgggg_assemble_res=hgggg_assemble_res+Dcoeff(k)*Dint(k)
c        write(6,*) 'k,Dcoeff(k),Dint(k)',k,Dcoeff(k),Dint(k)
      enddo
c        write(6,*)
      do k=1,cmax
        hgggg_assemble_res=hgggg_assemble_res+Ccoeff(k)*Cint(k)
c        write(6,*) 'k,Ccoeff(k),Cint(k)',k,Ccoeff(k),Cint(k)
      enddo
c        write(6,*)
      do k=1,bmax
        hgggg_assemble_res=hgggg_assemble_res+Bcoeff(k)*Bint(k)
c        write(6,*) 'k,Bcoeff(k),Bint(k)',k,Bcoeff(k),Bint(k)
      enddo
c        write(6,*)
c      write(6,*) 'Rat',Rat
      hgggg_assemble_res=hgggg_assemble_res+Rat

      return

