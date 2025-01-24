!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!


      subroutine reset_aem(ae_in)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'ewcouple.f'
      include 'zcouple_cms.f'
      include 'first.f'
      include 'constants.f'
      include 'mpicommon.f'
      real(kind=dp), intent(in) :: ae_in
      real(kind=dp) :: esq_old,ae_old

      esq_old=esq
      ae_old=esq/(fourpi)

      esq=ae_in*fourpi
      zesq=cmplx(esq,0d0,kind=dp)

      if(first) then
         first=.false.
         if (rank == 0) then
         write(6,*) '*********** Changed alpha_EM **********'
         write(6,88) '*  old value of esq :',esq_old,'     *'
         write(6,88) '*  new value of esq :',esq,'     *'
         write(6,88) '*  new 1/alpha_EM   :',fourpi/esq,'     *'
         write(6,*) '***************************************'
         endif
      endif

 88   format(1x,a21,f12.8,a6)
      return
      end

