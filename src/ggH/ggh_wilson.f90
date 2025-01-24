!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
module ggHwilson
    use types
    implicit none

    integer, save, public  :: expansionorder, Wilsonorder
!$omp threadprivate(Wilsonorder)

    public :: ggHexpand

    private

    contains

      function ggHexpand(expandWilson)
      implicit none
      include 'types.f'
      include 'nflav.f'
      include 'masses.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'kpart.f'
      integer, intent(in) :: expandWilson
      real(dp) :: ggHexpand
      
!      write(6,*) 'C1',ason2pi*11._dp
!      write(6,*) 'C2',ason2pi**2*(
!     &    1933._dp/18._dp - 67._dp/12._dp*real(nflav,dp)
!     &    -(19._dp/2._dp + 8._dp/3._dp*real(nflav,dp))*log(mt**2/musq))

      if (expandWilson < 1) then
        ggHexpand=1._dp
        return
      endif
      
      if (expandWilson == 1) then
        ggHexpand=1._dp+ason2pi*11._dp
        if ((expansionorder == 2) .and. (coeffonly)) ggHexpand=ason2pi*11._dp
        return
      endif
      
      if (expandWilson == 2) then
        ggHexpand=1._dp+ason2pi*11._dp+ason2pi**2*( &
          1933._dp/18._dp - 67._dp/12._dp*real(nflav,dp) &
          -(19._dp/2._dp + 8._dp/3._dp*real(nflav,dp))*log(mt**2/musq))
        if ((expansionorder == 2) .and. (coeffonly)) then
          ggHexpand=ason2pi*11._dp+ason2pi**2*( &
          1933._dp/18._dp - 67._dp/12._dp*real(nflav,dp) &
          -(19._dp/2._dp + 8._dp/3._dp*real(nflav,dp))*log(mt**2/musq))
        endif
        return
      endif
      
      write(6,*) 'Unexpected parameter: expandWilson = ',expandWilson
      stop
      
      return
      end


end module
