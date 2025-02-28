!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function etmiss(p,etvec)
      implicit none
      include 'types.f'
      real(dp):: etmiss
      include 'mxpart.f'
      integer:: j,k
      logical:: is_neutrino,is_darkmatter
      real(dp):: etvec(4),p(mxpart,4)

      do k=1,4
        etvec(k)=0._dp
      enddo

      do j=1,mxpart
        if (is_neutrino(j) .or. is_darkmatter(j)) then
          do k=1,4
            etvec(k)=etvec(k)+p(j,k)
          enddo
        endif
      enddo

      etmiss=sqrt(etvec(1)**2+etvec(2)**2)

      return
      end

