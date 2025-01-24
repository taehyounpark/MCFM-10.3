!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine msquarer(AB,BA,res)
c---Author: R.K. Ellis
c---Squares colour-ordered amplitude
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::AB(2,2,2,2,2),BA(2,2,2,2,2)
      real(dp)::res
      integer::j1,j2,jq,jl1,jl2
      res=0._dp
      do j1=1,2
      do j2=1,2
      do jq=1,2
      do jl1=1,2
      do jl2=1,2
      res=res
     & +cdabs(AB(j1,j2,jq,jl1,jl2))**2+cdabs(BA(j1,j2,jq,jl1,jl2))**2
     & -cdabs(AB(j1,j2,jq,jl1,jl2)+BA(j1,j2,jq,jl1,jl2))**2/xn**2
       enddo
       enddo
       enddo
       enddo
       enddo

       return
       end
