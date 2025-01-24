!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function DOTKS(I,J)
      implicit none
      include 'types.f'
      integer:: I,J
      real(dp):: DOTKS
      real(dp):: PLAB(4,10)
      COMMON/MOM/PLAB
!$omp threadprivate(/MOM/)

         DOTKS=PLAB(4,I)*PLAB(4,J)-PLAB(3,I)*PLAB(3,J)
     &        -PLAB(2,I)*PLAB(2,J)-PLAB(1,I)*PLAB(1,J)
      RETURN
      END

