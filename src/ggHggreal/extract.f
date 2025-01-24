!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine extract(target_struc,target_struc_v)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell                                            *
c     February, 2005.                                                  *
c                                                                      *
c     Extracts the matrix element colour structure from the common     *
c     include file 'msq_struc' and returns it in the arrays            *
c     target_struc (lowest order) and target_struc_v (lowest order     *
c     contracted with a vector)                                        *
c                                                                      *
c     For efficiency, this is only done for indices (-1,0,+1)          *
c                                                                      *
c***********************************************************************
      include 'nf.f'
      include 'msq_struc.f'
      integer::istruc,j,k
      real(dp)::target_struc(8,-1:1,-1:1),
     &               target_struc_v(8,-1:1,-1:1)

      do j=-1,1
      do k=-1,1
      do istruc=1,8
        target_struc(istruc,j,k)=msq_struc(istruc,j,k)
c--- dummy for now: needs to be filled once the appropriate
c--- modifications to the gvec routine have been made
        target_struc_v(istruc,j,k)=msq_strucv(istruc,j,k)
      enddo
      enddo
      enddo

      return
      end

