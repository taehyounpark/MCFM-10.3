!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function dot(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: dot
      include 'mxpart.f'
      integer:: i,j
      real(dp):: p(mxpart,4),dotpr
!     dot=p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3)
      dot=dotpr(p(i,:),p(j,:))
      return
      end

      function dotvec(pi,pj)
      implicit none
      include 'types.f'
      real(dp) :: dotvec,dotpr
      real(dp), intent(in) :: pi(4), pj(4)

!      dotvec=pi(4)*pj(4)-pi(1)*pj(1)-pi(2)*pj(2)-pi(3)*pj(3)
      dotvec=dotpr(pi(:),pj(:))
      return
      end function

      pure function massvec(p)
        implicit none
        include 'types.f'

        real(dp) :: massvec
        real(dp), intent(in) :: p(4)

        massvec=p(4)*p(4)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)
      end


      function dotpr(p,q)
      use types
      implicit none
      real(dp):: dotpr
!     returns the dot product of the two four vectors p(4) and q(4)

      real(dp):: p(4),q(4),x1,x2,x3,sinth

      dotpr=p(4)*q(4)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)

      if ((p(4) == 0d0) .or. (q(4) == 0d0)) return 
!     because at least one of the vectors is not light-like, 
!     so the rescue is inappropriate

      if (abs(dotpr/(p(4)*q(4))) > 1.e-6_dp) return
!     because the rescaled dotproduct is large enough

! Check to ensure vectors are massless before using rescue code
      if (abs(p(4)**2-p(1)**2-p(2)**2-p(3)**2) > 1.e-8_dp) return
      if (abs(q(4)**2-q(1)**2-q(2)**2-q(3)**2) > 1.e-8_dp) return

! If the dotpr product is very small then replace the above calculation
! with one using the cross product -- only true for massless vectors!
      x1=p(2)*q(3)-q(2)*p(3)
      x2=q(1)*p(3)-p(1)*q(3)
      x3=p(1)*q(2)-q(1)*p(2)
      sinth=sqrt(x1**2+x2**2+x3**2)
     & /sqrt(p(1)**2+p(2)**2+p(3)**2)
     & /sqrt(q(1)**2+q(2)**2+q(3)**2)

      dotpr=2*p(4)*q(4)*sin(0.5_dp*asin(sinth))**2

! Protect against numerical problems from exactly-vanishing invariants
      if (dotpr == 0._dp) dotpr = 1.e-32_dp

      return
      end function dotpr

