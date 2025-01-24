!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      integer,intent(in):: N
      real(dp), intent(in)::p(mxpart,4)
      complex(dp), intent(out) :: za(mxpart,mxpart),zb(mxpart,mxpart)
      real(dp):: rt(mxpart)
      complex(dp):: c23(mxpart),f(mxpart)
      integer:: i,j

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=cmplx(0,0,kind=dp)
         zb(j,j)=za(j,j)

c-----positive energy case
         if (p(j,4) > 0) then
            rt(j)=sqrt(p(j,4)+p(j,1))
            c23(j)=cmplx(p(j,3),-p(j,2),kind=dp)
            f(j)=cmplx(1,0,kind=dp)
         else
c-----negative energy case
            rt(j)=sqrt(-p(j,4)-p(j,1))
            c23(j)=cmplx(-p(j,3),p(j,2),kind=dp)
            f(j)=cmplx(0,1,kind=dp)
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=2._dp*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))

c Disabled this check August 2018 since it can cause poor cancellation
c Reinstated with a much tighter check for 4f single top calculation
         if (abs(s(i,j))<1.e-10_dp) then
           zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
         else
           zb(i,j)=-s(i,j)/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

