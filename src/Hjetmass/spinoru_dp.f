!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine spinoru_dp(N,p,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      complex(dp), intent(out) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      real(dp) :: s(mxpart,mxpart)
      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: N

      call spinoru_dp_s(N,p,za,zb,s)
      return
      end


      subroutine spinoru_dp_s(N,p,za,zb,s)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      double precision s(mxpart,mxpart)
      double precision p(mxpart,4),rt(mxpart)
      double complex c23(mxpart),f(mxpart)
      integer i,j,N

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

c-----positive energy case
         if (p(j,4)  >  0d0) then
            rt(j)=sqrt(p(j,4)+p(j,1))
            c23(j)=cmplx(p(j,3),-p(j,2),kind=dp)
            f(j)=cone
         else
c-----negative energy case
            rt(j)=sqrt(-p(j,4)-p(j,1))
            c23(j)=cmplx(-p(j,3),p(j,2),kind=dp)
            f(j)=im
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*cmplx(rt(j)/rt(i),kind=dp)
     &    -c23(j)*cmplx(rt(i)/rt(j),kind=dp))

         if (abs(s(i,j)) < 1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
         else
         zb(i,j)=-cmplx(s(i,j),kind=dp)/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
