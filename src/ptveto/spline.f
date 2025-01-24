!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)  
      implicit none
      include 'types.f'      
      INTEGER:: n,i,k  
      integer,PARAMETER::NMAX=500
      REAL(dp):: yp1,ypn,x(n),y(n),y2(n)  
      REAL(dp):: p,qn,sig,un,u(NMAX)  
      if (yp1 > 0.99d30) then  
         y2(1)=0._dp  
         u(1)=0._dp  
      else  
         y2(1)=-0.5_dp  
         u(1)=(3._dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
      endif  
      do i=2,n-1  
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
         p=sig*y2(i-1)+2._dp  
         y2(i)=(sig-1._dp)/p  
         u(i)=(6._dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &  /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p  
	enddo
       if (ypn > 0.99d30) then  
         qn=0.  
         un=0.  
       else  
         qn=0.5_dp  
         un=(3._dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
       endif  
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._dp)  
       do k=n-1,1,-1  
         y2(k)=y2(k)*y2(k+1)+u(k)  
       enddo
       return  
       end
  
