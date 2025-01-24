!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      SUBROUTINE splint(xa,ya,y2a,n,x,y)  
      implicit none
      include 'types.f'
      INTEGER::n  
      REAL(dp):: x,y,xa(n),y2a(n),ya(n)  
      INTEGER::k,khi,klo  
      REAL(dp):: a,b,h  
      klo=1  
      khi=n  
 1    if (khi-klo > 1) then  
        k=(khi+klo)/2  
        if(xa(k) > x)then  
          khi=k  
        else  
          klo=k  
        endif  
      goto 1  
      endif  
      h=xa(khi)-xa(klo)  
      if (h == 0.) stop 'bad xa input in splint'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)
     & +((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6._dp  
      return  
      end
