C=============================================================================
C---  HPLs  of Rank 1  
C=============================================================================
c --- main forking function
      double  complex function HPL1(n1, x)
      implicit none 
      integer n1
      double complex x,ris
      double complex HPL1at0,HPL1at1,HPL1atm1
      double complex HPL1ar1,HPL1arm1,HPL1ar0
      double complex HPL1else
      double precision rad1,radm1,rad0
      
      rad0 = 0.025d0
      rad1 = 0.01d0
      radm1 = 0.025d0      

      if (abs(n1).gt.1) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL1:"
         print*, "Index",n1," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL1at0(n1)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL1at1(n1)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL1atm1(n1)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL1ar1(n1,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL1arm1(n1,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL1ar0(n1,x)
      else 
         ris = HPL1else(n1,x)
      endif

      HPL1=ris 
      return
      end
c ------------------------------------------------
      double complex function HPL1at0(n1)
      implicit none
      integer n1
      double complex ris

      ris = dcmplx(0d0,0d0)

      if (n1.eq.0) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1at0=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1at1(n1)
      implicit none
      integer n1
      double complex ris
      double precision ll2

      ll2 = dlog(2d0)
      ris = dcmplx(0d0,0d0)

      if(n1.ne.1) then
         select case (n1)
         case(-1)
            ris = ll2
         case(0)
            ris = 0d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1at1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1atm1(n1)
      implicit none
      integer n1
      double complex ris,myi
      double precision ll2,pi

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)

      if(n1.ne.-1) then
         select case (n1)
         case(0)
            ris = myi*pi
         case(1)
            ris = -ll2
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1atm1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1ar1(n1,x)
      implicit none
      integer n1,bcflag
      double complex x,ris,zp,llzp
      double precision pi,ll2,xre

      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)
      bcflag = 0
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
         case(-1)            !-1

            zp = 1d0-x

            ris = -((zp)/2d0) - (zp**2)/8d0 - (zp**3)/2
     &4d0 - (zp**4)/64d0 - (zp**5)/160d0 - (zp**6)/384d0 + ll
     &2

         case(0)            !0

            zp = 1d0-x

            ris = -zp - (zp**2)/2d0 - (zp**3)/3d0 - (zp
     &**4)/4d0 - (zp**5)/5d0 - (zp**6)/6d0

         case(1)            !1

            zp = 1d0-x
            llzp = log(zp)

            ris = -llzp
c End of expansions around x = +1
      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c            
         else if (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      
      HPL1ar1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1arm1(n1,x)
      implicit none
      integer n1,s,szp,bcflag
      double complex x,ris,zp,llzp,myi
      double precision pi,ll2,xre
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)
      myi = dcmplx(0d0,1d0)
      
      ris = dcmplx(0d0,0d0)
      bcflag = 0
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
      case(-1)                  !-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = llzp
         
      case(0)                   !0
         
         zp = x+1d0
         szp = s(zp)
         
         ris = myi*pi*szp - zp - (zp**2)/2d0 - (zp**
     &        3)/3d0 - (zp**4)/4d0 - (zp**5)/5d0 - (zp**6)/6d0 - (zp*
     &        *7)/7d0 - (zp**8)/8d0 - (zp**9)/9d0
         
      case(1)                   !1
         
         zp = x+1d0
         
         ris = (zp)/2d0 + (zp**2)/8d0 + (zp**3)/24d0
     &        + (zp**4)/64d0 + (zp**5)/160d0 + (zp**6)/384d0 + (zp**
     &        7)/896d0 + (zp**8)/2048d0 + (zp**9)/4608d0 - ll2
c     End of expansions around x = -1
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      HPL1arm1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1ar0(n1,x)
      implicit none
      integer n1,bcflag
      double complex x,ris,llx
      double precision pi,ll2,xre
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
      case(-1)                  !-1
         
         
         ris = x - (x**2)/2d0 + (x**3)/3d0 - (x**4)/
     &        4d0 + (x**5)/5d0 - (x**6)/6d0 + (x**7)/7d0 - (x**8)/8d0
     &        + (x**9)/9d0 - (x**10)/10d0
         
      case(0)                   !0
         
         llx = log(x)
         
         ris = llx
         
      case(1)                   !1
         
         
         ris = x + (x**2)/2d0 + (x**3)/3d0 + (x**4)/
     &        4d0 + (x**5)/5d0 + (x**6)/6d0 + (x**7)/7d0 + (x**8)/8d0
     &        + (x**9)/9d0 + (x**10)/10d0
c     End of expansions around x = 0
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      HPL1ar0=ris
      return
      end function
c ------------------------------------------------
      double  complex function HPL1else(n1, x)
      implicit none
      double complex x, ris
      integer n1,bcflag
      double precision xre
      
      bcflag = 0
      ris = dcmplx(0d0,0d0)
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

      select case(n1)
      case(-1)
         ris =log(1.0d0 + x)
      case(0)
         ris=log(x)
      case(1)
         ris=-log(1.0d0 - x)
      end select

c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        

      HPL1else=ris 
      return
      end
c ------------------------------------------------
c --- Real part of HPL1     
      double precision function HPL1real(n1,xr,xi)
      implicit none
      double precision xr,xi
      integer n1
      double complex x,HPL1
      x=dcmplx(xr,xi)
      HPL1real = dreal(HPL1(n1,x))
      return
      end

c --- Imaginary part of HPL1     
      double precision function HPL1im(n1,xr,xi)
      implicit none
      double precision xr,xi
      integer n1
      double complex x,HPL1
      x=dcmplx(xr,xi)
      HPL1im = dimag(HPL1(n1,x))
      return
      end
