C=============================================================================
C---  HPLs  of Rank 2  
C=============================================================================
c --- main forking function
      double  complex function HPL2(n1,n2, x)
      implicit none 
      integer n1,n2
      double complex x,ris
      double complex HPL2at0,HPL2at1,HPL2atm1
      double complex HPL2ar1,HPL2arm1,HPL2ar0
      double complex HPL2else
      double precision rad1,radm1,rad0
      
      rad0 = 0.025d0
      rad1 = 0.01d0
      radm1 = 0.025d0      

      if ((abs(n1).gt.1).or.(abs(n2).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL2:"
         print*, "Indices",n1,n2," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL2at0(n1,n2)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL2at1(n1,n2)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL2atm1(n1,n2)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL2ar1(n1,n2,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL2arm1(n1,n2,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL2ar0(n1,n2,x)
      else 
         ris = HPL2else(n1,n2,x)
      endif

      HPL2=ris 
      return
      end
c ------------------------------------------------
      double complex function HPL2at0(n1, n2)
      implicit none
      integer n1,n2,j
      double complex ris
    
      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if (j.eq.5) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL2at0=ris
      return
      end function

c --- Real part of HPL2     
      double precision function HPL2real(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2real = dreal(HPL2(n1,n2,x))
      return
      end

c --- Imaginary part of HPL2     
      double precision function HPL2im(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2im = dimag(HPL2(n1,n2,x))
      return
      end
