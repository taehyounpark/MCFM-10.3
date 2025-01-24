C=============================================================================
C---  HPLs  of Rank 3  
C=============================================================================
c--- main forking function
      double complex function HPL3(n1, n2, n3, x)
      implicit none
      integer n1,n2,n3
      double complex x,ris
      double complex HPL3at0,HPL3at1,HPL3atm1
      double complex HPL3ar1,HPL3arm1,HPL3ar0,HPL3else
      double precision rad1,radm1,rad0

      rad1 = 0.01d0
      radm1 = 0.025d0
      rad0 = 0.025d0
      
      if ((abs(n1).gt.1).or.(abs(n2).gt.1).or.(abs(n3).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL3:"
         print*, "Indices",n1,n2,n3," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL3at0(n1,n2,n3)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL3at1(n1,n2,n3)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL3atm1(n1,n2,n3)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL3ar1(n1,n2,n3,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL3arm1(n1,n2,n3,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL3ar0(n1,n2,n3,x)
      else 
         ris = HPL3else(n1,n2,n3,x)
      endif
      HPL3=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL3at0(n1, n2, n3)
      implicit none
      integer n1,n2,n3,j
      double complex ris
    
      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

      if (j.eq.14) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3at0=ris
      return
      end function
c ------------------------------------------------
c --- Real part of HPL3     
      double precision function HPL3real(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3real = dreal(HPL3(n1,n2,n3,x))
      return
      end

c --- Imaginary part of HPL3     
      double precision function HPL3im(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3im = dimag(HPL3(n1,n2,n3,x))
      return
      end
