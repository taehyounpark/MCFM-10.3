      function pvfndd(n,x,iep)
C----Implementation of DD Eq. 4.11
        use mod_qcdloop_c
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRonshellcutoff.f'
      integer j,n,infty
      complex(dp):: pvfndd,xm1,x,cone
      real(dp):: iep
      parameter(cone=(1d0,0d0),infty=16) ! number of terms in sum
      
      xm1=x-cone
      if (abs(x) .lt. 10d0) then
        if (abs(x-cone) .lt. onshellcutoff) then
          pvfndd=0d0
        else
          pvfndd=(cone-x**(n+1))*(cln(xm1,iep)-cln(x,iep))
        endif
        do j=0,n
          pvfndd=pvfndd-x**(n-j)/dfloat(j+1)
        enddo
      else
        pvfndd=cln(cone-cone/x,iep)
        do j=n+1,n+infty
          pvfndd=pvfndd+x**(n-j)/dfloat(j+1)
        enddo
      endif
           
      return
      end
