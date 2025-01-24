      subroutine runCF_00iii(i1,i2,i3,f,Gr,Shat5,N0)
C---  Expression for rearrangement of extension of Eq. 5.70
C---  Calculates C00iii
C---  Small terms of order f(i)*Cijkl,Gr(i,j)*Cijklm
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np),den
      complex(dp):: Shat5(np,z4max,-2:0)
       
      do ep=-2,0
      if     ((i1 .eq. i2) .and. (i1 .eq. i3)) then
        den=8d0
        k=i1
      elseif (i1 .eq. i2) then
        den=6d0
        k=i1
      elseif (i1 .eq. i3) then
        den=6d0
        k=i1
      elseif (i2 .eq. i3) then
        den=6d0
        k=i2
      else
        den=4d0
        k=i1
      endif      
      Cv(czziii(z3(i1,i2,i3))+N0,ep)=
     . (Shat5(k,z4(k,i1,i2,i3),ep)
     . -f(k)*Cv(ciiii(z4(k,i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Cv(ciiiii(z5(1,k,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiii(z5(2,k,i1,i2,i3))+N0,ep))/den

      enddo
      
      return
      end
