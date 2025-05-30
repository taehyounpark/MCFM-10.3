      subroutine runCF_00iiii(i1,i2,i3,i4,f,Gr,Shat6,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,k,i1,i2,i3,i4,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np),den
      complex(dp):: Shat6(np,z5max,-2:0)
       
      do ep=-2,0
      if     ((i1 .eq. i2) .and. (i1 .eq. i3) .and. (i1 .eq. i4)) then
        den=10d0
        k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i3)) then
        den=8d0
        k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i4)) then
        den=8d0
        k=i1
      elseif ((i2 .eq. i3) .and. (i2 .eq. i4)) then
        den=8d0
        k=i2
      elseif ((i1 .eq. i2) .or. (i1 .eq. i3) .or. (i1 .eq. i4)) then
        den=6d0
        k=i1
      elseif ((i2 .eq. i3) .or. (i2 .eq. i4)) then
        den=6d0
        k=i2
      elseif (i3 .eq. i4) then
        den=6d0
        k=i3
      else
        den=4d0
        k=i1
      endif      
      Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (Shat6(k,z5(k,i1,i2,i3,i4),ep)
     . -f(k)*Cv(ciiiii(z5(k,i1,i2,i3,i4))+N0,ep)
     . -Gr(k,1)*Cv(ciiiiii(z6(1,k,i1,i2,i3,i4))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiiii(z6(2,k,i1,i2,i3,i4))+N0,ep))/den

      enddo
      
      return
      end
