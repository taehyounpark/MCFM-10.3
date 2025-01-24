      subroutine pvextAtensor(m1s,FA0,FA1,FA2)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvAnames.f'
      include 'lib/TensorReduction/Include/pvextAv.f'
      include 'lib/TensorReduction/Include/TRydef.f'
      include 'lib/TensorReduction/Include/TRmetric.f'
      complex(dp):: FA0(-2:0),FA1(y1max,-2:0),FA2(y2max,-2:0)
      real(dp)::m1s
      integer n1,n2,A0i,pvextAcache
      logical,save:: first=.true.
!$omp threadprivate(first)
      if (first) then
      first=.false.
      call pvarraysetup
      endif

      A0i=pvextAcache(m1s)

      FA0(:)=Av(A0i+aa0,:)

      do n1=1,4
      FA1(n1,:)=czip
      enddo

      do n1=1,4
      do n2=n1,4 
      FA2(y2(n1,n2),:)=g(n1,n2)*Av(A0i+aa00,:)
      enddo
      enddo

      return
      end

