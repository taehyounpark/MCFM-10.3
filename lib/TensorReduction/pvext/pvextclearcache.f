      subroutine pvextclearcache
      implicit none
      include 'lib/TensorReduction/Include/TRextclear.f'
      include 'lib/TensorReduction/Include/pvRespectmaxcindex.f'
      include 'lib/TensorReduction/Include/TRbadpoint.f'
      integer j
      do j=1,5
      clear(j)=.true.
      enddo
      pvRespectmaxcindex=.true.
      pvbadpoint=.false.
      return
      end
