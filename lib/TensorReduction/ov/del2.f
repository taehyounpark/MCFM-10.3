      function del2(p1,p2,p3,p4)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      real(dp):: del2,ovdot,p1(4),p2(4),p3(4),p4(4)
      del2=ovdot(p1,p3)*ovdot(p2,p4)-ovdot(p1,p4)*ovdot(p2,p3)
      return
      end
      
