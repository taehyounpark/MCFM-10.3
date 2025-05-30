      subroutine pvvcalc(p2,p3,p4,p5,v2,v3,v4,v5)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      real(dp):: p2(4),p3(4),p4(4),p5(4),v2(4),v3(4),v4(4),v5(4)
      real(dp):: v2Dp2,v3Dp3,v4Dp4,v5Dp5
      integer j
      v2(1)=+(p3(2)*p4(3)*p5(4)+p3(4)*p4(2)*p5(3)+p3(3)*p4(4)*p5(2)
     .       -p3(2)*p4(4)*p5(3)-p3(3)*p4(2)*p5(4)-p3(4)*p4(3)*p5(2))
      v2(2)=-(p3(3)*p4(4)*p5(1)+p3(1)*p4(3)*p5(4)+p3(4)*p4(1)*p5(3)
     .       -p3(3)*p4(1)*p5(4)-p3(4)*p4(3)*p5(1)-p3(1)*p4(4)*p5(3))
      v2(3)=+(p3(4)*p4(1)*p5(2)+p3(2)*p4(4)*p5(1)+p3(1)*p4(2)*p5(4)
     .       -p3(4)*p4(2)*p5(1)-p3(1)*p4(4)*p5(2)-p3(2)*p4(1)*p5(4))
      v2(4)=+(p3(1)*p4(2)*p5(3)+p3(3)*p4(1)*p5(2)+p3(2)*p4(3)*p5(1)
     .       -p3(1)*p4(3)*p5(2)-p3(2)*p4(1)*p5(3)-p3(3)*p4(2)*p5(1))

      v3(1)=+(p2(2)*p4(3)*p5(4)+p2(4)*p4(2)*p5(3)+p2(3)*p4(4)*p5(2)
     .       -p2(2)*p4(4)*p5(3)-p2(3)*p4(2)*p5(4)-p2(4)*p4(3)*p5(2))
      v3(2)=-(p2(3)*p4(4)*p5(1)+p2(1)*p4(3)*p5(4)+p2(4)*p4(1)*p5(3)
     .       -p2(3)*p4(1)*p5(4)-p2(4)*p4(3)*p5(1)-p2(1)*p4(4)*p5(3))
      v3(3)=+(p2(4)*p4(1)*p5(2)+p2(2)*p4(4)*p5(1)+p2(1)*p4(2)*p5(4)
     .       -p2(4)*p4(2)*p5(1)-p2(1)*p4(4)*p5(2)-p2(2)*p4(1)*p5(4))
      v3(4)=+(p2(1)*p4(2)*p5(3)+p2(3)*p4(1)*p5(2)+p2(2)*p4(3)*p5(1)
     .       -p2(1)*p4(3)*p5(2)-p2(2)*p4(1)*p5(3)-p2(3)*p4(2)*p5(1))

      v4(1)=+(p2(2)*p3(3)*p5(4)+p2(4)*p3(2)*p5(3)+p2(3)*p3(4)*p5(2)
     .       -p2(2)*p3(4)*p5(3)-p2(3)*p3(2)*p5(4)-p2(4)*p3(3)*p5(2))
      v4(2)=-(p2(3)*p3(4)*p5(1)+p2(1)*p3(3)*p5(4)+p2(4)*p3(1)*p5(3)
     .       -p2(3)*p3(1)*p5(4)-p2(4)*p3(3)*p5(1)-p2(1)*p3(4)*p5(3))
      v4(3)=+(p2(4)*p3(1)*p5(2)+p2(2)*p3(4)*p5(1)+p2(1)*p3(2)*p5(4)
     .       -p2(4)*p3(2)*p5(1)-p2(1)*p3(4)*p5(2)-p2(2)*p3(1)*p5(4))
      v4(4)=+(p2(1)*p3(2)*p5(3)+p2(3)*p3(1)*p5(2)+p2(2)*p3(3)*p5(1)
     .       -p2(1)*p3(3)*p5(2)-p2(2)*p3(1)*p5(3)-p2(3)*p3(2)*p5(1))

      v5(1)=+(p2(2)*p3(3)*p4(4)+p2(4)*p3(2)*p4(3)+p2(3)*p3(4)*p4(2)
     .       -p2(2)*p3(4)*p4(3)-p2(3)*p3(2)*p4(4)-p2(4)*p3(3)*p4(2))
      v5(2)=-(p2(3)*p3(4)*p4(1)+p2(1)*p3(3)*p4(4)+p2(4)*p3(1)*p4(3)
     .       -p2(3)*p3(1)*p4(4)-p2(4)*p3(3)*p4(1)-p2(1)*p3(4)*p4(3))
      v5(3)=+(p2(4)*p3(1)*p4(2)+p2(2)*p3(4)*p4(1)+p2(1)*p3(2)*p4(4)
     .       -p2(4)*p3(2)*p4(1)-p2(1)*p3(4)*p4(2)-p2(2)*p3(1)*p4(4))
      v5(4)=+(p2(1)*p3(2)*p4(3)+p2(3)*p3(1)*p4(2)+p2(2)*p3(3)*p4(1)
     .       -p2(1)*p3(3)*p4(2)-p2(2)*p3(1)*p4(3)-p2(3)*p3(2)*p4(1))

      v2Dp2=v2(4)*p2(4)-v2(1)*p2(1)-v2(2)*p2(2)-v2(3)*p2(3)
      v3Dp3=v3(4)*p3(4)-v3(1)*p3(1)-v3(2)*p3(2)-v3(3)*p3(3)
      v4Dp4=v4(4)*p4(4)-v4(1)*p4(1)-v4(2)*p4(2)-v4(3)*p4(3)
      v5Dp5=v5(4)*p5(4)-v5(1)*p5(1)-v5(2)*p5(2)-v5(3)*p5(3)
      do j=1,4
      v2(j)=v2(j)/v2Dp2
      v3(j)=v3(j)/v3Dp3
      v4(j)=v4(j)/v4Dp4
      v5(j)=v5(j)/v5Dp5
      enddo
      return
      end
