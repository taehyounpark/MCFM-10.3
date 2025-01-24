!----- kinematic file for WpWmg production

      xscale=10d0

      p3true(1)=7._dp
      p3true(2)=2._dp
      p3true(3)=-6._dp
      p3true(4)=+3._dp

      p4true(1)=7._dp
      p4true(2)=-2._dp
      p4true(3)=+6._dp
      p4true(4)=-3._dp

      p1true(1)=-9._dp
      p1true(2)=-8._dp
      p1true(3)=-1._dp
      p1true(4)=-4._dp

      p2true(1)=-15._dp
      p2true(2)=+11._dp
      p2true(3)=+2._dp
      p2true(4)=+10._dp

      p7true(1)=-21._dp
      p7true(2)=-13._dp
      p7true(3)=+4._dp
      p7true(4)=-16._dp

      p5true(1)=5._dp
      p5true(2)=-3._dp
      p5true(3)=0._dp
      p5true(4)=4._dp

      k56(:)=
     & -p1true(:)-p2true(:)-p3true(:)-p4true(:)-p7true(:)
      k56sq=k56(1)**2-k56(2)**2-k56(3)**2-k56(4)**2
      k56Dp5=k56(1)*p5true(1)-k56(2)*p5true(2)
     & -k56(3)*p5true(3)-k56(4)*p5true(4)
      alpha=k56sq/(2._dp*k56Dp5)
      p5true(:)=alpha*p5true(:)
      p6true(:)=k56(:)-p5true(:)

!      write(6,*)
!      write(6,*) p1true(1)**2-p1true(2)**2-p1true(3)**2-p1true(4)**2
!      write(6,*) p2true(1)**2-p2true(2)**2-p2true(3)**2-p2true(4)**2
!      write(6,*) p3true(1)**2-p3true(2)**2-p3true(3)**2-p3true(4)**2
!      write(6,*) p4true(1)**2-p4true(2)**2-p4true(3)**2-p4true(4)**2
!      write(6,*) p5true(1)**2-p5true(2)**2-p5true(3)**2-p5true(4)**2
!      write(6,*) p6true(1)**2-p6true(2)**2-p6true(3)**2-p6true(4)**2
!      write(6,*) p7true(1)**2-p7true(2)**2-p7true(3)**2-p7true(4)**2
!      write(6,*) p1true(:)+p2true(:)+p3true(:)+p4true(:)+p5true(:)
!     & +p6true(:)+p7true(:)

      p1true(:)=p1true(:)/xscale
      p2true(:)=p2true(:)/xscale
      p3true(:)=p3true(:)/xscale
      p4true(:)=p4true(:)/xscale
      p5true(:)=p5true(:)/xscale
      p6true(:)=p6true(:)/xscale
      p7true(:)=p7true(:)/xscale

!      write(6,4) p1true(1),',',p1true(2),',',p1true(3),',',p1true(4),' '
!      write(6,4) p2true(1),',',p2true(2),',',p2true(3),',',p2true(4),' '
!      write(6,4) p3true(1),',',p3true(2),',',p3true(3),',',p3true(4),' '
!      write(6,4) p4true(1),',',p4true(2),',',p4true(3),',',p4true(4),' '
!      write(6,4) p5true(1),',',p5true(2),',',p5true(3),',',p5true(4),' '
!      write(6,4) p6true(1),',',p6true(2),',',p6true(3),',',p6true(4),' '
!      write(6,4) p7true(1),',',p7true(2),',',p7true(3),',',p7true(4),' '
!      pause

 4    format(4(f20.15,A))
c--- now form the momenta that will be used in the Kirill routines
c--- NOTE that u and ubar need to be swapped
      do nu=1,4
      p2(nu)=p2true(nu)
      p1(nu)=p1true(nu)
      p3(nu)=p3true(nu)+p4true(nu)
      p4(nu)=p5true(nu)+p6true(nu)
      p5(nu)=p7true(nu)
      enddo


