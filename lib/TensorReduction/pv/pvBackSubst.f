      subroutine pvBackSubst(A,n,p,b)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      integer n,j,p(n)
      integer,parameter::nmax=4
      complex(dp):: A(n,n),b(n,-2:0)
      complex(dp)::xm2(nmax),xm1(nmax),xm0(nmax)
      if (n .gt. nmax) then
      write(6,*) 'Error in BackSubst, n .gt. nmax'
      stop
      endif
 
      do j=1,n
      xm2(j)=b(j,-2)
      xm1(j)=b(j,-1)
      xm0(j)=b(j, 0)
      enddo

      call XLUBackSubst(A, n, p, xm2)
      call XLUBackSubst(A, n, p, xm1)
      call XLUBackSubst(A, n, p, xm0)

      do j=1,n
      b(j,-2)=xm2(j)
      b(j,-1)=xm1(j)
      b(j, 0)=xm0(j)
      enddo

      return
      end
