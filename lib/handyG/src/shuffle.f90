
MODULE shuffle
  use globals
  use ieps
  implicit none

CONTAINS
  
  FUNCTION append_to_each_row(a, m) result(res)
    ! appends element a to each row of m
    type(inum) :: a, m(:,:)
    integer :: i
    type(inum) :: res(size(m,1),size(m,2)+1)
    do i=1,size(m,1)
      res(i,:) = [a,m(i,:)]
    end do
  END FUNCTION append_to_each_row

  FUNCTION stack_matrices_vertically(m1, m2) result(res)
    ! appends to matrix m1 the rows of matrix m2
    type(inum) :: m1(:,:), m2(:,:)
    type(inum) :: res(size(m1,1)+size(m2,1), size(m1,2))
    res(1:size(m1,1), :) = m1
    res(size(m1,1)+1:size(res,1),:) = m2 
  END FUNCTION stack_matrices_vertically

  RECURSIVE FUNCTION shuffle_product(v1, v2) result(res)
    type(inum) :: v1(:), v2(:)
    integer :: i
    type(inum) :: res(product((/(i,i=1,size(v1)+size(v2))/))/  &
      (product((/(i,i=1,size(v1))/))*product((/(i,i=1,size(v2))/))), & 
      size(v1) + size(v2))
    type(inum) :: alpha, beta, w1(size(v1)-1), w2(size(v2)-1)

    res = izero
    if(size(v1) == 0) then 
      res(1,:) = v2
      return
    else if(size(v2) == 0) then 
      res(1,:) = v1
      return
    end if

    alpha = v1(1)
    beta = v2(1)
    w1 = v1(2:size(v1))
    w2 = v2(2:size(v2))

    res = stack_matrices_vertically( &
      append_to_each_row(alpha, shuffle_product(w1, v2)), & 
      append_to_each_row(beta, shuffle_product(v1, w2)) )

  END FUNCTION shuffle_product


  FUNCTION shuffle_with_zero(a) result(res)
    ! rows of result are shuffles of a with 0
    type(inum) :: a(:)
    type(inum) :: res(size(a)+1,size(a)+1)
    integer :: i,j, N
    N = size(a)+1
    do i = 1,N
      ! i is the index of the row
      ! j is the index of the zero
      j  = N+1-i
      res(i,j) = izero
      res(i,1:j-1) = a(1:j-1)
      res(i,j+1:N) = a(j:)
    end do
  END FUNCTION shuffle_with_zero


END MODULE shuffle

! PROGRAM test
!   use utils
!   use shuffle
!   implicit none

!   type(inum) :: v1(3), v2(2)
!   integer :: amount_shuffles

!   v1 = cmplx((/1,2,3/))
!   v2 = cmplx((/4,5/))

!   call print_matrix(shuffle_product(v1,v2))

! END PROGRAM test

