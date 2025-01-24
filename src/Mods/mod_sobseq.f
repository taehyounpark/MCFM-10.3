       !> Module containing a state type for implementing parallel
       !> Sobol series generators with strided and skip-ahead generation.
       !> Note that this uses the gray code implementation, so
       !> Generated numbers are shuffled compared to the original series.

       ! original code from https://github.com/Exteris/sobseq
       ! authors: Daan van Vugt daanvanvugt@gmail.com
       !    and Koen Beljaars k.p.beljaars@tue.nl

       ! modified to use 64 bit integers to allow for 2^63 instead of 2^31 points

cCopyright (c) 2016 Daan van Vugt, Koen Beljaars
cCopyright (C) 2019-2022, respective authors of MCFM.

cPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
c(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify,merge,
cpublish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do
cso, subject to the following conditions:

cThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

cTHE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
cMERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
cFOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
cCONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

       module mod_sobseq
         use types
         implicit none
         private

         integer(kind=8), parameter :: N_M = 63 ! Generate at most 2^63 points
         ! Problems with sign bit when generating N_M = 32

         !> Type containing the state of a sobol sequence
         type, public :: sobol_state
             integer(kind=8) :: v(N_M)   !< Direction numbers
             integer(kind=8) :: i = 1    !< Current number
             integer(kind=8) :: x = 0    !< Current value
             integer(kind=8) :: stride=0 !< Skip 2^this many values when generating
         contains
             procedure, public :: initialize !< Initialize direction numbers
             procedure, public :: skip_ahead !< Skip ahead to a specific position and return this value
             procedure, public :: next         !< Generate the next value in the sequence
             procedure, public :: next_strided !< Generate the next value in the sequence (strided version)
         end type sobol_state
       contains


       !> Initialize the direction numbers using a primitive polynomial
       subroutine initialize(state, s, a, m_in, stride)
         implicit none
         class (sobol_state), intent(inout) :: state
         integer(kind=8), intent(in) :: s !< Number of direction numbers / Mathematical polynomial basis of degree s
         integer(kind=8), intent(in) :: a !< Coefficients of primitive polynomial
         integer(kind=8), intent(in), dimension(s) :: m_in !< First direction numbers
         integer(kind=8), intent(in), optional :: stride

         integer(kind=8), dimension(N_M) :: m
         integer(kind=8) :: k, i, tmp

         m(1:s) = m_in

         do k=s+1, N_M
           tmp=ieor(2**s * m(k-s), m(k-s))
           do i=1,s-1
             tmp = ieor(m(k-i) * 2**i * ai(a, s-i), tmp)
           end do
           m(k) = tmp
         end do

         do k=1, N_M
           state%v(k) = (2**(N_M-k))*m(k)
         end do

         state%x = 0
         state%i = 0
         state%stride = 0
         if (present(stride)) state%stride = stride
       end subroutine initialize


       !> Generate a value at a specific position i
       function skip_ahead(state, i) result(output)
         implicit none
         class (sobol_state), intent(inout) :: state
         integer(kind=8), intent(in) :: i
         real(dp)    :: output
         integer(kind=8) :: g ! Gray code representation of i
         integer(kind=8) :: j, tmp

         g = ieor(i,i/2)
         state%x = 0
         state%i = i

         tmp = ai(g,1_8) * state%v(1)
         do j=2, N_M
           tmp = ieor(tmp,ai(g,j) * state%v(j))
         end do
         output = real(tmp) * 2.d0**(-N_M)
         state%x = tmp

       end function skip_ahead


       !> Generate the next value in a series
       !> And update the state function
       function next(state)
         implicit none
         class (sobol_state), intent(inout) :: state
         real(dp) :: next
         state%x = ieor(state%x, state%v(i4_bit_lo0(state%i)))
         state%i = state%i + 1
         next = real(state%x) * 2.d0**(-N_M)
       end function next

       !> Generate the next value in a series
       !> And update the state function
       function next_strided(state)
         implicit none
         class (sobol_state), intent(inout) :: state
         real(dp) :: next_strided
         state%x = ieor(state%x, ieor(state%v(state%stride), state%v(
     &             i4_bit_lo0(ior(state%i, 2**state%stride - 1)))))
         state%i = state%i + 2**state%stride
         next_strided = real(state%x) * 2.d0**(-N_M)
       end function next_strided


       !> Returns the value of the bit at position i (1-based index)
       function ai(a,i)
       implicit none
       integer(kind=8), intent(in) :: a, i
       integer(kind=8) :: ai

       if (btest(a,i-1)) then
         ai = 1
       else
         ai = 0
       end if
       end function ai

       !> Return the position of the lowest (rightmost) 0 bit in a int(4)
       function i4_bit_lo0(num)
         implicit none

         integer(kind=8), intent(in) :: num
         integer(kind=8) :: i4_bit_lo0

         do i4_bit_lo0=1,bit_size(num)
           if (.not. btest(num,i4_bit_lo0-1)) return
         enddo
       end function i4_bit_lo0

       end module mod_sobseq
