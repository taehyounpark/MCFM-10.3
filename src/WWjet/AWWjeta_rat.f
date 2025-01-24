!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjeta_rat(j1,j2,j5,j6,j4,j3,j7,za,zb,coeff,Alo)
c---- Fills appropriate entries in coeff with the rational parts
c---- present at leading and subleading color

c--- MCFM notation
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c--- correspoding to Tania notation
c---   u(6) ubar(1) nu(3) e+(2) e-(5) nubar(4) g(7)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7,i
      complex(dp):: Alo,rat(20),rat_lc,rat_slc

c--- compute all rational contributions
      call rat_t1(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(1))
      call rat_t2(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(2))
      call rat_t4(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(4))
      call rat_t5(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(5))
      call rat_t6(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(6))
      call rat_t7(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(7))
      call rat_t8(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(8))
      call rat_t9(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(9))
      call rat_t11(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(11))
      call rat_t12(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(12))
      call rat_t13(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(13))
      call rat_t14(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(14))
      call rat_t15(j1,j2,j3,j4,j5,j6,j7,za,zb,rat(15))
      rat(3)=czip
      rat(10)=czip
c      do j=1,15
c      write(6,*) j,rat(j)
c      enddo
c      pause
      rat_slc=czip
      do i=1,15
      rat_slc=rat_slc+rat(i)
      enddo
      rat_lc=-rat(11)+rat(13)+rat(14)+rat(4)+rat(5)
     &       -rat( 1)+rat( 6)+rat(15)-rat(2)+rat(7)

c      write(6,*) 'rat_lc',rat_lc
c      write(6,*) 'rat_slc',rat_slc
c      pause
c--- relationship between mu^2 terms and rational pieces
      rat_slc=rat_slc*(-0.5_dp)
      rat_lc=rat_lc*(-0.5_dp)

      coeff(0,irat)=rat_lc
      coeff(0,iratsl)=rat_slc

c      write(6,*) 'rat_slc/Alo',rat_slc/Alo
c      write(6,*) 'rat_lc/Alo',rat_lc/Alo
c      pause

c      write(6,*) 'old rat13',rat(13)
c      call rat_t13(j1,j2,j3,j4,j5,j6,j7,rat(13))
c      write(6,*) 'new rat13',rat(13)
c      write(6,*)
c      write(6,*) 'old rat15',rat(15)
c      call rat_t15(j1,j2,j3,j4,j5,j6,j7,rat(15))
c      write(6,*) 'new rat15',rat(15)
c      pause

c--- debugging output
      if (1  ==  2) then
        write(6,*) 'rat_t1    ',rat(1)    ,cdabs(rat(1)    )
c        write(6,*) 'rat_t1/Alo',rat(1)/Alo,cdabs(rat(1)/Alo)
c        write(6,*)
        write(6,*) 'rat_t2    ',rat(2)    ,cdabs(rat(2)    )
c        write(6,*) 'rat_t2/Alo',rat(2)/Alo,cdabs(rat(2)/Alo)
c        write(6,*)
        write(6,*) 'rat_t4    ',rat(4)    ,cdabs(rat(4)    )
c        write(6,*) 'rat_t4/Alo',rat(4)/Alo,cdabs(rat(4)/Alo)
c        write(6,*)
        write(6,*) 'rat_t5    ',rat(5)    ,cdabs(rat(5)    )
c        write(6,*) 'rat_t5/Alo',rat(5)/Alo,cdabs(rat(5)/Alo)
c        write(6,*)
        write(6,*) 'rat_t6    ',rat(6)    ,cdabs(rat(6)    )
c        write(6,*) 'rat_t6/Alo',rat(6)/Alo,cdabs(rat(6)/Alo)
c        write(6,*)
        write(6,*) 'rat_t7    ',rat(7)    ,cdabs(rat(7)    )
c        write(6,*) 'rat_t7/Alo',rat(7)/Alo,cdabs(rat(7)/Alo)
c        write(6,*)
        write(6,*) 'rat_t8    ',rat(8)    ,cdabs(rat(8)    )
c        write(6,*) 'rat_t8/Alo',rat(8)/Alo,cdabs(rat(8)/Alo)
c        write(6,*)
        write(6,*) 'rat_t9    ',rat(9)    ,cdabs(rat(9)    )
c        write(6,*) 'rat_t9/Alo',rat(9)/Alo,cdabs(rat(9)/Alo)
c        write(6,*)
        write(6,*) 'rat_t11    ',rat(11)    ,cdabs(rat(11)    )
c        write(6,*) 'rat_t11/Alo',rat(11)/Alo,cdabs(rat(11)/Alo)
c        write(6,*)
        write(6,*) 'rat_t12    ',rat(12)    ,cdabs(rat(12)    )
c        write(6,*) 'rat_t12/Alo',rat(12)/Alo,cdabs(rat(12)/Alo)
c        write(6,*)
        write(6,*) 'rat_t13    ',rat(13)    ,cdabs(rat(13)    )
c        write(6,*) 'rat_t13/Alo',rat(13)/Alo,cdabs(rat(13)/Alo)
c        write(6,*)
        write(6,*) 'rat_t14    ',rat(14)    ,cdabs(rat(14)    )
c        write(6,*) 'rat_t14/Alo',rat(14)/Alo,cdabs(rat(14)/Alo)
c        write(6,*)
        write(6,*) 'rat_t15    ',rat(15)    ,cdabs(rat(15)    )
c        write(6,*) 'rat_t15/Alo',rat(15)/Alo,cdabs(rat(15)/Alo)
        write(6,*)
      endif

      return
      end

