!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqbgen(e1,p2,p3,e4,za,zb,zab,zba,ppmp,amp)
      implicit none
      include 'types.f'
c-----Authors: John Campbell and Keith Ellis, November 2011
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      complex(dp):: amp(2),ppmp
      integer:: e1,p2,p3,e4,e1p,e4p
c---we are calculating the amplitudes mmpp and mpmp
c---from our base amplitudes which are ppmp

c--- perform the swaps 5 <-> 7 and 6 <-> 8
      e1p=12-e1
      e4p=14-e4

c      write(6,*) 'e1',e1
c      write(6,*) 'e4',e4
c---mpmp
      amp(1)=+ppmp(e1p,p2,p3,e4,za,zb,zab,zba)
c      amp(2)=ppmp(e1p,p3,p2,e4,za,zb,zab,zba)
c      write(6,*) 'amp(1)',amp(1),abs(amp(1))
c      write(6,*) 'amp(2)',amp(2),abs(amp(2))
c---mmpp
c--- JC: added minus sign here based on form of tree level amplitude
c---     (complex conjugation = (-1) * [za <-> zb])
      amp(2)=-ppmp(e1,p2,p3,e4p,zb,za,zba,zab)
c      write(6,*) 'amp(2)',amp(2),abs(amp(2))
c      write(6,*) 'rat',abs(amp(1))/abs(amp(2))
c      pause
      return
      end
