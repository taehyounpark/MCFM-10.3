!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine amp_WWg_v(p,j1,j2,j3,j4,j5,j6,j7,Qid,amp0a,amp0b,amp0sr,amp1a,amp1b,amp1ax,amp1sr)
c--- master routine for assembling amplitude for
c---  ub(i1) + u(i2) -> nu(i3) + e+(i4) + e-(i5) + nu~(i6) + g(i7)
c---
c--- amplitude is returned as amp1a(h1,h2) and amp1b(h1,h2) with h1,h2
c--- the helicities of the quark and gluon respectively
c--- ampa refers to diagrams with both W's attached to the quark line
c--- ampb refers to diagrams with Z->WW
c--- routine also returns the corresponding LO amplidues amp0a,amp0b
c---
c--- Flavor of quark line is passed in as 'Qid'
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'srdiags.f'
      integer j1,j2,j3,j4,j5,j6,j7,Qid
      real(dp):: p(mxpart,4)

      complex(dp):: amp0a(2,2),amp0b(2,2),amp1a(2,2),amp1b(2,2),
     & Alom,Alop,A7vm,A7vp,amp1ax(2,2),amp0sr(2,2),amp1sr(2,2),
     & Alomp,Alomm,Alopp,Alopm,A7vmp,A7vmm,A7vpp,A7vpm

      amp0a(:,:)=czip
      amp0b(:,:)=czip
      amp1a(:,:)=czip
      amp1b(:,:)=czip

c--- a amplitude: left-handed quark line only
      call AWWjeta(p,j1,j2,j3,j4,j5,j6,j7,Alop,Alom,A7vp,A7vm)
      amp0a(1,2)=Alop
      amp0a(1,1)=Alom
      amp1a(1,2)=A7vp
      amp1a(1,1)=A7vm
c      write(6,*) 'A7vp_a,|A7p_a|',A7vp_a,cdabs(A7vp_a)
c      write(6,*) 'A7vm_a,|A7m_a|',A7vm_a,cdabs(A7vm_a)

c--- b amplitude: both left-handed and right-handed quarks
      call AWWjetb(j1,j2,j3,j4,j5,j6,j7,Alomp,Alomm,Alopp,Alopm,A7vmp,A7vmm,A7vpp,A7vpm)
      amp0b(1,2)=Alomp
      amp0b(1,1)=Alomm
      amp1b(1,2)=A7vmp
      amp1b(1,1)=A7vmm

      amp0b(2,2)=Alopp
      amp0b(2,1)=Alopm
      amp1b(2,2)=A7vpp
      amp1b(2,1)=A7vpm

      if (srdiags) then
c--- b amplitude: both left-handed and right-handed quarks
c note that this routine correctly handles u and d contributions without any j3,j4,j5,j6 permuation
c so we do some work to undo that
        if (((j3==3).and.(j4==4).and.(j5==5).and.(j6==6)) .or.
     &      ((j3==6).and.(j4==5).and.(j5==4).and.(j6==3))) then
          call AWWjetsr(j1,j2,3,4,5,6,j7,Qid,Alomp,Alomm,Alopp,Alopm,A7vmp,A7vmm,A7vpp,A7vpm)
          call AWWjetq_sr(j1,j2,3,4,5,6,j7,amp1ax)
        else
          write(6,*) 'Unexpected j3,j4,j5,j6 when calling AWWjetsr in amp_WWg_v.f:',j3,j4,j5,j6
          stop
        endif
        amp0sr(1,2)=Alomp
        amp0sr(1,1)=Alomm
        amp0sr(2,2)=Alopp
        amp0sr(2,1)=Alopm

        amp1sr(1,2)=A7vmp
        amp1sr(1,1)=A7vmm
        amp1sr(2,2)=A7vpp
        amp1sr(2,1)=A7vpm

c Add axial contribution since all propagators have already been applied
        amp1sr(:,:)=amp1sr(:,:)+amp1ax(:,:)

      else
        amp0sr(:,:)=czip
        amp1sr(:,:)=czip
      endif

c--- For efficiency these are now filled via a separate call to amp_wwg_v_qloop
c---
c---! DEBUG: qloop only
c---!      amp1a(:,:)=czip
c---!      amp1b(:,:)=czip
c---
c---c--- contribution of closed fermion loops; add to "a" since no further
c---c--- coupling factors are required, overall factor the same
c---c--- note that axial contribution is singled out since it should always
c---c--- enter with the same sign, and not depend on whether u or d is attached
c---c--- (changes sign here between 3456 and 6543, then flipped back in qqb_wwg_v)
c---!      call AWWjetq(p,j1,j2,j3,j4,j5,j6,j7,A7q,amp1ax)
c---!      amp1a(:,:)=amp1a(:,:)*0d0+A7q(:,:)

      return
      end

