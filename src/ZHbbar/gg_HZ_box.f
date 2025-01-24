!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c----- fills amplitudes for gg=>HZ process
c===== these are the pieces which come from box topologies
      subroutine gg_HZ_box(p,amp,mt2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),mt2
      complex(dp):: amp(2,2,2)
      complex(dp):: ggHZ_pp_box
      complex(dp):: ggHZ_mp_box
      external ggHZ_pp_box
      external ggHZ_mp_box
      complex(dp):: prop_34



c----- debug
c----- KC phase space point check
c      p(1,4)=  -3.0000000000000000_dp
c      p(1,1)=   2.1213203435596424_dp
c      p(1,2)=   1.0606601717798212_dp
c      p(1,3)=   1.8371173070873839_dp
c      p(2,4)=  -3.0000000000000000_dp
c      p(2,1)=  -1.3016512173526746_dp
c      p(2,2)=   2.6346884688170253_dp
c      p(2,3)= -0.60342421284441206_dp
c      p(3,4)=  0.85714285714285710_dp
c      p(3,1)= -0.31578947368421051_dp
c      p(3,2)=  0.79685060448070799_dp
c      p(3,3)=   0.0000000000000000_dp
c      p(4,4)=  -1.5000000000000000_dp
c      p(4,1)= -0.88167787843870971_dp
c      p(4,2)=   1.1895840590472970_dp
c      p(4,3)= -0.23986222114450001_dp


c      scale=one
c      musq=one
c      mt=0.4255266775_dp
c----- end debug

      call spinoru(4,p,za,zb)
c===== z propogator
      prop_34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)

c------- fill amplitudes which come from box diagrams

c------ left handed lepton coupling (- sign from line reversal)
      amp(2,2,1)=ggHZ_pp_box(1,2,3,4,za,zb,mt2)
      amp(1,2,1)=ggHZ_mp_box(1,2,3,4,za,zb,mt2)
      amp(2,1,1)=-ggHZ_mp_box(1,2,4,3,zb,za,mt2)
      amp(1,1,1)=-ggHZ_pp_box(1,2,4,3,zb,za,mt2)

c------- right handed lepton coupling
      amp(2,2,2)=ggHZ_pp_box(1,2,4,3,za,zb,mt2)
      amp(1,2,2)=ggHZ_mp_box(1,2,4,3,za,zb,mt2)
      amp(2,1,2)=-ggHZ_mp_box(1,2,3,4,zb,za,mt2)
      amp(1,1,2)=-ggHZ_pp_box(1,2,3,4,zb,za,mt2)


c----- debug, all checked against numeric D-dim code CW 2/10/15
c      write(6,*) '******** L ********* '
c      write(6,*) '2,2,1',amp(2,2,1)*im
c      write(6,*) '1,2,1',amp(1,2,1)*im
c      write(6,*) '2,1,1',amp(2,1,1)*im
c      write(6,*) '1,1,1',amp(1,1,1)*im
c      write(6,*) '******************** '
c      write(6,*)
c      write(6,*) '******** R ********* '
c      write(6,*) '2,2,2',amp(2,2,2)*im
c      write(6,*) '1,2,2',amp(1,2,2)*im
c      write(6,*) '2,1,2',amp(2,1,2)*im
c      write(6,*) '1,1,2',amp(1,1,2)*im
c      write(6,*) '******************** '
c      pause

c---- dress with propagators
      amp(:,:,:)=prop_34*amp(:,:,:)

      return
      end

