!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function virt5_VH(ip,za,zb,includeaxial)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: virt5_VH
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1999.                                                      *
c   Given za and zb calculate the                                      *
c   the interference of the amplitude for the process                  *
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L/R(5)                         *
c   at one loop with the corresponding lowest order amplitude          *
c   summed over the polarizations of the emitted gluon                 *
c   Virtual terms are in units of
c   (as/4/pi) (4 pi)^ep Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)
c
c   Vector boson coupling to closed fermion loop through an axial
c   coupling only included if includeaxial = .true.
c   (required for ZH processes, absent for WH)
c
c***********************************************************************
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      complex(dp):: virt5ax
      common/virt5ax/virt5ax
      integer:: ip(5)
      logical includeaxial
      complex(dp):: A5LOm,A5NLOm,A5LOp,A5NLOp,A5axp,A5axm
!$omp threadprivate(/virt5ax/)

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLO_VH(ip(1),ip(2),ip(3),ip(4),ip(5),za,zb,includeaxial,A5LOm,A5NLOm,A5axm)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLO_VH(ip(2),ip(1),ip(4),ip(3),ip(5),zb,za,includeaxial,A5LOp,A5NLOp,A5axp)

      virt5ax=czip
      if (includeaxial) then
      virt5ax=
     & +ason2pi*(conjg(A5LOp)*A5axp+conjg(A5LOm)*A5axm)
      endif

      virt5_VH=
     & +ason2pi*(real(conjg(A5LOp)*A5NLOp,dp)
     &          +real(conjg(A5LOm)*A5NLOm,dp))

      return
      end

