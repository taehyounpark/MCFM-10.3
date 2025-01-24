!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ZWWcurr_ab(j1,jx,j2,jy,za,zb)
      implicit none
      include 'types.f'
c--- this function implements the replacement of a term
c--- of the form <j1 jx> [j2 jy] in the V+1 jet amplitudes
c--- with the expression corresponding to replacing the
c--- current giving |jx> |jy] with the Z->WW current

c--- note that jx,jy are actually unused here (except for
c--- returning the "usual" value for checking)

c--- the labels for the momenta in the Z->WW process are
c--- given the prefix "p" and are passed in via common block
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'srdiags.f'
      include 'kprocess.f'
      integer j1,jx,j2,jy
      integer p1,p2,p3,p4,p5,p6,p7
      complex(dp):: zab2,ZWWcurr_ab
      common/momWWZ/p1,p2,p3,p4,p5,p6,p7
!$omp threadprivate(/momWWZ/)

c--- statement function
      zab2(j1,jx,jy,j2)=za(j1,jx)*zb(jx,j2)+za(j1,jy)*zb(jy,j2)

c--- expression we are replacing: za(j1,jx)*zb(j2,jy)

      ZWWcurr_ab=za(j1,p3)*zb(j2,p4)*zab2(p5,p3,p4,p6)
     &          -za(j1,p5)*zb(j2,p6)*zab2(p3,p5,p6,p4)
     &          -za(p3,p5)*zb(p6,p4)
     &           *(zab2(j1,p5,p6,j2)-zab2(j1,p3,p4,j2))*0.5d0

c This fix is appropriate for WWjet only when not including all diagrams
      if ((srdiags .eqv. .false.) .or. (kcase==kWZ_jet)) then
        ZWWcurr_ab=ZWWcurr_ab
     &          +0.5d0*(zab2(j1,p3,p4,j2)+zab2(j1,p5,p6,j2))
     &          /(s(p3,p4)+s(p3,p5)+s(p3,p6)+s(p4,p5)+s(p4,p6)+s(p5,p6))
     &          *za(p3,p5)*zb(p6,p4)*(s(p5,p6)-s(p3,p4))
      endif

      ZWWcurr_ab=ZWWcurr_ab/(s(p3,p4)*s(p5,p6))

      return
      end

