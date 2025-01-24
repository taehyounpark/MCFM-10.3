!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function cdlogwrap(x)
      implicit none
      include 'types.f'
      real(dp):: x
      complex(dp)::cdlogwrap
      cdlogwrap = cdlog(cmplx(x,0._dp,kind=dp))
      end function cdlogwrap

      function dilogc(x)
      implicit none
      include 'types.f'
      real(dp):: x
      complex(dp):: dilogc,cli2
      dilogc = cli2(cmplx(x,0d0,kind=dp))
      end function dilogc


c!! This function is not used except for debugging

c      function virtme_gq_new(sman, tman, uman, mtex)
c      implicit none
c      include 'types.f'
c      real(dp):: virtme_gq_new,sman, tman, uman
c      integer, intent(in) :: mtex
c      include 'constants.f'
c      include 'masses.f'
c      include 'scale.f'
c      include 'epinv.f'
c      include 'qcdcouple.f'
c      include 'ewcouple.f'
c      complex(dp):: cdlogwrap
c      complex(dp):: cli2
c      real(dp):: kgluon, kquark, nl, pcut, mH
c      complex(dp):: LogMums, LogMumt, LogMumu, LogMumH
c      real(dp):: LogMuMtop
c      complex(dp):: finite
c      complex(dp):: ASMh
c      complex(dp):: alomt0, alomt2, alomt4, alomt6
c      complex(dp):: amt0, amt2, amt4, amt6
c      complex(dp):: anlo
c      complex(dp):: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
c      complex(dp):: t11,t12,t13,t14,t15,t16,t17,t18
c      complex(dp):: t19,t20,t21,t22,t23,t24,t25,t26
c      complex(dp):: t27,t28,t29,t30,t31,t32,t33,t34
c      complex(dp):: t35,t36,t37,t38,t39,t40,t41,t42
c      complex(dp):: t43,t44,t45,t46,t47,t48,t49,t50
c      complex(dp):: t51,t52,t53,t54,t55,t56,t57,t58
c      complex(dp):: t59,t60,t61,t62,t63,t64,t65,t66
c      complex(dp):: t67,t68,t69

c      mH = hmass

c      nl = 5d0
c      pcut = 1d0
c      kgluon = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
c     &        - ca*dlog(pcut)**2
c     &        + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
c      kquark = cf*(7d0/2d0 - pisqo6)
c     &        - cf*dlog(pcut)**2
c     &        + 3d0/2d0*cf*(pcut - 1 - dlog(pcut))


c      LogMums = cdlogwrap(-musq/sman)
c      LogMumt = cdlogwrap(-musq/tman)
c      LogMumu = cdlogwrap(-musq/uman)
c      LogMumH = cdlogwrap(-musq/mH**2)

c      LogMuMtop = log(musq/mt**2)

c      finite = 0d0

c      alomt0 = 16/(3.*sman)
c      alomt2 = 0.8 + (14*tman)/(45.*sman) + (14*uman)/(45.*sman)
c      alomt4 =
c     &   (16*sman)/105. + (4*tman)/35. + (2*tman**2)/(63.*sman) +
c     &   (4*uman)/35. + (4*tman*uman)/(63.*sman) +
c     &   (2*uman**2)/(63.*sman)
c      alomt6 =
c     &   (2*sman**2)/63. + (11*sman*tman)/315. + (2*tman**2)/105. +
c     &   (13*tman**3)/(3150.*sman) + (11*sman*uman)/315. +
c     &   (4*tman*uman)/105. + (13*tman**2*uman)/(1050.*sman) +
c     &   (2*uman**2)/105. + (13*tman*uman**2)/(1050.*sman) +
c     &   (13*uman**3)/(3150.*sman)

c      anlo = 0d0

c       if (mtex >= 0) then
c         include 'src/Hjetmass/qag/a_mt0.f'
c         anlo = amt0
c       end if

c       if (mtex >= 2) then
c         include 'src/Hjetmass/qag/a_mt2.f'
c         anlo = anlo + amt2/mt**2
c       end if

c       if (mtex >= 4) then
c         include 'src/Hjetmass/qag/a_mt4.f'
c         anlo = anlo + amt4/mt**4
c       end if

c       if (mtex >= 6) then
c         include 'src/Hjetmass/qag/a_mt6.f'
c         anlo = anlo + amt6/mt**6
c       end if

c       finite = 0d0

c       finite = conjg(ASMh(sman,tman+uman,mt**2)*mt**2*alomt0) * anlo

c       finite = finite * sman*(tman**2 + uman**2) * 9d0/256d0

c      virtme_gq_new = 2d0*dreal(finite) * as/2/pi *
c     &           (as/3/pi)**2 / vevsq *
c     &          4*pi*as*4d0/(4d0*3d0*8d0)

c      end function virtme_gq_new

