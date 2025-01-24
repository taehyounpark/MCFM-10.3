!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c ---------------------------------------------------------------------- C
c     Born for quark induced dijets processes. It differentiates         C
c     from channel to channel, namely, each piece of the function        C
c     returns a matrix element of one channel interferes with another,   C
c     which is similar with the functions in "small.f".                  C
c ---------------------------------------------------------------------- C

      function sxs(ss,tt,uu)
c     function sxs denotes s-channel interferes with s-channel, and similar
c     notations are used in the following little functions
c     the s-channel process is defined as qa_ii_jj
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: ss,tt,uu,ssq,tsq,usq,sxs

      ssq = ss**2
      tsq = tt**2
      usq = uu**2

c      sxs = 2._dp*V*(tsq + usq)/ssq
c --- strip off the color structure
      sxs = + V*(tsq + usq)/ssq

      end function sxs



      function sxt(ss,tt,uu)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: ss,tt,uu,usq,sxt

      usq = uu**2

c      sxt = -2._dp*V/xn*usq/(ss*tt)
c --- strip off the color structure; served as mixed born instead of qcd born
c --- refer to routine 'dijet_qqb_born_mix' in 'dijet_qqb_mix.f'
c      sxt = V*usq/(ss*tt)
c      sxt =  - V*uu**2/ss/tt
c --- add '-' due to single fermion-loop in s- and t-channel interference
      sxt =  + V*uu**2/ss/tt

      end function sxt


c -- list the crossing symmetries at tree level
c     1. qa_ii_jj = sxs(ss,tt,uu) ~ (tsq + usq)/ssq
c     2. aq_ii_jj = sxs(ss,uu,tt) ~ (tsq + usq)/ssq
c     3. qq_ij_ij = sxs(tt,uu,ss) ~ (ssq + usq)/tsq
c     4. aa_ij_ij = sxs(tt,uu,ss) ~ (ssq + usq)/tsq
c     5. qa_ij_ij = sxs(tt,ss,uu) ~ (ssq + usq)/tsq
c     6. aq_ij_ij = sxs(tt,ss,uu) ~ (ssq + usq)/tsq
c     7. qa_ii_ii = sxs(ss,tt,uu) + sxs(tt,ss,uu) + 2*sxt(ss,tt,uu)
c                 ~ smallb(uu,ss,tt)
c     8. aq_ii_ii = sxs(ss,tt,uu) + sxs(tt,uu,ss) + 2*sxt(ss,tt,uu)
c                 ~ smallb(uu,ss,tt)
c     9. qq_ii_ii = sxs(tt,uu,ss) + sxs(uu,tt,ss) + 2*sxt(tt,uu,ss)
c                 ~ smallb(ss,tt,uu)
c     10.aa_ii_ii = sxs(tt,uu,ss) + sxs(uu,tt,ss) + 2*sxt(tt,uu,ss)
c                 ~ smallb(ss,tt,uu)

c     since we know the vertex correction is proportional to the Born,
c     and the correction is contained in the so-called form factor, we
c     can use the Born and the form factor information to get the vertex
c     correction.
