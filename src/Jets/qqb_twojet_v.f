!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_twojet_v(p,msq)
c*******************************************************************************
c                                                                              *
c   Author: J. Campbell, September 2014                                        *
c   Virtual matrix elements for two-jet production                             *
c                                                                              *
c*******************************************************************************
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'blha.f'
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  QjQk_QjQk,AjAk_AjAk,AjQk_AjQk,QjAk_QjAk,
     &  QjQj_QjQj,AjAj_AjAj,AjQj_AjQj,QjAj_QjAj,
     &  AjQj_AkQk,QjAj_QkAk,
     &  GQ_GQ,GA_GA,QG_QG,AG_AG,GG_GG,GG_QA,QA_GG,AQ_GG,ss,tt,uu,
     &  virta,virtb,virtc,virtd

      scheme='tH-V'

      call dotem(4,p,s)
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      fac=gsq**2*ason2pi

c--- four quarks
      if (useblha == 0 .or. (blhatype == 4)) then
      QjQk_QjQk=fac*aveqq*virta(ss,tt,uu)
      AjAk_AjAk=QjQk_QjQk
      endif

      if (useblha == 0) then
      QjAj_QkAk=fac*aveqq*virta(tt,ss,uu)
      QjAk_QjAk=fac*aveqq*virta(uu,tt,ss)
      AjQj_AkQk=QjAj_QkAk
      AjQk_AjQk=QjAk_QjAk
      endif

      if (useblha == 0 .or. (blhatype == 3)) then
      QjQj_QjQj=fac*half*aveqq*virtb(ss,tt,uu)
      AjAj_AjAj=QjQj_QjQj
      endif

      if (useblha == 0) then
      QjAj_QjAj=fac*aveqq*virtb(uu,tt,ss)
      AjQj_AjQj=QjAj_QjAj
      endif

c--- two quarks, two gluons
      if (useblha == 0 .or. (blhatype == 2)) then
      GG_QA=fac*avegg*virtc(ss,uu,tt)
      QA_GG=half*aveqq/avegg*GG_QA
      AQ_GG=half*aveqq/avegg*GG_QA
      endif

      if (useblha == 0) then
      QG_QG=-fac*aveqg*virtc(tt,ss,uu)
      AG_AG=-fac*aveqg*virtc(tt,uu,ss)
      GQ_GQ=QG_QG
      GA_GA=AG_AG
      endif

c--- four gluons
      if (useblha == 0 .or. (blhatype == 1)) then
      GG_GG=fac*avegg*virtd(ss,tt,uu)
      endif

c--- initalize matrix elements to zero
      msq(:,:)=0d0

      nflav=nf

      do j=-nflav,nflav
      do k=-nflav,nflav

c-- QQ
      if ((j  >  0) .and. (k  >  0)) then
          if (j  ==  k) then
            msq(j,k)=QjQj_QjQj
          else
            msq(j,k)=QjQk_QjQk
          endif

c-- QA
      elseif ((j  >  0) .and. (k  <  0)) then
          if (j  ==  -k) then
            msq(j,k)=QjAj_QjAj+dfloat(nflav-1)*QjAj_QkAk+QA_GG
            if (useblha == 1) then
              if (blhatype == 2) then
                msq(j,k)=QA_GG
              elseif (blhatype == 3) then
                msq(j,k)=QjAj_QjAj
              elseif (blhatype == 4) then
                msq(j,k)=QjAj_QkAk
              endif
              if (p(3,4)*p(4,4) < 0._dp) msq(j,k)=-msq(j,k)
            endif
          else
            msq(j,k)=QjAk_QjAk
          endif

c-- AA
      elseif ((j  <  0) .and. (k  <  0)) then
          if (j  ==  k) then
            msq(j,k)=AjAj_AjAj
          else
            msq(j,k)=AjAk_AjAk
          endif

c-- AQ
      elseif ((j  >  0) .and. (k  <  0)) then
          if (j  ==  -k) then
            msq(j,k)=AjQj_AjQj+dfloat(nflav-1)*AjQj_AkQk+AQ_GG
            if (blhatype == 2) then
              msq(j,k)=AQ_GG
            elseif (blhatype == 3) then
              msq(j,k)=AjQj_AjQj
            elseif (blhatype == 4) then
              msq(j,k)=AjQj_AkQk
            endif
          else
            msq(j,k)=AjQk_AjQk
          endif

c-- QG
      elseif ((j  >  0) .and. (k  ==  0)) then
            msq(j,k)=QG_QG

c-- AG
      elseif ((j  <  0) .and. (k  ==  0)) then
            msq(j,k)=AG_AG

c-- GQ
      elseif ((j  ==  0) .and. (k  >  0)) then
            msq(j,k)=GQ_GQ

c-- GA
      elseif ((j  ==  0) .and. (k  <  0)) then
            msq(j,k)=GA_GA

c-- GG
      elseif ((j  ==  0) .and. (k  ==  0)) then
            msq(j,k)=half*GG_GG+dfloat(nflav)*GG_QA
            if (blhatype == 1) then
              msq(j,k)=half*GG_GG
            else
              msq(j,k)=GG_QA
            endif

      endif

      enddo
      enddo

      return
      end
