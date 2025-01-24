!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_2jnogg(p,msq)
      implicit none
      include 'types.f'
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(ss,uu,tt)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)



      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)


      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0._dp


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then
            msq(j,k)=qaii_ii+real(nf-1,dp)*qaii_jj+qa_gg
          else
            msq(j,k)=qaij_ij
          endif

c--aa
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=aqii_ii+real(nf-1,dp)*aqii_jj+aq_gg
          else
            msq(j,k)=aqij_ij
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_qg
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_ag
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gq
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_ga
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end

      subroutine qqb_2jnoggswap(pin,msq)
      implicit none
      include 'types.f'
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: pin(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc

      p=0._dp

      do j=1,4
         p(1,j)=pin(1,j)
         p(2,j)=pin(2,j)
         p(4,j)=pin(3,j)
         p(3,j)=pin(4,j)
      enddo



      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(ss,uu,tt)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half


      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)


      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0._dp


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then
            msq(j,k)=qaij_ij+0._dp*(real(nf-1,dp)*qaii_jj+qa_gg)
          else
            msq(j,k)=qaij_ij
          endif

c--aa
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=aqij_ij+0._dp*(real(nf-1,dp)*aqii_jj+aq_gg)
          else
            msq(j,k)=aqij_ij
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_qg
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_ag
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gq
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_ga
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end


      subroutine qqb_2j_t(p,msq)
      implicit none
      include 'types.f'
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
c----- This routine contains t-channel type diagrams which will contain an intial-final
c---- photon singularity other terms are set to zero

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,
     &  qqii_ii,aaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,ss,tt,uu,
     &  smalla,smallb,smallc
c    & ,gq_qg,qaii_jj,aqii_ii,qaii_ii,
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)


      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
c      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
c      aqij_ij=fac*aveqq*smalla(uu,tt,ss)
      aqij_ij=fac*aveqq*smalla(ss,uu,tt)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
c      aqii_ii=fac*aveqq*smallb(tt,uu,ss)
c      qaii_ii=fac*aveqq*smallb(uu,tt,ss)



      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
c      gq_qg=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

c      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)


      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0._dp


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then

            msq(j,k)=qaij_ij !--- want only t-channel scattering
          else
            msq(j,k)=qaij_ij
          endif

c--aa
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=aqij_ij+0._dp*(real(nf-1,dp)*aqii_jj+aq_gg) !-- want only t-channel scattering
          else
            msq(j,k)=aqij_ij
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_qg
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_ag
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gq
 !        msq(j,k)=gq_qg
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_ga
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end



      subroutine qqb_2j_s(p,msq)
      implicit none
      include 'types.f'
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)



      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
c      gq_qg=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)


      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0._dp


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then
            msq(j,k)=0._dp*qaii_ii+qaii_jj+0._dp*qa_gg  !-- want 1 s-channel scattering
          else
            msq(j,k)=qaij_ij*0._dp
          endif

c--aa
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=0._dp*aqii_ii+aqii_jj+0._dp*aq_gg !-- want 1 s-channel scattering
          else
            msq(j,k)=aqij_ij*0._dp
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_qg
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_ag
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gq
 !           pause
 !        msq(j,k)=gq_qg
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_ga
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end


      subroutine qqb_2j_sswap(pin,msq)
      implicit none
      include 'types.f'
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: pin(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc

      p=0._dp

      do j=1,4
         p(1,j)=pin(1,j)
         p(2,j)=pin(2,j)
         p(4,j)=pin(3,j)
         p(3,j)=pin(4,j)
      enddo




      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)



      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(ss,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)



      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)


      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0._dp


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then
            msq(j,k)=0._dp*qaii_ii+qaii_jj+0._dp*qa_gg  !-- want 1 s-channel scattering
          else
            msq(j,k)=qaij_ij*0._dp
          endif

c--aa
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=0._dp*aqii_ii+aqii_jj+0._dp*aq_gg !-- want 1 s-channel scattering
          else
            msq(j,k)=aqij_ij*0._dp
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_qg
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_ag
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gq
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_ga
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end


