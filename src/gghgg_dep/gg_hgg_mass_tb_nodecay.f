!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hgg_mass_tb_nodecay(p,iglue1,iglue2,msq)
      use double_precision
      use sprod_dp
      use spinor
      use hgggg_mass_tb_generic
      use haqgg_mass_tb_generic
      use haqaq_mass_tb_generic
      use gghgg_dep_params
      implicit none

c---Matrix element squared averaged over initial colors and spins

c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6)

c   Including the effect of the finite top-quark and bottom-quark masses

      include 'masses.f'
      include 'zprods_com.f'
      include 'nflav.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'yukawas.f'
      include 'blha.f'

      real(dp), intent(in) :: p(mxpart,4)
      integer, intent(in) :: iglue1, iglue2
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

      integer:: j,k
      real(dp):: pglue(mxpart,4)
      real(dp):: Hgggg,Hqagg,Haqgg,Hgqqg,Hgaag,Hqgqg,Hagag,Hggqa,
     & Hqrqr,Hqqqq,Habab,Haaaa,Hbqbq,
     & Hqarb,Hqaqa,Hqbqb,Haqbr,Haqaq

c--- Set up momenta for partons
      pglue(1,:)=p(1,:)
      pglue(2,:)=p(2,:)
      pglue(3,:)=p(iglue1,:)
      pglue(4,:)=p(iglue2,:)

c---fill spinor products
      call spinoru(4,pglue,za,zb)

c--- load double precision values of parameters
      as_dp=as
      vevsq_dp=vevsq
      mt_dp=mt
      mt_yuk_dp=mt_yuk
      mb_dp=mb
      mb_yuk_dp=mb_yuk

c--- four gluon terms
      if ((useblha==0).or.(blhatype==1)) then
      call hgggg_mass_tb(za,zb,Hgggg)
      endif

c--- four quark terms
      if ((useblha==0).or.(blhatype>=7)) then
      call haqaq_mass_tb(1,3,2,4,za,zb,Hqrqr,Hqqqq)
      Habab=Hqrqr
      Haaaa=Hqqqq
      Hqbqb=Hqrqr
      Hbqbq=Hqrqr
      endif

      if (useblha==0) then
      call haqaq_mass_tb(1,2,4,3,za,zb,Hqarb,Hqaqa)
      Haqbr=Hqarb
      Haqaq=Hqaqa
      endif

c--- two quark two gluon terms
      if (useblha==0) then
      call haqgg_mass_tb(1,2,3,4,za,zb,Hqagg)
      Haqgg=Hqagg
      endif

      if (useblha==0) then
      call haqgg_mass_tb(1,3,2,4,za,zb,Hqgqg)
      Hagag=Hqgqg
      endif

      if (useblha==0) then
      call haqgg_mass_tb(2,3,1,4,za,zb,Hgqqg)
      Hgaag=Hgqqg
      endif

      if ((useblha==0).or.(blhatype==2)) then
      call haqgg_mass_tb(4,3,1,2,za,zb,Hggqa)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if ((j>0).and.(k>0)) then
        if (j==k) then
          msq(j,k)=half*aveqq*Hqqqq
        else
          msq(j,k)=aveqq*Hqrqr
        endif
      endif

      if ((j<0).and.(k<0)) then
        if (j==k) then
          msq(j,k)=half*aveqq*Haaaa
        else
          msq(j,k)=aveqq*Habab
        endif
      endif

      if ((j>0).and.(k<0)) then
        if (j==-k) then
          msq(j,k)=aveqq*(half*Hqagg+Hqaqa+real(nflav-1,dp)*Hqarb)
          if (useblha==1) then
             if ((blhafl(5)==0).and.(blhafl(6)==0)) then
                msq(j,k)=aveqq*half*Hqagg
             elseif ((blhafl(5)==blhafl(1)).and.
     &               (blhafl(6)==blhafl(2))) then
                msq(j,k)=aveqq*Hqaqa
             else
                msq(j,k)=aveqq*Hqarb
             endif
          endif
        else
          msq(j,k)=aveqq*Hqbqb
        endif
      endif

      if ((j<0).and.(k>0)) then
        if (j==-k) then
          msq(j,k)=aveqq*(half*Haqgg+Haqaq+real(nflav-1,dp)*Haqbr)
          if (useblha==1) then
             if ((blhafl(5)==0).and.(blhafl(6)==0)) then
                msq(j,k)=aveqq*half*Haqgg
             elseif ((blhafl(5)==blhafl(1)).and.
     &               (blhafl(6)==blhafl(2))) then
                msq(j,k)=aveqq*Haqaq
             else
                msq(j,k)=aveqq*Haqbr
             endif
          endif
        else
          msq(j,k)=aveqq*Hbqbq
        endif
      endif

      if ((j>0).and.(k==0)) then
        msq(j,0)=aveqg*Hqgqg
      endif

      if ((j<0).and.(k==0)) then
        msq(j,0)=aveqg*Hagag
      endif

      if ((j==0).and.(k>0)) then
        msq(0,k)=aveqg*Hgqqg
      endif

      if ((j==0).and.(k<0)) then
        msq(0,k)=aveqg*Hgaag
      endif

      if ((j==0).and.(k==0)) then
        msq(0,0)=avegg*(half*Hgggg+real(nflav,dp)*Hggqa)
        if (useblha==1) then
           if ((blhafl(5)==0).and.(blhafl(6)==0)) then
              msq(j,k)=avegg*half*Hgggg
           else
              msq(j,k)=avegg*Hggqa
           endif
        endif
      endif

      enddo
      enddo

      return
      end


