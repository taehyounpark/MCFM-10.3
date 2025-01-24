!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use fillpenttobox_generic
      use hgggg_integralfill_generic
      use hgggg_pppp_generic
      use hgggg_pppm_generic
      use hgggg_pmpm_generic
      use hgggg_ppmm_generic
      use hgggg_assemble_generic
      use gghgg_dep_params
      implicit none
      include 'Inc/zprods_decl.f'
      include 'Inc/hgggglabels.f'
      include 'Inc/IntResults.f'
      include 'Inc/Cred.f'
      include 'Inc/I5to4.f'
      integer:: j,c1,c2,c3,c4,pppm,ppmm,h1,h2,h3,h4
      real(dp)::mtsq,mt,mt_yuk,as,vevsq
      real(dp):: msq
      complex(dp) :: Dcoeff(dmax),Ccoeff(cmax),Bcoeff(bmax),Rat,
     & amp(2,2,2,2,3)

      as=real(as_dp,kind=dp)
      vevsq=real(vevsq_dp,kind=dp)
      mt=real(mt_dp,kind=dp)
      mt_yuk=real(mt_yuk_dp,kind=dp)
      mtsq=mt**2

      call fillpenttobox(mtsq,Cred,I5to4)

c Loop over non-cyclic color orderings
      do j=1,3
        if (j == 1) then
          c1=1; c2=2; c3=3; c4=4
        elseif (j == 2) then
          c1=1; c2=3; c3=4; c4=2
        elseif (j == 3) then
          c1=1; c2=4; c3=2; c4=3
        endif

c Permutation (1,2,3,4)
        call hgggg_integralfill(c1,c2,c3,c4,mtsq,Dint,Cint,Bint)
        call hgggg_pppp(c1,c2,c3,c4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(2,2,2,2,j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(1,1,1,1,j)=hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

        call hgggg_pppm(c1,c2,c3,c4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(pppm(1,c4),pppm(2,c4),pppm(3,c4),pppm(4,c4),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-pppm(1,c4),3-pppm(2,c4),3-pppm(3,c4),3-pppm(4,c4),j)=
     &   hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

        call hgggg_pmpm(c1,c2,c3,c4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(ppmm(1,c2,c4),ppmm(2,c2,c4),ppmm(3,c2,c4),ppmm(4,c2,c4),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-ppmm(1,c2,c4),3-ppmm(2,c2,c4),3-ppmm(3,c2,c4),3-ppmm(4,c2,c4),j)=
     &   hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

        call hgggg_ppmm(c1,c2,c3,c4,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(ppmm(1,c3,c4),ppmm(2,c3,c4),ppmm(3,c3,c4),ppmm(4,c3,c4),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-ppmm(1,c3,c4),3-ppmm(2,c3,c4),3-ppmm(3,c3,c4),3-ppmm(4,c3,c4),j)=
     &   hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

c Cyclic permutation (4,1,2,3)
        call hgggg_integralfill(c4,c1,c2,c3,mtsq,Dint,Cint,Bint)
        call hgggg_pppm(c4,c1,c2,c3,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(pppm(1,c3),pppm(2,c3),pppm(3,c3),pppm(4,c3),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-pppm(1,c3),3-pppm(2,c3),3-pppm(3,c3),3-pppm(4,c3),j)=
     &   hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

        call hgggg_ppmm(c4,c1,c2,c3,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(ppmm(1,c2,c3),ppmm(2,c2,c3),ppmm(3,c2,c3),ppmm(4,c2,c3),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-ppmm(1,c2,c3),3-ppmm(2,c2,c3),3-ppmm(3,c2,c3),3-ppmm(4,c2,c3),j)=
     &   hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

c Cyclic permutation (3,4,1,2)
        call hgggg_integralfill(c3,c4,c1,c2,mtsq,Dint,Cint,Bint)
        call hgggg_pppm(c3,c4,c1,c2,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(pppm(1,c2),pppm(2,c2),pppm(3,c2),pppm(4,c2),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-pppm(1,c2),3-pppm(2,c2),3-pppm(3,c2),3-pppm(4,c2),j)
     &   =hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)

c Cyclic permutation (2,3,4,1)
        call hgggg_integralfill(c2,c3,c4,c1,mtsq,Dint,Cint,Bint)
        call hgggg_pppm(c2,c3,c4,c1,mtsq,za,zb,Dcoeff,Ccoeff,Bcoeff,Rat,Cred,I5to4)
        amp(pppm(1,c1),pppm(2,c1),pppm(3,c1),pppm(4,c1),j)=hgggg_assemble(Dcoeff,Ccoeff,Bcoeff,Rat,Dint,Cint,Bint)
        amp(3-pppm(1,c1),3-pppm(2,c1),3-pppm(3,c1),3-pppm(4,c1),j)
     &   =hgggg_assemble(conjg(Dcoeff),conjg(Ccoeff),conjg(Bcoeff),conjg(Rat),Dint,Cint,Bint)
      enddo

c Add helicities
      msq=zero
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
c DEBUG: removed subleading color contribution
c      msq=msq
c     & +two*xn**2*(abs(amp(h1,h2,h3,h4,1))**2+abs(amp(h1,h2,h3,h4,2))**2+abs(amp(h1,h2,h3,h4,3))**2)
c Full result
      msq=msq
     & +two*xn**2*(abs(amp(h1,h2,h3,h4,1))**2+abs(amp(h1,h2,h3,h4,2))**2+abs(amp(h1,h2,h3,h4,3))**2)
     & -four*(xn**2-three)/xn**2*abs(amp(h1,h2,h3,h4,1)+amp(h1,h2,h3,h4,2)+amp(h1,h2,h3,h4,3))**2
c      write(6,*) 'h1,h2,h3,h4,abs(amp(h1,h2,h3,h4,1))',h1,h2,h3,h4,abs(amp(h1,h2,h3,h4,3))
      enddo
      enddo
      enddo
      enddo
      !write(6,*) '|H1234 ++++|',abs(amp(2,2,2,2,1))
      !write(6,*) '|H1234 +++-|',abs(amp(2,2,2,1,1))
      !write(6,*) '|H1234 +-+-|',abs(amp(2,1,2,1,1))
      !write(6,*) '|H1234 ++--|',abs(amp(2,2,1,1,1))

c Apply overall factor
      msq=msq*V*(as**2*mt*mt_yuk)**2/vevsq

      return


