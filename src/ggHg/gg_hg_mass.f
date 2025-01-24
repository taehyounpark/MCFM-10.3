!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gg_hg_mass(p,ig,msq)
c Note: added from original MCFM routine named (confusingly) qqb_higgs.f
c Notable changes:
c     (1) Gluon momentum is in position ig here
c     (2) Automatic sum over t, b and c loops
c     (3) Conversion to complex mass scheme
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> H + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)

c--all momenta incoming

c--- Matrix elements are taken from:
c--- R.~K.~Ellis, I.~Hinchliffe, M.~Soldate and J.~J.~van der Bij,
c--- %``Higgs Decay To Tau+ Tau-: A Possible Signature Of Intermediate
c--- % Mass Higgs Bosons At The SSC,''
c--- Nucl.\ Phys.\ B {\bf 297}, 221 (1988).
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'hdecaymode.f'
      include 'sprods_decl.f'
      include 'blha.f'
      include 'yukawas.f'
      integer:: ig,j,k,nmin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      real(dp):: ehsvm3_tbc,ehsvm4_tbc,origmbsq,s34
      real(dp):: msqgamgam
      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      msq(:,:)=zip

      call dotem(ig,p,s)

      if (hdecaymode == 'none') then
         hdecay = 1._dp
      else
c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
          write(6,*) 'Unimplemented process in gg_hg_mass'
          stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      endif

      nmin=4
      if (mc_yuk == 0) nmin=5
      if (mb_yuk == 0) nmin=6

      origmbsq=mbsq
      if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype == 1))) then
      gg=+avegg*ehsvm3_tbc(nmin,s(1,2),s(1,ig),s(2,ig))*hdecay
      endif
      if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype == 2))) then
      qq=+aveqq*ehsvm4_tbc(nmin,s(1,2),s(1,ig),s(2,ig))*hdecay
      endif
      if ((useblha == 0) .or. ((useblha == 1) .and. (blhatype > 2))) then
      qg=-aveqg*ehsvm4_tbc(nmin,s(1,ig),s(1,2),s(2,ig))*hdecay
      gq=-aveqg*ehsvm4_tbc(nmin,s(2,ig),s(1,ig),s(1,2))*hdecay
      endif
      mbsq=origmbsq

      do j=-nf,nf
      do k=-nf,nf

      if ((j== 0) .or. (k==0)) then
           if ((j== 0) .and. (k==0)) then
                msq(j,k)=gg
           elseif ((j==0).and.(k /= 0)) then
                msq(j,k)=gq
           elseif ((j /= 0).and.(k==0)) then
                msq(j,k)=qg
           endif
      elseif ((j==-k).and. (j /= 0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end


      function ehsvm3_tbc(nmin,s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ehsvm3_tbc
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'yukawas.f'
      integer n,nmin
c---Matrix element squared Eqn 2.2 of EHSV
      complex(dp):: ehsva2,ehsva4,ampa2_stu,ampa2_ust,ampa2_tus,ampa4_stu
      real(dp):: hmass2,s,t,u,yuk_rescale
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      hmass2=s+t+u
      if (approx) then
        ehsvm3_tbc=gwsq/pi*as**3*xn*V/9._dp*(
     &          hmass2**4+s**4+t**4+u**4)/s/t/u/wmass**2
      else
        ampa2_stu=czip
        ampa2_ust=czip
        ampa2_tus=czip
        ampa4_stu=czip
        do n=nmin,6
          if (n == 4) then
            mbsq=mc**2
            yuk_rescale=mc_yuk/mc
          endif
          if (n == 5) then
            mbsq=mb**2
            yuk_rescale=mb_yuk/mb
          endif
          if (n == 6) then
            mbsq=mt**2
            yuk_rescale=mt_yuk/mt
          endif
          ampa2_stu=ampa2_stu+ehsva2(s,t,u)*yuk_rescale
          ampa2_ust=ampa2_ust+ehsva2(u,s,t)*yuk_rescale
          ampa2_tus=ampa2_tus+ehsva2(t,u,s)*yuk_rescale
          ampa4_stu=ampa4_stu+ehsva4(s,t,u)*yuk_rescale
        enddo
c        ehsvm3_tbc=
c     &   abs(ehsva2(s,t,u))**2+abs(ehsva2(u,s,t))**2+abs(ehsva2(t,u,s))**2
c     &   +abs(ehsva4(s,t,u))**2
        ehsvm3_tbc=
     &   abs(ampa2_stu)**2+abs(ampa2_ust)**2+abs(ampa2_tus)**2+abs(ampa4_stu)**2
        ehsvm3_tbc=4._dp/vevsq/pi*as**3*xn*V*hmass2**4/(s*t*u)*ehsvm3_tbc
      endif

      return
      end


      function ehsvm4_tbc(nmin,s,t,u)
      implicit none
      include 'types.f'
      real(dp):: ehsvm4_tbc
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'yukawas.f'
      integer n,nmin
c---Matrix element squared Eqn 2.6 of EHSV
      complex(dp):: ehsva5,amp
      real(dp):: hmass2,s,t,u,yuk_rescale
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      hmass2=s+t+u
      if (approx) then
        ehsvm4_tbc=gwsq/pi*as**3*V/18._dp*(u**2+t**2)/s/wmass**2
      else
        amp=czip
        do n=nmin,6
          if (n == 4) then
            mbsq=mc**2
            yuk_rescale=mc_yuk/mc
          endif
          if (n == 5) then
            mbsq=mb**2
            yuk_rescale=mb_yuk/mb
          endif
          if (n == 6) then
            mbsq=mt**2
            yuk_rescale=mt_yuk/mt
          endif
          amp=amp+ehsva5(s,t,u)*yuk_rescale
        enddo
        ehsvm4_tbc=abs(amp)**2
        ehsvm4_tbc=1._dp/vevsq/pi*as**3*V/2._dp*(u**2+t**2)/(s)
     &   *hmass2**2/(u+t)**2*ehsvm4_tbc
      endif

      return
      end



