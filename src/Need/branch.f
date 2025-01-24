!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function branch_zqq()
          !! returns lowest order branching ratio for Y -> q qbar
          !! for quark flavors as in Inc/nf.f
          implicit none
          include 'types.f'
          include 'cplx.h'
          include 'constants.f'
          include 'nf.f'
          include 'masses.f'
          include 'ewcouple.f'
          include 'ewinput.f'
          include 'zcouple.f'
          include 'zcouple_cms.f'
          include 'zerowidth.f'

          real(dp) :: branch_zqq

          real(dp) :: facz,zwidth_cms,pwidth_d, pwidth_u
          complex(dp) :: zzmass2

          if (zerowidth) then
                  zwidth_cms = 0d0
          else
                  zwidth_cms = zwidth
          endif


          if (ewscheme  ==  4) then
                zzmass2=cplx2(zmass**2,-zmass*zwidth_cms)
                facz=abs(zesq/4._dp*zzmass2**0.5/(6._dp*pi))
                pwidth_d=n_light_down*facz*real(zL(1)**2+zR(1)**2)
                pwidth_u=n_light_up*facz*real(zL(2)**2+zR(2)**2)
          else
                facz=esq/4._dp*zmass/(6._dp*pi)
                pwidth_d=n_light_down*facz*(L(1)**2+R(1)**2)
                pwidth_u=n_light_up*facz*(L(2)**2+R(2)**2)
          endif
          branch_zqq = ca*(pwidth_d + pwidth_u) / zwidth
      end function

      subroutine branch(brwen,brzee,brznn,brtau,brtop,brcharm)
      implicit none
      include 'types.f'

c     Returns the lowest order branching ratios for
c     1) W   --> e nu
c     2) Z   --> e e
c     3) Z   --> (nu nubar) x 3
c     4) tau --> e nu nubar
c     5) t   --> b W
c     6) c   --> s W (with Vcs omitted)
      include 'constants.f'
      include 'nf.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewinput.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'zerowidth.f'

      real(dp):: facz,facw,factau,faccharm,
     & pwidth_e,pwidth_n,zwidth_cms,wwidth_cms
c      real(dp):: pwidth_u,pwidth_d,width
      real(dp):: brwen,brzee,brznn,brtau,brtop,brcharm
      complex(dp) :: zzmass2,zwmass2

      if (zerowidth) then
          zwidth_cms = 0d0
          wwidth_cms = 0d0
      else
          zwidth_cms = zwidth
          wwidth_cms = wwidth
      endif

      if (ewscheme == 4) then
              zzmass2=cplx2(zmass**2,-zmass*zwidth_cms)
              zwmass2=cplx2(wmass**2,-wmass*wwidth_cms)
              gwsq=abs(zesq/zxw)
              facz=abs(zesq/4._dp*zzmass2**0.5/(6._dp*pi))
              facw=abs(gwsq/8._dp*zwmass2**0.5/(6._dp*pi))
              factau=gwsq**2/32._dp/wmass**4*mtau**5/192._dp/pi**3
              pwidth_e=facz*(abs(zle)**2+abs(zre)**2)
              pwidth_n=facz*(abs(zln)**2)*3._dp
      else
        facz=esq/4._dp*zmass/(6._dp*pi)
        facw=gwsq/8._dp*wmass/(6._dp*pi)
        factau=gwsq**2/32._dp/wmass**4*mtau**5/192._dp/pi**3
        pwidth_e=facz*(le**2+re**2)
        pwidth_n=facz*(ln**2)*3._dp
      endif
c      xwsq=(wmass/mt)**2
c      xbsq=(mb/mt)**2
c      write(6,*) '(mb/mt)**2',xbsq
c      root=sqrt((1._dp+xbsq-xwsq)**2-4._dp*xbsq)
c      factop=(gw/wmass)**2*mt**3/(64._dp*pi)(1._dp-xwsq)**2*(1._dp+2._dp*xwsq)
c      factop=(gw/wmass)**2*mt**3/(64._dp*pi)*root
c     & *(1._dp+xwsq-2._dp*xwsq**2-2._dp*xbsq+xbsq*xwsq+xbsq**2)

      faccharm=(gwsq/wmass**2)**2/32._dp*mc**5/192._dp/pi**3
c      pwidth_d=3*facz*(L(1)**2+R(1)**2)
c      pwidth_u=3*facz*(L(2)**2+R(2)**2)
c calculated zwidth=3*pwidth_.e+2_dp*pwidth_u+3*pwidth_e+3*pwidth_n

      brzee=pwidth_e/zwidth
      brznn=pwidth_n/zwidth
      brwen=facw/wwidth
      brtau=factau/tauwidth
c      brtop=factop/twidth
      brcharm=faccharm

      brtop=1._dp

      return
      end
