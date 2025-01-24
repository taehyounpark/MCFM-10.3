!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c This routine is hard-wired for s > 0, t < 0, u < 0
      subroutine hjetfill(s,t,u,virtgg,virtqa,virtaq,virtqg,virtgq)
	use ggHwilson
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'b0.f'
      real(dp):: virtgg,virtqa,virtaq,virtqg,virtgq,
     & logg,loqa,loaq,loqg,logq,ddilog,Li2s,Li2t,Li2u,
     & lnm,lns,lnt,lnu,ln2t,ln2u,mhsq,s,t,u,xlf,subuv,Delta

      mhsq=s+t+u

      xlf=real(nf,dp)
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)
      lns=log(s/mhsq)
      lnt=log(-t/mhsq)
      lnu=log(-u/mhsq)
      lnm=log(musq/mhsq)
      ln2t=log((mhsq-t)/mhsq)
      ln2u=log((mhsq-u)/mhsq)

      logg=+V*xn*(mhsq**4+s**4+t**4+u**4)/(s*t*u)
      loqa=xn*cf/s*(t**2+u**2)
      loaq=loqa
      logq=-xn*cf/u*(s**2+t**2)
      loqg=-xn*cf/t*(s**2+u**2)

c--- UV counterterm in MS bar scheme.
      subuv=-epinv*b0
      if (expansionorder > 0) then
! we are handling Wilson coefficient expansion separately
        Delta=0._dp
      else
        Delta=11._dp
      endif
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
c--- Three powers of as in Born --> 3
      subuv=3._dp*subuv+Delta

      virtgg=-3._dp*epinv*epinv2*xn*logg
     & +epinv*xn*logg*(lns+lnt+lnu-3._dp*lnm )
     & +xn*logg
     & *(2._dp*(Li2t+Li2u+Li2s)
     & +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     & +0.5_dp*(lns**2-lnt**2-lnu**2)-1.5_dp*lnm**2
     & +2._dp*(lnu*ln2u+lnt*ln2t)+4._dp/3._dp*pisq)
     & +V*xn*(xn-xlf)/3._dp*mhsq*(1._dp+mhsq/s+mhsq/t+mhsq/u)
     & +subuv*logg

      virtqa=+(-2._dp*xn+1._dp/xn)*loqa*epinv*epinv2
     & -2._dp/3._dp*xlf*epinv*loqa
     & +epinv*xn*loqa*(13._dp/6._dp-2._dp*lnm+lnt+lnu)
     & +epinv/xn*loqa*(1.5_dp-lns+lnm)
     & +loqa*xlf*(-10._dp/9._dp+2._dp/3._dp*lns-2._dp/3._dp*lnm)
     & +xn*loqa* (40._dp/9._dp+Li2t+Li2u+2._dp*Li2s-13._dp/6._dp*(lns-lnm)
     & +(lnm-lns)*(lnt+lnu)+lns**2-lnm**2-0.5_dp*lnt**2-0.5_dp*lnu**2
     & +lnt*ln2t+lnu*ln2u)
     & +loqa/xn
     & *(4._dp-Li2t-Li2u-1.5_dp*(lns-lnm)+0.5_dp*(lns-lnm)**2
     & +lnt*lnu-lnt*ln2t-lnu*ln2u)
     & -4._dp/3._dp*pi**2/xn*loqa
     & -0.25_dp*(xn**3-1._dp/xn)*(t+u)
     & +subuv*loqa

      virtaq=(-2._dp*xn+1._dp/xn)*loaq*epinv*epinv2
     & -2._dp/3._dp*xlf*epinv*loaq
     & +epinv*xn*loaq*(13._dp/6._dp-2._dp*lnm+lnu+lnt)
     & +epinv/xn*loaq*(1.5_dp-lns+lnm)
     & +loaq*xlf*(-10._dp/9._dp+2._dp/3._dp*lns-2._dp/3._dp*lnm)
     & +xn*loaq* (40._dp/9._dp+Li2u+Li2t+2._dp*Li2s-13._dp/6._dp*(lns-lnm)
     & +(lnm-lns)*(lnu+lnt)+lns**2-lnm**2-0.5_dp*lnu**2-0.5_dp*lnt**2
     & +lnu*ln2u+lnt*ln2t)
     & +loaq/xn
     & *(4._dp-Li2u-Li2t-1.5_dp*(lns-lnm)+0.5_dp*(lns-lnm)**2
     & +lnu*lnt-lnu*ln2u-lnt*ln2t)
     & -4._dp/3._dp*pi**2/xn*loaq
     & -0.25_dp*(xn**3-1._dp/xn)*(u+t)
     & +subuv*loaq


      virtgq=(-2._dp*xn+1._dp/xn)*epinv*epinv2*logq
     & -2._dp/3._dp*xlf*epinv*logq
     & +epinv*xn*logq*(13._dp/6._dp+lns-2._dp*lnm+lnt)
     & +epinv/xn*logq*(3._dp/2._dp+lnm-lnu)
     & +logq*xlf*(-10._dp/9._dp-2._dp/3._dp*lnm+2._dp/3._dp*lnu)
     & +xn*logq*(40._dp/9._dp+Li2t+2._dp*Li2u+Li2s
     & +lns*lnm-lns*lnu-13._dp/6._dp*(lnu-lnm)
     & +lnm*lnt-lnm**2-lnt*lnu-0.5_dp*lnt**2
     & +2._dp*lnu*ln2u+lnt*ln2t)
     & +logq/xn*(4._dp-Li2t-Li2s+lns*lnt+0.5_dp*lnu**2-0.5_dp*lns**2
     & -lnm*lnu+0.5_dp*lnm**2-lnt*ln2t-1.5_dp*(lnu-lnm))
     & +4._dp/3._dp*pi**2*xn*logq
     & +0.25_dp*(xn**3-1._dp/xn)*(t+s)
     & +subuv*logq

      virtqg=(-2._dp*xn+1._dp/xn)*epinv*epinv2*loqg
     & -2._dp/3._dp*xlf*epinv*loqg
     & +epinv*xn*loqg*(13._dp/6._dp+lns-2._dp*lnm+lnu)
     & +epinv/xn*loqg*(3._dp/2._dp+lnm-lnt)
     & +loqg*xlf*(-10._dp/9._dp-2._dp/3._dp*lnm+2._dp/3._dp*lnt)
     & +xn*loqg*(40._dp/9._dp+Li2u+2._dp*Li2t+Li2s
     & +lns*lnm-lns*lnt-13._dp/6._dp*(lnt-lnm)
     & +lnm*lnu-lnm**2-lnu*lnt-0.5_dp*lnu**2
     & +2._dp*lnt*ln2t+lnu*ln2u)
     & +loqg/xn*(4._dp-Li2u-Li2s+lns*lnu+0.5_dp*lnt**2-0.5_dp*lns**2
     & -lnm*lnt+0.5_dp*lnm**2-lnu*ln2u-1.5_dp*(lnt-lnm))
     & +4._dp/3._dp*pi**2*xn*loqg
     & +0.25_dp*(xn**3-1._dp/xn)*(u+s)
     & +subuv*loqg

      return
      end


c This routine is meant to be valid for all s, t, u
      subroutine hjetfillgg(s,t,u,virtgg)
	use ggHwilson
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'b0.f'
      real(dp):: virtgg,logg,ddilog,Li2s,Li2t,Li2u,
     & lnm,mhsq,s,t,u,xlf,subuv,Delta
      complex(dp):: lns,lnt,lnu,ln2t,ln2u,lnrat

      mhsq=s+t+u

      xlf=real(nf,dp)
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)

      lns=lnrat(s,mhsq)
      lnt=lnrat(-t,mhsq)
      lnu=lnrat(-u,mhsq)
      lnm=log(musq/mhsq)
      ln2t=lnrat((mhsq-t),mhsq)
      ln2u=lnrat((mhsq-u),mhsq)

      logg=+V*xn*(mhsq**4+s**4+t**4+u**4)/(s*t*u)

c--- UV counterterm in MS bar scheme.
      subuv=-epinv*b0
      if (expansionorder > 0) then
! we are handling Wilson coefficient expansion separately
        Delta=0._dp
      else
        Delta=11._dp
      endif
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
c--- Three powers of as in Born --> 3
      subuv=3._dp*subuv+Delta

      virtgg=real(-3._dp*epinv*epinv2*xn*logg
     & +epinv*xn*logg*(lns+lnt+lnu-3._dp*lnm)
     & +xn*logg
     & *(2._dp*(Li2t+Li2u+Li2s)
     & +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
     & +0.5_dp*(lns**2-lnt**2-lnu**2)-1.5_dp*lnm**2
     & +2._dp*(lnu*ln2u+lnt*ln2t)+4._dp/3._dp*pisq)
     & +V*xn*(xn-xlf)/3._dp*mhsq*(1._dp+mhsq/s+mhsq/t+mhsq/u)
     & +subuv*logg)

c explicit pi^2 terms from analytic continuation not captured above
      if ((s > 0) .and. (t > 0) .and. (u > 0)) then
        virtgg=virtgg-logg*pisq*2._dp*xn
      endif

      return
      end


c This routine is meant to be valid for all s, t, u
      subroutine hjetfillqa(s,t,u,virtqa)
	use ggHwilson
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'b0.f'
      real(dp):: virtqa,loqa,ddilog,Li2s,Li2t,Li2u,
     & lnm,mhsq,s,t,u,xlf,subuv,Delta
      complex(dp):: lns,lnt,lnu,ln2t,ln2u,lnrat

      mhsq=s+t+u

      xlf=real(nf,dp)
      Li2t=ddilog(t/mhsq)
      Li2u=ddilog(u/mhsq)
      Li2s=ddilog((s-mhsq)/s)

      lns=lnrat(s,mhsq)
      lnt=lnrat(-t,mhsq)
      lnu=lnrat(-u,mhsq)
      lnm=log(musq/mhsq)
      ln2t=lnrat((mhsq-t),mhsq)
      ln2u=lnrat((mhsq-u),mhsq)

      loqa=xn*cf/s*(t**2+u**2)

c--- UV counterterm in MS bar scheme.
      subuv=-epinv*b0
      if (expansionorder > 0) then
! we are handling Wilson coefficient expansion separately
        Delta=0._dp
      else
        Delta=11._dp
      endif
c--- See C.R.Schmidt, PLB (413) 391, eq. (16),(17)
c--- Factor of ason2pi included in gg_hg_v.f
c--- Three powers of as in Born --> 3
      subuv=3._dp*subuv+Delta

      virtqa=real((-2._dp*xn+1._dp/xn)*loqa*epinv*epinv2
     & -2._dp/3._dp*xlf*epinv*loqa
     & +epinv*xn*loqa*(13._dp/6._dp-2._dp*lnm+lnt+lnu)
     & +epinv/xn*loqa*(1.5_dp-lns+lnm)
     & +loqa*xlf*(-10._dp/9._dp+2._dp/3._dp*lns-2._dp/3._dp*lnm)
     & +xn*loqa
     & *(40._dp/9._dp+Li2t+Li2u+2._dp*Li2s-13._dp/6._dp*(lns-lnm)
     & +(lnm-lns)*(lnt+lnu)+lns**2-lnm**2-0.5_dp*lnt**2-0.5_dp*lnu**2
     & +lnt*ln2t+lnu*ln2u)
     & +loqa/xn
     & *(4._dp-Li2t-Li2u-1.5_dp*(lns-lnm)+0.5_dp*(lns-lnm)**2
     & +lnt*lnu-lnt*ln2t-lnu*ln2u)
     & -4._dp/3._dp*pi**2/xn*loqa
     & -0.25_dp*(xn**3-1._dp/xn)*(t+u)
     & +subuv*loqa)

c explicit pi^2 terms from analytic continuation not captured above
      if (s < 0._dp) then
        virtqa=virtqa+loqa*(pisq*10._dp/3._dp)
        virtqa=-virtqa
      endif

      if ((s > 0) .and. (t > 0) .and. (u > 0)) then
        virtqa=virtqa+loqa*pisq*2._dp/3._dp
      endif

      return
      end
