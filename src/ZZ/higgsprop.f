      function higgsprop(s)
      implicit none
      include 'types.f'
      complex(dp):: higgsprop
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      real(dp):: s,mhbarsq,mhbar,gammahbar
      complex(dp):: cfac
      save mhbarsq,cfac

      if (CPscheme) then
c--- complex pole scheme propagator      
        if (first) then
          mhbarsq=hmass**2+hwidth**2
          mhbar=sqrt(mhbarsq)
          gammahbar=mhbar/hmass*hwidth
          cfac=cplx2(one,gammahbar/mhbar)
          first=.false.
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'*  Using complex pole scheme for Higgs propagator  *'
        write(6,99) mhbar,gammahbar
        write(6,*)'****************************************************'
        write(6,*)
        endif
        higgsprop=cfac/(s*cfac-cplx2(mhbarsq,zip))
      else
c--- Breit Wigner propagator      
        higgsprop=one/cplx2(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end

c--- 
c---  One-loop correction to the Higgs propagator
c---      
      function sigmah(s,c6,wfr)
      use loopI2_generic
      implicit none
      include 'types.f'
      complex(dp):: sigmah

      include 'constants.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      ! include 'qlfirst.f'
      real(dp):: s, c6, wfr, MH2, dB0h, dZh, dmhsq
      ! complex(dp):: qlI2

      ! if (qlfirst) then
      !   qlfirst=.false. 
      !   call qlinit
      ! endif

      MH2=hmass**2

c--- Wavefunction renormalisation
      dB0h = (-9 + 2*Sqrt(3.d0)*Pi)/(9*MH2)
      dZh  = -wfr*dB0h  
c--- Mass renormalisation
      dmhsq = loopI2(MH2,MH2,MH2,1D0,0)
c--- Renormalised correction
      sigmah = (9*c6*(2.d0 + c6)*MH2**2*
     & (-dZh*(MH2 - s)))/
     &     (32.d0*Pi**2*vevsq)

      sigmah = sigmah + (9*c6*(2.d0 + c6)*MH2**2*
     & (loopI2(s,MH2,MH2,1D0,0) - dmhsq))/
     & (32.d0*Pi**2*vevsq)

      sigmah= sigmah/cplx2(s-MH2,hmass*hwidth) 


      end

c--- 
c---  One-loop correction to the Higgs propagator due to marginal Higgs portal
c---      
c      function sigmahx(s,cx,mx)
c      implicit none
c      include 'types.f'
c      complex(dp):: sigmahx
c      real(dp):: MH2 
c      include 'constants.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'ewcouple.f'
c      include 'qlfirst.f'
c      real(dp):: s, cx, mx
c      complex(dp):: qlI2
c
c     if (qlfirst) then
c        qlfirst=.false. 
c        call qlinit
c      endif
c
c      MH2=hmass**2
c
c--- Renormalised correction without wavefunction correction
c
c      sigmahx = cx**2*vevsq* 
c     &     (loopI2(s,mx**2,mx**2,1D0,0) -
c     &      dreal(loopI2(MH2,mx**2,mx**2,1D0,0)))/
c     & (8.0d0*Pi**2)
c
c      sigmahx = sigmahx/cplx2(s-MH2,hmass*hwidth) 
c
c      end

c--- 
c---  One-loop correction to the Higgs propagator due to fermion Higgs portal
c---      
      function sigmahx(s,cx,mx)
      use loopI2_generic
      implicit none
      include 'types.f'
      complex(dp):: sigmahx
c      complex(dp):: sigmahx, sigmahxcheck, Bderivative
      real(dp):: MH2 
      include 'constants.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      ! include 'qlfirst.f'
      real(dp):: s, cx, mx, f2, mu2, Nd, CHBox
      complex(dp):: qlI2

      ! if (qlfirst) then
      !   qlfirst=.false. 
      !   call qlinit
      ! endif

      MH2=hmass**2

c--- Renormalised correction without wavefunction correction

c      sigmahx = 3d0*cx**2*vevsq* 
c     &     ((s-4*mx**2)*loopI2(s,mx**2,mx**2,s,0) -
c     &      (MH2-4*mx**2)*dreal(loopI2(MH2,mx**2,mx**2,MH2,0)))/
c     & (8.0d0*Pi**2)/(9d0*vevsq)

c--- Renormalised correction with wavefunction correction 
c--- My implementation

c      sigmahx = 3d0*cx**2*vevsq*
c     &     ((s-4*mx**2)*
c     & 		(loopI2(s,mx**2,mx**2,1d0,0)-dreal(loopI2(MH2,mx**2,mx**2,1d0,0))) -
c     &      (s-MH2)*(2d0*mx**2*dreal(loopI2(MH2,mx**2,mx**2,mx**2,0))/MH2-1d0))/
c     & (8d0*Pi**2)/(9d0*vevsq)

c--- Konstantin's implementation
c--- This agrees with mine

c      Bderivative = ((-1d0)/MH2)+((2d0*mx**2)/(MH2*(MH2-4d0*mx**2)))*
c     & (dreal(loopI2(MH2,mx**2,mx**2,mx**2,0))-2d0)
c
c      sigmahxcheck = (3d0*cx**2*vevsq*((s-4*mx**2)*
c     & ((loopI2(s,mx**2,mx**2,1d0,0)-dreal(loopI2(MH2,mx**2,mx**2,1d0,0))))
c     & +(4*mx**2-MH2)*(s-MH2)*Bderivative))/(9d0*8d0*Pi**2*vevsq)      
c
c      write(6,99) sigmahxcheck-sigmahx      

c--- Renormalised correction with inclusion of wave function 
c--- and higher-dimensional operator correction   
c--- My implementation

      Nd = 3d0
      f2 = 3d0*vevsq
      mu2 = 16d0*Pi**2*f2

c      Nd = 1d0
c      f2 = 1000d0**2
c      mu2 = 1000d0**2

      CHBox = 1d0/(2d0*f2)

      sigmahx = -Nd*cx**2*vevsq*
     &     ((s-4*mx**2)*loopI2(s,mx**2,mx**2,mx**2,0) -
     &      (MH2-4*mx**2)*dreal(loopI2(MH2,mx**2,mx**2,mx**2,0)) + 
     &      (s-MH2)*log(mu2/(mx**2)))/
     & (8d0*Pi**2)/f2 - 2d0*vevsq*CHBox*(s-MH2)		

      sigmahx = sigmahx/cplx2(s-MH2,hmass*hwidth) 
      
      return

c   99 format(' * Check = ',f9.4x)

      end      

     

























