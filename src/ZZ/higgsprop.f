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
c---  UH 17/03/25: Checked that this agrees with the MCFM7 and MCFM8 versions
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
      real(dp):: s, c6, wfr, MH2, dB0h, B0h, dZh, dmhsq
      ! complex(dp):: qlI2

      ! if (qlfirst) then
      !   qlfirst=.false.
      !   call qlinit
      ! endif

      MH2=hmass**2

c      c6=10D0
c      s=300D0**2
      
c--- Wavefunction renormalisation
      dB0h = (-9D0 + 2D0*Sqrt(3D0)*Pi)/(9D0*MH2)
      dZh  = -wfr*dB0h  
c--- Loop integral
      B0h = loopI2(s,MH2,MH2,1D0,0)
c--- Mass renormalisation
      dmhsq = loopI2(MH2,MH2,MH2,1D0,0)
c--- Renormalised correction
      sigmah = (9D0*c6*(2D0 + c6)*MH2**2*(B0h - dmhsq + dZh*(-MH2 + s)))/(32D0*Pi**2D0*vevsq)

c      write(*,*) dreal(sigmah)

      sigmah= sigmah/cplx2(s-MH2,hmass*hwidth) 

      end

c--- 
c---  One-loop correction to the Higgs propagator due to marginal Higgs portal
c---      
      function sigmahx(s,cx,mx)
      use loopI2_generic
      implicit none
      include 'types.f'
      complex(dp):: sigmahx
      real(dp):: MH2 
      include 'constants.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      ! include 'qlfirst.f'
      real(dp):: s, cx, mx
      ! complex(dp):: qlI2

      ! if (qlfirst) then
      ! qlfirst=.false. 
      ! call qlinit
      ! endif

      MH2=hmass**2

c--- Renormalised correction without wavefunction correction

      sigmahx = cx**2d0*vevsq* 
     &     (loopI2(s,mx**2,mx**2,1D0,0) -
     &      dreal(loopI2(MH2,mx**2,mx**2,1D0,0)))/
     & (8.0d0*Pi**2d0)

      sigmahx = sigmahx/cplx2(s-MH2,hmass*hwidth) 

      end

























