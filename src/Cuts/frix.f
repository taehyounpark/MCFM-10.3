!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- Implementation of "Frixione" or "smooth cone" isolation. according to
c--- the procedure described in S. Frixione,hep-ph/9801442

c---  C. Williams July 2011
c---  C. Williams July 2015 -- extended to allow for multiple partons
c---                           in the cone (can happen at NNLO)
c---  T. Neumann 2020: Allow for fixed cone energy
c---  T. Neumann 2021: Hybrid isolation scheme

      subroutine frix(p,passed,j,isub)
c----- p -momentum array passed
c----- j - photon identification in p i.e. p(j,nu) = photon(nu)
c----- isub whether we are working with a dipole or not
c----- paramaters eps, delta_0 and n_pow are taken from input.DAT
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'frag.f'
      include 'npart.f'
      include 'first.f'
      include 'mpicommon.f'
      real(dp):: p(mxpart,4),R,pref,ret_ET
      integer:: j,isub,k
      logical:: passed,is_hadronic,in_cone_n
      integer:: i
      real(dp):: vsmall
      parameter(vsmall=1.e-10_dp)
      real(dp):: ET_had,pt
      integer:: itag

      passed=.true.

      if(first) then
        first=.false.
c----- check for non-zero parameters, if zero exit with warning
        if((epsilon_h<vsmall).or.(cone_ang<vsmall)) then
!$omp master
          if (rank == 0) then
          write(6,*)
          write(6,*) '************** Frixione Isolation    ***************'
          write(6,*) '*   Read zero parameters, not isolating            *'
          write(6,*) '*   WARNING: this may be unsafe in general         *'
          write(6,99)'*  eps_phot = ',epsilon_h,' delta_0 = ',cone_ang, '*'
          write(6,97)'*  n = ',n_pow,'                                   *'
          write(6,*) '****************************************************'
          endif
!$omp end master
          return
        endif


!$omp master
        if (rank == 0) then
        write(6,*)
        write(6,*) '************** Frixione Isolation    ***************'
        write(6,*) '*                                                  *'
        write(6,99)'*  eps_phot = ',epsilon_h,', delta_0 = ',cone_ang,  '*'
        write(6,97)'*  n = ',n_pow,'                                   *'
        write(6,*) '****************************************************'
        endif
!$omp end master
      endif

 99   format (1x,a14,f5.3,a12,f5.3,a16)
c 97   format (1x,a7,i1,a44)
 97   format (1x,a7,f5.2,a40)

c----- Cycle over final state particles, if it is hadronic and inside
c----- photon cone then check its energy passes frixione requirement
c----- else fail and return

      if (.not. fixed_coneenergy) then
          pref = ret_ET(p,j)*epsilon_h
      else
          pref=epsilon_h
      endif
      ET_had=0._dp

c===== this section is altered now to allow for an additional parton
c===== < R_ij inside the cone too, which can happen at NNLO
      itag=0
      do i=3,2+npart-isub
c===== reset ET_had for each initial hadron
         ET_had =0._dp
c======first thing, find a hadron inside isolation cone
         if(is_hadronic(i).and.(R(p,i,j)<cone_ang)) then
            ET_had=ET_had+ret_ET(p,i)
c========= now we need to sum over the additional hadronic momenta looking for
c========= friends inside the cone
            do k=3,2+npart-isub
c========= is k nearer the photon then i?
               if(is_hadronic(k).and.(k /= i)) then
               if(R(p,k,j) < R(p,i,j)) then
                  ET_had=ET_had+ret_ET(p,k)
c                  itag=2
               endif
               endif
c====== include this hadron in total energy
            enddo
c========= now check total energy is isolated against
            passed=in_cone_n(R(p,i,j),ET_had,pref)

            if(itag==2) then
               call writeout(p)
               write(6,*) 'photon = ',j
               write(6,*) 'pt(photon) = ',pt(j,p)
               write(6,*) 'R(5,j) = ',R(p,5,j)
               if(isub==0) write(6,*) 'R(6,j)= ',R(p,6,j)
               write(6,*) 'PT(5) = ',pt(5,p)
               write(6,*) 'PT(6) = ',pt(6,p)
               write(6,*) 'RET_ET(5)  = ',ret_ET(p,5)
               if(isub==0) write(6,*) 'RET_ET(6)  = ', ret_ET(p,6)

               write(6,*) 'total ET_had in cone',ET_had
               !write(6,*) 'iso condition',pref*(1._dp-cos(R(p,i,j)))
               write(6,*) 'passed = ',passed
               error stop 'in frix.f'
            endif
            if(passed.eqv..false.) return
         endif
      enddo

      return
      end

      function ret_ET(p,j)
      implicit none
      include 'types.f'
      real(dp):: ret_ET

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
      integer:: j
      real(dp):: ptsq

      ptsq=p(j,1)**2+p(j,2)**2
      ret_ET=p(j,4)*sqrt(ptsq)/(sqrt(ptsq+p(j,3)**2))

      if(ptsq==0) ret_ET=0._dp
      return
      end

      function in_cone_n(Rij,Ejet,pref)
          use types
          implicit none
          include 'frag.f'
          logical:: in_cone_n

          real(dp), intent(in) :: Rij,Ejet,pref

          real(dp) :: chi_inner, chi_outer, chi

          if (photiso_hybrid) then
              chi_inner = (1._dp-cos(Rij))**n_pow * ((1._dp-cos(R_inner))**(-n_pow))
              chi_outer = 1._dp
              if (Rij < R_inner) then
                  chi = chi_inner
              elseif (Rij < cone_ang) then
                  chi = chi_outer
              else
                  in_cone_n = .false.
                  return
              endif
          else
              chi = (1._dp-cos(Rij))**n_pow * ((1._dp-cos(cone_ang))**(-n_pow))
          endif

          if ( Ejet < pref*chi) then
             in_cone_n = .true.
          else
             in_cone_n = .false.
          endif
      end function













