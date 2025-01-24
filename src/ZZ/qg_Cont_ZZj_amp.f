!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qg_Cont_ZZj_amp(i1,i2,i7,za,zb,ampu,ampd)
c====== amplitudes for q(qb)+qb(q)+g+ZZ
c==== These amplitudes include virtual photon contributions.
c==== ordering is  0 --> q(-p1)+qb(-p2)+e-(p3)+e+(p4)+mu-(p5)+mu(p6)+g(p7)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'srdiags.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'zcouple_cms.f'
      integer:: i1,i2,i7,hg,hq,h34,h56,k
      reaL(dp):: q34,q56,s127
      complex(dp):: v34(2),v56(2)
      complex(dp):: ampu(2,2,2,2),ampd(2,2,2,2),
     & aq12(2,2,2,2),aq34(2,2,2,2),aq56(2,2,2,2),
     & prop34,prop56,prop127,rescale,fac

      ampu(:,:,:,:)=czip
      ampd(:,:,:,:)=czip
      v34(1)=zl1
      v34(2)=zr1
      q34=q1
      v56(1)=zl2
      v56(2)=zr2
      q56=q2

c===== everything but color in prefactor
      fac=im*four*rt2*abs(zesq)**2*sqrt(gsq)

c=== for now no neutrinos or 4 fermion interferences
      rescale=one

      srdiags=.true.
c==== propagators
      s127=s(i1,i2)+s(i1,i7)+s(i2,i7)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)
      prop34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)

c==== amplitude
      call zzgamp(i1,i2,3,4,5,6,i7,za,zb,aq12,aq34,aq56)
c==== dress with appropriate propagators and charge structure

      do hq=1,2
         do h34=1,2
            do h56=1,2
               do hg=1,2

c===== k=2 Up type amplitudes

      k=2
      if (hq == 1) then
      ampu(hq,hg,h34,h56)=(prop56*v56(h56)*zL(k)+q56*q(k))
     & *(prop34*v34(h34)*zL(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      ampu(hq,hg,h34,h56)=ampu(hq,hg,h34,h56)
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*zL(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale
     & *(prop127*v56(h56)*zL(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      ampu(hq,hg,h34,h56)=(prop56*v56(h56)*zR(k)+q56*q(k))
     & *(prop34*v34(h34)*zR(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      ampu(hq,hg,h34,h56)=ampu(hq,hg,h34,h56)
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*zR(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale
     & *(prop127*v56(h56)*zR(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      endif

c===== k=1 Down type amplitudes

      k=1
      if (hq == 1) then
      ampd(hq,hg,h34,h56)=(prop56*v56(h56)*zL(k)+q56*q(k))
     & *(prop34*v34(h34)*zL(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      ampd(hq,hg,h34,h56)=ampd(hq,hg,h34,h56)
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*zL(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale
     & *(prop127*v56(h56)*zL(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      ampd(hq,hg,h34,h56)=(prop56*v56(h56)*zR(k)+q56*q(k))
     & *(prop34*v34(h34)*zR(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      ampd(hq,hg,h34,h56)=ampd(hq,hg,h34,h56)
     & +(prop56*v34(h34)*v56(h56)+q34*q56)
     & *(prop127*v34(h34)*zR(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale
     & *(prop127*v56(h56)*zR(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      endif

c---- dress amplitudes with appropriate factors (everything but color (ave and V))

      ampu(hq,hg,h34,h56)=fac*ampu(hq,hg,h34,h56)
      ampd(hq,hg,h34,h56)=fac*ampd(hq,hg,h34,h56)

      enddo
      enddo
      enddo
      enddo


      return
      end




