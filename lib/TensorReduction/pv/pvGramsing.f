      logical function pvGramsing(G,n)
c--- JC: added comparison to onshellcutoff for each element of the
c---     matrix passed in; if the element is smaller, set to zero
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/TRonshellcutoff.f'
      include 'lib/TensorReduction/Include/TRscale.f'
      integer:: nmax,n,j,k
      complex(dp):: G(n,n)
      real(dp)::preci,wmax,wmin
      parameter(nmax=4,preci=1d-7)
c--- Regular PV reduction fails checks at the C4 level at approx.
c---  10^-5 level of precision; lower tensors correspond to
c---  much smaller values of preci.      
      real(dp)::Ga(nmax,nmax),V(nmax,nmax),w(nmax)
C--- logical function which return true if the
C--- Gram matrix is singular
      if (n .gt. nmax) then
      write(6,*) 'Error in pvGramsing, n .gt. nmax'
      stop
      endif

      do j=1,n
      do k=1,n
        Ga(j,k)=dble(G(j,k))
c--- the next line improves convergence in pvdsvdcmp
        if (abs(Ga(j,k)/musq) .lt. onshellcutoff) Ga(j,k)=zero
      enddo
      enddo
            
      call pvdsvdcmp(Ga,n,n,nmax,nmax,w,v)

      wmax=zero
      do j=1,n
      if (w(j) .gt. wmax)wmax=w(j)
      enddo

      wmin=preci*wmax
      pvGramsing=.false.

      do j=1,n
c        write(6,*) 'wj ',w(j)/wmax
        if (w(j) .lt. wmin) then
          pvGramsing=.true.
          return
        endif
      enddo

      end
