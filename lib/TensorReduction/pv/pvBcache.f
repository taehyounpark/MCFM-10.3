      integer function pvBcache(p1sq,m1sq,m2sq)
      implicit none
C---p1sq is the square of the momentum
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvBnames.f'
      include 'lib/TensorReduction/Include/TRclear.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/TRonshellcutoff.f'
      real(dp):: para(Pbb),p1sq,m1sq,m2sq
      integer jtable,j,Ntrue
      real(dp),save:: tableB(Pbb,Nbmax)      
      integer,save:: Nstore=0
!$omp threadprivate(tableB,Nstore)

      if (clear(2)) then
      clear(2)=.false.
      Nstore=0
      endif

      if (Nstore .gt. NBmax) then
      print * 
      print *, 'pvBcache: Nstore .gt. Nbmax'
      print *, 'pvBcache:Nstore,Nbmax',Nstore,Nbmax
      print *, 'Either adjust Nbmax in Bnames.f and recompile'
      print *, 'or call clearcache to clear the cache.'
      stop
      endif
      para(1)=p1sq
      para(2)=m1sq
      para(3)=m2sq
C if parameter set is found set pvBcache equal to the starting
C value
      if (Nstore .eq. 0) go to 20
      do jtable=1,Nstore
      Ntrue=0
        do j=1,Pbb
        if (abs(para(j)-tableB(j,jtable)) .lt. 1d-8) Ntrue=Ntrue+1 
        enddo
      if (Ntrue .eq. Pbb) then
      pvBcache=(jtable-1)*Nbb
      return
      endif
      enddo

C    if parameter set is not found we have to calculate
 20   pvBcache=Nstore*Nbb
      Nstore=Nstore+1
      do j=1,Pbb
        if(abs(para(j)) .lt. onshellcutoff) para(j)=zero
      enddo
      do j=1,Pbb
      tableB(j,Nstore)=para(j)
      enddo
c      call pvBfill(p1sq,m1sq,m2sq,pvBcache)
      call pvBfill(para(1),para(2),para(3),pvBcache)
      return
      end 
