      integer:: nplot,maxhisto
      parameter(maxhisto=200,nplot=4*maxhisto)
      character(len=3):: linlog(nplot)
      character(len=8):: titlearray(nplot)
      common/topd/titlearray,linlog
      integer:: nextnplot
      common/plotindex/nextnplot
!$omp threadprivate(/plotindex/)
