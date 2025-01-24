! Fortran
      integer,parameter::nlegborn=4,nlegreal=nlegborn+1
      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1+1)
      integer,parameter:: maxprocborn=31,maxprocreal=176
      integer maxalr
      parameter (maxalr=maxprocreal*nlegreal*(nlegreal-1)/2)
