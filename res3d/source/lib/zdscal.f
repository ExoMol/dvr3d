      subroutine zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c
      double complex zx(1)
      double precision da
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increments not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increments equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end
