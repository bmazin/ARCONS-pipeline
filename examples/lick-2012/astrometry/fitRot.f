      program fitFromTxt
      integer ip
      double precision x1, y1, x2, y2
      double precision xye(2,100), xym(2,100), coeffs(6)
      double precision xz, yz, xs, ys, perp, orient
      integer j
      integer np
      integer itype
      double precision xc, yc
      double precision a,b,c,d,e,f
      print *, "Hello"
      np = 0
      do ip = 1, 100
        read(*,*,end=100) x1, y1, x2, y2
        print *, "x1=",x1
        xym(1,ip) = x1
        xym(2,ip) = y1
        xye(1,ip) = x2
        xye(2,ip) = y2
        np = np + 1
      enddo
 100  print *," DONE with np=",np
      itype = 4
      call sla_FITXY(itype, np, xye, xym, coeffs, j)
      print *," j=",j
      do ip=1,6
         print *," ip=",ip," coeffs=",coeffs(ip)
      enddo
      call sla_DCMPF(coeffs, xz, yz, xs, ys, perp, orient)
      print *, "orient=",orient
      print *, "perp=",perp
      print *, "xs=", xs
      print *, "ys=", ys
      print *, "xz=", xz
      print *, "yz=", yz

      a = coeffs(1)
      b = coeffs(2)
      c = coeffs(3)
      d = coeffs(4)
      e = coeffs(5)
      f = coeffs(6)
      xc = (a + (c*d)/(1-f))/(1 - b - (c*e)/(1-f))
      yc = (d+e*xc)/(1-f)
      print *, "xc=",xc," yc=",yc
      end program fitFromTxt
