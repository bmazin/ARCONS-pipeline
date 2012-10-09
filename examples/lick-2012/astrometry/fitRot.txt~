      program fitFromTxt
      integer ip
      double precision x, y, rah, ram, ras, decd, decm, decs
      double precision xye(2,100), xym(2,100), coeffs(6)
      double precision xz, yz, xs, ys, perp, orient
      integer j
      integer np
      integer itype

      print *, "Hello"
      np = 0
      do ip = 1, 100
        read(*,*,end=100) x, y, ray, ram, ras, decd, decm, decs
        print *, "x=",x
        xym(1,ip) = x
        xym(2,ip) = y
        xye(1,ip) = 15.0*(rah+ram/60.0+ras/3600.0)
        xye(2,ip) = decd + decm/60.0 + decs/3600.0
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
      end program fitFromTxt
