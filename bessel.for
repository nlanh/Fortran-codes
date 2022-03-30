      program Bessel

      implicit real*8 (a-h,o-z)

      open(unit=1, file='besj.dat')


      print*, 'Enter n'
      read*, n
      print*, 'Enter xmax'
      read*, xmax
      print*, 'Enter step'
      read*, step

      ni = int((xmax+0.00000001)/step)

      do ix=0,ni
        besj=0
        x=ix*step
      do ns=0,100
        fac1=1
        if (ns.ne.0) then
            fac1=1
            do j=1,ns
                fac1=fac1*j
            end do
        end if
        fac2=1
        if ((n+ns).ne.0) then
            fac2=1
            do j=1,n+ns
                fac2=fac2*j
            end do
        end if
        besj = besj + (x/2)**(n+2*ns)*(-1)**ns/(fac1*fac2)
        end do
        write(1,'((f8.2),(f20.4))') x, besj

      end do

      end program
