program Bessel

implicit none

real*16 x,besj,xmax,step,fb,gb
real a,g1,g2,xfrho,xgrho,drho,rhomax,xf,xg,c

integer*16 i,j,n,s,nmax,ni,fac1,fac2,nrho
integer krho 

open(unit=1, file='besj.dat')
open(unit=2, file='frho.dat')
open(unit=3, file='grho.dat')

print*, 'Enter n (order of bessel function, hint: 0)'
read*, n
print*, 'Enter nmax (infinity index for series, hint: 20)'
read*, nmax
print*, 'Enter xmax (maximum range of x for output, hint: 10)'
read*, xmax
print*, 'Enter step (step of x, hint: 0.01)'
read*, step

print*, 'Calculation of frho and grho (yes = 1, no = 0)'
read*, krho
if (krho.eq.1) then
	print*, 'Enter a (hint: 1)'
	read*, a 
	print*, 'Enter c for grho (hint: 1)'
	read*, c
	print*, 'Enter gamma1 for frho (hint: 2.4048)'
	read*, g1
	print*, 'Enter gamma2 for grho (hint: 5.5201)'
	read*, g2
	xfrho = g1/a
	xgrho = g2/a 
	print*, 'Enter rhomax of rho (hint: 1)'
	read*, rhomax
	print*, 'Enter step of rho (hint: 0.01)'
	read*, drho
	nrho = int((rhomax+0.00000001)/drho)
end if 


ni = int((xmax+0.00000001)/step)


do i=0,ni
	besj=0
	x=i*step
	do s=0,nmax
		fac1=1
		if (s.ne.0) then
		    fac1=1
			do j=1,s
				fac1=fac1*j
			end do
		end if
		fac2=1
		if ((n+s).ne.0) then
			fac2=1
			do j=1,n+s
				fac2=fac2*j
			end do
		end if
		besj = besj + (x/2)**(n+2*s)*(-1)**s/(fac1*fac2)
	end do
	write(1,'((f8.2),(f20.4))') x, besj
end do


do i=0,nrho
	fb=0
	gb=0
	xf=i*xfrho*drho
	xg=i*xgrho*drho
	do s=0,nmax
		fac1=1
		if (s.ne.0) then
		    fac1=1
			do j=1,s
				fac1=fac1*j
			end do
		end if
		fac2=1
		if ((n+s).ne.0) then
			fac2=1
			do j=1,n+s
				fac2=fac2*j
			end do
		end if
		fb = fb + (xf/2)**(n+2*s)*(-1)**s/(fac1*fac2)
		gb = gb + (xg/2)**(n+2*s)*(-1)**s/(fac1*fac2)
	end do
	write(2,'((f8.2),2(f20.4))') i*drho, fb, -0.01*a*fb
	write(3,'((f8.2),2(f20.4))') i*drho, gb, -0.1*xgrho*c*gb
end do

end program
