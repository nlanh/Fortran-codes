program Bessel

implicit none

real*16 x,besj,xmax,step,fb,gb,tmax,ztrho
real a,g1,g2,xfrho,xgrho,drho,rhomax,xf,xg,c,dit,cta

integer*16 j,n,s,nmax,fac1,fac2
integer irho,ix,ni,krho,nrho,nit,jt

open(unit=1, file='besj.dat')
open(unit=2, file='frho.dat')
open(unit=3, file='grho.dat')
open(unit=4, file='ztrho.dat')

print*, 'Enter n (order of bessel function, hint: 0)'
read*, n
print*, 'Enter nmax (infinity index for series, hint: 20)'
read*, nmax
print*, 'Enter xmax (maximum range of x for output, hint: 10)'
read*, xmax
print*, 'Enter step (step of x, hint: 0.01)'
read*, step

print*, 'Calculation of frho, grho, zthro (yes = 1, no = 0)'
read*, krho
if (krho.eq.1) then
	print*, 'Enter a (hint: 1)'
	read*, a 
	print*, 'Enter c for grho (hint: 0.1)'
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
	print*, 'Enter tmax of t (hint: 10)'
	read*, tmax
	print*, 'Enter step of t (hint: 0.1)'
	read*, dit
	nrho = int((rhomax+0.00000001)/drho)
	nit = int((tmax+0.00000001)/dit)
end if 

ni = int((xmax+0.00000001)/step)

do ix=0,ni
	besj=0
	x=ix*step
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


do irho=0,nrho
	fb=0
	gb=0
	xf=irho*xfrho*drho
	xg=irho*xgrho*drho
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
	write(2,'((f8.2),2(f20.4))') irho*drho, fb, -0.01*a*fb
	write(3,'((f8.2),2(f20.4))') irho*drho, gb, -0.1*xgrho*c*gb
	do jt=1,nit
		cta=c*jt*dit/a
		ztrho=-0.01*a*fb*cos(g1*cta)-0.1*gb*sin(g2*cta)
		write(4,'(2(f15.4),(f20.8))')jt*dit,irho*drho,ztrho
	end do
end do

end program
