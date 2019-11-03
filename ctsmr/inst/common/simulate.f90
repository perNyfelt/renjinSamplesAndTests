subroutine simulate(x, y, x0, uo, xm, t, nt, dt)

implicit none

include 'global.h'

! arguments
integer :: nt
double precision :: dt
double precision, dimension(n) :: x0
double precision, dimension(nparam) :: xm
double precision, dimension(nt) :: t


double precision, dimension(n, nt) :: x
double precision, dimension(s, nt) :: y

double precision, dimension(m, nt) :: uo

! variables
double precision, dimension(n) :: f, xi, dw, diff
double precision, dimension(s) :: h, e
double precision :: tc, sqrtdt
integer :: i, j, nsteps

double precision :: jacob1(n,n), jacob2(s,s), sig(n, n), sm(s,s)

double precision, dimension(m) :: uoc, alpha


! Initialise the RNG seed
call rndstart()

! if linear FAIL

! Assume non linear model for now

xi = x0

sqrtdt = sqrt(dt)

tc = t(1)

do i=1,(nt-1)

   nsteps = INT((t(i+1) - t(i)) / dt)

   alpha = (uo(:,i+1) - uo(:,i)) / REAL(nsteps)

   do j = 1,nsteps

      ! constant

      ! linear interpolation
      uoc = alpha * REAL(j-1) + uo(:,i)

      call fvecxj(xm,xi,uoc,tc,JACOB1,f,NPARAM,N,M)
      call sigmat(xm, uoc, tc, sig, NPARAM, N, M)

      call normrnd(dw, n)

      ! multiply sqrt(dt) * sig * rnorm(n)
      call dgemv('N',n,n,sqrtdt,sig,n,dw,1,0.0d0,diff,1)

      xi = xi + f * dt + diff
      tc = tc + dt
   end do

   x(:,i) = xi;

   ! y = h(xi)
   call HVECXJ(XM,xi,uoc,tc,JACOB2,y(:,i),NPARAM,N,M,S)

   ! noise
   call normrnd(e, s)
   call smat(xm,uoc,tc,SM,NPARAM,N,M,S)
   
   ! y = 1.0 * S * e + y
!  DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
   call dgemv('N',s,s,1.0d0,SM,s,e,1,1.0d0,y(:,i),1)
   
end do

! Finalise the RNG seed
call rndend()


end subroutine simulate
