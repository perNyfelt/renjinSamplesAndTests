module covariance
implicit none

integer, private, save :: curpos

contains

subroutine clear_cov_counter
   curpos = 0
end subroutine clear_cov_counter

subroutine store_covariance(wk,Nwk,vec,mat,N,full)
   implicit none
   integer :: i,j,N,full,Nwk
   double precision :: wk(Nwk), mat(N*(N+1)/2), vec(N)

   wk((curpos+1) : (curpos + N)) = vec

   curpos = curpos + N

   if (full == 0) then
      i=1
      do j=1,N
         wk(CURPOS+J) = DSQRT(mat(i))
         i=i+N-j+1
      end do
      curpos = curpos + N
   else
      wk((curpos + 1) : (curpos + (N*(N+1)/2))) = mat
      curpos = curpos + (N*(N+1)/2)
   end if

end subroutine store_covariance

end module covariance
