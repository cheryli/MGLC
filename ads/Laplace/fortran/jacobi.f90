! /*  data mesh arrangement
!  *  ny
!  *  ^
!  *  | (1,3)  
!  *  | (1,2)   
!  *  | (1,1)  (2,1)  (3,1)
!  *  ---------------------> nx
!  * 
!     for fortran language, data fill row first, (1,1),(2,1) ...
! */


module jacobi
    implicit none
    public :: calc_next
    public :: swap
contains
    function calc_next(A, A_new, nx, ny)
        implicit none
        integer, parameter          :: fp_kind=kind(1.0d0)
        integer,intent(in)          :: nx, ny
        real(fp_kind),intent(inout) :: A(0:nx-1,0:ny-1)
        real(fp_kind),intent(inout) :: A_new(0:nx-1,0:ny-1)
        real(fp_kind)               :: error, calc_next
        integer                     :: i, j

        error = 0.0d0

        do j=1,ny-2
            do i=1,nx-2
              A_new(i,j) = 0.25d0 * ( A(i+1,j) + A(i-1,j) + &
                                           A(i,j-1) + A(i,j+1) )
              error = max( error, abs(A_new(i,j)-A(i,j)) )
            end do
          end do        

        calc_next = error
    end function calc_next

    subroutine swap(A, Anew, nx, ny)
        implicit none
        integer, parameter        :: fp_kind=kind(1.0d0)
        real(fp_kind),intent(out) :: A(0:nx-1,0:ny-1)
        real(fp_kind),intent(in)  :: Anew(0:nx-1,0:ny-1)
        integer,intent(in)        :: nx, ny
        integer                   :: i, j
  
        do j=1,ny-2
          do i=1,nx-2
            A(i,j) = Anew(i,j)
          end do
        end do
      end subroutine swap
end module jacobi