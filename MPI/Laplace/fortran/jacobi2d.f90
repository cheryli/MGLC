! /* This code simulates a 2-dimensional laplace equation
!  * using jacobi interation
!  * accelated by mpi
!  * author: btli(2021)
!  */

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

program main
    implicit none
    integer, parameter :: nx = 1000, ny = 1000, itc_max = 1000
    real(8), parameter :: tolerance = 1.0e-5
    integer :: i, j, itc
    real(8) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1), f(0:nx+1, 0:ny+1), A_p(0:nx+1, 0:ny+1)
    real(8) :: error, error_max, check_diff
    real(8) :: start_time, end_time



    call init(A, A_new, f, nx, ny)

    call CPU_TIME(start_time)

    !! jacobi iteration
    itc = 0
    error_max = 1.0
    A_p = A
    do while(error_max > tolerance .AND. itc <= itc_max)
        itc =  itc + 1

        call jacobi(A, A_new, f, nx, ny)
        call swap(A, A_new)

        if(MOD(itc, 100) == 0) then
        error = check_diff(A, A_p, nx, ny)
        write(*, *) itc, error_max
        endif
    enddo

    call CPU_TIME(end_time)
    write(*,*) "Total run time =", end_time-start_time, "s."
    

end program main



subroutine init(A, A_new, f, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1), f(0:nx+1, 0:ny+1)
    integer :: i, j

    ! initialize
    do j = 0, ny+1
        do i = 0, nx+1
            A(i, j) = 0.0d0
            A_new(i, j) = 0.0d0
            f(i, j) = 0.0d0
        enddo
    enddo

    ! boundary condition (top boundary = 1.0, others = 0.0)
    do i = 0, nx+1 
        A(i, ny+1) = 1.0d0
        A_new(i, ny+1) = 1.0d0
    enddo

end subroutine init


subroutine jacobi(A, A_new, f, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1), f(0:nx+1, 0:ny+1)
    integer :: i, j

    do j = 1, ny
        do i = 1, nx
            A_new(i, j) = 0.25 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) + f(i, j)) 
        enddo
    enddo

end subroutine jacobi



function check_diff(A, A_p, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: A(0:nx+1, 0:ny+1)
    real(8), intent(inout) :: A_p(0:nx+1, 0:ny+1)
    integer :: i, j
    real(8) :: error, check_diff

    error = 0.0d0
    do j = 1, ny
        do i = 1, nx
            error = max( error, abs(A_p(i,j)-A(i,j)) )
        enddo
    enddo

    A_p = A
    check_diff = error

end function check_diff


subroutine swap(A, A_new, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(inout) :: A(0:nx+1, 0:ny+1), A_new(0:nx+1, 0:ny+1)
    integer :: i, j

    do j = 1, ny
        do i = 1, nx
            A(i, j) = A_new(i, j)
        enddo
    enddo

end subroutine swap