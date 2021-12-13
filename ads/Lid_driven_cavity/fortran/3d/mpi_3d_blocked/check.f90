subroutine check()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    real(8) :: error1, error2
    real(8) :: total_error1, total_error2

    error1 = 0.0d0
    error2 = 0.0d0

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                error1 = error1 + (u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k))+(v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k))
                error2 = error2 + u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
            enddo
        enddo
    enddo

    up = u
    vp = v
    wp = w

    call MPI_Barrier(MPI_COMM_WORLD, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    errorU = sqrt(total_error1)/sqrt(total_error2)

    if (rank == 0) then
        write(*,*) itc,' ',errorU
    endif

    return
end subroutine check
