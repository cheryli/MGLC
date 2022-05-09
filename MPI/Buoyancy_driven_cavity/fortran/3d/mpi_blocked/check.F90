subroutine check()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    

    do k=1,nz
        do j=1,ny
            do i=1,nx
                error1 = error1+(u(i,j,k)-up(i,j,k))*(u(i,j,k)-up(i,j,k))+(v(i,j,k)-vp(i,j,k))*(v(i,j,k)-vp(i,j,k))+(w(i,j,k)-wp(i,j,k))*(w(i,j,k)-wp(i,j,k))
                error2 = error2+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                
                error5 = error5+dABS( T(i,j,k)-Tp(i,j,k) )
                error6 = error6+dABS( T(i,j,k) )
                
                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
                wp(i,j,k) = w(i,j,k)
                Tp(i,j,k) = T(i,j,k)
            enddo
        enddo
    enddo

    call MPI_Barrier(comm3d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    return
end subroutine check