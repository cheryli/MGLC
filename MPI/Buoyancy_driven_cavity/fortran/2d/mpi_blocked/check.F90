#ifdef steadyFlow
    subroutine check()
    use mpi
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: error1, error2, error5, error6
    real(8) :: total_error1, total_error2, total_error5, total_error6

    error1 = 0.0d0
    error2 = 0.0d0

    error5 = 0.0d0
    error6 = 0.0d0
    

    do j=1,ny
        do i=1,nx
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                error5 = error5+dABS( T(i,j)-Tp(i,j) )
                error6 = error6+dABS( T(i,j) )
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j)
                Tp(i,j) = T(i,j)
        enddo
    enddo

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error5, total_error5, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error6, total_error6, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)
    errorT = total_error5 / total_error6

    if (rank == 0) then
        open(unit=01,file='convergence.log',status='unknown',position='append')
            write(01,*) itc,' ',errorU,' ',errorT
        close(01)
        
        write(*,*) itc,' ',errorU,' ',errorT
    endif

    
    return
    end subroutine check
#endif