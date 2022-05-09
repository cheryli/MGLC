program main  
    use mpi
    use omp_lib   
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads
    real(8) :: start_time, end_time
    integer :: name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    dims(0) = 1
    call MPI_Dims_create(num_process, 2, dims, rc)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, rc)
    if(rank == 0) then
        write(*,*) "dimens is x*y = ", dims(0), "x", dims(1)
    endif

    ! get my new rank in decomposition
    call MPI_Comm_rank(comm2d, rank2d, rc)
    ! write(*,*) "process ", rank2d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    call MPI_Cart_get(comm2d, 2, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    ! write(*,*) "coords, nx, ny, is, js = ", coords(0), coords(1), nx, ny, i_start_global, j_start_global
    ! write(*,*) "nx*ny = ", nx, ny

    ! get the neighbors
    call MPI_Cart_shift(comm2d, 0, 1, nbr_left, nbr_right, rc)
    call MPI_Cart_shift(comm2d, 1, 1, nbr_bottom, nbr_top, rc)
    ! write(*,*) "I'm process ", rank2d, "My neighbor lrbt is", nbr_left, nbr_right, nbr_bottom, nbr_top
    call MPI_Cart_find_corners()

    ! construct the datatype for the exchange 
    ! row_x exchange in y direction(non-contiguous) -- top and bottom
    call MPI_Type_vector(nx, 1, 9, MPI_REAL8, f_row_x, rc)
    call MPI_Type_commit(f_row_x, rc)
    ! column_y exchange in x direction(non-contiguous) -- left and right
    call MPI_Type_vector(ny, 1, 9 * (nx+2), MPI_REAL8, f_column_y, rc)
    call MPI_Type_commit(f_column_y, rc)


    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------
#ifdef _OPENMP
    call OMP_set_num_threads(4)
    write(*,*) "Start OpenMP......"
    myMaxThreads = OMP_get_max_threads()
    write(*,*) "Running threads =",myMaxThreads
#endif
    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------

    call initial()
    call output()
    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( (errorU.GT.eps).AND.(itc.LE.itc_max) )

        call collision()

        call message_passing_sendrecv()

        call streaming()

        call bounceback()

        call bounceback_particle()

        call macro()

        call calForce()

        itc = itc+1

        call MPI_Barrier(comm2d, rc)
        ! write(*,*) itc
        
        if(MOD(itc,1000).EQ.0) then
            call check()
        endif
        !if(MOD(itc,2000).EQ.0) then
        !    call output_Tecplot()
        !endif

        ! call updateCenter()
        
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()

    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif
    
    if (rank == 0) then
        write(*,*) "Time (CPU) = ",finish-start, "s"
#ifdef _OPENMP
        write(*,*) "Time (OMP) = ", finish2-start2, "s"
#endif
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
    endif
    
    itc = itc+1
    call output()

    call freeall()

    call MPI_Finalize(rc)
    
    stop
end program main


subroutine decompose_1d(total_n, local_n, rank, num_process, i_start_global)
    implicit none
    integer, intent(in) :: total_n, rank, num_process
    integer, intent(out) :: local_n, i_start_global

    local_n = total_n / num_process

    if (rank < MOD(total_n, num_process)) then
        local_n = local_n + 1
    endif

    if (local_n > total_n / num_process) then ! --- 5 5 '5' 4 4 4
        i_start_global = local_n * rank
    else                    ! --- 5 5 5 4 '4' 4
        i_start_global = local_n * rank + mod(total_n, num_process)
    endif

end subroutine decompose_1d


subroutine MPI_Cart_find_corners()
    use mpi
    use commondata
    implicit none

    call MPI_Cart_shift_2d(1, 1, cnr_top_right)
    call MPI_Cart_shift_2d(1, -1, cnr_bottom_right)
    call MPI_Cart_shift_2d(-1, 1, cnr_top_left)
    call MPI_Cart_shift_2d(-1, -1, cnr_bottom_left)

end subroutine MPI_Cart_find_corners


subroutine MPI_Cart_shift_2d(idx0, idx1, corner_rank)
    use mpi
    use commondata
    implicit none
    integer, intent(in) :: idx0, idx1
    integer, intent(out) :: corner_rank 
    integer :: new_coords(0:1)

    new_coords(0) = coords(0) + idx0
    new_coords(1) = coords(1) + idx1
    if (new_coords(0) < 0 .OR. new_coords(0) > dims(0)-1) then  
        ! beyond the left/right boundary
        if (.not. periods(0)) then
            corner_rank = MPI_PROC_NULL
            return
        endif
    endif
    if (new_coords(1) < 0 .OR. new_coords(1) > dims(1)-1) then
        ! beyond the top/bottom boundary
        if (.not. periods(1)) then
            corner_rank = MPI_PROC_NULL
            return
        endif
    endif

    new_coords(0) = mod(new_coords(0) + dims(0), dims(0))
    new_coords(1) = mod(new_coords(1) + dims(1), dims(1))
    call MPI_Cart_rank(comm2d, new_coords, corner_rank, rc)

end subroutine MPI_Cart_shift_2d