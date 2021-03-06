!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


! block arrangement is same as the data arrangement

program main
    use mpi
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start_time, end_time
    integer :: name_len, stride
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)


    !!! decomposition the domain 
    call MPI_Dims_create(num_process, 3, dims, rc)
    call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, .true., comm3d, rc)
    if(rank == 0) then
        write(*,*) "Using ", num_process, "processers."
        write(*,*) "dimens is x*y*z = ", dims(0), "x", dims(1), "x", dims(2)
    endif

    ! get my new rank in decomposition
    call MPI_Comm_rank(comm3d, rank3d, rc)
    ! write(*,*) "process ", rank3d, " of total ", num_process, "is alive."

    ! determain sub-domain size
    call MPI_Cart_get(comm3d, 3, dims, periods, coords, rc)
    call decompose_1d(total_nx, nx, coords(0), dims(0))
    call decompose_1d(total_ny, ny, coords(1), dims(1))
    call decompose_1d(total_nz, nz, coords(2), dims(2))
    ! write(*,*) "coords = ", coords(1), coords(2), coords(3), "nx*ny = ", nx, ny, nz

    ! get the neighbors
    call MPI_Cart_shift(comm3d, 0, 1, nbr_surface(2), nbr_surface(1), rc)
    call MPI_Cart_shift(comm3d, 1, 1, nbr_surface(4), nbr_surface(3), rc)
    call MPI_Cart_shift(comm3d, 2, 1, nbr_surface(6), nbr_surface(5), rc)
    ! write(*,*) "I'm process ", rank3d, "My neighbor surfaces are", nbr_surface(1), nbr_surface(2), nbr_surface(3), nbr_surface(4), nbr_surface(5), nbr_surface(6)
    call MPI_Cart_find_corners()


    ! construct the datatype for the exchange 
    ! line data parallel x direction (don't exchange in x direction)
    call MPI_Type_vector(nx, 1, 19, MPI_REAL8, line_x, rc)
    call MPI_Type_commit(line_x, rc)
    ! line data parallel y direction (don't exchange in y direction)
    call MPI_Type_vector(ny, 1, 19 * (nx+2), MPI_REAL8, line_y, rc)
    call MPI_Type_commit(line_y, rc)
    ! line data parallel z direction (don't exchange in z direction)
    call MPI_Type_vector(nz, 1, 19 * (nx+2) * (ny+2), MPI_REAL8, line_z, rc)
    call MPI_Type_commit(line_z, rc)

    ! surface data exchange in x direction -- line_y * nz, stride = 19*(nx+2)*(ny+2)
    stride = 19 * (nx+2) * (ny+2) 
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), line_y, surface_x, rc)
    call MPI_Type_commit(surface_x, rc)
    ! surface data exchange in y direction -- line_x * nz, stride = 19*(nx+2)*(ny+2)
    stride = 19 * (nx+2) * (ny+2)
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), line_x, surface_y, rc)
    call MPI_Type_commit(surface_y, rc)
    ! surface data exchange in z direction -- line_x * ny, stride = 19*(nx+2)
    stride = 19 * (nx+2)
    call MPI_Type_create_hvector(ny, 1, stride * sizeof(start), line_x, surface_z, rc)
    call MPI_Type_commit(surface_z, rc)


    call initial()
    
    ! call output_field()
    ! call output_velocity()

    call CPU_TIME(start)

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()
    
    do while((errorU.GT.eps).AND.(itc.LE.itc_max))

        itc = itc+1

        call collision()

        call message_passing_sendrecv()

        call streaming()

        call bounceback()

        call macro()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()
    
    call CPU_TIME(finish)

    if (rank == 0) then
        write(*,*) "Time (CPU) = ", real(finish-start), "s"
        write(*,*) "Wall time = ", real(end_time - start_time), "s"
    endif

    itc = itc+1
    ! call output_field()
    ! call output_velocity()
    
    if (rank == 0) then
        write(*,*) "Deallocate Array..."
    endif

    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(wp)
    deallocate(f)
    deallocate(f_post)

     if (rank == 0) then
        write(*,*) "    "
    
        write(*,*) "Successfully: DNS completed!"
    endif

    call MPI_Finalize(rc)

end program main


subroutine decompose_1d(total_n, local_n, rank, num_process)
    implicit none
    integer, intent(in) :: total_n, rank, num_process
    integer, intent(out) :: local_n

    local_n = total_n / num_process

    if (rank < MOD(total_n, num_process)) then
        local_n = local_n + 1
    endif

end subroutine decompose_1d


subroutine MPI_Cart_find_corners()
    use mpi
    use commondata
    implicit none
    integer :: alpha

    ! name the neighbor line by the discrete velocity
    do alpha = 7, 18
        call MPI_Cart_shift_3d(ex(alpha), ey(alpha), ez(alpha), nbr_line(alpha))
    enddo

end subroutine MPI_Cart_find_corners


subroutine MPI_Cart_shift_3d(idx0, idx1, idx2, corner_rank)
    use mpi
    use commondata
    implicit none
    integer, intent(in) :: idx0, idx1, idx2
    integer, intent(out) :: corner_rank 
    integer :: new_coords(0:2)

    new_coords(0) = coords(0) + idx0
    new_coords(1) = coords(1) + idx1
    new_coords(2) = coords(2) + idx2

    if (new_coords(0) < 0 .OR. new_coords(0) > dims(0)-1) then  
        ! beyond the dims_0 boundary
        if (periods(0) == .false.) then
            corner_rank = MPI_PROC_NULL
            return
        else
            new_coords(0) = mod(new_coords(0) + dims(0), dims(0))
        endif
    else if (new_coords(1) < 0 .OR. new_coords(1) > dims(1)-1) then
        ! beyond the dims_1 boundary
        if (periods(1) == .false.) then
            corner_rank = MPI_PROC_NULL
            return
        else
            new_coords(1) = mod(new_coords(1) + dims(1), dims(1))
        endif
    else if (new_coords(2) < 0 .OR. new_coords(2) > dims(2)-1) then
        ! beyond the dims_2 boundary
        if (periods(2) == .false.) then
            corner_rank = MPI_PROC_NULL
            return
        else
            new_coords(2) = mod(new_coords(2) + dims(2), dims(2))
        endif
    endif

    call MPI_Cart_rank(comm3d, new_coords, corner_rank, rc)

end subroutine MPI_Cart_shift_3d