!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model


! block arrangement is same as the data arrangement

program main
    use mpi
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start_time, end_time
    integer :: name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    if (rank == 0) then
        write(*,*) "Using ", num_process, "processers."
        write(*,*) "nx_block = ", nx_block
        write(*,*) "ny_block = ", ny_block
        write(*,*) "nz_block = ", nz_block
    endif

    if (nx_block * ny_block * nz_block /= num_process) then
        if (rank == 0) then
            write(*,*) "Please make sure nx_block * ny_block * nz_block = num_process."
        endif
        stop
    endif

    ! block index --- block start from index 1 but rank start from rank 0
    block_x = mod(rank, nx_block) + 1
    block_y = mod(rank, nx_block * ny_block) / nx_block + 1
    block_z = rank / (nx_block * ny_block) + 1

    write(*,"('Rank ', I3, ' has block index = (', I3, ',', I3, ',', I3, ').')") rank, block_x, block_y, block_z
    
    !!! zone divide
    nx = total_nx / nx_block
    ny = total_ny / ny_block
    nz = total_nz / nz_block
    ! handle non-divisible part
    if(block_x <= mod(total_nx, nx_block)) then
        nx = nx + 1
    endif
    if(block_y <= mod(total_ny, ny_block)) then
        ny = ny + 1
    endif
    if(block_z <= mod(total_nz, nz_block)) then
        nz = nz + 1
    endif

    write(*,"('Rank ', I3, ' has nx*ny*nz = ', I3, 'x', I3, 'x', I3)") rank, nx, ny, nz

    call initial()
    
    call output_field()
    call output_velocity()

    call CPU_TIME(start)

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()
    
    do while((errorU.GT.eps).AND.(itc.LE.itc_max))

        itc = itc+1

        call collision()

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
    call output_field()
    call output_velocity()
    
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