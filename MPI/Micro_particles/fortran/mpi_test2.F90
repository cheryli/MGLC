program main
    use mpi
    implicit none
    integer :: name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
    integer(KIND=MPI_ADDRESS_KIND) :: lower_bound, extent
    real(8) :: type_real_8
    real(8), allocatable :: array(:,:)
    integer :: rc, rank, num_process
    integer :: i, j, nx
    real(8) :: force(3), total_force(3)

    force = 0.0d0
    total_force = 0.0d0
    
    call MPI_Init(rc)
    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)


    if (rank == 0) then
        force(1) = 1.4d0
        force(2) = 2.5d0
        force(3) = 3.6d0
    endif

    if (rank == 1) then
        force(1) = 4.1d0
        force(2) = 5.2d0
        force(3) = 6.3d0
    endif

    call MPI_ALLreduce(force, total_force, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

    if (rank == 0) then
        write(*,*) total_force(1), total_force(2), total_force(3)
    endif

    call MPI_Finalize(rc)

end program main