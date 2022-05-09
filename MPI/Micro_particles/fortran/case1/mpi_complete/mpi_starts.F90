subroutine mpi_starts()
    use mpi 
    use commondata
    implicit none
    integer :: name_len
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
    integer(KIND=MPI_ADDRESS_KIND) :: lower_bound, extent
    real(8) :: type_real_8

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    ! dims(0) = 1
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
    ! write(*,*) "I'm process ", rank2d, "My cooords are", coords(0), coords(1), "My cnr 1234 is", cnr_top_left, cnr_top_right, cnr_bottom_left, cnr_bottom_right

    ! construct the datatype for the exchange 
    ! fp_row_x exchange in y direction(non-contiguous) -- top and bottom
    call MPI_Type_vector(nx, 1, 9, MPI_REAL8, fp_row_x, rc)
    call MPI_Type_commit(fp_row_x, rc)
    ! fp_column_y exchange in x direction(non-contiguous) -- left and right
    call MPI_Type_vector(ny, 1, 9 * (nx+4), MPI_REAL8, fp_column_y, rc)
    call MPI_Type_commit(fp_column_y, rc)

    ! particles data
    lower_bound = 0
    call MPI_Type_contiguous(9, MPI_REAL8, type_f, rc)
    call MPI_Type_commit(type_f, rc)

    extent = 9 * sizeof(type_real_8)
    call MPI_Type_create_resized(type_f, lower_bound, extent, type_f_x, rc)     ! f_{1-9} or fp_{1-9} in row
    call MPI_Type_commit(type_f_x, rc)

    call MPI_Type_create_resized(MPI_REAL8, lower_bound, extent, type_fi_x, rc) ! f_i or fp_i in row
    call MPI_Type_commit(type_fi_x, rc)

    extent = 9 * (nx + 6) * sizeof(type_real_8)
    call MPI_Type_create_resized(type_f, lower_bound, extent, type_f_y, rc)     ! f_{1-9} in column
    call MPI_Type_commit(type_f_y, rc)

    call MPI_Type_create_resized(MPI_REAL8, lower_bound, extent, type_fi_y, rc) ! f_i in column
    call MPI_Type_commit(type_fi_y, rc)

    extent = 9 * (nx + 4) * sizeof(type_real_8)
    call MPI_Type_create_resized(type_f, lower_bound, extent, type_fp_y, rc)     ! fp_{1-9} in column
    call MPI_Type_commit(type_fp_y, rc)

    call MPI_Type_create_resized(MPI_REAL8, lower_bound, extent, type_fpi_y, rc) ! fp_i in column
    call MPI_Type_commit(type_fpi_y, rc)

    ! square data
    call MPI_Type_vector(3, 3, nx+6, type_f, type_square3_f, rc)      ! f in square 3
    call MPI_Type_commit(type_square3_f, rc)

    call MPI_Type_vector(2, 2, nx+6, type_f, type_square2_f, rc)      ! f in square 2
    call MPI_Type_commit(type_square2_f, rc)

    call MPI_Type_vector(2, 2, nx+4, type_f, type_square2_fp, rc)      ! fp in square 2
    call MPI_Type_commit(type_square2_fp, rc)

    call MPI_Type_vector(2, 2, nx+4, type_fi_x, type_square2_fpi, rc)  ! fp_i in square 2
    call MPI_Type_commit(type_square2_fpi, rc)

end subroutine mpi_starts




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
