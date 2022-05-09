module commondata
    implicit none
    
    integer, parameter :: total_nx = 65, total_ny = 65, total_nz = 65
    integer :: nx, ny, nz
    real(8), parameter :: Reynolds=1000.0d0
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: U0=0.1d0
    real(8), parameter :: tauf=U0*dble(total_nx)/Reynolds*3.0d0+0.5d0
    
    integer :: itc
    integer, parameter :: itc_max=INT(50000000)
    
    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    
    real(8) :: xp(0:total_nx+1), yp(0:total_ny+1), zp(0:total_nz+1)
    real(8), allocatable :: u(:, :, :), v(:, :, :), w(:, :, :)
    real(8), allocatable :: rho(:, :, :)
    real(8), allocatable :: up(:, :, :), vp(:, :, :), wp(:, :, :)
    
    real(8), allocatable :: f(:, :, :, :), f_post(:, :, :, :)
    
    integer :: index
    ! real(8), parameter :: omega(0:18) = (/ 1.0d0/3.0d0, & 
    !         1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, &
    !         1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
    !         1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
    real(8), parameter :: omega(0:18) = (/ 1.0d0/3.0d0, &
                                        (1.0d0/18.0d0, index = 1, 6), &
                                        (1.0d0/36.0d0, index = 7, 18) /) 
    integer, parameter :: ex(0:18) = (/ 0,  &
                                        1, -1,  0,  0,  0,  0, &
                                        1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0 /)
    integer, parameter :: ey(0:18) = (/ 0, &
                                        0,  0,  1, -1,  0,  0, &
                                        1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1 /)  
    integer, parameter :: ez(0:18) = (/ 0, &
                                        0,  0,  0,  0,  1, -1, &
                                        0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    

    ! mpi data
    integer :: rc, rank, num_process 
    integer :: dims(0:2) = (/ 0, 0, 0 /), coords(0:2)
    logical :: periods(0:2) = (/ .false., .false., .false. /)
    integer :: comm3d, rank3d
    integer :: nbr_surface(1:6)
    integer :: nbr_line(7:18) 
    integer :: surface_x, surface_y, surface_z  ! surface data perpendicular to axis x, y, z
    integer :: line_x, line_y, line_z       ! line data parallel to axis x, y, z

end module commondata


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


    !!! ---- Decomposition the domain
    ! dims(0) = 2
    ! dims(1) = 2
    ! dims(2) = 3
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
    
    call output()
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
    call output()
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


subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(8) :: un(0:18)
    real(8) :: us2

    itc = 0
    errorU = 100.0d0

    if (rank == 0) then
        write(*,*) "nx=",nx,", ny=",ny, ",  nz = ",nz
        write(*,*) "Reynolds=",real(Reynolds)
        write(*,*) "U0=",real(U0),",    tauf=",real(tauf)
    endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i=1,total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j=1,total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(total_nz+1) = dble(total_nz)
    do k=1,total_nz
        zp(k) = dble(k)-0.5d0
    enddo

    
    allocate (u(nx, ny, nz))
    allocate (v(nx, ny, nz))
    allocate (w(nx, ny, nz))
    allocate (rho(nx, ny, nz))
    allocate (up(nx, ny, nz))
    allocate (vp(nx, ny, nz))
    allocate (wp(nx, ny, nz))
    
    allocate (f(0:18, nx, ny, nz))
    allocate (f_post(0:18, 0:nx+1, 0:ny+1, 0:nz+1))

    rho = rho0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    up = 0.0d0
    vp = 0.0d0
    wp = 0.0d0

    ! top bpundary
    if (coords(2) == dims(2) - 1) then
        do j = 1, ny
            do i = 1, nx
                u(i, j, nz) = U0
            enddo
        enddo
    endif

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                us2 = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                do alpha = 0, 18
                    un(alpha) = u(i,j,k)*ex(alpha) + v(i,j,k)*ey(alpha) + w(i,j,k)*ez(alpha)
                    f(alpha,i,j,k) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine initial


subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !---------------------------
    real(8) :: s(0:18)
    real(8) :: m(0:18)
    real(8) :: m_post(0:18)
    real(8) :: meq(0:18)
    !---------------------------
    ! real(8) :: us2, un, feq
    !---------------------------

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
    ! loop body, re-indented for readability 
        
    m(0) = f(0,i,j,k) &
        + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(1) = -30.0d0 * f(0,i,j,k) &
        - 11.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + 8.0d0 * (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k)) &
        + 8.0d0 * (f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(2) = 12.0d0 * f(0,i,j,k) &
        - 4.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(3) = f(1,i,j,k) - f(2,i,j,k) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)
    m(4) = -4.0d0 * (f(1,i,j,k) - f(2,i,j,k)) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)       
    m(5) = f(3,i,j,k) - f(4,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(6) = -4.0d0 * (f(3,i,j,k) - f(4,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(7) = f(5,i,j,k) - f(6,i,j,k) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(8) = -4.0d0 * (f(5,i,j,k) - f(6,i,j,k)) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(9) = 2.0d0 * (f(1,i,j,k) + f(2,i,j,k)) - f(3,i,j,k) - f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(10) = -4.0d0 * (f(1,i,j,k) + f(2,i,j,k)) + 2.0d0 * (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(11) = f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(12) = -2.0d0 * (f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(13) = f(7,i,j,k) - f(8,i,j,k) - f(9,i,j,k) + f(10,i,j,k)
    m(14) = f(15,i,j,k) - f(16,i,j,k) - f(17,i,j,k) + f(18,i,j,k)
    m(15) = f(11,i,j,k) - f(12,i,j,k) - f(13,i,j,k) + f(14,i,j,k)
    m(16) = f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) - f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) + f(14,i,j,k)
    m(17) = -f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(18) = f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) - f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
        
        
        meq(0) = rho(i,j,k) 
        meq(1) = rho(i,j,k) * (-11.0d0 + 19.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(2) = rho(i,j,k) * (3.0d0 - 11.0d0/2.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(3) = rho(i,j,k) * u(i,j,k)
        meq(4) = -2.0d0/3.0d0 * rho(i,j,k) * u(i,j,k)
        meq(5) = rho(i,j,k) * v(i,j,k)
        meq(6) = -2.0d0/3.0d0 * rho(i,j,k) * v(i,j,k)
        meq(7) = rho(i,j,k) * w(i,j,k)
        meq(8) = -2.0d0/3.0d0 * rho(i,j,k) * w(i,j,k)
        meq(9) = rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(10) = -1.0d0/2.0d0 * rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(11) = rho(i,j,k) * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(12) = -1.0d0/2.0d0 * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(13) = rho(i,j,k) * u(i,j,k) * v(i,j,k)
        meq(14) = rho(i,j,k) * v(i,j,k) * w(i,j,k)
        meq(15) = rho(i,j,k) * u(i,j,k) * w(i,j,k)
        meq(16) = 0.0d0
        meq(17) = 0.0d0
        meq(18) = 0.0d0


        s(0) = 0.0d0    !! s_{\rho}
        s(1) = Snu      !! s_{e}
        s(2) = Snu      !! s_{epsilon}
        s(3) = 0.0d0    !! s_{j}
        s(4) = Sq       !! s_{q}
        s(5) = 0.0d0    !! s_{j}
        s(6) = Sq       !! s_{q}
        s(7) = 0.0d0    !! s_{j}
        s(8) = Sq       !! s_{q}
        s(9) = Snu      !! s_{\nu}
        s(10) = Snu     !! s_{\pi}
        s(11) = Snu     !! s_{\nu}
        s(12) = Snu     !! s_{\pi}
        s(13) = Snu     !! s_{nu}
        s(14) = Snu     !! s_{nu}
        s(15) = Snu     !! s_{nu}
        s(16) = Sq      !! s_{m}
        s(17) = Sq      !! s_{m}
        s(18) = Sq      !! s_{m}
        
        do alpha=0,18
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
        enddo
 
    f_post(0,i,j,k) = m_post(0)/19.0d0 - 5.0d0/399.0d0*m_post(1) + m_post(2)/21.0d0

    ! ------------
    f_post(1,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(3)/10.0d0 &
        - m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(3)/10.0d0 &
        + m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0
    
    ! ------------
    f_post(3,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(5)/10.0d0 &
        - m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(4,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(5)/10.0d0 &
        + m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(5,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(7)/10.0d0 &
        - m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(7)/10.0d0 &
        + m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    ! ------------
    f_post(7,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 + m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(8,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 - m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(9,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 + m_post(16)/8.0d0 + m_post(17)/8.0d0

    f_post(10,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 - m_post(16)/8.0d0 + m_post(17)/8.0d0

    ! -----------
    f_post(11,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 - m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(12,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 + m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(13,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 - m_post(16)/8.0d0 - m_post(18)/8.0d0
    
    f_post(14,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 + m_post(16)/8.0d0 - m_post(18)/8.0d0

    ! -----------
    f_post(15,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 + m_post(17)/8.0d0 - m_post(18)/8.0d0
    
    f_post(16,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 - m_post(17)/8.0d0 - m_post(18)/8.0d0

    f_post(17,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 + m_post(17)/8.0d0 + m_post(18)/8.0d0

    f_post(18,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 - m_post(17)/8.0d0 + m_post(18)/8.0d0


                ! us2 = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                ! do alpha = 0, 18
                !     un = u(i,j,k)*ex(alpha) + v(i,j,k)*ey(alpha) + w(i,j,k)*ez(alpha)
                !     feq = rho(i,j,k) * omega(alpha) * (1.0d0 + 3.0d0*un + 4.5d0*un*un - 1.5d0*us2)
                    
                !     f_post(alpha,i,j,k) = f(alpha,i,j,k) - Snu*(f(alpha,i,j,k) - feq)
                ! enddo
    
            enddo
        enddo
    enddo

    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp   
    integer :: alpha
    
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha=0,18
                    ip = i - ex(alpha)
                    jp = j - ey(alpha)
                    kp = k - ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

    if (coords(0) == 0) then
        do k = 1, nz
            do j = 1, ny
                ! (i = 1)
                f(1,1,j,k) = f_post(2,1,j,k)
                f(7,1,j,k) = f_post(10,1,j,k)
                f(9,1,j,k) = f_post(8,1,j,k)
                f(11,1,j,k) = f_post(14,1,j,k)
                f(13,1,j,k) = f_post(12,1,j,k)
            enddo
        enddo
    endif

    if (coords(0) == dims(0)-1) then
        do k = 1, nz
            do j = 1, ny
                !(i = nx)
                f(2,nx,j,k) = f_post(1,nx,j,k)
                f(10,nx,j,k) = f_post(7,nx,j,k)
                f(8,nx,j,k) = f_post(9,nx,j,k)
                f(14,nx,j,k) = f_post(11,nx,j,k)
                f(12,nx,j,k) = f_post(13,nx,j,k)
            enddo
        enddo
    endif


    if (coords(1) == 0) then
        do k = 1, nz
            do i = 1, nx
                !(j = 1)
                f(3,i,1,k) = f_post(4,i,1,k)
                f(7,i,1,k) = f_post(10,i,1,k)
                f(8,i,1,k) = f_post(9,i,1,k)
                f(15,i,1,k) = f_post(18,i,1,k)
                f(17,i,1,k) = f_post(16,i,1,k)
            enddo
        enddo
    endif

    if (coords(1) == dims(1)-1) then
        do k = 1, nz
            do i = 1, nx
                !(j = ny)
                f(4,i,ny,k) = f_post(3,i,ny,k)
                f(10,i,ny,k) = f_post(7,i,ny,k)
                f(9,i,ny,k) = f_post(8,i,ny,k)
                f(18,i,ny,k) = f_post(15,i,ny,k)
                f(16,i,ny,k) = f_post(17,i,ny,k)
            enddo
        enddo
    endif

    if (coords(2) == 0) then
        do j = 1, ny
            do i = 1, nx
                ! (k = 1)
                f(5,i,j,1) = f_post(6,i,j,1)
                f(11,i,j,1) = f_post(14,i,j,1)
                f(12,i,j,1) = f_post(13,i,j,1)
                f(15,i,j,1) = f_post(18,i,j,1)
                f(16,i,j,1) = f_post(17,i,j,1)
            enddo
        enddo
    endif

    if (coords(2) == dims(2)-1) then
        do j = 1, ny
            do i = 1, nx
                ! (k = nz)
                f(6,i,j,nz) = f_post(5,i,j,nz)
                f(14,i,j,nz) = f_post(11,i,j,nz) - rho(i,j,nz) / 6.0d0 * (U0)
                f(13,i,j,nz) = f_post(12,i,j,nz) - rho(i,j,nz) / 6.0d0 * (-U0)
                f(18,i,j,nz) = f_post(15,i,j,nz)
                f(17,i,j,nz) = f_post(16,i,j,nz)
            enddo
        enddo
    endif
    
    return
end subroutine bounceback



subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k, alpha

    rho = 0.0d0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 18
                    rho(i,j,k) = rho(i,j,k) + f(alpha,i,j,k)
                    u(i,j,k) = u(i,j,k) + f(alpha,i,j,k) * ex(alpha)
                    v(i,j,k) = v(i,j,k) + f(alpha,i,j,k) * ey(alpha)
                    w(i,j,k) = w(i,j,k) + f(alpha,i,j,k) * ez(alpha)
                enddo

                u(i,j,k) = u(i,j,k)/rho(i,j,k)
                v(i,j,k) = v(i,j,k)/rho(i,j,k)
                w(i,j,k) = w(i,j,k)/rho(i,j,k)
            enddo
        enddo
    enddo

    return
end subroutine macro


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

    call MPI_Barrier(comm3d, rc)

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm3d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm3d, rc)

    errorU = sqrt(total_error1)/sqrt(total_error2)

    if (rank3d == 0) then
        write(*,*) itc,' ',errorU
    endif

    return
end subroutine check



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




subroutine message_passing_sendrecv()
    use mpi
    use commondata
    implicit none
    integer :: tag(5), idx

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)s
    tag = (/ 1, 7, 9, 11, 13 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), nx, 1, 1), 1, surface_x, nbr_surface(1), tag(idx), &
                    f_post(tag(idx), 0, 1, 1), 1, surface_x, nbr_surface(2), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (i--)
    tag = (/ 2, 8, 10 ,12, 14 /) 
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_x, nbr_surface(2), tag(idx), &
                    f_post(tag(idx), nx+1, 1, 1), 1, surface_x, nbr_surface(1), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    tag = (/ 3, 7, 8, 15, 17 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, ny, 1), 1, surface_y, nbr_surface(3), tag(idx), &
                    f_post(tag(idx), 1, 0, 1), 1, surface_y, nbr_surface(4), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (j--)
    tag = (/ 4, 9, 10, 16, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_y, nbr_surface(4), tag(idx), &
                    f_post(tag(idx), 1, ny+1, 1), 1, surface_y, nbr_surface(3), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    tag = (/ 5, 11, 12, 15, 16 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, nz), 1, surface_z, nbr_surface(5), tag(idx), &
                    f_post(tag(idx), 1, 1, 0), 1, surface_z, nbr_surface(6), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (k--)
    tag = (/ 6, 13, 14, 17, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, surface_z, nbr_surface(6), tag(idx), &
                    f_post(tag(idx), 1, 1, nz+1), 1, surface_z, nbr_surface(5), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! /* ----------------------------- line data ----------------------------- */
    ! ------------ exchange message perpendicular to z ----------------
    ! message passing to (i++, j++) (7)
    call MPI_Sendrecv(f_post(7, nx, ny, 1), 1, line_z, nbr_line(7), 7, &
                    f_post(7, 0, 0, 1), 1, line_z, nbr_line(10), 7, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j--) (10)
    call MPI_Sendrecv(f_post(10, 1, 1, 1), 1, line_z, nbr_line(10), 10, &
                    f_post(10, nx+1, ny+1, 1), 1, line_z, nbr_line(7), 10, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, j--) (9)
    call MPI_Sendrecv(f_post(9, nx, 1, 1), 1, line_z, nbr_line(9), 9, &
                    f_post(9, 0, ny+1, 1), 1, line_z, nbr_line(8), 9, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j++) (8)
    call MPI_Sendrecv(f_post(8, 1, ny, 1), 1, line_z, nbr_line(8), 8, &
                    f_post(8, nx+1, 0, 1), 1, line_z, nbr_line(9), 8, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to y ----------------
    ! message passing to (i++, k++) (11)
    call MPI_Sendrecv(f_post(11, nx, 1, nz), 1, line_y, nbr_line(11), 11, &
                    f_post(11, 0, 1, 0), 1, line_y, nbr_line(14), 11, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k--) (14)
    call MPI_Sendrecv(f_post(14, 1, 1, 1), 1, line_y, nbr_line(14), 14, &
                    f_post(14, nx+1, 1, nz+1), 1, line_y, nbr_line(11), 14, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, k--) (13)
    call MPI_Sendrecv(f_post(13, nx, 1, 1), 1, line_y, nbr_line(13), 13, &
                    f_post(13, 0, 1, nz+1), 1, line_y, nbr_line(12), 13, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k++) (12)
    call MPI_Sendrecv(f_post(12, 1, 1, nz), 1, line_y, nbr_line(12), 12, &
                    f_post(12, nx+1, 1, 0), 1, line_y, nbr_line(13), 12, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to x ----------------
    ! message passing to (j++, k++) (15)
    call MPI_Sendrecv(f_post(15, 1, ny, nz), 1, line_x, nbr_line(15), 15, &
                    f_post(15, 1, 0, 0), 1, line_x, nbr_line(18), 15, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k--) (18)
    call MPI_Sendrecv(f_post(18, 1, 1, 1), 1, line_x, nbr_line(18), 18, &
                    f_post(18, 1, ny+1, nz+1), 1, line_x, nbr_line(15), 18, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j++, k--) (17)
    call MPI_Sendrecv(f_post(17, 1, ny, 1), 1, line_x, nbr_line(17), 17, &
                    f_post(17, 1, 0, nz+1), 1, line_x, nbr_line(16), 17, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k++) (16)
    call MPI_Sendrecv(f_post(16, 1, 1, nz), 1, line_x, nbr_line(16), 16, &
                    f_post(16, 1, ny+1, 0), 1, line_x, nbr_line(17), 16, &
                    comm3d, MPI_STATUS_IGNORE, rc)
    

    ! don't need to echange point data


end subroutine message_passing_sendrecv


subroutine output()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    integer :: p_rank, num(0:2), new_coords(0:2), dx = 0, dy = 0, dz = 0
    real(8), allocatable :: total_u(:, :, :), total_v(:, :, :), total_w(:, :, :), total_rho(:, :, :)
    real(8), allocatable :: tmp_u(:, :, :), tmp_v(:, :, :), tmp_w(:, :, :), tmp_rho(:, :, :)
    
    ! use rank 0 to receive data and output the results

    if (rank3d > 0) then  !!! ----  rank != 0 send data
        ! collect the rank information
        num(0) = nx
        num(1) = ny
        num(2) = nz
        ! send to rank 0
        call MPI_Send(num, 3, MPI_INTEGER, 0, 0, comm3d, rc)    ! rank information
        call MPI_Send(u, nx*ny*nz, MPI_REAL8, 0, 1, comm3d, rc)
        call MPI_Send(v, nx*ny*nz, MPI_REAL8, 0, 2, comm3d, rc)
        call MPI_Send(w, nx*ny*nz, MPI_REAL8, 0, 3, comm3d, rc)
        call MPI_Send(rho, nx*ny*nz, MPI_REAL8, 0, 4, comm3d, rc)
    else    
        !!! ---- rank 0 collect data
        ! allocate array
        allocate(total_u(total_nx, total_ny, total_nz))
        allocate(total_v(total_nx, total_ny, total_nz))
        allocate(total_w(total_nx, total_ny, total_nz))
        allocate(total_rho(total_nx, total_ny, total_nz))

        ! determine the origin
        if (nx > total_nx / dims(0)) then ! --- 5 5 '5' 4 4 4
            dx = nx * coords(0)
        else                    ! --- 5 5 5 4 '4' 4
            dx = nx * coords(0) + mod(total_nx, dims(0))
        endif

        if (ny > total_ny / dims(1)) then ! --- 5 5 '5' 4 4 4
            dy = ny * coords(1)
        else                    ! --- 5 5 5 4 '4' 4
            dy = ny * coords(1) + mod(total_ny, dims(1))
        endif

        if (nz > total_nz / dims(2)) then
            dz = nz * coords(2)
        else
            dz = nz * coords(2) + mod(total_nz, dims(2))
        endif

        ! collect data from rank 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    total_u(dx + i, dy + j, dz + k) = u(i, j, k)
                    total_v(dx + i, dy + j, dz + k) = v(i, j, k)
                    total_w(dx + i, dy + j, dz + k) = w(i, j, k)
                    total_rho(dx + i, dy + j, dz + k) = rho(i, j, k)
                enddo
            enddo
        enddo

        ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) * dims(2) - 1

            call MPI_Cart_coords(comm3d, p_rank, 3, new_coords, rc)

            ! receive the rank information
            call MPI_Recv(num, 3, MPI_INTEGER, p_rank, 0, comm3d, MPI_STATUS_IGNORE, rc)

            ! creat buffer
            allocate(tmp_u(num(0), num(1), num(2)))
            allocate(tmp_v(num(0), num(1), num(2)))
            allocate(tmp_w(num(0), num(1), num(2)))
            allocate(tmp_rho(num(0), num(1), num(2)))

            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 1, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 2, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_w, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 3, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 4, comm3d, MPI_STATUS_IGNORE, rc)


            ! determine the origin
            if (num(0) > total_nx / dims(0)) then ! --- 5 5 '5' 4 4 4
                dx = num(0) * new_coords(0)
            else                    ! --- 5 5 5 4 '4' 4
                dx = num(0) * new_coords(0) + mod(total_nx, dims(0))
            endif

            if (num(1) > total_ny / dims(1)) then ! --- 5 5 '5' 4 4 4
                dy = num(1) * new_coords(1)
            else                    ! --- 5 5 5 4 '4' 4
                dy = num(1) * new_coords(1) + mod(total_ny, dims(1))
            endif

            if (num(2) > total_nz / dims(2)) then
                dz = num(2) * new_coords(2)
            else
                dz = num(2) * new_coords(2) + mod(total_nz, dims(2))
            endif

            ! assign data
            do k = 1, num(2)
                do j = 1, num(1)
                    do i = 1, num(0)
                        total_u(i + dx, j + dy, k + dz) = tmp_u(i, j, k)
                        total_v(i + dx, j + dy, k + dz) = tmp_v(i, j, k)
                        total_w(i + dx, j + dy, k + dz) = tmp_w(i, j, k)
                        total_rho(i + dx, j + dy, k + dz) = tmp_rho(i, j, k)
                    enddo
                enddo
            enddo

            ! de-allocate buffer arrays
            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_w)
            deallocate(tmp_rho)
        enddo

        ! after collect total_* data, then output
        ! call output_ASCII(xp, yp, zp, total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        ! call output_binary(total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        call output_Tecplot(xp, yp, zp, total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        call getVelocity(xp, yp, zp, total_u, total_v, total_w, total_nx, total_ny, total_nz, U0, itc)

        ! de-allocate total arrays
        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_w)
        deallocate(total_rho)
    endif

end subroutine output



subroutine output_ASCII(xp, yp, zp, u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=77,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
        write(77,*) 'TITLE="Lid Driven Cavity"'
        write(77,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "Pressure" '
        write(77,101) nx, ny, nz
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(77,100) xp(i), yp(j), zp(k), u(i,j,k), v(i,j,k), w(i,j,k), rho(i,j,k)/3.0d0
                enddo
            enddo
        enddo
100     format(1x,3(e11.4,' '),10(e13.6,' '))
101     format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'K=',1x,i5,1x,'F=POINT')
    close(77)

    return
end subroutine output_ASCII


!!!c--------------------------------
!!!c--------------------------------
subroutine output_Tecplot(xp, yp, zp, u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    
    write(B2,'(i9.9)') itc

    open(41,file='MRTcavity-'//B2//'.plt',form='binary')
    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"

    !c--Integer value of 1
    write(41) 1

    Title="MyFirst"
    call dumpstring(title)

    !c-- Number of variables in this data file (here 5 variables)
    write(41) 7

    !c-- Variable names.
    V1 = 'X'
    call dumpstring(V1)
    V2 = 'Y'
    call dumpstring(V2)
    V3 = 'Z'
    call dumpstring(V3)
    V4 = 'U'
    call dumpstring(V4)
    V5 = 'V'
    call dumpstring(V5)
    V6 = 'W'
    call dumpstring(V6)
    V7 = 'Pressure'
    call dumpstring(V7)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName="ZONE 001"
    call dumpstring(zoneName)

    !---------Zone Color
    write(41) -1

    !---------ZoneType
    write(41) 0

    !---------DataPacking 0=Block, 1=Point
    write(41) 1

    !---------Specify Var Location. 0 = Do not specify, all data
    !---------is located at the nodes. 1 = Specify
    write(41) 0

    !---------Number of user defined face neighbor connections
    ! (value >= 0)
    write(41) 0

    !---------IMax,JMax,KMax
    write(41) nx
    write(41) ny
    write(41) nz

    !-----------1=Auxiliary name/value pair to follow
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker

    !----zone ------------------------------------------------------------
    write(41) zoneMarker

    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(zp(k))
                write(41) real(u(i,j,k))
                write(41) real(v(i,j,k))
                write(41) real(w(i,j,k))
                write(41) real(rho(i,j,k)/3.0d0)
            enddo
        enddo
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
    end subroutine output_Tecplot
!!!c--------------------------------
    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer :: stringLength
    integer :: ii
    integer :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring
!!!c--------------------------------
!!!c-------------------------------- 


subroutine getVelocity(xp, yp, zp, u, v, w, nx, ny, nz, U0, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz)
    real(8), intent(in) :: U0
    integer :: i, j, k
    integer :: nxHalf, nyHalf, nzHalf
    character(len=9) :: B2

    write(B2,'(i9.9)') itc

    nxHalf = (nx - 1) / 2 + 1
    nyHalf = (ny - 1) / 2 + 1
    nzHalf = (nz - 1) / 2 + 1

    open(unit=02,file='./u-z_'//B2//'.dat',status='unknown')
    do k=1,nz
        write(02,*) u(nxHalf, nyHalf, k)/U0, zp(k)/dble(nz)
    enddo
    close(02)

    open(unit=03,file='./x-w_'//B2//'.dat',status='unknown')
    do i=1,nx
        write(03,*) xp(i)/dble(nx), w(i, nyHalf, nzHalf)/U0
    enddo
    close(03)

    return
end subroutine getVelocity


subroutine output_binary(u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)
    
    open(unit=01,file='MRTcavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((v(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((rho(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    close(01)

    return
end subroutine output_binary