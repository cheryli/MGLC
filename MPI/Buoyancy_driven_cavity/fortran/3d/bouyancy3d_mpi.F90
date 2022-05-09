!!!    This program sloves Buoyancy Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with MRT-LBE model

!!==velocity B.C.==
#define noslipWalls
!!==velocity B.C.==

!~~temperature B.C. (for RB convection)~~
!#define RBconvection
!#define BackFrontWallsAdiabatic
!#define LeftRightWallsAdiabatic
!#define TopBottomPlatesConstT
!~~temperature B.C.~~

!~~temperature B.C. (for cavity flow benchmark)~~
#define benchmarkCavity
#define BackFrontWallsAdiabatic
#define LeftRightWallsConstT
#define TopBottomPlatesAdiabatic
!~~temperature B.C.~~

module commondata
    implicit none
    
    integer, parameter :: loadInitField=0
    integer, parameter :: flowReversal=0  !! value 0: do NOT simulate reversal; value 1: do simulate reversal

    integer, parameter :: total_nx = 51, total_ny = total_nx, total_nz = total_nx
    integer :: nx, ny, nz
    integer, parameter :: nxHalf = (total_nx-1)/2+1, nyHalf = (total_ny-1)/2+1, nzHalf = (total_nz-1)/2+1
    
    real(kind=8), parameter :: Rayleigh=1e6
    real(kind=8), parameter :: Prandtl=0.71d0
    real(kind=8), parameter :: Mach=0.1d0

    real(kind=8), parameter :: tauf=0.5d0+Mach*dble(total_nz)*DSQRT(3.0d0*Prandtl/Rayleigh)
    real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
    real(kind=8), parameter :: viscosity=(tauf-0.5d0)/3.0d0
    real(kind=8), parameter :: diffusivity=viscosity/Prandtl
    
    real(kind=8), parameter :: Ekman=0.001d0
    real(kind=8), parameter :: omegaRatating=viscosity/2.0d0/Ekman/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: paraA=42.0d0*dsqrt(3.0d0)*diffusivity-6.0d0
    real(kind=8), parameter :: gBeta1=Rayleigh*viscosity*diffusivity/dble(total_nz)
    real(kind=8), parameter :: gBeta=gBeta1/dble(total_nz*total_nz)
    
    real(kind=8), parameter :: timeUnit=dsqrt(dble(total_nz)/gBeta)  !!dble(ny*ny)/diffusivity
    integer, parameter :: dimensionlessTimeMax=1000
    integer, parameter :: flowReversalTime=20000
    integer :: itc
    integer, parameter :: itc_max=INT(dimensionlessTimeMax*timeUnit)
    
    real(kind=8), parameter :: epsU=1e-8
    real(kind=8), parameter :: epsT=1e-8
    real(kind=8) :: errorU, errorT
    
    real(kind=8) :: xp(0:total_nx+1), yp(0:total_ny+1), zp(0:total_nz+1)
    real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), T(:,:,:)
    real(kind=8), allocatable :: rho(:,:,:)
    real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), Tp(:,:,:)

    real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)
    real(kind=8), allocatable :: g(:,:,:,:), g_post(:,:,:,:)
    
    integer :: ex(0:18), ey(0:18), ez(0:18)
    data ex/0, 1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0/
    data ey/0, 0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1/
    data ez/0, 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1/
    
    real(kind=8), allocatable :: Fx(:,:,:), Fy(:,:,:), Fz(:,:,:)
    
    real(kind=8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    real(kind=8), parameter :: Qd=3.0d0-dsqrt(3.0d0), Qnu=4.0d0*dsqrt(3.0d0)-6.0d0
    
    integer :: fileNum
    integer :: dimensionlessTime
    integer :: statisticallyStationaryState, statisticallyStationaryTime
    integer :: statisticallyConvergeState
    real(kind=8) :: NuVolAvg(0:dimensionlessTimeMax), ReVolAvg(0:dimensionlessTimeMax)
    real(kind=8) :: NuVolAvg_mean(0:dimensionlessTimeMax), ReVolAvg_mean(0:dimensionlessTimeMax)
    

    ! mpi data
    integer :: rc, rank, num_process 
    integer :: dims(0:2) = (/ 0, 0, 0 /), coords(0:2)
    logical :: periods(0:2) = (/ .false., .false., .false. /)
    integer :: comm3d, rank3d
    integer :: nbr_surface(1:6)
    integer :: nbr_line(7:18) 
    integer :: f_surface_x, f_surface_y, f_surface_z  ! surface data perpendicular to axis x, y, z
    integer :: f_line_x, f_line_y, f_line_z       ! line data parallel to axis x, y, z
    integer :: g_surface_x, g_surface_y, g_surface_z  
    integer :: g_line_x, g_line_y, g_line_z
    integer :: i_start_global, j_start_global, k_start_global

end module commondata





program main
    use omp_lib
    use mpi    
    use commondata
    implicit none
    real(kind=8) :: start, finish
    real(kind=8) :: start2, finish2
    real(8) :: start_time, end_time
    integer :: myMaxThreads
    integer :: name_len, stride
    character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

#ifdef _OPENMP
    write(*,*) "Starting OpenMP >>>>>>"
    call OMP_set_num_threads(24)
    myMaxThreads = OMP_get_max_threads()
    write(*,100) "|----------Max Running threads:",myMaxThreads,"----------|"
100 format(1X,A,I5)
#endif

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
    call decompose_1d(total_nx, nx, coords(0), dims(0), i_start_global)
    call decompose_1d(total_ny, ny, coords(1), dims(1), j_start_global)
    call decompose_1d(total_nz, nz, coords(2), dims(2), k_start_global)
    ! write(*,*) "coords = ", coords(1), coords(2), coords(3), "nx*ny = ", nx, ny, nz

    ! get the neighbors
    call MPI_Cart_shift(comm3d, 0, 1, nbr_surface(2), nbr_surface(1), rc)
    call MPI_Cart_shift(comm3d, 1, 1, nbr_surface(4), nbr_surface(3), rc)
    call MPI_Cart_shift(comm3d, 2, 1, nbr_surface(6), nbr_surface(5), rc)
    ! write(*,*) "I'm process ", rank3d, "My neighbor surfaces are", nbr_surface(1), nbr_surface(2), nbr_surface(3), nbr_surface(4), nbr_surface(5), nbr_surface(6)
    call MPI_Cart_find_corners()


    ! construct the datatype for the exchange  f
    ! line data parallel x direction (don't exchange in x direction)
    call MPI_Type_vector(nx, 1, 19, MPI_REAL8, f_line_x, rc)
    call MPI_Type_commit(f_line_x, rc)
    ! line data parallel y direction (don't exchange in y direction)
    call MPI_Type_vector(ny, 1, 19 * (nx+2), MPI_REAL8, f_line_y, rc)
    call MPI_Type_commit(f_line_y, rc)
    ! line data parallel z direction (don't exchange in z direction)
    call MPI_Type_vector(nz, 1, 19 * (nx+2) * (ny+2), MPI_REAL8, f_line_z, rc)
    call MPI_Type_commit(f_line_z, rc)

    ! surface data exchange in x direction -- line_y * nz, stride = 19*(nx+2)*(ny+2)
    stride = 19 * (nx+2) * (ny+2) 
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), f_line_y, f_surface_x, rc)
    call MPI_Type_commit(f_surface_x, rc)
    ! surface data exchange in y direction -- line_x * nz, stride = 19*(nx+2)*(ny+2)
    stride = 19 * (nx+2) * (ny+2)
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), f_line_x, f_surface_y, rc)
    call MPI_Type_commit(f_surface_y, rc)
    ! surface data exchange in z direction -- line_x * ny, stride = 19*(nx+2)
    stride = 19 * (nx+2)
    call MPI_Type_create_hvector(ny, 1, stride * sizeof(start), f_line_x, f_surface_z, rc)
    call MPI_Type_commit(f_surface_z, rc)


    ! construct the datatype for the exchange  g
    ! line data parallel x direction (don't exchange in x direction)
    call MPI_Type_vector(nx, 1, 7, MPI_REAL8, g_line_x, rc)
    call MPI_Type_commit(g_line_x, rc)
    ! line data parallel y direction (don't exchange in y direction)
    call MPI_Type_vector(ny, 1, 7 * (nx+2), MPI_REAL8, g_line_y, rc)
    call MPI_Type_commit(g_line_y, rc)
    ! line data parallel z direction (don't exchange in z direction)
    call MPI_Type_vector(nz, 1, 7 * (nx+2) * (ny+2), MPI_REAL8, g_line_z, rc)
    call MPI_Type_commit(g_line_z, rc)

    ! surface data exchange in x direction -- line_y * nz, stride = 19*(nx+2)*(ny+2)
    stride = 7 * (nx+2) * (ny+2) 
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), g_line_y, g_surface_x, rc)
    call MPI_Type_commit(g_surface_x, rc)
    ! surface data exchange in y direction -- line_x * nz, stride = 19*(nx+2)*(ny+2)
    stride = 7 * (nx+2) * (ny+2)
    call MPI_Type_create_hvector(nz, 1, stride * sizeof(start), g_line_x, g_surface_y, rc)
    call MPI_Type_commit(g_surface_y, rc)
    ! surface data exchange in z direction -- line_x * ny, stride = 19*(nx+2)
    stride = 7 * (nx+2)
    call MPI_Type_create_hvector(ny, 1, stride * sizeof(start), g_line_x, g_surface_z, rc)
    call MPI_Type_commit(g_surface_z, rc)


    call initial()

    call output()

    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( ((errorU.GT.epsU).OR.(errorT.GT.epsT)).AND.(itc.LE.itc_max).AND.(statisticallyConvergeState.EQ.0) )

        itc = itc+1

        call collision()

        call f_message_passing_sendrecv()

        call streaming()

        call bounceback()

        call collisionT()

        call g_message_passing_sendrecv()

        call streamingT()

        call bouncebackT()

        call macro()

        call macroT()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

        ! ! timer test
        ! if (mod(itc, 6000) == 0) then
        !     exit
        ! endif

! #ifdef benchmarkCavity
!         if(MOD(itc,1000000).EQ.0) then
!             call backupData()
!         endif
! #endif
        
! #ifdef RBconvection
!         if( MOD(itc,int(timeUnit)).EQ.0 ) then
            
!             call calNuRe()
            
!             if(statisticallyStationaryState.EQ.1) then

!                 ! call output_binary()

!             endif
            
!         endif
! #endif
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()
    
    
    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif

    if (rank == 0) then
        write(*,*) "---------------------------------------------"
        write(*,*) "Time (CPU) = ", real(finish-start), "s"
#ifdef _OPENMP
        write(*,*) "Time (OMP) = ", real(finish2-start2), "s"
#endif
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
        write(*,*) "---------------------------------------------"
    endif
    ! write(*,*) "Backup Distribution Function......"
    ! call backupData()

    call output()

    if (rank == 0) then
        write(*,*) "Deallocate Array..."
    endif

    deallocate(f)
    deallocate(g)
    deallocate(f_post)
    deallocate(g_post)
    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(T)
    deallocate(up)
    deallocate(vp)
    deallocate(wp)
    deallocate(Tp)
    deallocate(rho)
    deallocate(Fx)
    deallocate(Fy)
    deallocate(Fz)

    if (rank == 0) then
        write(*,*) "    "
    
        write(*,*) "Successfully: DNS completed!"
    endif

    call MPI_Finalize(rc)

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



subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: un(0:18), unT(0:6)
    real(kind=8) :: us2
    real(kind=8) :: omega(0:18), omegaT(0:6)
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    
    if (rank == 0) then
        write(*,*) 'Mesh:',nx,ny,nz
        write(*,*) 'Rayleigh=',real(Rayleigh), ', Prandtl=',real(Prandtl), ', Mach=',real(Mach)
        write(*,*) "Ekman=",real(Ekman)
        write(*,*) "   "
        write(*,*) 'tauf=',real(tauf)
        write(*,*) "paraA=",real(paraA)
        write(*,*) "viscosity=",real(viscosity), ", diffusivity=",real(diffusivity)
        write(*,*) "omegaRatating=", real(omegaRatating)
        write(*,*) "itc_max=",itc_max
        write(*,*) "Ouput will begin at", int(timeUnit)
        write(*,*) "Output interval is", int(timeUnit)
        write(*,*) "Time unit: Sqrt(L0/(gBeta*DeltaT))=",dsqrt(dble(Ny)/gBeta)
        write(*,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT)=",dsqrt(gBeta*dble(Ny))
        write(*,*) "    "
    
#ifdef RBconvection
    write(*,*) "    "
    write(*,*) "I am RBconvection"
    if(flowReversal.EQ.1) then
        write(*,*) "I am also flow reversal"
    write(*,*) "    "
    elseif(flowReversal.EQ.0) then
        write(*,*) "I am NOT flow reversal"
    write(*,*) "    "
    else
        write(*,*) "Error: Please check flowReversal setting!"
        stop
    endif
#endif

#ifdef benchmarkCavity
    write(*,*) "I am benchmarkCavity"
    write(*,*) "    "
#endif
    endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i = 1, total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j = 1, total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(total_nz+1) = dble(nz)
    do k = 1, total_nz
        zp(k) = dble(k)-0.5d0
    enddo
    
    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (T(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    allocate (up(nx,ny,nz))
    allocate (vp(nx,ny,nz))
    allocate (wp(nx,ny,nz))
    allocate (Tp(nx,ny,nz))
    
    allocate (f(0:18,nx,ny,nz))
    allocate (f_post(0:18,0:nx+1,0:ny+1,0:nz+1))
    allocate (g(0:6,nx,ny,nz))
    allocate (g_post(0:6,0:nx+1,0:ny+1,0:nz+1))
    
    allocate (Fx(nx,ny,nz))
    allocate (Fy(nx,ny,nz))
    allocate (Fz(nx,ny,nz))
    
    rho = 1.0d0
    
    omega(0) = 1.0d0/3.0d0
    do alpha=1,6
        omega(alpha) = 1.0d0/18.0d0
    enddo
    do alpha=7,18
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/7.0d0
    do alpha=1,6
        omegaT(alpha) = (paraA+6.0d0)/42.0d0
    enddo

    
    if(loadInitField.EQ.0) then 
        u = 0.0d0
        v = 0.0d0
        w = 0.0d0
        T = 0.0d0
        
    if (rank == 0) then
        write(*,*) "Initial field is set exactly"
#ifdef noslipWalls
    write(*,*) "Velocity B.C. for left/right walls are: ===No-slip wall==="
#endif
#ifdef TopBottomPlatesConstT
    write(*,*) "Temperature B.C. for bottom plate is:==Hot bottom plate==="
    write(*,*) "Temperature B.C. for top plate is:====Cold top plate==="
#endif
#ifdef LeftRightWallsConstT
    write(*,*) "Temperature B.C. for left wall is:===Hot left wall==="
    write(*,*) "Temperature B.C. for right wall is:==Cold right wall==="
#endif
#ifdef LeftRightWallsAdiabatic
    write(*,*) "Temperature B.C. for left/right walls are:===Adiabatic walls==="
#endif
#ifdef TopBottomPlatesAdiabatic
    write(*,*) "Temperature B.C. for top/bottom plates are:===Adiabatic plates==="
#endif

#ifdef BackFrontWallsAdiabatic
    write(*,*) "Temperature B.C. for back/front walls are:===Adiabatic walls==="
#endif    
    endif


#ifdef LeftRightWallsConstT
    ! left wall y=1
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                T(i,1,k) = Thot
            enddo
        enddo
    endif

    ! right wall y = ny
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                T(i,ny,k) = Tcold
            enddo
        enddo
    endif

#endif


#ifdef TopBottomPlatesConstT
    ! bottom wall k = 1
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                T(i,j,1) = Thot
            enddo
        enddo
    endif


    ! top wall k = nz
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                T(i,j,nz) = Tcold
            enddo
        enddo
    endif
#endif

        !$omp parallel do default(none) shared(f,g,u,v,w,T,ex,ey,ez,omega,omegaT,rho) private(i,j,k,alpha,us2,un,unT)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                    do alpha=0,18
                        un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        f(alpha,i,j,k) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                    enddo
                    do alpha=0,6
                        unT(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        g(alpha,i,j,k) = omegaT(alpha)*T(i,j,k)*(1.0d0+21.0d0/(6.0d0+paraA)*unT(alpha))
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

    ! elseif(loadInitField.EQ.1) then
    !     write(*,*) "Load initial field from previous simulation"
    !     write(*,*) "Start to read raw data >>>>>>>>>>>>"
    !     open(unit=01,file='./reload/backupFile-1000.bin',form="unformatted",access="sequential",status='old')
    !     read(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) ((((f(alpha,i,j,k),alpha=0,18),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) ((((g(alpha,i,j,k),alpha=0,6),i=1,nx),j=1,ny),k=1,nz)
    !     close(01)
    !     write(*,*) "Read raw data!"
    else
        write(*,*) "Error: initial field is not properly set"
    endif
    
    up = 0.0d0
    vp = 0.0d0
    wp = 0.0d0
    Tp = 0.0d0
    
    f_post = 0.0d0
    g_post = 0.0d0
        
    fileNum = 0
    dimensionlessTime = 0
    
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 0
    statisticallyConvergeState = 0

    return
end subroutine initial


subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,w,Fx,Fy,Fz,T) private(i,j,k,alpha,s,m,m_post,meq,fSource) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
    !--------------------------------------------------------------------------------------------------------------------
    !---m0    
    m(0) =f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m1
    m(1) = -30.0d0*f(0,i,j,k)-11.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
    +8.0d0*( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m2
    m(2) = 12.0d0*f(0,i,j,k)-4.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m3
    m(3) = f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m4
    m(4) = -4.0d0*(f(1,i,j,k)-f(2,i,j,k))+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k) &
                    +f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m5
    m(5) = f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m6
    m(6) = -4.0d0*(f(3,i,j,k)-f(4,i,j,k))+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k) &
                +f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m7
    m(7) = f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m8
    m(8) = -4.0d0*(f(5,i,j,k)-f(6,i,j,k))+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k) &
            +f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m9
    m(9) = 2.0d0*(f(1,i,j,k)+f(2,i,j,k))-f(3,i,j,k)-f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m10
     m(10) = -4.0d0*(f(1,i,j,k)+f(2,i,j,k))+2.0d0*(f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k)) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m11
    m(11) = f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m12
     m(12) = -2.0d0*(f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k))+( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k) )-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m13
    m(13) = f(7,i,j,k)-f(8,i,j,k)-f(9,i,j,k)+f(10,i,j,k)
    !---m14
    m(14) = f(15,i,j,k)-f(16,i,j,k)-f(17,i,j,k)+f(18,i,j,k)
    !---m15
    m(15) = f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m16
    m(16) = f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)-f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m17
    m(17) = -f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m18
    m(18) = f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)-f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !--------------------------------------------------------------------------------------------------------------------

    meq(0) = rho(i,j,k)
    meq(1) = -11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(2) = 3.0d0*rho(i,j,k)-11.0d0/2.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(3) = rho(i,j,k)*u(i,j,k)
    meq(4) = -2.0d0/3.0d0*meq(3)
    meq(5) = rho(i,j,k)*v(i,j,k)
    meq(6) = -2.0d0/3.0d0*meq(5)
    meq(7) = rho(i,j,k)*w(i,j,k)
    meq(8) = -2.0d0/3.0d0*meq(7)
    meq(9) = rho(i,j,k)*(2.0d0*u(i,j,k)*u(i,j,k)-v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(10) = -0.5d0*meq(9)
    meq(11) = rho(i,j,k)*(v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(12) = -0.5d0*meq(11)
    meq(13) = rho(i,j,k)*(u(i,j,k)*v(i,j,k))
    meq(14) = rho(i,j,k)*(v(i,j,k)*w(i,j,k))
    meq(15) = rho(i,j,k)*(w(i,j,k)*u(i,j,k))
    meq(16) = 0.0d0
    meq(17) = 0.0d0
    meq(18) = 0.0d0

            s(0) = 0.0d0 
            s(1) = Snu  !!!s_{e}
            s(2) = Snu   !!! s_{epsilon}
            s(3) = 0.0d0 
            s(4) = Sq   !!! s_{q}
            s(5) = 0.0d0 
            s(6) = Sq   !!! s_{q}
            s(7) = 0.0d0 
            s(8) = Sq   !!! s_{q}
            s(9) = Snu !!! s_{nu}
            s(10) = Snu   !!! s_{pi}
            s(11) = Snu   !!! s_{nu}
            s(12) = Snu !!! s_{pi}
            s(13) = Snu !!! s_{nu}
            s(14) = Snu   !!! s_{nu}
            s(15) = Snu   !!! s_{nu}
            s(16) = Sq   !!! s_{m}
            s(17) = Sq   !!! s_{m}
            s(18) = Sq   !!! s_{m}

            Fx(i,j,k) = -2.0d0*rho(i,j,k)*v(i,j,k)*omegaRatating
            Fy(i,j,k) = 2.0d0*rho(i,j,k)*u(i,j,k)*omegaRatating
            Fz(i,j,k) = rho(i,j,k)*gBeta*(T(i,j,k)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = 38.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(2) = -11.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(3) = Fx(i,j,k)
            fSource(4) = -2.0d0/3.0d0*Fx(i,j,k)
            fSource(5) = Fy(i,j,k)
            fSource(6) = -2.0d0/3.0d0*Fy(i,j,k)
            fSource(7) = Fz(i,j,k)
            fSource(8) = -2.0d0/3.0d0*Fz(i,j,k)
            fSource(9) = 4.0d0*u(i,j,k)*Fx(i,j,k)-2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(10) = -2.0d0*u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k) 
            fSource(11) = 2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(12) = -v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k)
            fSource(13) = u(i,j,k)*Fy(i,j,k)+v(i,j,k)*Fx(i,j,k)
            fSource(14) = v(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fy(i,j,k)
            fSource(15) = u(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fx(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+(1.0d0-0.5d0*s(alpha))*fSource(alpha)
            enddo

    f_post(0,i,j,k) = m_post(0)/19.0d0-5.0d0/399.0d0*m_post(1)+m_post(2)/21.0d0

    f_post(1,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(3,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(4,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(5,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

!---
    f_post(7,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0+( m_post(16)-m_post(17) )*0.125d0

    f_post(8,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0-( m_post(16)+m_post(17) )*0.125d0

    f_post(9,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0+( m_post(16)+m_post(17) )*0.125d0

    f_post(10,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0-( m_post(16)-m_post(17) )*0.125d0

!---
    f_post(11,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)-0.1250d0*( m_post(16)-m_post(18) )

    f_post(12,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)+0.125d0*( m_post(16)+m_post(18) )

    f_post(13,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)-0.125d0*( m_post(16)+m_post(18) )

    f_post(14,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)+0.125d0*( m_post(16)-m_post(18) )

!---
    f_post(15,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)+0.125d0*( m_post(17)-m_post(18) )

    f_post(16,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)-0.125d0*( m_post(17)+m_post(18) )

    f_post(17,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)+0.125d0*( m_post(17)+m_post(18) )

    f_post(18,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)-0.125d0*( m_post(17)-m_post(18) )

            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,18
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef noslipWalls
    !Back plane (i = 1)
    if (coords(0) == 0) then
        do k=1,nz
            do j=1,ny
                f(1,1,j,k) = f_post(2,1,j,k)
                f(7,1,j,k) = f_post(10,1,j,k)
                f(9,1,j,k) = f_post(8,1,j,k)
                f(11,1,j,k) = f_post(14,1,j,k)
                f(13,1,j,k) = f_post(12,1,j,k)
            enddo
        enddo
    endif

    ! !Front plane (i=nx)
    if (coords(0) == dims(0) - 1) then
        do k=1,nz
            do j=1,ny
                f(2,nx,j,k) = f_post(1,nx,j,k)
                f(8,nx,j,k) = f_post(9,nx,j,k)
                f(10,nx,j,k) = f_post(7,nx,j,k)
                f(12,nx,j,k) = f_post(13,nx,j,k)
                f(14,nx,j,k) = f_post(11,nx,j,k)
            enddo
        enddo
    endif

    ! Left plane (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                f(3,i,1,k) = f_post(4,i,1,k)
                f(7,i,1,k) = f_post(10,i,1,k)
                f(8,i,1,k) = f_post(9,i,1,k)
                f(15,i,1,k) = f_post(18,i,1,k)
                f(17,i,1,k) = f_post(16,i,1,k)
            enddo
        enddo
    endif

    ! Right plane (j=ny)
    if  (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                !Right plane (j=ny)
                f(4,i,ny,k) = f_post(3,i,ny,k)
                f(9,i,ny,k) = f_post(8,i,ny,k)
                f(10,i,ny,k) = f_post(7,i,ny,k)
                f(16,i,ny,k) = f_post(17,i,ny,k)
                f(18,i,ny,k) = f_post(15,i,ny,k)
            enddo
        enddo
    endif
#endif

    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                f(5,i,j,1) = f_post(6,i,j,1)
                f(11,i,j,1) = f_post(14,i,j,1)
                f(12,i,j,1) = f_post(13,i,j,1)
                f(15,i,j,1) = f_post(18,i,j,1)
                f(16,i,j,1) = f_post(17,i,j,1)
            enddo
        enddo
    endif


    !Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                f(6,i,j,nz) = f_post(5,i,j,nz)
                f(13,i,j,nz) = f_post(12,i,j,nz)
                f(14,i,j,nz) = f_post(11,i,j,nz)
                f(17,i,j,nz) = f_post(16,i,j,nz)
                f(18,i,j,nz) = f_post(15,i,j,nz)
            enddo
        enddo
    endif

    return
end subroutine bounceback
    

subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(f,rho,u,v,w,Fx,Fy,Fz) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k) = f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
            
                u(i,j,k) = ( f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)+0.5d0*Fx(i,j,k) )/rho(i,j,k)
            
                v(i,j,k) = ( f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fy(i,j,k) )/rho(i,j,k)
            
                w(i,j,k) = ( f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fz(i,j,k) )/rho(i,j,k)
            
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macro



subroutine collisionT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !------------------------
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: q(0:6)

    !$omp parallel do default(none) shared(g,g_post,u,v,w,T) private(i,j,k,alpha,n,neq,q,n_post) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
            
    n(0) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(1) = g(1,i,j,k)-g(2,i,j,k)
    n(2) = g(3,i,j,k)-g(4,i,j,k)
    n(3) = g(5,i,j,k)-g(6,i,j,k)
    n(4) = -6.0d0*g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(5) = 2.0d0*g(1,i,j,k)+2.0d0*g(2,i,j,k)-g(3,i,j,k)-g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
    n(6) = g(3,i,j,k)+g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
        
            neq(0) = T(i,j,k)
            neq(1) = T(i,j,k)*u(i,j,k)
            neq(2) = T(i,j,k)*v(i,j,k)
            neq(3) = T(i,j,k)*w(i,j,k)
            neq(4) = T(i,j,k)*paraA
            neq(5) = 0.0d0
            neq(6) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qd
            q(4) = Qnu
            q(5) = Qnu
            q(6) = Qnu
        
            do alpha=0,6
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j,k) = n_post(0)/7.0d0-n_post(4)/7.0d0
    g_post(1,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(2,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(3,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(4,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(5,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
    g_post(6,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
        
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collisionT


subroutine streamingT()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(g,g_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,6
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    g(alpha,i,j,k) = g_post(alpha,ip,jp,kp)
                    
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine streamingT


subroutine bouncebackT()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef TopBottomPlatesAdiabatic
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                g(5,i,j,1) = g_post(6,i,j,1)
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                g(6,i,j,nz) = g_post(5,i,j,nz)
            enddo
        enddo
    endif
#endif

#ifdef TopBottomPlatesConstT
    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                g(5,i,j,1) = -g_post(6,i,j,1)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                g(6,i,j,nz) = -g_post(5,i,j,nz)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif


#ifdef LeftRightWallsConstT
    ! Left side (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                g(3,i,1,k) = -g_post(4,i,1,k)+(6.0d0+paraA)/21.0d0*Thot
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                g(4,i,ny,k) = -g_post(3,i,ny,k)+(6.0d0+paraA)/21.0d0*Tcold
            enddo
        enddo
    endif
#endif

#ifdef LeftRightWallsAdiabatic
    ! Left side (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                g(3,i,1,k) = g_post(4,i,1,k)
            enddo
        enddo
    endif

    ! Right side (j=ny)
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                g(4,i,ny,k) = g_post(3,i,ny,k)
            enddo
        enddo
    endif
#endif
    
    
#ifdef BackFrontWallsAdiabatic
    ! Back side (i=1)
    if (coords(0) == 0) then
        do k=1,nz
            do j=1,ny
                g(1,1,j,k) = g_post(2,1,j,k)
            enddo
        enddo
    endif

    !Front side (i=nx)
    if (coords(0) == dims(0) - 1) then
        do k=1,nz
            do j=1,ny          
                g(2,nx,j,k) = g_post(1,nx,j,k)
            enddo
        enddo
    endif

#endif

    return
end subroutine bouncebackT




subroutine macroT()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(g,T) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                T(i,j,k) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macroT


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





subroutine f_message_passing_sendrecv()
    use mpi
    use commondata
    implicit none
    integer :: tag(5), idx

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)
    tag = (/ 1, 7, 9, 11, 13 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), nx, 1, 1), 1, f_surface_x, nbr_surface(1), tag(idx), &
                    f_post(tag(idx), 0, 1, 1), 1, f_surface_x, nbr_surface(2), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (i--)
    tag = (/ 2, 8, 10 ,12, 14 /) 
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, f_surface_x, nbr_surface(2), tag(idx), &
                    f_post(tag(idx), nx+1, 1, 1), 1, f_surface_x, nbr_surface(1), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    tag = (/ 3, 7, 8, 15, 17 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, ny, 1), 1, f_surface_y, nbr_surface(3), tag(idx), &
                    f_post(tag(idx), 1, 0, 1), 1, f_surface_y, nbr_surface(4), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (j--)
    tag = (/ 4, 9, 10, 16, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, f_surface_y, nbr_surface(4), tag(idx), &
                    f_post(tag(idx), 1, ny+1, 1), 1, f_surface_y, nbr_surface(3), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    tag = (/ 5, 11, 12, 15, 16 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, nz), 1, f_surface_z, nbr_surface(5), tag(idx), &
                    f_post(tag(idx), 1, 1, 0), 1, f_surface_z, nbr_surface(6), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! message passing to (k--)
    tag = (/ 6, 13, 14, 17, 18 /)
    do idx = 1, 5
        call MPI_Sendrecv(f_post(tag(idx), 1, 1, 1), 1, f_surface_z, nbr_surface(6), tag(idx), &
                    f_post(tag(idx), 1, 1, nz+1), 1, f_surface_z, nbr_surface(5), tag(idx), &
                    comm3d, MPI_STATUS_IGNORE, rc)
    enddo

    ! /* ----------------------------- line data ----------------------------- */
    ! ------------ exchange message perpendicular to z ----------------
    ! message passing to (i++, j++) (7)
    call MPI_Sendrecv(f_post(7, nx, ny, 1), 1, f_line_z, nbr_line(7), 7, &
                    f_post(7, 0, 0, 1), 1, f_line_z, nbr_line(10), 7, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j--) (10)
    call MPI_Sendrecv(f_post(10, 1, 1, 1), 1, f_line_z, nbr_line(10), 10, &
                    f_post(10, nx+1, ny+1, 1), 1, f_line_z, nbr_line(7), 10, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, j--) (9)
    call MPI_Sendrecv(f_post(9, nx, 1, 1), 1, f_line_z, nbr_line(9), 9, &
                    f_post(9, 0, ny+1, 1), 1, f_line_z, nbr_line(8), 9, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, j++) (8)
    call MPI_Sendrecv(f_post(8, 1, ny, 1), 1, f_line_z, nbr_line(8), 8, &
                    f_post(8, nx+1, 0, 1), 1, f_line_z, nbr_line(9), 8, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to y ----------------
    ! message passing to (i++, k++) (11)
    call MPI_Sendrecv(f_post(11, nx, 1, nz), 1, f_line_y, nbr_line(11), 11, &
                    f_post(11, 0, 1, 0), 1, f_line_y, nbr_line(14), 11, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k--) (14)
    call MPI_Sendrecv(f_post(14, 1, 1, 1), 1, f_line_y, nbr_line(14), 14, &
                    f_post(14, nx+1, 1, nz+1), 1, f_line_y, nbr_line(11), 14, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i++, k--) (13)
    call MPI_Sendrecv(f_post(13, nx, 1, 1), 1, f_line_y, nbr_line(13), 13, &
                    f_post(13, 0, 1, nz+1), 1, f_line_y, nbr_line(12), 13, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--, k++) (12)
    call MPI_Sendrecv(f_post(12, 1, 1, nz), 1, f_line_y, nbr_line(12), 12, &
                    f_post(12, nx+1, 1, 0), 1, f_line_y, nbr_line(13), 12, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message perpendicular to x ----------------
    ! message passing to (j++, k++) (15)
    call MPI_Sendrecv(f_post(15, 1, ny, nz), 1, f_line_x, nbr_line(15), 15, &
                    f_post(15, 1, 0, 0), 1, f_line_x, nbr_line(18), 15, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k--) (18)
    call MPI_Sendrecv(f_post(18, 1, 1, 1), 1, f_line_x, nbr_line(18), 18, &
                    f_post(18, 1, ny+1, nz+1), 1, f_line_x, nbr_line(15), 18, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j++, k--) (17)
    call MPI_Sendrecv(f_post(17, 1, ny, 1), 1, f_line_x, nbr_line(17), 17, &
                    f_post(17, 1, 0, nz+1), 1, f_line_x, nbr_line(16), 17, &
                    comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (j--, k++) (16)
    call MPI_Sendrecv(f_post(16, 1, 1, nz), 1, f_line_x, nbr_line(16), 16, &
                    f_post(16, 1, ny+1, 0), 1, f_line_x, nbr_line(17), 16, &
                    comm3d, MPI_STATUS_IGNORE, rc)
    

    ! don't need to echange point data


end subroutine f_message_passing_sendrecv





subroutine g_message_passing_sendrecv()
    use mpi
    use commondata
    implicit none
    integer :: tag(5), idx

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)
    call MPI_Sendrecv(g_post(1, nx, 1, 1), 1, g_surface_x, nbr_surface(1), 1, &
                g_post(1, 0, 1, 1), 1, g_surface_x, nbr_surface(2), 1, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--)
    call MPI_Sendrecv(g_post(2, 1, 1, 1), 1, g_surface_x, nbr_surface(2), 2, &
                g_post(2, nx+1, 1, 1), 1, g_surface_x, nbr_surface(1), 2, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    call MPI_Sendrecv(g_post(3, 1, ny, 1), 1, g_surface_y, nbr_surface(3), 3, &
                g_post(3, 1, 0, 1), 1, g_surface_y, nbr_surface(4), 3, &
                comm3d, MPI_STATUS_IGNORE, rc)


    ! message passing to (j--)
    call MPI_Sendrecv(g_post(4, 1, 1, 1), 1, g_surface_y, nbr_surface(4), 4, &
                g_post(4, 1, ny+1, 1), 1, g_surface_y, nbr_surface(3), 4, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    call MPI_Sendrecv(g_post(5, 1, 1, nz), 1, g_surface_z, nbr_surface(5), 5, &
                g_post(5, 1, 1, 0), 1, g_surface_z, nbr_surface(6), 5, &
                comm3d, MPI_STATUS_IGNORE, rc)


    ! message passing to (k--)
    call MPI_Sendrecv(g_post(6, 1, 1, 1), 1, g_surface_z, nbr_surface(6), 6, &
                g_post(6, 1, 1, nz+1), 1, g_surface_z, nbr_surface(5), 6, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! don't need to echange line data and point data


end subroutine g_message_passing_sendrecv




subroutine output()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    integer :: p_rank, num(0:5), new_coords(0:2), dx = 0, dy = 0, dz = 0
    real(8), allocatable :: total_u(:, :, :), total_v(:, :, :), total_w(:, :, :), total_rho(:, :, :), total_T(:, :, :)
    real(8), allocatable :: tmp_u(:, :, :), tmp_v(:, :, :), tmp_w(:, :, :), tmp_rho(:, :, :), tmp_T(:, :, :)
    
    ! use rank 0 to receive data and output the results

    if (rank3d > 0) then  !!! ----  rank != 0 send data
        ! collect the rank information
        num(0) = nx
        num(1) = ny
        num(2) = nz
        num(3) = i_start_global
        num(4) = j_start_global
        num(5) = k_start_global
        ! send to rank 0
        call MPI_Send(num, 6, MPI_INTEGER, 0, 0, comm3d, rc)    ! rank information
        call MPI_Send(u, nx*ny*nz, MPI_REAL8, 0, 1, comm3d, rc)
        call MPI_Send(v, nx*ny*nz, MPI_REAL8, 0, 2, comm3d, rc)
        call MPI_Send(w, nx*ny*nz, MPI_REAL8, 0, 3, comm3d, rc)
        call MPI_Send(rho, nx*ny*nz, MPI_REAL8, 0, 4, comm3d, rc)
        call MPI_Send(T, nx*ny*nz, MPI_REAL8, 0, 5, comm3d, rc)
    else    
        !!! ---- rank 0 collect data
        ! allocate array
        allocate(total_u(total_nx, total_ny, total_nz))
        allocate(total_v(total_nx, total_ny, total_nz))
        allocate(total_w(total_nx, total_ny, total_nz))
        allocate(total_rho(total_nx, total_ny, total_nz))
        allocate(total_T(total_nx, total_ny, total_nz))

        ! determine the origin
        dx = i_start_global
        dy = j_start_global
        dz = k_start_global

        ! collect data from rank 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    total_u(dx + i, dy + j, dz + k) = u(i, j, k)
                    total_v(dx + i, dy + j, dz + k) = v(i, j, k)
                    total_w(dx + i, dy + j, dz + k) = w(i, j, k)
                    total_rho(dx + i, dy + j, dz + k) = rho(i, j, k)
                    total_T(dx + i, dy + j, dz + k) = T(i, j, k)
                enddo
            enddo
        enddo

        ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) * dims(2) - 1

            call MPI_Cart_coords(comm3d, p_rank, 3, new_coords, rc)

            ! receive the block size and origion
            call MPI_Recv(num, 6, MPI_INTEGER, p_rank, 0, comm3d, MPI_STATUS_IGNORE, rc)

            ! creat buffer
            allocate(tmp_u(num(0), num(1), num(2)))
            allocate(tmp_v(num(0), num(1), num(2)))
            allocate(tmp_w(num(0), num(1), num(2)))
            allocate(tmp_rho(num(0), num(1), num(2)))
            allocate(tmp_T(num(0), num(1), num(2)))

            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 1, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 2, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_w, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 3, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 4, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_T, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 5, comm3d, MPI_STATUS_IGNORE, rc)


            ! determine the origin
            dx = num(3)
            dy = num(4)
            dz = num(5)

            ! assign data
            do k = 1, num(2)
                do j = 1, num(1)
                    do i = 1, num(0)
                        total_u(i + dx, j + dy, k + dz) = tmp_u(i, j, k)
                        total_v(i + dx, j + dy, k + dz) = tmp_v(i, j, k)
                        total_w(i + dx, j + dy, k + dz) = tmp_w(i, j, k)
                        total_rho(i + dx, j + dy, k + dz) = tmp_rho(i, j, k)
                        total_T(i + dx, j + dy, k + dz) = tmp_T(i, j, k)
                    enddo
                enddo
            enddo

            ! de-allocate buffer arrays
            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_w)
            deallocate(tmp_rho)
            deallocate(tmp_T)
        enddo

        ! after collect total_* data, then output
        ! call output_ASCII(xp, yp, zp, total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
        call output_binary(total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
        call output_Tecplot(xp, yp, zp, total_u, total_v, total_w, total_rho, total_T, total_nx, total_ny, total_nz, itc)
        ! call out_Velocity_Nu(total_u, total_v, total_w, total_T, total_nx, total_ny, total_nz, diffusivity, lengthUnit)

        ! de-allocate total arrays
        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_w)
        deallocate(total_rho)
        deallocate(total_T)
    endif

end subroutine output



subroutine output_binary(u, v, w, rho, T, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz), T(nx, ny, nz)
    integer :: i, j, k
    character(len=100) :: filename

#ifdef RBconvection
    fileNum = fileNum+1
    write(filename,*) fileNum
    filename = adjustl(filename)
#endif
#ifdef benchmarkCavity
    write(filename,*) itc
    filename = adjustl(filename)
#endif

    open(unit=01,file='buoyancyCavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    write(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(01)

    return
end subroutine output_binary




subroutine output_Tecplot(xp, yp, zp, u, v, w, rho, T, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz), T(nx, ny, nz)
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    character(len=100) :: filename
  
#ifdef RBconvection
    fileNum = fileNum+1
    write(filename,*) fileNum
    filename = adjustl(filename)
#endif
#ifdef benchmarkCavity
    write(filename,*) itc
    filename = adjustl(filename)
#endif
    open(41,file='buoyancyCavity-'//trim(filename)//'.plt',form='binary')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--Integer value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file (here 5 variables)
    write(41) 7

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='Z'
    call dumpstring(V3)
    
    V4='U'
    call dumpstring(V4)
    V5='V'
    call dumpstring(V5)
    V6='W'
    call dumpstring(V6)
    
    V7='T'
    call dumpstring(V7)

    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName='ZONE 001'
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
                write(41) real(T(i,j,k))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
end subroutine output_Tecplot


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





! subroutine backupData()
!     use commondata
!     implicit none
!     integer :: i, j, k, alpha
!     character(len=100) :: filename

!     write(filename,*) itc
!     filename = adjustl(filename)

!     open(unit=01,file='backupFile-'//trim(filename)//'.bin',form="unformatted",access="sequential")
!     write(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!     write(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!     write(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!     write(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!     write(01) ((((f(alpha,i,j,k),alpha=0,18),i=1,nx),j=1,ny),k=1,nz)
!     write(01) ((((g(alpha,i,j,k),alpha=0,6),i=1,nx),j=1,ny),k=1,nz)
!     close(01)

!     return
! end subroutine backupData
