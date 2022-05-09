module commondata
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)

    integer :: nx, ny
    integer, parameter :: total_nx=201, total_ny=101
    integer, parameter :: nxHalf=(total_nx-1)/2+1, nyHalf=(total_ny-1)/2+1
    integer, parameter :: nxFourth=(total_nx-1)/4+1, nyFourth=(total_ny-1)/4+1
    integer :: itc
    integer, parameter :: itc_max=INT(5000)
    integer, parameter :: cNumMax=1

    real(8) :: xCenter(cNumMax), yCenter(cNumMax)
    real(8) :: xCenterOld(cNumMax), yCenterOld(cNumMax)
    
    real(8) :: rationalOmega(cNumMax)
    real(8) :: rationalOmegaOld(cNumMax)
    
    real(8) :: Uc(cNumMax), Vc(cNumMax)
    real(8) :: UcOld(cNumMax), VcOld(cNumMax)
    
    real(8), parameter :: radius0=25.25d0/2.0d0 
    real(8) :: radius(cNumMax)
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: rhoSolid=2.0d0
    real(8) :: rhoAvg
    
    real(8), parameter :: U0=0.02d0 
    real(8), parameter :: Uwall=0.1d0
    real(8), parameter :: viscosity=1.0d0/9.0d0
    real(8), parameter :: tauf=3.0d0*viscosity+0.5d0
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    real(8), parameter :: shearRate=2.0d0*Uwall/dble( total_ny)
    real(8), parameter :: Reynolds=shearRate*(2.0d0*radius0)**2.0d0/viscosity

    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    real(8), allocatable :: X(:), Y(:)
    real(8), allocatable :: u(:, :), v(:, :), rho(:, :)
    real(8), allocatable :: up(:, :), vp(:, :)
    real(8), allocatable :: f(:, :, :), f_post(:, :, :)
    integer, allocatable :: obst(:, :)
    integer, allocatable :: obstNew(:, :)

    real(8) :: omega(0:8)
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
    integer :: r(1:8)
    data r/3, 4, 1, 2, 7, 8, 5, 6/
    
    real(8) :: wallTotalForceX(cNumMax), wallTotalForceY(cNumMax)
    real(8) :: totalTorque(cNumMax)

    ! mpi data
    integer :: rc, rank, num_process
    integer :: dims(0:1) = (0, 0), coords(0:1)
    logical :: periods(0:1)
    data periods/.true., .false./
    integer :: comm2d, rank2d
    integer :: f_row_x, f_column_y, type_f
    integer :: type_fi_x, type_f_x, type_fi_y, type_f_y
    integer :: nbr_left, nbr_right, nbr_top, nbr_bottom
    integer :: cnr_top_left, cnr_top_right, cnr_bottom_left, cnr_bottom_right
    integer :: i_start_global, j_start_global
    
end module commondata

#define movingFrame
! #define stationaryFrame
    
#define linear
! #define quadratic

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
    integer(KIND=MPI_ADDRESS_KIND) :: lower_bound, extent

    call MPI_Init(rc)

    call MPI_Comm_size(MPI_COMM_WORLD, num_process, rc)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, rc)
    call MPI_Get_processor_name(processor_name, name_len, rc)

    !!! decomposition the domain 
    dims(1) = 1
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
    ! f_column_y exchange in x direction(non-contiguous) -- left and right
    call MPI_Type_vector(ny, 1, 9 * (nx+4), MPI_REAL8, f_column_y, rc)
    call MPI_Type_commit(f_column_y, rc)

    ! particles data
    call MPI_Type_contiguous(9, MPI_REAL8, type_f, rc)
    call MPI_Type_commit(type_f, rc)
    extent = 9 * sizeof(start)
    call MPI_Type_create_resized(MPI_REAL8, 0, extent, type_fi_x, rc) ! f_i in row
    call MPI_Type_create_resized(type_f, 0, extent, type_f_x, rc)     ! f_{1-9} in row
    call MPI_Type_commit(type_fi_x, rc)
    call MPI_Type_commit(type_f_x, rc)

    extent = 9 * (nx + 4) * sizeof(start)
    call MPI_Type_create_resized(MPI_REAL8, 0, extent, type_fi_y, rc) ! f_i in column
    call MPI_Type_create_resized(type_f, 0, extent, type_f_y, rc)     ! f_{1-9} in column
    call MPI_Type_commit(type_fi_y, rc)
    call MPI_Type_commit(type_f_y, rc)


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

        ! call message_passing_sendrecv()
        call send_all_fp()

        call streaming()

        call bounceback()

        call bounceback_particle()

        call send_all_f()

        call macro()

        call calForce()
    

        ! call message_particle()
        ! call send_all()

        itc = itc+1
        
        if(MOD(itc,1000).EQ.0) then
            call check()
            ! write(*, *) "rhoavg = ", rhoAvg
        endif

        call updateCenter()
        
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


subroutine initial()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
    real(8) :: us2
    integer :: cNum
    real(8) :: dx, dy, dt

if(rank == 0) then

#ifdef linear
    write(*,*) "    "
    write(*,*) "Linear interpolation!"
    write(*,*) "    "
#endif
#ifdef quadratic
    write(*,*) "    "
    write(*,*) "Quadratic interpolation!"
    write(*,*) "    "
#endif

#ifdef stationaryFrame
    write(*,*) "I am stationaryFrame"
#endif
#ifdef movingFrame
    write(*,*) "I am movingFrame"
#endif
    write(*,*) "Reynolds=",real(Reynolds),",   tauf=",real(tauf)
    write(*,*) "U0 (cylinder)=",real(U0)
    write(*,*) "Uwall=",real(Uwall)
    write(*,*) "    "

endif

    itc = 0
    errorU = 100.0d0

    dx = 1.0d0
    dy = 1.0d0
    dt = 1.0d0
    
    radius(1) = dble(radius0) 
    if (rank == 0) then
        do cNum=1,cNumMax
            write(*,*) "diameter=",real(2.0d0*radius(cNum))
        enddo
    endif

    allocate(X(total_nx))
    allocate(Y(total_ny))
    ! define mesh
    if (rank == 0) then
        do i=1,total_nx
            X(i) = (i-1)*dx
        enddo
        do j=1,total_ny
            Y(j) = (j-1)*dy
        enddo
    endif

    allocate(u(nx,ny))
    allocate(v(nx,ny))
    allocate(rho(nx,ny))
    allocate(up(nx,ny))
    allocate(vp(nx,ny))

    allocate(obst(0:nx+1, 0:ny+1))
    allocate(obstNew(0:nx+1, 0:ny+1))

    allocate(f(0:8, -1:nx+2, -1:ny+2))
    allocate(f_post(0:8, -1:nx+2, -1:ny+2))

    
    xCenter(1) = 30.0d0 
    yCenter(1) = 54.0d0 
#ifdef stationaryFrame
    Uc(1) = U0
    Vc(1) = 0.0d0
#endif
#ifdef movingFrame
    Uc(1) = 0.0d0
    Vc(1) = 0.0d0
#endif
    rationalOmega(1) = 0.0d0

    obst = 0
    obstNew = 0
    rho = rho0

    do j=0,ny+1
        do i=0,nx+1
            do cNum=1,cNumMax
                ! golbal coordinate
                if( ((i+i_start_global-xCenter(cNum))**2.0d0 + (j+j_start_global-yCenter(cNum))**2.0d0) &
                                .LE.radius(cNum)**2.0d0 ) then 
                    obst(i,j) = 1 ! solid node                  
                endif
                ! if( ((i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then 
                !     obst(i,j) = 1 ! solid node
                !     rho(i,j) = rhoSolid
                ! endif
            enddo
        enddo
    enddo

    do j=1,ny
        do i=1,nx
            if( obst(i,j) == 1 ) then 
                rho(i,j) = rhoSolid
            endif
        enddo
    enddo




#ifdef stationaryFrame
    u = U0
    ! bottom wall(j = 1)
    if (coords(1) == 0) then
        do i=1,nx
            u(i,1) = -Uwall
        enddo
    endif
    ! upper wall(j = ny)
    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            u(i,ny) = Uwall
        enddo
    endif
#endif

#ifdef movingFrame
    u = 0.0d0
    ! bottom wall(j = 1)
    if (coords(1) == 0) then
        do i=1,nx
            u(i,1) = -Uwall - U0
        enddo
    endif
    ! upper wall(j = ny)
    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            u(i,ny) = Uwall - U0
        enddo
    endif
#endif
    v = 0.0d0

    up = u
    vp = v

    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo

    f = 0.0d0
    !$omp parallel do default(none) shared(f,u,v,ex,ey,omega,obst) private(i,j,alpha,us2,un)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
                do alpha=0,8
                    un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                    f(alpha,i,j) = omega(alpha)*( 1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2 )
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine initial



subroutine collision()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: m(0:8), m_post(0:8), meq(0:8)
    real(8) :: s(0:8)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,obst) private(i,j,alpha,s,m,m_post,meq) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
    m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -meq(3)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -meq(5)
            meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
            meq(8) = rho(i,j)*( u(i,j)*v(i,j) )

            s(0) = 0.0d0      !!s_{\rho}
            s(1) = Snu !!s_{e}
            s(2) = Snu !!s_{\epsilon}
            s(3) = 0.0d0      !!s_{j} 
            s(4) = Sq !!s_{q}
            s(5) = 0.0d0      !!s_{j}
            s(6) = Sq       !!s_{q}
            s(7) = Snu !!s_{\nu}
            s(8) = Snu       !!s_{\nu}

            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0 
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0

            endif
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, alpha
    integer :: ip, jp

    ! !$omp parallel do default(none) shared(ex,ey,f,f_post,obst) private(i,j,alpha,ip,jp)
    ! do j=1,ny
    !     do i=1,nx
    !         if(obst(i,j).EQ.0) then
    !             do alpha=0,8
    !                 ip = i+ex(alpha)
    !                 jp = j+ey(alpha)

                    ! if(ip.EQ.0) ip = nx
                    ! if(ip.EQ.nx+1) ip = 1

    !                 f(alpha,ip,jp) = f_post(alpha,i,j)
    !             enddo
    !         endif
    !     enddo
    ! enddo
    ! !$omp end parallel do

    do j=1,ny
        do i=1,nx
            do alpha=0,8
                ip = i-ex(alpha)
                jp = j-ey(alpha)

                ! only have one subdomain in x direction (peroid)
                if (dims(0) == 1) then
                    if(ip.EQ.0) ip = nx
                    if(ip.EQ.nx+1) ip = 1
                endif

                if(obst(ip,jp).EQ.0) then
                    f(alpha,i,j) = f_post(alpha,ip,jp)
                endif
            enddo
        enddo
    enddo

    return
end subroutine streaming



subroutine bounceback_particle()
    use mpi
    use commondata
    implicit none
    integer :: i, j, alpha
    integer :: ip, jp
    real(8) :: q
    integer :: cNum
    integer :: myFlag
    integer :: fluidNum, total_fluidNum = 0
    real(8) :: temp1, temp2
    real(8) :: x0, y0
    real(8) :: total_rho = 0.0d0

    rhoAvg = 0.0d0
    fluidNum = 0
    !$omp parallel do default(none) shared(rho,obst) private(i,j) reduction(+:rhoAvg,fluidNum)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                rhoAvg = rhoAvg+rho(i,j)
                fluidNum = fluidNum+1
            endif
        enddo
    enddo
    !$omp end parallel do

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(rhoAvg, total_rho, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(fluidNum, total_fluidNum, 1, MPI_INTEGER, MPI_SUM, comm2d, rc)

    rhoAvg = total_rho / dble(total_fluidNum)

    !$omp parallel do default(none) &
    !$omp shared(ex,ey,r,omega,obst,xCenter,yCenter,rationalOmega,radius,f,f_post,Uc,Vc,rhoAvg) &
    !$omp private(i,j,alpha,ip,jp,myFlag,cNum,x0,y0,q,temp1,temp2) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then

                do alpha=0,8
                    ip = i+ex(alpha)
                    jp = j+ey(alpha)

                    if(obst(ip,jp).EQ.1) then
                        
                        myFlag = 0
                        
                        do cNum=1,cNumMax
                            if( ((ip + i_start_global-xCenter(cNum))**2.0d0+(jp + j_start_global - yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i + i_start_global),dble(j + j_start_global),alpha,x0,y0,q) 

                                ! add the contribution of rotational velocity: v = r * \omega
                                temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)   
                                temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)

!---------------------------------------------------------------------------------------------------------------------
#ifdef linear        
if(q.LT.0.5d0) then
    f(r(alpha),i,j) = 2.0d0*q*f_post(alpha,i,j) &
                     +(1.0d0-2.0d0*q)*f_post(alpha,i-ex(alpha),j-ey(alpha)) &
                     +6.0d0*omega(alpha)*rhoAvg*(ex(r(alpha))*(Uc(cNum)+temp1)+ey(r(alpha))*(Vc(cNum)+temp2))                                                    
elseif(q.GE.0.5d0) then
    f(r(alpha),i,j) = 0.5d0/q*f_post(alpha,i,j) &
                     +(1.0d0-0.50d0/q)*f_post(r(alpha),i,j) &
                     +3.0d0*omega(alpha)*rhoAvg/q*(ex(r(alpha))*(Uc(cNum)+temp1)+ey(r(alpha))*(Vc(cNum)+temp2))
endif
#endif
!---------------------------------------------------------------------------------------------------------------------
#ifdef quadratic
if(q.LT.0.5d0) then                     
    f(r(alpha),i,j) = q*(1.0d0+2.0d0*q)*f_post(alpha,i,j) &
                     +(1.0d0-4.0d0*q*q)*f_post(alpha,i-ex(alpha),j-ey(alpha)) &
                     -q*(1.0d0-2.0d0*q)*f_post(alpha,i-2*ex(alpha),j-2*ey(alpha))    &
                     +6.0d0*omega(alpha)*rhoAvg*(ex(r(alpha))*(Uc(cNum)+temp1)+ey(r(alpha))*(Vc(cNum)+temp2))                                
elseif(q.GE.0.5d0) then                     
    f(r(alpha),i,j) = f_post(alpha,i,j)/q/(1.0d0+2.0d0*q) &
                     +f_post(r(alpha),i,j)*(2.0d0*q-1.0d0)/q &
                     -f_post(r(alpha),i-ex(alpha),j-ey(alpha))*(2.0d0*q-1.0d0)/(2.0d0*q+1.0d0) &
                     +6.0d0*omega(alpha)*rhoAvg/q/(1.0d0+2.0d0*q)*(ex(r(alpha))*(Uc(cNum)+temp1)+ey(r(alpha))*(Vc(cNum)+temp2))
endif
#endif
!--------------------------------------------------------------------------------------------------------------------

                            endif
                        enddo

                        if(myFlag.EQ.0) then
                            write(*,*) "    "
                            write(*,*) "Did not find the center owning the boundary points!"
                            write(*,*) "    "
                            stop
                        endif

                    endif

                enddo

            endif
        enddo
    enddo
    !$omp end parallel do
end subroutine bounceback_particle





subroutine calQ(cNum,i,j,alpha,x0,y0,q)
    use commondata
    implicit none
    real(8) :: i, j
    real(8) :: q
    real(8) :: x0, y0
    integer :: cNum
    integer :: alpha
    real(8), parameter :: epsRadius=1e-9
    real(8) :: qTemp

    q = 0.5d0
    qTemp = 0.5d0
    x0 = i+qTemp*dble(ex(alpha))
    y0 = j+qTemp*dble(ey(alpha))
    do while( dabs(dsqrt( (x0-xCenter(cNum))**2.0d0+(y0-yCenter(cNum))**2.0d0 )-radius(cNum)).GE.epsRadius) 
        if(dsqrt((x0-xCenter(cNum))**2.0d0+(y0-yCenter(cNum))**2.0d0).GT.radius(cNum)) then
            qTemp = qTemp/2.0d0
            x0 = x0+qTemp*dble(ex(alpha))
            y0 = y0+qTemp*dble(ey(alpha))
            q = q+qTemp
        elseif(dsqrt((x0-xCenter(cNum))**2.0d0+(y0-yCenter(cNum))**2.0d0).LT.radius(cNum)) then
            qTemp = qTemp/2.0d0
            x0 = x0-qTemp*dble(ex(alpha))
            y0 = y0-qTemp*dble(ey(alpha))
            q = q-qTemp
        else
            write(*,*) "error calQ!"
            stop
        endif
    enddo

    if( (q.GT.1.0d0).OR.(q.LT.0.0d0) ) then
        write(*,*) "error q!"
        write(*,*) "q =",q
        write(*,*) 'i=',i,' j=',j,' alpha=',alpha
        stop
    endif

    return
end subroutine calQ


subroutine calQ2(cNum,i,j,alpha,x0,y0,q)
    use commondata
    implicit none
    real(8) :: i, j
    real(8) :: q
    real(8) :: x0, y0
    integer :: cNum
    integer :: alpha
    real(8), parameter :: epsRadius=1e-9
    real(8) :: qTemp
    real(8) :: coe_a, coe_b, coe_c, delta

    x0 = i
    y0 = j

    coe_a = dble(ex(alpha)) * dble(ex(alpha)) + dble(ey(alpha)) * dble(ey(alpha))

    coe_b = 2.0d0 * (dble(ex(alpha)) * (x0 - xCenter(cNum)) + dble(ey(alpha)) * (y0 - yCenter(cNum))) 
    
    coe_c = (x0 - xCenter(cNum)) ** 2 + (y0 - yCenter(cNum)) ** 2 - radius(cNum) ** 2

    delta = coe_b**2 - 4.0d0 * coe_a * coe_c

    if( delta < 0 ) then
        write(*,*) "error delta of q!"
        write(*,*) "q =",q
        write(*,*) 'i=',i,' j=',j,' alpha=',alpha
        stop
    endif

    q = (-coe_b - dsqrt(delta)) / (2.0d0 * coe_a)

    if( (q.GT.1.0d0).OR.(q.LT.0.0d0) ) then
        write(*,*) "error q!"
        write(*,*) "q =",q
        write(*,*) 'i=',i,' j=',j,' alpha=',alpha
        stop
    endif

    return

end subroutine


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j
    
#ifdef movingFrame
    ! bottom side (j = 1)   
    if (coords(1) == 0) then
        do i=1,nx
            !Bottom side
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1) + (-Uwall - U0)/6.0d0
            f(6,i,1) = f_post(8,i,1) - (-Uwall - U0)/6.0d0
        enddo
    endif

    ! Top side (j = ny)
    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            !Top side
            f(4,i,ny) = f_post(2,i,ny)
            f(7,i,ny) = f_post(5,i,ny) - (Uwall - U0)/6.0d0
            f(8,i,ny) = f_post(6,i,ny) + (Uwall - U0)/6.0d0
        enddo
    endif
#endif

#ifdef stationaryFrame
    ! bottom side (j = 1)   
    if (coords(1) == 0) then
        do i=1,nx
            !Bottom side
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)+(-Uwall)/6.0d0
            f(6,i,1) = f_post(8,i,1)-(-Uwall)/6.0d0
        enddo
    endif

    ! Top side (j = ny)
    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            !Top side
            f(4,i,ny) = f_post(2,i,ny)
            f(7,i,ny) = f_post(5,i,ny)-(Uwall)/6.0d0
            f(8,i,ny) = f_post(6,i,ny)+(Uwall)/6.0d0
        enddo
    endif
#endif

    return
end subroutine bounceback


subroutine macro()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v,obst) private(i,j) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0)  then
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j) )/rho(i,j)
                v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j) )/rho(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macro



subroutine calForce()
    use mpi
    use commondata
    implicit none
    integer :: alpha
    integer :: i, j
    integer :: ip, jp
    integer :: fluidNum
    real(8) :: temp1, temp2
    real(8) :: q
    real(8) :: x0, y0
    integer :: myFlag
    integer :: cNum
    real(8) :: tempForceX, tempForceY, tempTorque
    real(8) :: force_x(cNumMax), force_y(cNumMax), torque(cNumMax)

    do cNum=1,cNumMax
        force_x(cNumMax) = 0.0d0
        force_y(cNumMax) = 0.0d0
        torque(cNumMax) = 0.0d0
    enddo
        
    !$omp parallel do default(none) &
    !$omp shared(ex,ey,r,obst,xCenter,yCenter,rationalOmega,radius,f,f_post,Uc,Vc,rhoAvg) &
    !$omp private(i,j,alpha,ip,jp,myFlag,cNum,x0,y0,q,temp1,temp2,tempForceX,tempForceY,tempTorque) &
    !$omp reduction(+:wallTotalForceX,wallTotalForceY,totalTorque)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
        
                do alpha=1,8 !neighboring nodes, i.e., those are fluid nodes which is near solid wall
                    ip = i+ex(alpha)
                    jp = j+ey(alpha)
                        
                    if(obst(ip,jp).EQ.1) then
                                
                        myFlag = 0
                            
                        do cNum=1,cNumMax
                            if( ((ip + i_start_global-xCenter(cNum))**2.0d0+(jp + j_start_global - yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i + i_start_global),dble(j + j_start_global),alpha,x0,y0,q) 
                                    
                                temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)
                                temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)

    !---------No. 2: Wen JCP 2014----------------------------------------
    tempForceX = (ex(alpha)-Uc(cNum)-temp1)*f_post(alpha,i,j)-(ex(r(alpha))-Uc(cNum)-temp1)*f(r(alpha),i,j)
    tempForceY = (ey(alpha)-Vc(cNum)-temp2)*f_post(alpha,i,j)-(ey(r(alpha))-Vc(cNum)-temp2)*f(r(alpha),i,j)
    tempTorque = (x0-xCenter(cNum))*tempForceY-(y0-yCenter(cNum))*tempForceX 
                                force_x(cNum) = force_x(cNum) + tempForceX 
                                force_y(cNum) = force_y(cNum) + tempForceY
                                torque(cNum) = torque(cNum) + tempTorque
                            endif
                        enddo
                                
                        if(myFlag.EQ.0) then
                            write(*,*) "    "
                            write(*,*) "Did not find the center owning the boundary points! (when calculating force)"
                            write(*,*) "    "
                            stop
                        endif

                    endif

                enddo
                
            endif
        enddo
    enddo
    !$omp end parallel do


    call MPI_ALLreduce(force_x, wallTotalForceX, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(force_y, wallTotalForceY, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(torque, totalTorque, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)


if (rank2d == 0) then    
    open(unit=01,file='force.dat',status='unknown',position='append')
    write(01,*) itc, wallTotalForceX(1), wallTotalForceY(1) 
    close(01)
    
    open(unit=01,file='torque.dat',status='unknown',position='append')
    write(01,*) itc, totalTorque(1)
    close(01)
endif

    do cNum=1,cNumMax
        xCenterOld(cNum) = xCenter(cNum)
        yCenterOld(cNum) = yCenter(cNum)
        
        UcOld(cNum) = Uc(cNUm)
        VcOld(cNum) = Vc(cNum)
        
        rationalOmegaOld(cNum) = rationalOmega(cNum)
    enddo
    
    return
end subroutine calForce



subroutine check()
    use mpi
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2
    real(8) :: total_error1, total_error2

    error1 = 0.0d0
    error2 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,obst) private(i,j) reduction(+:error1,error2)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j) 
            endif
        enddo
    enddo
    !$omp end parallel do

    call MPI_ALLreduce(error1, total_error1, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(error2, total_error2, 1, MPI_REAL8, MPI_SUM, comm2d, rc)

    
    errorU = dsqrt(total_error1)/dsqrt(total_error2)

    if (rank == 0) then
        write(*,*) itc,' ',errorU
    endif
    return
end subroutine check



subroutine updateCenter()
    use mpi
    use commondata
    implicit none
    integer :: i, j
    integer :: fluidNum
    real(8) :: un(0:8)
    real(8) :: us2
    integer :: alpha
    integer :: cNum
    integer :: myFlag
    real(8) :: tempNormal
    real(8) :: outNormal
    integer :: ec
    real(8) :: m(0:8)
    real(8) :: total_rho = 0.0d0
    integer :: total_fluidNum = 0

    do cNum=1,cNumMax
#ifdef movingFrame
        Uc(cNum) = 0.0d0
        Vc(cNum) = 0.0d0
#endif
#ifdef stationaryFrame
        Uc(cNum) = U0
        Vc(cNum) = 0.0d0
#endif

        rationalOmega(cNum) = 0.0d0
        
        xCenter(cNum) = xCenterOld(cNum)+UcOld(cNum)
        yCenter(cNum) = yCenterOld(cNum)+VcOld(cNum)
    enddo

    !$omp parallel do default(none) shared(rho,u,v,obstNew,xCenter,yCenter,radius) private(i,j,cNum)
    do j=0,ny+1
        do i=0,nx+1
            obstNew(i,j) = 0 !re-initial
            do cNum=1,cNumMax
                if( ((i + i_start_global - xCenter(cNum))**2.0d0+(j + j_start_global - yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then
                    obstNew(i,j) = 1 ! solid node
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do
    do j=1,ny
        do i=1,nx
            if( obstNew(i,j) == 1) then
                rho(i,j) = rhoSolid
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0
            endif
        enddo
    enddo

    !--------------------------------------------
    rhoAvg = 0.0d0
    fluidNum = 0
    !$omp parallel do default(none) shared(rho,obstNew) private(i,j) reduction(+:rhoAvg,fluidNum)
    do j=1,ny
        do i=1,nx
            if(obstNew(i,j).EQ.0) then
                rhoAvg = rhoAvg+rho(i,j)
                fluidNum = fluidNum+1
            endif
        enddo
    enddo
    !$omp end parallel do
    
    call MPI_ALLreduce(rhoAvg, total_rho, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(fluidNum, total_fluidNum, 1, MPI_INTEGER, MPI_SUM, comm2d, rc)
    rhoAvg = total_rho / dble(total_fluidNum)
    !--------------------------------------------

    !$omp parallel do default(none) &
    !$omp shared(obst,obstNew,xCenter,yCenter,xCenterOld,yCenterOld,Uc,Vc,radius,f,rho,u,v,ex,ey,rhoAvg) &
    !$omp private(i,j,myFlag,cNum,outNormal,ec,alpha,tempNormal,m) 
    do j=1,ny
        do i=1,nx
            if( (obst(i,j).EQ.1).AND.(obstNew(i,j).EQ.0) ) then ! a solid node turns to new birth fluid node
                myFlag = 0
                do cNum=1,cNumMax
                    if( ((i+i_start_global-xCenterOld(cNum))**2.0d0+(j+j_start_global-yCenterOld(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then
                        myFlag = 1

!--------determine extrapolating direction-----------------------------
                        outNormal = 0.0d0
                        ec = 0
                        do alpha=1,8
                            tempNormal = ( ((i+i_start_global-xCenter(cNum))*ex(alpha)+(j+j_start_global-yCenter(cNum))*ey(alpha)) )/ &
                                        dsqrt( (i+i_start_global-xCenter(cNum))**2.0d0+(j+j_start_global-yCenter(cNum))**2.0d0 )
                            if( tempNormal.GT.outNormal ) then
                                outNormal =tempNormal
                                ec = alpha
                            endif
                        enddo
                        if(ec.EQ.0) then
                            write(*,*) "error determine extrapolaing direction"
                            stop
                        endif
                        ! write(*,*) "ec==",ec
!--------determine extrapolating direction-----------------------------

!---refilling using extrpolating scheme--------------------------------
! #ifdef linear
            do alpha=0,8
                f(alpha,i,j) = 2.0d0*f(alpha,i+ex(ec),j+ey(ec))-f(alpha,i+2*ex(ec),j+2*ey(ec))
            enddo
! #endif
! #ifdef quadratic
!             do alpha=0,8
!                 f(alpha,i,j) = 3.0d0*f(alpha,i+ex(ec),j+ey(ec))-3.0d0*f(alpha,i+2*ex(ec),j+2*ey(ec))+f(alpha,i+3*ex(ec),j+3*ey(ec))
!             enddo
! #endif
!---refilling using extrpolating scheme--------------------------------

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*f(5,i,j)+2.0d0*f(6,i,j)+2.0d0*f(7,i,j)+2.0d0*f(8,i,j)
    m(2) = 4.0d0*f(0,i,j)-2.0d0*f(1,i,j)-2.0d0*f(2,i,j)-2.0d0*f(3,i,j)-2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = rhoAvg*Uc(1) 
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = rhoAvg*Vc(1) 
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

    f(0,i,j) = ( m(0)-m(1)+m(2) )/9.0d0
    f(1,i,j) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0+m(3)/6.0d0-m(4)/6.0d0 &
                    +m(7)*0.25d0
    f(2,i,j) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
                    +m(5)/6.0d0-m(6)/6.0d0-m(7)*0.25d0 
    f(3,i,j) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0-m(3)/6.0d0+m(4)/6.0d0 &
                    +m(7)*0.25d0
    f(4,i,j) = m(0)/9.0d0-m(1)/36.0d0-m(2)/18.0d0 &
                    -m(5)/6.0d0+m(6)/6.0d0-m(7)*0.25d0
    f(5,i,j) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
                    +m(5)/6.0d0+m(6)/12.0d0+m(8)*0.25d0
    f(6,i,j) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
                    +m(5)/6.0d0+m(6)/12.0d0-m(8)*0.25d0
    f(7,i,j) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0-m(3)/6.0d0-m(4)/12.0d0 &
                    -m(5)/6.0d0-m(6)/12.0d0+m(8)*0.25d0
    f(8,i,j) = m(0)/9.0d0+m(1)/18.0d0+m(2)/36.0d0+m(3)/6.0d0+m(4)/12.0d0 &
                    -m(5)/6.0d0-m(6)/12.0d0-m(8)*0.25d0
                    
    rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j) )/rho(i,j)
    v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j) )/rho(i,j)

!---refilling using extrpolating scheme--------------------------------
                    endif
                enddo
                if(myFlag.EQ.0) then
                    write(*,*) 'error'
                    stop
                endif
            endif
        enddo
    enddo
    !$omp end parallel do
    
    obst = obstNew

    return
end subroutine updateCenter



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




subroutine message_passing_sendrecv()
    use mpi
    use commondata
    implicit none

    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Sendrecv(f_post(2, 1, ny), 1, f_row_x, nbr_top, 2, &
                    f_post(2, 1, 0), 1, f_row_x, nbr_bottom, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, 1, ny), 1, f_row_x, nbr_top, 5, &
                    f_post(5, 1, 0), 1, f_row_x, nbr_bottom, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(6, 1, ny), 1, f_row_x, nbr_top, 6, &
                    f_post(6, 1, 0), 1, f_row_x, nbr_bottom, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to bottom(j--)
    call MPI_Sendrecv(f_post(4, 1, 1), 1, f_row_x, nbr_bottom, 4, &
                    f_post(4, 1, ny+1), 1, f_row_x, nbr_top, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(7, 1, 1), 1, f_row_x, nbr_bottom, 7, &
                    f_post(7, 1, ny+1), 1, f_row_x, nbr_top, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, 1, 1), 1, f_row_x, nbr_bottom, 8, &
                    f_post(8, 1, ny+1), 1, f_row_x, nbr_top, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f_post(1, nx, 1), 1, f_column_y, nbr_right, 1, &
                    f_post(1, 0, 1), 1, f_column_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(5, nx, 1), 1, f_column_y, nbr_right, 5, &
                    f_post(5, 0, 1), 1, f_column_y, nbr_left, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(8, nx, 1), 1, f_column_y, nbr_right, 8, &
                    f_post(8, 0, 1), 1, f_column_y, nbr_left, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)
    ! message passing to left(i--)
    call MPI_Sendrecv(f_post(3, 1, 1), 1, f_column_y, nbr_left, 3, &
                    f_post(3, nx+1, 1), 1, f_column_y, nbr_right, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)       

    call MPI_Sendrecv(f_post(6, 1, 1), 1, f_column_y, nbr_left, 6, &
                    f_post(6, nx+1, 1), 1, f_column_y, nbr_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    call MPI_Sendrecv(f_post(7, 1, 1), 1, f_column_y, nbr_left, 7, &
                    f_post(7, nx+1, 1), 1, f_column_y, nbr_right, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)     

    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Sendrecv(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, &
                    f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Sendrecv(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, &
                    f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Sendrecv(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, &
                    f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Sendrecv(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, &
                    f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)

end subroutine message_passing_sendrecv


subroutine freeall()
    use commondata
    implicit none

    deallocate(X)
    deallocate(Y)

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)

    deallocate(f)
    deallocate(f_post)

    deallocate(obst)
    deallocate(obstNew)

end subroutine freeall


subroutine output()
    use mpi
    use commondata
    integer :: i, j
    integer :: p_rank, num(0:3) ,dx = 0, dy = 0
    real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :)
    real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :)
    ! real(8), allocatable :: stream(:, :), vorticity(:, :)
    ! integer :: total_obst(:, :), tmp_obst(:, :)

    if (rank2d > 0) then
        ! rank != 0 send data
        num(0) = nx
        num(1) = ny
        num(2) = i_start_global
        num(3) = j_start_global
        ! send to rank 0
        call MPI_Send(num, 4, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size and origion
        call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
        call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
        call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
        ! call MPI_Send(T, nx*ny, MPI_REAL8, 0, 4, comm2d, rc)
    else
        ! rank 0 collect data
        allocate(total_u(total_nx, total_ny))
        allocate(total_v(total_nx, total_ny))
        allocate(total_rho(total_nx, total_ny))
        ! allocate(total_T(total_nx, total_ny))

        dx = i_start_global
        dy = j_start_global

        ! collect data from rank 0
        do j = 1, ny
            do i = 1, nx
                total_u(dx + i, dy + j) = u(i, j)
                total_v(dx + i, dy + j) = v(i, j)
                total_rho(dx + i, dy + j) = rho(i, j)
                ! total_T(dx + i, dy + j) = T(i, j)
            enddo
        enddo

         ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) - 1
            ! receive the block size and origion
            call MPI_Recv(num, 4, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)
            ! creat buffer
            allocate(tmp_u(num(0), num(1)))
            allocate(tmp_v(num(0), num(1)))
            allocate(tmp_rho(num(0), num(1)))
            ! allocate(tmp_T(num(0), num(1)))
            
            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)
            ! call MPI_Recv(tmp_T, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 4, comm2d, MPI_STATUS_IGNORE, rc)

            ! determine the origin
            dx = num(2)
            dy = num(3)

            ! assign data
            do j = 1, num(1)
                do i = 1, num(0)
                    total_u(dx + i, dy + j) = tmp_u(i, j)
                    total_v(dx + i, dy + j) = tmp_v(i, j)
                    total_rho(dx + i, dy + j) = tmp_rho(i, j)
                    ! total_T(dx + i, dy + j) = tmp_T(i, j)
                enddo
            enddo

            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_rho)
            ! deallocate(tmp_T)
        enddo


        ! allocate(stream(total_nx, total_ny))
        ! allocate(vorticity(total_nx, total_ny))

        ! call compute_stream_vorticity(stream, vorticity, total_u, total_v, total_nx, total_ny)
        call output_Tecplot(X, Y, total_u, total_v, total_rho, total_nx, total_ny, itc)
        ! call output_binary(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        ! call output_ASCII(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        ! call out_Velocity_Nu(total_u, total_v, total_T, total_nx, total_ny, diffusivity, lengthUnit)

        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_rho)
        ! deallocate(total_T)
        ! deallocate(stream)
        ! deallocate(vorticity)

    endif


end subroutine output




subroutine output_Tecplot(xp, yp, u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(nx), yp(ny)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7,V8
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName

    
    write(B2,'(i9.9)') itc
    open(unit=41,file='movingCylinder-'//B2//'.plt',form='binary')

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

    !c-- Number of variables in this data file

    write(41) 5

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='rho'
    call dumpstring(V5)
    ! V6='T'
    ! call dumpstring(V6)
    ! V7='stream'
    ! call dumpstring(V7)
    ! V8='vorticity'
    ! call dumpstring(V8)


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
    write(41) kmax

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
    ! write(41) 1
    ! write(41) 1
    ! write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(u(i,j))
                write(41) real(v(i,j))
                write(41) real(rho(i,j))
                ! write(41) real(T(i,j))
                ! write(41) real(stream(i,j))
                ! write(41) real(vorticity(i,j))
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
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
end subroutine dumpstring


subroutine message_particle()
    use mpi
    use commondata
    implicit none
    integer :: cNum
    integer :: flag(4)  ! four corners
    integer :: js, i, j, k, idx(2), n_idx(2)
    real(8) :: lf_b, rg_b, tp_b, bt_b, len_x, len_y
    integer :: req(cNumMax * 16), cnt = 0

    idx(1) = 1
    idx(2) = 2

    ! assume effect zone is a square
    do cNum=1,cNumMax
        ! boundary of the zone
        lf_b = ceiling(xCenter(cNum) - radius(cNum) - 2)
        rg_b = floor(xCenter(cNum) + radius(cNum) + 2)
        bt_b = ceiling(yCenter(cNum) - radius(cNum) - 2)
        tp_b = floor(yCenter(cNum) + radius(cNum) + 2)
        len_x = rg_b - lf_b + 1
        len_y = tp_b - bt_b + 1

        if (i_start_global - 1 <= rg_b .AND. &
            i_start_global + nx + 1 >= lf_b .AND. &
            j_start_global - 1 <= tp_b .AND. &
            j_start_global + ny + 1 >= bt_b) then   ! particle effect sub_domain 

                ! left 1 line inside effect zone
                i = i_start_global + 1
                if (i <= rg_b .AND. i >= lf_b) then
                    ! send 1st line to left neighbor
                    call MPI_Isend(f(0, 1, bt_b), len_y, type_f_y, nbr_left, 1, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(0, 1, bt_b), len_y, type_f_y, nbr_left, 2, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! left 2 line inside effect zone
                i = i_start_global + 2
                if (i <= rg_b .AND. i >= lf_b) then
                    ! send 2nd line to left neighbor
                    call MPI_Isend(f(0, 2, bt_b), len_y, type_f_y, nbr_left, 3, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(0, 2, bt_b), len_y, type_f_y, nbr_left, 4, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! right nx line inside effect zone
                i = i_start_global + nx
                if (i <= rg_b .AND. i >= lf_b) then
                    ! send nx line to right neighbor
                    call MPI_Isend(f(0, nx, bt_b), len_y, type_f_y, nbr_right, 1, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(0, nx, bt_b), len_y, type_f_y, nbr_right, 2, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! right nx-1 line inside effect zone
                i = i_start_global + nx - 1
                if (i <= rg_b .AND. i >= lf_b) then
                    ! send nx-1 line to right neighbor
                    call MPI_Isend(f(0, nx-1, bt_b), len_y, type_f_y, nbr_right, 3, comm2d, req(cnt+1), rc)
                    call MPI_Isend(f_post(0, nx-1, bt_b), len_y, type_f_y, nbr_right, 4, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! left -1 line inside effect zone
                i = i_start_global - 1
                if (i >= lf_b .AND. i <= rg_b) then
                    ! receive -1 line from left neighbor
                    call MPI_Irecv(f(0, -1, bt_b), len_y, type_f_y, nbr_left, 1, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(0, -1, bt_b), len_y, type_f_y, nbr_left, 2, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! left -2 line inside effect zone
                i = i_start_global - 2
                if (i >= lf_b .AND. i <= rg_b) then
                    ! receive -2 line from left neighbor
                    call MPI_Irecv(f(0, -2, bt_b), len_y, type_f_y, nbr_left, 3, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(0, -2, bt_b), len_y, type_f_y, nbr_left, 4, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! right nx+1 line inside effect zone
                i = i_start_global + nx + 1
                if (i >= lf_b .AND. i <= rg_b) then
                    ! receive nx+1 line from right neighbor
                    call MPI_Irecv(f(0, nx+1, bt_b), len_y, type_f_y, nbr_left, 1, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(0, nx+1, bt_b), len_y, type_f_y, nbr_left, 2, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif

                ! right nx+2 line inside effect zone
                i = i_start_global + nx +2
                if (i >= lf_b .AND. i <= rg_b) then
                    ! receive nx+2 line from right neighbor
                    call MPI_Irecv(f(0, nx+2, bt_b), len_y, type_f_y, nbr_left, 3, comm2d, req(cnt+1), rc)
                    call MPI_Irecv(f_post(0, nx+2, bt_b), len_y, type_f_y, nbr_left, 4, comm2d, req(cnt+2), rc)
                    cnt = cnt + 2
                endif
        endif
    enddo

    call MPI_Waitall(cnt, req, MPI_STATUSES_IGNORE, rc)

end subroutine message_particle




subroutine send_all_f()
    use mpi
    use commondata
    implicit none

        ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f(0, nx, -1), ny+4, type_f_y, nbr_right, 2, &
                    f(0, 0, -1), ny+4, type_f_y, nbr_left, 2, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f(0, nx-1, -1), ny+4, type_f_y, nbr_right, 4, &
                    f(0, -1, -1), ny+4, type_f_y, nbr_left, 4, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to left(i--) 
    call MPI_Sendrecv(f(0, 1, -1), ny+4, type_f_y, nbr_left, 6, &
                    f(0, nx+1, -1), ny+4, type_f_y, nbr_right, 6, &
                    comm2d, MPI_STATUS_IGNORE, rc)          
                
    call MPI_Sendrecv(f(0, 2, -1), ny+4, type_f_y, nbr_left, 8, &
                    f(0, nx+2, -1), ny+4, type_f_y, nbr_right, 8, &
                    comm2d, MPI_STATUS_IGNORE, rc)   

end subroutine send_all_f


subroutine send_all_fp()
    use mpi
    use commondata
    implicit none

        ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Sendrecv(f_post(0, nx, -1), ny+4, type_f_y, nbr_right, 1, &
                    f_post(0, 0, -1), ny+4, type_f_y, nbr_left, 1, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    call MPI_Sendrecv(f_post(0, nx-1, -1), ny+4, type_f_y, nbr_right, 3, &
                    f_post(0, -1, -1), ny+4, type_f_y, nbr_left, 3, &
                    comm2d, MPI_STATUS_IGNORE, rc)

    ! message passing to left(i--)
    call MPI_Sendrecv(f_post(0, 1, -1), ny+4, type_f_y, nbr_left, 5, &
                    f_post(0, nx+1, -1), ny+4, type_f_y, nbr_right, 5, &
                    comm2d, MPI_STATUS_IGNORE, rc)    

    call MPI_Sendrecv(f_post(0, 2, -1), ny+4, type_f_y, nbr_left, 7, &
                    f_post(0, nx+2, -1), ny+4, type_f_y, nbr_right, 7, &
                    comm2d, MPI_STATUS_IGNORE, rc)     
end subroutine send_all_fp