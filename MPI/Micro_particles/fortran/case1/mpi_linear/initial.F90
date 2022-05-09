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

    allocate(f(0:8, 1:nx, 1:ny))
    allocate(f_post(0:8, 0:nx+1, 0:ny+1))

    
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
    do j=1,ny
        do i=1,nx
            do cNum=1,cNumMax
                ! golbal coordinate
                if( ((i+i_start_global-xCenter(cNum))**2.0d0 + (j+j_start_global-yCenter(cNum))**2.0d0) &
                                .LE.radius(cNum)**2.0d0 ) then 
                    
                    rho(i,j) = rhoSolid
                endif
                ! if( ((i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then 
                !     obst(i,j) = 1 ! solid node
                !     rho(i,j) = rhoSolid
                ! endif
            enddo
        enddo
    enddo

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