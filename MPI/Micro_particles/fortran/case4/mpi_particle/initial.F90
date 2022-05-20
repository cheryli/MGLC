subroutine initial()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
    real(8) :: us2
    integer :: cNum
    real(8) :: tmpx, tmpy, rdx, rdy
    integer :: particle_num, particle_start
    integer, allocatable :: seed(:)
    integer :: n

    
if(rank == 0) then

    write(*,*) "l0=",real(l0),"cm"
    write(*,*) "t0=",real(t0),"s"
    write(*,*) "Mesh:",nx,ny
    write(*,*) "    "
    write(*,*) "tMax=",real(tMax),"s"
    write(*,*) "itcMax=",itc_max
    write(*,*) "    "
    write(*,*) "tauf=",real(tauf)
    write(*,*) "gravity=",real(gravity)
    write(*,*) "    "
    write(*,*) "thresholdWall=",real(thresholdWall)
    write(*,*) "stiffWall=",real(stiffWall)
    write(*,*) "    "
    write(*,*) "thresholdParticle=",real(thresholdParticle)
    write(*,*) "stiffParticle=",real(stiffParticle)
    write(*,*) "    "

endif
    
    itc = 0
    errorU = 100.0d0

    ! define mesh

    do i=1,total_nx
        X(i) = dble(i-1)
    enddo
    do j=1,total_ny
        Y(j) = dble(j-1)
    enddo

    tmpx = 25.0d0
    tmpy = 25.0d0

    ! repeated pseudo-random number for debug 
    call random_seed(size = n)
    allocate(seed(n))
    do i = 1, n
      seed(i) = 1
    enddo
    call random_seed(put=seed)

    do cNum = 1, cNumMax
        call random_number(rdx)
        call random_number(rdy)

        xCenter(cNum) = tmpx + (rdx - 0.5d0) * 20
        yCenter(cNum) = tmpy + (rdy - 0.5d0) * 20

        tmpx = tmpx + 50.0d0

        if (tmpx > 200.0d0) then
            tmpx = 25.0d0
            tmpy = tmpy + 50.0d0
        endif
    enddo

    deallocate(seed)

    ! xCenter(1) = 101.01d0
    ! yCenter(1) = 720.00d0
    ! xCenter(2) = 101.0d0
    ! yCenter(2) = 680.0d0

    xCenterOld = xCenter
    yCenterOld = yCenter

    radius = dble(radius0)

    if(rank == 0) then
        do cNum=1,cNumMax
            write(*,*) "I am particle ", cNum
            write(*,*) "xCenter =", real(xCenter(cNum)), ",    yCenter =", real(yCenter(cNum))
            write(*,*) "diameter=",real(2.0d0*radius(cNum))
            write(*,*) "    "
        enddo
    endif

    Uc = 0.0d0
    Vc = 0.0d0
    rationalOmega = 0.0d0

    UcOld = 0.0d0
    VcOld = 0.0d0
    rationalOmegaOld = 0.0d0

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


    u = 0.0d0
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
    f_post = 0.0d0
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

    ! ghost points
    do j=-1,ny+2
        do alpha=0,8
            f(alpha,0,j) = omega(alpha)*rho0 
            f(alpha,-1,j) = omega(alpha)*rho0

            f(alpha,nx+1,j) = omega(alpha)*rho0 
            f(alpha,nx+2,j) = omega(alpha)*rho0 
            
            f_post(alpha,0,j) = omega(alpha)*rho0 
            f_post(alpha,-1,j) = omega(alpha)*rho0

            f_post(alpha,nx+1,j) = omega(alpha)*rho0 
            f_post(alpha,nx+2,j) = omega(alpha)*rho0 
        enddo
    enddo
    
    do i=-1,nx+2
        do alpha=0,8
            f(alpha,i,0) = omega(alpha)*rho0 
            f(alpha,i,-1) = omega(alpha)*rho0

            f(alpha,i,ny+1) = omega(alpha)*rho0 
            f(alpha,i,ny+2) = omega(alpha)*rho0 
            
            f_post(alpha,i,0) = omega(alpha)*rho0 
            f_post(alpha,i,-1) = omega(alpha)*rho0

            f_post(alpha,i,ny+1) = omega(alpha)*rho0 
            f_post(alpha,i,ny+2) = omega(alpha)*rho0 
        enddo
    enddo
    

    ! determine particle masks
    local_mask = 0

    do cNum = 1, cNumMax
        if (MOD(cNum-1, num_process) == rank) then
            local_mask(cNum) = 1
        endif
    enddo



    return
end subroutine initial