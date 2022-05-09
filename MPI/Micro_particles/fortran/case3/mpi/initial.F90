subroutine initial()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
    real(8) :: us2
    integer :: cNum

    
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
    write(*,*) "wallStiff=",real(wallStiff)
    write(*,*) "thresholdR=",real(thresholdR)
    write(*,*) "    "

endif
    
    itc = 0
    errorU = 100.0d0

    radius(1) = dble(radius0)
    if(rank == 0) then
        do cNum=1,cNumMax
            write(*,*) "diameter=",real(2.0d0*radius(cNum))
        enddo
    endif

    ! define mesh
    if (rank == 0) then
        do i=1,total_nx
            X(i) = dble(i-1)
        enddo
        do j=1,total_ny
            Y(j) = dble(j-1)
        enddo
    endif

    xCenter(1) = 101.0d0
    yCenter(1) = 401.0d0
    Uc(1) = 0.0d0
    Vc(1) = 0.0d0
    rationalOmega(1) = 0.0d0

    obst = 0
    obstNew = 0
    do j=1,ny
        do i=1,nx
            rho(i,j) = rho0
            do cNum=1,cNumMax
                ! golbal coordinate
                if( ((i+i_start_global-xCenter(cNum))**2.0d0+(j+j_start_global-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then ! circular cylinder
                    obst(i,j) = 1 ! solid node
                    rho(i,j) = rhoSolid
                endif
            enddo
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

    return
end subroutine initial