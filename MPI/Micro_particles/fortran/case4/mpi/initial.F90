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

    xCenter(1) = 101.01d0
    yCenter(1) = 720.00d0
    xCenter(2) = 101.0d0
    yCenter(2) = 680.0d0

    xCenterOld(1) = xCenter(1)
    xCenterOld(2) = xCenter(2)
    yCenterOld(1) = yCenter(1)
    yCenterOld(2) = yCenter(2)

    radius(1) = dble(radius0)
    radius(2) = dble(radius0)
    if(rank == 0) then
        do cNum=1,cNumMax
            write(*,*) "I am particle ", cNum
            write(*,*) "xCenter =", real(xCenter(cNum)), ",    yCenter =", real(yCenter(cNum))
            write(*,*) "diameter=",real(2.0d0*radius(cNum))
            write(*,*) "    "
        enddo
    endif

    Uc(1) = 0.0d0
    Vc(1) = 0.0d0
    rationalOmega(1) = 0.0d0
    Uc(2) = 0.0d0
    Vc(2) = 0.0d0
    rationalOmega(2) = 0.0d0
    UcOld(1) = Uc(1)
    VcOld(1) = Vc(1)
    UcOld(2) = UcOld(2)
    VcOld(2) = VcOld(2)
    rationalOmegaOld(1) = rationalOmega(1)
    rationalOmegaOld(2) = rationalOmega(2)

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

    return
end subroutine initial