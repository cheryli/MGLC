subroutine initial()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
    real(8) :: us2
    integer :: cNum

    itc = 0
    errorU = 100.0d0

    write(*,*) "l0=",real(l0),"cm"
    write(*,*) "t0=",real(t0),"s"
    write(*,*) "Mesh:",nx,ny
    write(*,*) "    "
    write(*,*) "tMax=",real(tMax),"s"
    write(*,*) "itcMax=",itc_max
    write(*,*) "    "
    write(*,*) "tauf=",real(tauf)
    radius(1) = dble(radius0)
    do cNum=1,cNumMax
        write(*,*) "diameter=",real(2.0d0*radius(cNum))
    enddo
    write(*,*) "gravity=",real(gravity)
    write(*,*) "    "
    write(*,*) "wallStiff=",real(wallStiff)
    write(*,*) "thresholdR=",real(thresholdR)
    write(*,*) "    "
    
    do i=1,nx
        X(i) = dble(i-1)
    enddo
    do j=1,ny
        Y(j) = dble(j-1)
    enddo

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
                if( ((i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then ! circular cylinder
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

    return
end subroutine initial