
module commondata
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)

    integer, parameter :: nx=201,ny=101
    integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
    integer, parameter :: nxFourth=(nx-1)/4+1, nyFourth=(ny-1)/4+1
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
    real(8) :: dx, dy, dt
    real(8), parameter :: rho0=1.0d0
    real(8), parameter :: rhoSolid=2.0d0
    real(8) :: rhoAvg
    
    real(8), parameter :: U0=0.02d0 
    real(8), parameter :: Uwall=0.1d0
    real(8), parameter :: viscosity=1.0d0/9.0d0
    real(8), parameter :: tauf=3.0d0*viscosity+0.5d0
    real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
    real(8), parameter :: shearRate=2.0d0*Uwall/dble(ny)
    real(8), parameter :: Reynolds=shearRate*(2.0d0*radius0)**2.0d0/viscosity

    real(8) :: errorU
    real(8), parameter :: eps=1e-6
    real(8) :: X(nx), Y(ny)
    real(8) :: u(nx,ny), v(nx,ny), rho(nx,ny)
    real(8) :: up(nx,ny), vp(nx,ny)
    real(8) :: f(0:8,nx,ny), f_post(0:8,0:nx+1,0:ny+1)
    real(8) :: omega(0:8)
    integer :: ex(0:8), ey(0:8)
    data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
    data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/
    integer :: obst(0:nx+1,0:ny+1)
    integer :: obstNew(0:nx+1,0:ny+1)
    integer :: r(1:8)
    data r/3, 4, 1, 2, 7, 8, 5, 6/
    
    real(8) :: wallTotalForceX(cNumMax), wallTotalForceY(cNumMax)
    real(8) :: totalTorque(cNumMax)
    
end module commondata


! #define movingFrame
#define stationaryFrame
    
#define linear
! #define quadratic

program main  
    use omp_lib   
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads

    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------
#ifdef _OPENMP
    call OMP_set_num_threads(4)
    write(*,*) "Start OpenMP......"
    myMaxThreads = OMP_get_max_threads()
    write(*,*) "Running threads =",myMaxThreads
#endif
    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------

    call initial()
    call output_Tecplot()

    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    do while( (errorU.GT.eps).AND.(itc.LE.itc_max) )

        call collision()

        call streaming()

        call macro()

        call calForce()

        itc = itc+1
        
        if(MOD(itc,1000).EQ.0) then
            call check()
        endif
        !if(MOD(itc,2000).EQ.0) then
        !    call output_Tecplot()
        !endif

        call updateCenter()
        
    enddo

    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif
    
    write(*,*) "Time (CPU) = ",finish-start, "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", finish2-start2, "s"
#endif

    itc = itc+1
    call output_Tecplot()
    
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

    itc = 0
    errorU = 100.0d0

    dx = 1.0d0
    dy = 1.0d0
    dt = 1.0d0
    write(*,*) "Reynolds=",real(Reynolds),",   tauf=",real(tauf)
    radius(1) = dble(radius0)
    do cNum=1,cNumMax
        write(*,*) "diameter=",real(2.0d0*radius(cNum))
    enddo
    write(*,*) "U0 (cylinder)=",real(U0)
    write(*,*) "Uwall=",real(Uwall)
    write(*,*) "    "

    do i=1,nx
        X(i) = (i-1)*dx
    enddo
    do j=1,ny
        Y(j) = (j-1)*dy
    enddo

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
    do j=1,ny
        do i=1,nx
            rho(i,j) = rho0
            do cNum=1,cNumMax
                if( ((i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then 
                    obst(i,j) = 1 ! solid node
                    rho(i,j) = rhoSolid
                endif
            enddo
        enddo
    enddo

#ifdef stationaryFrame
    u = U0
    do i=1,nx
        u(i,1) = -Uwall
        u(i,ny) = Uwall
    enddo
#endif

#ifdef movingFrame
    u = 0.0d0
    do i=1,nx
        u(i,1) = -Uwall - U0
        u(i,ny) = Uwall - U0
    enddo
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
    real(8) :: q
    integer :: cNum
    integer :: myFlag
    integer :: fluidNum
    real(8) :: temp1, temp2
    real(8) :: x0, y0

    ! !$omp parallel do default(none) shared(ex,ey,f,f_post,obst) private(i,j,alpha,ip,jp)
    ! do j=1,ny
    !     do i=1,nx
    !         if(obst(i,j).EQ.0) then
    !             do alpha=0,8
    !                 ip = i+ex(alpha)
    !                 jp = j+ey(alpha)

    !                 if(ip.EQ.0) ip = nx
    !                 if(ip.EQ.nx+1) ip = 1

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
                    if(ip.EQ.0) ip = nx
                    if(ip.EQ.nx+1) ip = 1

                if(obst(ip,jp).EQ.0) then
                    f(alpha,i,j) = f_post(alpha,ip,jp)
                endif
            enddo
        enddo
    enddo
    
    call bounceback()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    rhoAvg = rhoAvg/dble(fluidNum)
    
    ! do j=-1,ny+2
    !     do alpha=0,8
    !         f(alpha,0,j) = omega(alpha)*rho0 
    !         f(alpha,-1,j) = omega(alpha)*rho0
    !         f_post(alpha,0,j) = omega(alpha)*rho0 
    !         f_post(alpha,-1,j) = omega(alpha)*rho0

    !         f(alpha,nx+1,j) = omega(alpha)*rho0 
    !         f(alpha,nx+2,j) = omega(alpha)*rho0
    !         f_post(alpha,nx+1,j) = omega(alpha)*rho0 
    !         f_post(alpha,nx+2,j) = omega(alpha)*rho0             
    !     enddo
    ! enddo
    
    ! do i=-1,nx+2
    !     do alpha=0,8
    !         f(alpha,i,0) = omega(alpha)*rho0 
    !         f(alpha,i,-1) = omega(alpha)*rho0
    !         f_post(alpha,i,0) = omega(alpha)*rho0 
    !         f_post(alpha,i,-1) = omega(alpha)*rho0

    !         f(alpha,i,ny+1) = omega(alpha)*rho0 
    !         f(alpha,i,ny+2) = omega(alpha)*rho0
    !         f_post(alpha,i,ny+1) = omega(alpha)*rho0 
    !         f_post(alpha,i,ny+2) = omega(alpha)*rho0             
    !     enddo
    ! enddo

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
                            if( ((ip-xCenter(cNum))**2.0d0+(jp-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i),dble(j),alpha,x0,y0,q) 

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

    return
end subroutine streaming


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
    do i=1,nx
        !Top side
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny) - (Uwall - U0)/6.0d0
        f(8,i,ny) = f_post(6,i,ny) + (Uwall - U0)/6.0d0

        !Bottom side
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1) + (-Uwall - U0)/6.0d0
        f(6,i,1) = f_post(8,i,1) - (-Uwall - U0)/6.0d0
    enddo
#endif

#ifdef stationaryFrame
    do i=1,nx
        !Top side
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)-(Uwall)/6.0d0
        f(8,i,ny) = f_post(6,i,ny)+(Uwall)/6.0d0

        !Bottom side
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)+(-Uwall)/6.0d0
        f(6,i,1) = f_post(8,i,1)-(-Uwall)/6.0d0
    enddo
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

    do cNum=1,cNumMax
        wallTotalForceX(cNum) = 0.0d0
        wallTotalForceY(cNum) = 0.0d0
        totalTorque(cNum) = 0.0d0
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
                            if( ((ip-xCenter(cNum))**2.0d0+(jp-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i),dble(j),alpha,x0,y0,q) 
                                    
                                temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)
                                temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)

    !---------No. 2: Wen JCP 2014----------------------------------------
    tempForceX = (ex(alpha)-Uc(cNum)-temp1)*f_post(alpha,i,j)-(ex(r(alpha))-Uc(cNum)-temp1)*f(r(alpha),i,j)
    tempForceY = (ey(alpha)-Vc(cNum)-temp2)*f_post(alpha,i,j)-(ey(r(alpha))-Vc(cNum)-temp2)*f(r(alpha),i,j)
    tempTorque = (x0-xCenter(cNum))*tempForceY-(y0-yCenter(cNum))*tempForceX 
                           
                                wallTotalForceX(cNum) = wallTotalForceX(cNum)+tempForceX 
                                wallTotalForceY(cNum) = wallTotalForceY(cNum)+tempForceY 
                                totalTorque(cNum) = totalTorque(cNum)+tempTorque 
                                    
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
    
    open(unit=01,file='force.dat',status='unknown',position='append')
    write(01,*) itc, wallTotalForceX(1), wallTotalForceY(1) 
    close(01)
    
    open(unit=01,file='torque.dat',status='unknown',position='append')
    write(01,*) itc, totalTorque(1)
    close(01)

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
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2

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
    
    errorU = dsqrt(error1)/dsqrt(error2)

    write(*,*) itc,' ',errorU

    return
end subroutine check



subroutine updateCenter()
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
    do j=1,ny
        do i=1,nx
            obstNew(i,j) = 0 !re-initial
            do cNum=1,cNumMax
                if( ((i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then
                    obstNew(i,j) = 1 ! solid node
                    rho(i,j) = rhoSolid
                    u(i,j) = 0.0d0
                    v(i,j) = 0.0d0
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do

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
    rhoAvg = rhoAvg/dble(fluidNum)
    !--------------------------------------------

    !$omp parallel do default(none) &
    !$omp shared(obst,obstNew,xCenter,yCenter,xCenterOld,yCenterOld,Uc,Vc,radius,f,rho,u,v,ex,ey,rhoAvg) &
    !$omp private(i,j,myFlag,cNum,outNormal,ec,alpha,tempNormal,m) 
    do j=1,ny
        do i=1,nx
            if( (obst(i,j).EQ.1).AND.(obstNew(i,j).EQ.0) ) then ! a solid node turns to new birth fluid node
                myFlag = 0
                do cNum=1,cNumMax
                    if( ((i-xCenterOld(cNum))**2.0d0+(j-yCenterOld(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then
                        myFlag = 1

!--------determine extrapolating direction-----------------------------
                        outNormal = 0.0d0
                        ec = 0
                        do alpha=1,8
                            tempNormal = ( ((i-xCenter(cNum))*ex(alpha)+(j-yCenter(cNum))*ey(alpha)) )/ &
                                        dsqrt( (i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0 )
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
#ifdef linear
            do alpha=0,8
                f(alpha,i,j) = 2.0d0*f(alpha,i+ex(ec),j+ey(ec))-f(alpha,i+2*ex(ec),j+2*ey(ec))
            enddo
#endif
#ifdef quadratic
            do alpha=0,8
                f(alpha,i,j) = 3.0d0*f(alpha,i+ex(ec),j+ey(ec))-3.0d0*f(alpha,i+2*ex(ec),j+2*ey(ec))+f(alpha,i+3*ex(ec),j+3*ey(ec))
            enddo
#endif
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


subroutine output_Tecplot()
    use commondata
    implicit none
    integer :: i, j, k

    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6
    integer, parameter :: kmax=1
    character(len=40) :: zoneName

    write(B2,'(i9.9)') itc
    open(41,file='movingCylinder-'//B2//'.plt',form='binary')

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

    !c-- Number of variables in this data file (here 6 variables)
    write(41) 5

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    ! V3='Obst'
    ! call dumpstring(V3)
    V4='U'
    call dumpstring(V4)
    V5='V'
    call dumpstring(V5)
    V6='rho'
    call dumpstring(V6)
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

    !---------Specify Var Location. 0 = Donft specify, all data
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
    ! write(41) 1
    write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(X(i))
                write(41) real(Y(j))
                ! write(41) real(obst(i,j))
                write(41) real(U(i,j))
                write(41) real(V(i,j))
                write(41) real(rho(i,j))
            end do
        end do
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
