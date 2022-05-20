
    module commondata
        implicit none
        real(8), parameter :: Pi=4.0d0*datan(1.0d0)

        integer, parameter :: nx=121,ny=961
        integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
        integer, parameter :: nxFourth=(nx-1)/4+1, nyFourth=(ny-1)/4+1

        real(8), parameter :: viscosity=0.05d0
        real(8), parameter :: tauf=3.0d0*viscosity+0.5d0
        real(8), parameter :: Snu=1.0d0/tauf, Sq=8.0d0*(2.0d0*tauf-1.0d0)/(8.0d0*tauf-1.0d0)
        
        real(8), parameter :: blockageRatio=16.0d0/13.0d0
        real(8), parameter :: aspectRatio=2.0d0
        
        real(8), parameter :: diameter0=dble(nx)/blockageRatio
        real(8), parameter :: radius0=diameter0/2.0d0
        real(8), parameter :: charaLength=diameter0
        
        real(8), parameter :: l_phys=0.1d0 !! unit: cm !! length of major axis
        real(8), parameter :: nu_phys=0.01d0  !! unit: cm^2/s 

        real(8), parameter :: l0=l_phys/dble(diameter0)  
        real(8), parameter :: t0=viscosity*l0**2.0d0/nu_phys 
        real(8), parameter :: tMax=8.0d0 !! unit: s

        integer :: itc
        integer, parameter :: itc_max=INT(tMax/t0)
        integer, parameter :: cNumMax=1
        
        real(8) :: xCenter(cNumMax), yCenter(cNumMax)
        real(8) :: xCenterOld(cNumMax), yCenterOld(cNumMax) 
        
        real(8) :: theta(cNumMax)
        real(8) :: thetaOld(cNumMax)
        
        real(8) :: rationalOmega(cNumMax)
        real(8) :: rationalOmegaOld(cNumMax)
        
        real(8) :: Uc(cNumMax), Vc(cNumMax)
        real(8) :: UcOld(cNumMax), VcOld(cNumMax)
        
        real(8) :: radiusA(cNumMax), radiusB(cNumMax)
        
        real(8), parameter :: rho0=1.0d0
        real(8), parameter :: rhoSolid=1.1d0
        real(8) :: rhoAvg

        real(8) :: errorU
        real(8), parameter :: eps=1e-6
        real(8) :: X(nx), Y(ny)
        real(8) :: u(nx,ny), v(nx,ny), rho(nx,ny)
        real(8) :: up(nx,ny), vp(nx,ny)
        real(8) :: f(0:8,-1:nx+2,-1:ny+2), f_post(0:8,-1:nx+2,-1:ny+2)
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

        real(8), parameter :: gravity=980.0d0*t0**2.0d0/l0
        real(8), parameter :: Archimedes=dsqrt(gravity*charaLength**3/viscosity**2.0d0*(rhoSolid-rho0)/rho0)
        
        real(8), parameter :: thresholdR=1.0d0
        real(8), parameter :: wallStiff=0.25d0
        real(8), parameter :: thresholdRTemp=4.0d0
        integer :: shiftTimes
    end module commondata


    program main
    use omp_lib     
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads

    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------
#ifdef _OPENMP
    call OMP_set_num_threads(24)
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

    do while( (errorU.GT.eps).AND.(itc.LT.itc_max) )

        call collision()

        call streaming()

        call macro()

        call calForce()
        
        itc = itc+1

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
        !if(MOD(itc,50000).EQ.0) then
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

    itc = 0
    errorU = 100.0d0

    write(*,*) "l0=",real(l0),"cm"
    write(*,*) "t0=",real(t0),"s"
    write(*,*) "tMax=",real(tMax),"s"
    write(*,*) "Mesh:",nx,ny
    write(*,*) "itcMax=",itc_max
    write(*,*) "tauf=",real(tauf)
    write(*,*) "gravity=",real(gravity)
    write(*,*) "    "
    write(*,*) "Major axis=",real(l_phys),"cm"
    write(*,*) "aspect ratio=",real(aspectRatio)
    write(*,*) "|-----------------------------------------|"
    write(*,*) "  blockage ratio=",real(blockageRatio)
    write(*,*) "|-----------------------------------------|"
    write(*,*) "density ratio=",real(rhoSolid/rho0)
    write(*,*) "Archimedes=",real(Archimedes)
    write(*,*) "   "
    write(*,*) "thresholdR=",real(thresholdR)
    write(*,*) "wallStiff=",real(wallStiff)
    write(*,*) "    "
    
    do i=1,nx
        X(i) = dble(i-1)
    enddo
    do j=1,ny
        Y(j) = dble(j-1)
    enddo

    radiusA(1) = radius0
    radiusB(1) = radius0/aspectRatio
    xCenter(1) = dble(nxHalf)
    yCenter(1) = dble(nyHalf)
    !--wide channel: Pi/4.0d0; --narrow channel: Pi/3.0d0
    if(blockageRatio.GE.3.0d0) then
        theta(1) = Pi/4.0d0
    elseif(blockageRatio.LT.3.0d0) then
        theta(1) = Pi/3.0d0
    endif
    write(*,*) "theta(1)=",real(theta(1)/Pi*180.0d0),"degree"
    Uc(1) = 0.0d0
    Vc(1) = 0.0d0
    rationalOmega(1) = 0.0d0

    obst = 0
    obstNew = 0
    do j=1,ny
        do i=1,nx
            rho(i,j) = rho0
            do cNum=1,cNumMax
                
       if( (((i-xCenter(cNum))*dcos(theta(cNum))+(j-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
          +(-(i-xCenter(cNum))*dsin(theta(cNum))+(j-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then
        
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


    subroutine collision()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: m(0:8), m_post(0:8), meq(0:8)
    real(8) :: s(0:8)
    real(8) :: positionShift
    integer :: cNum
    
    positionShift = int(nyHalf-yCenter(1))
    if(abs(positionShift).GE.2) then
            write(*,*) "solid particle moving too fast!!"
            write(*,*) "positionShift=",real(positionShift)
            call output_Tecplot()
            stop
    endif
    if(positionShift.EQ.1) then
        do cNum=1,cNumMax
            yCenter(cNum) = yCenter(cNum)+1.0d0
        enddo
        shiftTimes = shiftTimes+1
        do j=ny-1,1,-1
            do i=1,nx
                do alpha=0,8
                    f(alpha,i,j+1) = f(alpha,i,j) 
                enddo
                u(i,j+1) = u(i,j)
                v(i,j+1) = v(i,j)
                rho(i,j+1) = rho(i,j)
                obst(i,j+1) = obst(i,j)
            enddo
        enddo
        j = 1
        do i=1,nx
            u(i,j) = 0.0d0
            v(i,j) = 0.0d0
            rho(i,j) = rho0
            obst(i,j) = 0
            do alpha=0,8
                f(alpha,i,j) = omega(alpha)*rho0
            enddo
        enddo
    elseif(positionShift.EQ.-1) then
        do cNum=1,cNumMax
            yCenter(cNum) = yCenter(cNum)-1.0d0
        enddo
        shiftTimes = shiftTimes-1
        do j=2,ny
            do i=1,nx
                do alpha=0,8
                    f(alpha,i,j-1) = f(alpha,i,j)
                enddo
                u(i,j-1) = u(i,j)
                v(i,j-1) = v(i,j)
                rho(i,j-1) = rho(i,j)
                obst(i,j-1) = obst(i,j)
            enddo
        enddo
        j = ny
        do i=1,nx
            u(i,j) = 0.0d0
            v(i,j) = 0.0d0
            rho(i,j) = rho0
            obst(i,j) = 0
            do alpha=0,8
                f(alpha,i,j) = omega(alpha)*rho0
            enddo
        enddo
    endif

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,obst) private(i,j,alpha,s,m,m_post,meq) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*f(5,i,j)+2.0d0*f(6,i,j)+2.0d0*f(7,i,j)+2.0d0*f(8,i,j)
    m(2) = 4.0d0*f(0,i,j)-2.0d0*f(1,i,j)-2.0d0*f(2,i,j)-2.0d0*f(3,i,j)-2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
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
    
    !$omp parallel do default(none) shared(ex,ey,f,f_post,obst) private(i,j,alpha,ip,jp)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=0,8
                    ip = i+ex(alpha)
                    jp = j+ey(alpha)

                    f(alpha,ip,jp) = f_post(alpha,i,j)
                
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do

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

    !$omp parallel do default(none) &
    !$omp shared(ex,ey,r,omega,obst,xCenter,yCenter,rationalOmega,radiusA,radiusB,theta,f,f_post,Uc,Vc,rhoAvg) &
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
                        
       if( (((ip-xCenter(cNum))*dcos(theta(cNum))+(jp-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
          +(-(ip-xCenter(cNum))*dsin(theta(cNum))+(jp-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then
        
                            myFlag = 1

                            call calQellipsoid(cNum,dble(i),dble(j),alpha,x0,y0,q) 

                            temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)
                            temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)
                                
!--------------------------------------------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------------------------------------------

        endif
                        enddo

                        if(myFlag.EQ.0) then
                                write(*,*) "    "
                                write(*,*) "Did not find the center owning the boundary points!"
                                write(*,*) "I am subroutine Streaming"
                                write(*,*) "i=",i,"     ,j=",j
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


    subroutine calQellipsoid(cNum,i,j,alpha,x0,y0,q)
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
        
    do while( dabs(((x0-xCenter(cNum))*dcos(theta(cNum))+(y0-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(x0-xCenter(cNum))*dsin(theta(cNum))+(y0-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0-1.0d0).GE.epsRadius )
        
        if( (((x0-xCenter(cNum))*dcos(theta(cNum))+(y0-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(x0-xCenter(cNum))*dsin(theta(cNum))+(y0-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).GT.1.0d0 ) then
            qTemp = qTemp/2.0d0
            x0 = x0+qTemp*dble(ex(alpha))
            y0 = y0+qTemp*dble(ey(alpha))
            q = q+qTemp
        elseif( (((x0-xCenter(cNum))*dcos(theta(cNum))+(y0-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(x0-xCenter(cNum))*dsin(theta(cNum))+(y0-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LT.1.0d0 ) then
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
    end subroutine calQellipsoid


    subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j
    real(8) :: Utypical

    do j=1,ny
        !Left side
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo

    do i=1,nx
        !Bottom side
        j = 1
        f(2,i,j) = f_post(4,i,j)  
        f(5,i,j) = f_post(7,i,j)  
        f(6,i,j) = f_post(8,i,j) 
    enddo
    
    Utypical = 0.0d0
    do i=1,nx
        j = ny-1
        Utypical = Utypical+v(i,j)
    enddo
    Utypical = Utypical/dble(nx)
    do i=1,nx
        ! Top Side
        j = ny
        f(4,i,j) = (f(4,i,j)+Utypical*f(4,i,j-1))/(Utypical+1.0d0)
        f(7,i,j) = (f(7,i,j)+Utypical*f(7,i,j-1))/(Utypical+1.0d0)
        f(8,i,j) = (f(8,i,j)+Utypical*f(8,i,j-1))/(Utypical+1.0d0)
    enddo

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
    
    real(8) :: dwLeft, dwRight
    real(8) :: Fwxij(cNumMax), Fwyij(cNumMax), Twij(cNumMax)
    real(8) :: forceScale
    
    real(8) :: i_minCoarse, i_maxCoarse
    integer :: j_minStart, j_minEnd
    integer :: j_maxStart, j_maxEnd
    real(8) :: i_min, j_min
    real(8) :: i_max, j_max
    integer :: temp_countX, temp_countY
    real(8) :: x_temp, y_temp
    real(8) :: tempEllip1, tempEllip2

    do cNum=1,cNumMax
        wallTotalForceX(cNum) = 0.0d0
        wallTotalForceY(cNum) = 0.0d0
        totalTorque(cNum) = 0.0d0
    enddo

    !$omp parallel do default(none) &
    !$omp shared(ex,ey,r,obst,xCenter,yCenter,rationalOmega,radiusA,radiusB,theta,f,f_post,Uc,Vc,rhoAvg) &
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

       if( (((ip-xCenter(cNum))*dcos(theta(cNum))+(jp-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(ip-xCenter(cNum))*dsin(theta(cNum))+(jp-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then
        
                            myFlag = 1

                            call calQellipsoid(cNum,dble(i),dble(j),alpha,x0,y0,q) 
                                
                            temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)
                            temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)

    !---------No. 2: Wen JCP 2014-----------------------------------------
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
    
    do cNum=1,cNumMax
!---------------particle-wall (spring force model)------------------
        forceScale = Pi*radiusA(cNum)*radiusB(cNum)*(rhoSolid-rhoAvg)*gravity/wallStiff

        Fwxij(cNum) = 0.0d0
        Fwyij(cNum) = 0.0d0
        Twij(cNum) = 0.0d0
            
        i_minCoarse = dble(nx)
        i_maxCoarse = 1.0d0
        do j=1,ny
            do i=1,nx
                if(obst(i,j).EQ.1) then
    !--Use with caution, should be extend to N-particles
    !!--if( (((i-xCenter(cNum))*dcos(theta(cNum))+(j-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
    !!--  +(-(i-xCenter(cNum))*dsin(theta(cNum))+(j-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then !! for the cNum-th particle
                    if(dble(i).LE.i_minCoarse) then
                        i_minCoarse = dble(i)
                    elseif(dble(i).GE.i_maxCoarse) then
                        i_maxCoarse = dble(i)
                    endif
    !!--endif
                endif
            enddo
        enddo
        i_min = i_minCoarse
        i_max = i_maxCoarse
        
        j_minStart = ny
        j_minEnd = 1
        j_maxStart = ny
        j_maxEnd = 1
        do j=1,ny
            i = int(i_minCoarse)
            if(obst(i,j).EQ.1) then
                if(j.LT.j_minStart) then
                    j_minStart = j
                elseif(j.GT.j_minEnd) then
                    j_minEnd = j
                endif
            endif
            
            i = int(i_maxCoarse)
            if(obst(i,j).EQ.1) then
                if(j.LT.j_maxStart) then
                    j_maxStart = j
                elseif(j.GT.j_maxEnd) then
                    j_maxEnd = j
                endif
            endif
        enddo
        
        dwLeft = i_min-1.0d0
        if(dwLeft.LT.thresholdRTemp) then
            do temp_countX=0,99
                x_temp = i_minCoarse-dble(temp_countX)/100.0d0
                do j=j_minStart,j_minEnd
                    do temp_countY=-99,99
                        y_temp = dble(j)+dble(temp_countY)/100.0d0
                        tempEllip1 = ( (x_temp-xCenter(cNum))*dcos(theta(cNum))+(y_temp-yCenter(cNum))*dsin(theta(cNum)) )**2.0d0/radiusA(cNum)**2.0d0 
                        tempEllip2 = (-(x_temp-xCenter(cNum))*dsin(theta(cNum))+(y_temp-yCenter(cNum))*dcos(theta(cNum)) )**2.0d0/radiusB(cNum)**2.0d0
                        if( (tempEllip2+tempEllip2).LE.1.0d0 ) then
                            if(x_temp.LT.i_min) then    
                                i_min = x_temp
                                j_min = y_temp
                            endif
                        endif
                    enddo
                enddo
            enddo
            dwLeft = i_min-1.0d0
        endif

        dwRight = dble(nx)-i_max
        if(dwRight.LT.thresholdRTemp) then
            do temp_countX=0,99
                x_temp = i_maxCoarse+dble(temp_countX)/100.0d0
                do j=j_maxStart,j_maxEnd
                    do temp_countY=-99,99
                        y_temp = dble(j)+dble(temp_countY)/100.0d0
                        tempEllip1 = ( (x_temp-xCenter(cNum))*dcos(theta(cNum))+(y_temp-yCenter(cNum))*dsin(theta(cNum)) )**2.0d0/radiusA(cNum)**2.0d0 
                        tempEllip2 = (-(x_temp-xCenter(cNum))*dsin(theta(cNum))+(y_temp-yCenter(cNum))*dcos(theta(cNum)) )**2.0d0/radiusB(cNum)**2.0d0
                        if( (tempEllip1+tempEllip2).LE.1.0d0 ) then
                            if(x_temp.GT.i_max) then    
                                i_max = x_temp
                                j_max = y_temp
                            endif
                        endif
                    enddo
                enddo
            enddo
            dwRight = dble(nx)-i_max 
        endif
            
        if(MOD(itc,2000).EQ.0) then
            open(unit=01,file='dw.dat',status='unknown',position='append')
            write(01,*) itc, real(dwLeft), real(dwRight)
            close(01)
        endif
            
        if(dwLeft.GE.thresholdR) then
            Fwxij(cNum) = Fwxij(cNum)
            Fwyij(cNum) = Fwyij(cNum)
            Twij(cNum) = Twij(cNum)
        elseif( (dwLeft.LT.thresholdR).AND.(dwLeft.GE.0.0d0) ) then
            Fwxij(cNum) = Fwxij(cNum)+forceScale*((dwLeft-thresholdR)/thresholdR)**2.0d0
            Fwyij(cNum) = Fwyij(cNum)
            Twij(cNum) = Twij(cNum)+forceScale*((dwLeft-thresholdR)/thresholdR)**2.0d0*(yCenter(cNum)-j_min)
            if(MOD(itc,1000).EQ.0) then
                write(*,*) "Repulsive force activated (near left wall)! itc=",itc
                !!call output_Tecplot()
            endif
        else
            write(*,*) "penetration left wall! cNum=",cNum
            write(*,*) "xCenter=",xCenter(cNum)
            write(*,*) "dwLeft=",dwLeft
            stop
        endif
            
        if(dwRight.GE.thresholdR) then
            Fwxij(cNum) = Fwxij(cNum)
            Fwyij(cNum) = Fwyij(cNum)
            Twij(cNum) = Twij(cNum)
        elseif( (dwRight.LT.thresholdR).AND.(dwRight.GE.0.0d0) ) then
            Fwxij(cNum) = Fwxij(cNum)-forceScale*((dwRight-thresholdR)/thresholdR)**2.0d0
            Fwyij(cNum) = Fwyij(cNum)
            Twij(cNum) = Twij(cNum)+(-forceScale*((dwRight-thresholdR)/thresholdR)**2.0d0)*(yCenter(cNum)-j_max)
            if(MOD(itc,1000).EQ.0) then
                write(*,*) "Repulsive force activated (near right wall)! itc=",itc
                !!call output_Tecplot()
            endif
        else
            write(*,*) "penetration right wall!"
            write(*,*) "xCenter=",xCenter(cNum)
            write(*,*) "dwRight=",dwRight
            stop
        endif
 !---------------particle-wall (spring force model)------------------
    enddo
    
!------------------------------------------------------------------------------------------------------------------------------
    
    do cNum=1,cNumMax        
        wallTotalForceX(cNum) = wallTotalForceX(cNum)+Fwxij(cNum)
        wallTotalForceY(cNum) = wallTotalForceY(cNum)-(rhoSolid-rhoAvg)*Pi*radiusA(cNum)*radiusB(cNum)*gravity+Fwyij(cNum)
        totalTorque(cNum) = totalTorque(cNum)+Twij(cNum)
        
        if(MOD(itc,200).EQ.0) then
            open(unit=01,file='springForce.dat',status='unknown',position='append')
            write(01,*) itc*viscosity/diameter0**2.0d0, Fwxij(cNum), wallTotalForceX(cNum)
            close(01)

            open(unit=01,file='springTorque.dat',status='unknown',position='append')
            write(01,*) itc*viscosity/diameter0**2.0d0, Twij(cNum), totalTorque(cNum)
            close(01)
        endif
        
        xCenterOld(cNum) = xCenter(cNum)
        yCenterOld(cNum) = yCenter(cNum)
        thetaOld(1) = theta(cNum)
        UcOld(cNum) = Uc(cNum)
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
    
    errorU = sqrt(error1)/sqrt(error2)

    write(*,*) itc,' ',errorU

    return
    end subroutine check


    subroutine output_Tecplot()
    use commondata
    implicit none
    integer :: i, j, k

    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5
    integer, parameter :: kmax=1
    character(len=40) :: zoneName

    write(B2,'(i9.9)') itc
    open(41,file='particleChannel-'//B2//'.plt', access='stream', form='unformatted')

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
    V3='obst'
    call dumpstring(V3)
    V4='U'
    call dumpstring(V4)
    V5='V'
    call dumpstring(V5)
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
                write(41) real(obst(i,j))
                write(41) real(U(i,j))
                write(41) real(V(i,j))
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
    real(8) :: ax, ay
    real(8) :: aOmega
    
    real(8) :: xp, yp
    real(8) :: temp1, temp2

    do cNum=1,cNumMax
        !--------------------- !--------------------- !--------------------- !--------------------- !---------------------|        
        ax = wallTotalForceX(cNum)/Pi/radiusA(cNum)/radiusB(cNum)/rhoSolid
        ay = wallTotalForceY(cNum)/Pi/radiusA(cNum)/radiusB(cNum)/rhoSolid
        
        aOmega = totalTorque(cNum)/(0.25d0*rhoSolid*Pi*radiusA(cNum)*radiusB(cNum)*(radiusA(cNum)**2.0d0+radiusB(cNum)**2.0d0))
        !--------------------- !--------------------- !--------------------- !--------------------- !---------------------|
        Uc(cNum) = UcOld(cNum)+ax
        Vc(cNum) = VcOld(cNum)+ay  
        
        rationalOmega(cNum) = rationalOmegaOld(cNum)+aOmega
        !--------------------- !--------------------- !--------------------- !--------------------- !---------------------|
        xCenter(cNum) = xCenterOld(cNum)+UcOld(cNum)+0.5d0*ax
        yCenter(cNum) = yCenterOld(cNum)+VcOld(cNum)+0.5d0*ay
        
        theta(cNum) = thetaOld(cNum)+rationalOmega(cNum)+0.5d0*aOmega
        !--------------------- !--------------------- !--------------------- !--------------------- !---------------------|
    enddo

    if(MOD(itc,200).EQ.0) then
        open(unit=01,file='trajectories.dat',status='unknown',position='append')
        write(01,*) (nyHalf-yCenter(1)+shiftTimes)/dble(nx), (xCenter(1)-0.5d0)/dble(nx)
        close(01)

        open(unit=01,file='orientations.dat',status='unknown',position='append')
        write(01,*) (nyHalf-yCenter(1)+shiftTimes)/dble(nx), theta(1)/Pi
        close(01)

        open(unit=01,file='torque.dat',status='unknown',position='append')
        write(01,*) itc, totalTorque(1)/((rhoSolid-1.0d0)*Pi*radiusA(1)*radiusB(1)*gravity*radiusA(1))
        close(01)
        
        open(unit=01,file='Ux-Uy.dat',status='unknown',position='append')
        write(01,*) itc*viscosity/diameter0**2.0d0, Uc(1)*diameter0/viscosity, Vc(1)*diameter0/viscosity
        close(01)
        
        open(unit=01,file='Reynolds.dat',status='unknown',position='append')
        write(01,*) itc*viscosity/diameter0**2.0d0, -Vc(1)*diameter0/viscosity
        close(01)
    endif

    obstNew = 0
    !$omp parallel do default(none) shared(rho,u,v,obstNew,xCenter,yCenter,radiusA,radiusB,theta) private(i,j,cNum)
    do j=1,ny
        do i=1,nx
            do cNum=1,cNumMax

       if( (((i-xCenter(cNum))*dcos(theta(cNum))+(j-yCenter(cNum))*dsin(theta(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(i-xCenter(cNum))*dsin(theta(cNum))+(j-yCenter(cNum))*dcos(theta(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then
        
                    obstNew(i,j) = 1 ! solid node
                    rho(i,j) = rhoSolid
                    u(i,j) = 0.0d0
                    v(i,j) = 0.0d0
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do

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

    !$omp parallel do default(none) &
    !$omp shared(obst,obstNew,xCenter,yCenter,xCenterOld,yCenterOld,thetaOld,Uc,Vc,radiusA,radiusB,theta,f,rho,u,v,ex,ey,rhoAvg,rationalOmega) &
    !$omp private(i,j,myFlag,cNum,outNormal,ec,alpha,tempNormal,m,xp,yp,temp1,temp2) 
    do j=1,ny
        do i=1,nx
            if( (obst(i,j).EQ.1).AND.(obstNew(i,j).EQ.0) ) then ! a solid node turns to new birth fluid node
                myFlag = 0
                do cNum=1,cNumMax
                    
       if( (((i-xCenterOld(cNum))*dcos(thetaOld(cNum))+(j-yCenterOld(cNum))*dsin(thetaOld(cNum)))**2.0d0/radiusA(cNum)**2.0d0 &          
        +(-(i-xCenterOld(cNum))*dsin(thetaOld(cNum))+(j-yCenterOld(cNum))*dcos(thetaOld(cNum)))**2.0d0/radiusB(cNum)**2.0d0).LE.1.0d0 ) then
        
                        myFlag = 1
                        
!--------determine extrapolating direction-----------------------------
                        outNormal = 0.0d0
                        ec = 0
                        do alpha=1,8
                            !!tempNormal = ( ((i-xCenter(cNum))*ex(alpha)+(j-yCenter(cNum))*ey(alpha)) )/ &
                            !!            dsqrt( (i-xCenter(cNum))**2.0d0+(j-yCenter(cNum))**2.0d0 )
                            xp = (i-xCenter(cNum))*dcos(thetaOld(cNum))+(j-yCenter(cNum))*dsin(thetaOld(cNum))
                            yp = -(i-xCenter(cNum))*dsin(thetaOld(cNum))+(j-yCenter(cNum))*dcos(thetaOld(cNum))
                            temp1 = (dcos(thetaOld(cNum))*xp/radiusA(cNum)**2.0d0-dsin(thetaOld(cNum))*yp/radiusB(cNum)**2.0d0)*ex(alpha) &
                                               +(dsin(thetaOld(cNum))*xp/radiusA(cNum)**2.0d0+dcos(thetaOld(cNum))*yp/radiusB(cNum)**2.0d0)*ey(alpha) 
                            temp2 = dsqrt( (dcos(thetaOld(cNum))*xp/radiusA(cNum)**2.0d0-dsin(thetaOld(cNum))*yp/radiusB(cNum)**2.0d0)**2.0d0+ &
                                                         (dsin(thetaOld(cNum))*xp/radiusA(cNum)**2.0d0+dcos(thetaOld(cNum))*yp/radiusB(cNum)**2.0d0)**2.0d0 )
                            tempNormal = temp1 / temp2
                                                
                            if( tempNormal.GT.outNormal ) then
                                outNormal =tempNormal
                                ec = alpha
                            endif
                        enddo
                        if(ec.EQ.0) then
                            write(*,*) "error determine extrapolaing direction"
                            write(*,*) "tempNormal=",tempNormal
                            write(*,*) "temp1=",temp1
                            write(*,*) "temp2=",temp2
                            write(*,*) "xp",xp
                            write(*,*) "yp",yp
                            write(*,*) "i, j",i,j
                            write(*,*) "xCenter,yCenter",xCenter(cNum), yCenter(cNum)
                            write(*,*) "thetaOld,theta",thetaOld(cNum),theta(cNum)
                            stop
                        endif
                        ! write(*,*) "ec==",ec
!--------determine extrapolating direction-----------------------------

!---refilling using extrpolating scheme--------------------------------
            do alpha=0,8
                f(alpha,i,j) = 3.0d0*f(alpha,i+ex(ec),j+ey(ec))-3.0d0*f(alpha,i+2*ex(ec),j+2*ey(ec))+f(alpha,i+3*ex(ec),j+3*ey(ec))
            enddo
!---refilling using extrpolating scheme--------------------------------

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
