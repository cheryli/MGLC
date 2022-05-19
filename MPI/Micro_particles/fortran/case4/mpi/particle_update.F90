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
    real(8) :: total_rho
    integer :: total_fluidNum
    real(8) :: ax, ay
    real(8) :: aOmega
    character(len=100) :: filename
    real(8) :: distance, distanceX, distanceY

    total_rho = 0.0d0
    total_fluidNum = 0

    do cNum=1,cNumMax
        
        ax = wallTotalForceX(cNum)/Pi/radius0**2.0d0/rhoSolid
        ay = wallTotalForceY(cNum)/Pi/radius0**2.0d0/rhoSolid
        
        aOmega = totalTorque(cNum)/(0.5d0*rhoSolid*Pi*radius(cNum)**4.0d0)
        
        !--------------------- !--------------------- !--------------------- !--------------------- !---------------------|
        
        Uc(cNum) = UcOld(cNum)+ax
        Vc(cNum) = VcOld(cNum)+ay  
        
        rationalOmega(cNum) = rationalOmegaOld(cNum)+aOmega
        
        xCenter(cNum) = xCenterOld(cNum)+UcOld(cNum)+0.5d0*ax
        yCenter(cNum) = yCenterOld(cNum)+VcOld(cNum)+0.5d0*ay
        
    enddo

    if (rank == 0) then
        do cNum=1,cNumMax
            write(filename,*) cNum
            filename = adjustl(filename)

            open(unit=01,file='x-'//trim(filename)//'.dat',status='unknown',position='append')
                write(01,*) itc*t0, (xCenter(cNum)-0.0d0)*l0
            close(01)

            open(unit=01,file='y-'//trim(filename)//'.dat',status='unknown',position='append')
                write(01,*) itc*t0, (yCenter(cNum)-0.0d0)*l0
            close(01)

            distanceX = xCenter(1)-xCenter(2)
            distanceY = yCenter(1)-yCenter(2)
            distance = dsqrt(distanceX**2.0d0+distanceY**2.0d0)-(radius(1)+radius(2))
            open(unit=01,file='distance.dat',status='unknown',position='append')
                write(01,100) itc*t0, distance*l0, distanceX*l0, distanceY*l0
            close(01)
100 format(2x, f13.6, f13.6, f13.6, f13.6)
        enddo
    endif


    !$omp parallel do default(none) shared(rho,u,v,obstNew,xCenter,yCenter,radius) private(i,j,cNum)
    do j=0,ny+1
        do i=0,nx+1
            obstNew(i,j) = 0 !re-initial
            do cNum=1,cNumMax
                if( ((i + i_start_global - xCenter(cNum))**2.0d0 &
                +(j + j_start_global - yCenter(cNum))**2.0d0) &
                .LE.radius(cNum)**2.0d0 ) then
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
    !$omp shared(obst,obstNew,xCenter,yCenter,xCenterOld,yCenterOld,Uc,Vc,radius,f,rho,u,v,ex,ey,rhoAvg,rationalOmega) &
    !$omp private(i,j,myFlag,cNum,outNormal,ec,alpha,tempNormal,m) 
    do j=1,ny
        do i=1,nx
            if( (obst(i,j).EQ.1).AND.(obstNew(i,j).EQ.0) ) then ! a solid node turns to new birth fluid node
                myFlag = 0
                do cNum=1,cNumMax
                    if( ((i + i_start_global - xCenterOld(cNum))**2.0d0 &
                    +(j + j_start_global - yCenterOld(cNum))**2.0d0) &
                    .LE.radius(cNum)**2.0d0 ) then
                        myFlag = 1
                        
!--------determine extrapolating direction-----------------------------
                        outNormal = 0.0d0
                        ec = 0
                        do alpha=1,8
                            tempNormal = ( ((i + i_start_global - xCenter(cNum))*ex(alpha) &
                            +(j + j_start_global - yCenter(cNum))*ey(alpha)) )/ &
                            dsqrt( (i + i_start_global - xCenter(cNum))**2.0d0 &
                            +(j + j_start_global - yCenter(cNum))**2.0d0 )
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
            do alpha=0,8
                f(alpha,i,j) = 3.0d0*f(alpha,i+ex(ec),j+ey(ec))-3.0d0*f(alpha,i+2*ex(ec),j+2*ey(ec))+f(alpha,i+3*ex(ec),j+3*ey(ec))
            enddo
!---refilling using extrpolating scheme--------------------------------

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*f(5,i,j)+2.0d0*f(6,i,j)+2.0d0*f(7,i,j)+2.0d0*f(8,i,j)
    m(2) = 4.0d0*f(0,i,j)-2.0d0*f(1,i,j)-2.0d0*f(2,i,j)-2.0d0*f(3,i,j)-2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = rhoAvg*( Uc(cNum)-(j + j_start_global - yCenter(cNum))*rationalOmega(cNum) ) !---------------------------------------------------
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = rhoAvg*( Vc(cNum)+(i + i_start_global - xCenter(cNum))*rationalOmega(cNum) ) !--------------------------------------------------
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