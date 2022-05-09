subroutine initial()
    use commondata
    use ioFolder
    implicit none
    integer(kind=4) :: i, j, start, end
    integer(kind=4) :: alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: us2
    character(len=100) :: reloadFileName
    real(kind=8) :: dev
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    fluidNumMax = total_nx*total_ny
    
if (rank == 0) then
    open(unit=00,file="SimulationSettings.txt",status='unknown',position='append')

    if(outputBinFile.EQ.1) then
        open(unit=01,file=trim(binFolderPrefix)//"-"//"readme",status="unknown")
        write(01,*) "binFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", binFolderPrefix
    endif 
    if(outputPltFile.EQ.1) then
        open(unit=01,file=trim(pltFolderPrefix)//"-"//"readme",status="unknown")
        write(01,*) "pltFile folder exist!"
        close(01)
        write(00,*) "Data will be stored in ", pltFolderPrefix
    endif
    
    if( (paraA.GE.1.0d0).OR.(paraA.LE.-4.0d0) ) then
        write(00,*) '----------------------------------'
        write(00,*) 'Error: paraA=', paraA
        write(00,*) '... Current Mach number is :', real(Mach)
        write(00,*) '... Please try to reduce Mach number!'
        write(00,*) '----------------------------------'
        stop
    endif

    if(flowReversal.EQ.1) then
        if(dimensionlessTimeMax.LE.flowReversalTimeMax) then
            write(00,*) "Error: in case of flor reversal,"
            write(00,*) "... please check dimensionlessTimeMax!"
            stop
        endif
    endif

    write(00,*)"-------------------------------------------------------------------------------"
    write(00,*) 'Mesh:',total_nx, total_ny
    write(00,*) 'Rayleigh=',real(Rayleigh), '; Prandtl =',real(Prandtl), '; Mach =',real(Mach)
    write(00,*) 'shearReynolds=', real(shearReynolds)
    ! if(shearReynolds .GT. 0.0d0) then
    !     write(00,*) "Richardson=", real(Rayleigh/(shearReynolds**2.0d0*Prandtl))
    !     write(00,*) "...UwallTopLeft=", real(UwallTopLeft), "; UwallTopRight=", real(UwallTopRight)
    !     write(00,*) "...UwallBottomLeft=", real(UwallBottomLeft), "; UwallBottomRight=", real(UwallBottomRight)
    !     write(00,*) "...UwallLeftTop=", real(UwallLeftTop), "; UwallLeftBottom=", real(UwallLeftBottom)
    !     write(00,*) "...UwallRightTop=", real(UwallRightTop), "; UwallRightBottom=", real(UwallRightBottom)
    ! endif
    write(00,*) "Length unit: L0 =", real(lengthUnit)
    write(00,*) "Time unit: Sqrt(L0/(gBeta*DeltaT)) =", real(timeUnit)
    write(00,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT) =", real(velocityUnit)
    write(00,*) "   "
    write(00,*) 'tauf=',real(tauf)
    write(00,*) "viscosity =",real(viscosity), "; diffusivity =",real(diffusivity)
    write(00,*) "outputFrequency", real(outputFrequency), "tf"
    write(00,*) "......................  or ",  int(outputFrequency*timeUnit), "in itc units"
    write(00,*) "dimensionlessTimeMax:"
    write(00,*) "...statistically stationary PLUS convergence = ", real(dimensionlessTimeMax*outputFrequency)
    if(flowReversal.EQ.1) then
        write(00,*) "flowReversalTimeMax:"
        write(00,*) "...statistically stationary PLUS convergence = ", real(flowReversalTimeMax*outputFrequency)
    endif
    write(00,*) "minAvgPeriod:"
    write(00,*) "...to obtain statistically convergence =", real(minAvgPeriod*outputFrequency)
    write(00,*) "backupInterval =", backupInterval, " free-fall time units"
    write(00,*) ".................... or ", int(backupInterval/outputFrequency)*int(outputFrequency*timeUnit), "itc units"
    if(loadInitField.EQ.1) then
        write(00,*) "reloadDimensionlessTime=", reloadDimensionlessTime
        write(00,*) "reloadStatisticallyConvergeTime=", reloadStatisticallyConvergeTime
    endif
    write(00,*) "itc_max =",itc_max
    write(00,*) "default epsU =", real(epsU),"; epsT =", real(epsT)
    write(00,*) "    "
    

    write(00,*) "I am regular geometry (Square Cavity)"
    if(flowGeometry.NE.squareGeo) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "... flowGeometry should be", squareGeo, "; while its actual value is ", flowGeometry
        stop
    endif
    if(dabs(lengthUnit-dble(total_ny)).GT.1e-6) then
        write(00,*) "Error: Please check the geometry setting of the cell!"
        write(00,*) "... lengthUnit should be equal to ny"
        stop
    endif

#ifdef RayleighBenardCell
    write(00,*) "I am Rayleigh Benard Cell"
#endif
#ifdef SideHeatedCell
    write(00,*) "I am Side Heated Cell"
#endif

#ifdef steadyFlow
        write(00,*) "I am steadyFlow"
#endif
#ifdef unsteadyFlow
        write(00,*) "I am unsteadyFlow"
#endif
    
    if(flowReversal.EQ.1) then
        write(00,*) "I am also flow reversal"
    elseif(flowReversal.EQ.0) then
        write(00,*) "I am NOT flow reversal"
    else
        write(00,*) "Error: Please check flowReversal setting!"
        stop
    endif

    write(00,*) "I am level cell"


    write(00,*) "Initial field is set exactly"
    if(reloadDimensionlessTime.NE.0) then
        write(00,*) "Error: since loadInitField.EQ.0, reloadDimensionlessTime should also be 0"
        stop
    endif
    
#ifdef VerticalWallsNoslip
    write(00,*) "Velocity B.C. for vertical walls are: ===No-slip wall==="
#endif
#ifdef VerticalWallsFreeslip
    write(00,*) "Velocity B.C. for vertical walls are: ===Free-slip wall==="
#endif
#ifdef VerticalWallsPeriodicalU
    write(00,*) "Velocity B.C. for vertical walls are: ===Periodical==="
#endif
#ifdef HorizontalWallsNoslip
    write(00,*) "Velocity B.C. for horizontal walls are: ===No-slip wall==="
#endif
#ifdef HorizontalWallsFreeslip
    write(00,*) "Velocity B.C. for horizontal walls are: ===Free-slip wall==="
#endif
#ifdef VerticalWallsConstT
    write(00,*) "Temperature B.C. for vertical walls are:===Hot/cold wall==="
#endif
#ifdef HorizontalWallsConstT
    write(00,*) "Temperature B.C. for horizontal walls are:===Hot/cold wall==="
#endif
#ifdef VerticalWallsPeriodicalT
write(00,*) "Temperature B.C. for vertical walls are:===Periodical wall==="
#endif

#ifdef VerticalWallsAdiabatic
write(00,*) "Temperature B.C. for vertical walls are:===Adiabatic wall==="
#endif

#ifdef HorizontalWallsAdiabatic
write(00,*) "Temperature B.C. for horizontal walls are:===Adiabatic wall==="
#endif

close(00)

endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i=1,total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j=1,total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (T(nx,ny))
    allocate (rho(nx,ny))
    
#ifdef steadyFlow
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (Tp(nx,ny))
#endif
    
    allocate (f(0:8,nx,ny))
    allocate (f_post(0:8,0:nx+1,0:ny+1))
    allocate (g(0:4,nx,ny))
    allocate (g_post(0:4,0:nx+1,0:ny+1))
    
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))
    
    rho = rho0
    
    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/5.0d0
    do alpha=1,4
        omegaT(alpha) = (paraA+4.0d0)/20.0d0
    enddo
    

    if(loadInitField.EQ.0) then 
    
        u = 0.0d0
        v = 0.0d0
        T = 0.0d0

        ! do i=1,nxHalf
        !     u(i,ny) = UwallTopLeft
        !     u(i,1)  = UwallBottomLeft
        ! enddo

        ! if (i_start_global + nx  <= nxHalf) then
        !     end = nx
        ! elseif (i_start_global >= nxHalf) then
        !     end = 0
        ! else
        !     end = nxHalf - i_start_global
        ! endif
        ! do i=1,end
        !     u(i,ny) = UwallTopLeft
        !     u(i,1)  = UwallBottomLeft
        ! enddo

        ! do i=nxHalf+1, nx
        !     u(i,ny) = UwallTopRight
        !     u(i,1) = UwallBottomRight
        ! enddo
        ! do j=2,nyHalf
        !     v(1,j) = UwallLeftBottom
        !     v(nx,j) = UwallRightBottom
        ! enddo
        ! do j=nyHalf+1, ny-1
        !     v(1,j) = UwallLeftTop
        !     v(nx,j) = UwallRightTop
        ! enddo


#ifdef VerticalWallsConstT
    do j=1,ny
        !~T(1,j) = Thot
        !~T(nx,j) = Tcold
        do i=1,nx
            ! T(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
            T(i, j) = dble(i_start_global + i - 1) / dble(total_nx - 1) * (Tcold - Thot) + Thot
        enddo
    enddo
#endif
#ifdef HorizontalWallsConstT
    do i=1,nx
        !~T(i,1) = Thot
        !~T(i,ny) = Tcold
        do j=1,ny
            ! T(i,j) = dble(j-1)/dble(ny-1)*(Tcold-Thot)+Thot
            T(i, j) = dble(j_start_global + j - 1) / dble(total_ny - 1) * (Tcold - Thot) + Thot
        enddo
    enddo
#endif

    f = 0.0d0
    g = 0.0d0
    !$omp parallel do default(none) shared(f,g,u,v,T,ex,ey,omega,omegaT,rho) private(i,j,alpha,us2,un)
    do j=1,ny
        do i=1,nx
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha,i,j) = rho(i,j)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
            enddo
            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(alpha,i,j) = T(i,j)*omegaT(alpha)*(1.0d0+10.0d0/(4.0d0+paraA)*un(alpha))
            enddo
        enddo
    enddo
    !$omp end parallel do

    ! elseif(loadInitField.EQ.1) then
    !     if(reloadDimensionlessTime.EQ.0) then
    !         write(00,*) "Error: since loadInitField.EQ.1, reloadDimensionlessTime should not be 0"
    !         stop
    !     endif
    !     write(00,*) "Load initial field from previous simulation >>>"
    !     write(reloadFileName, *) reloadbinFileNum
    !     reloadFileName = adjustl(reloadFileName)
    !     open(unit=01,file="../reloadFile/backupFile-"//trim(reloadFileName)//".bin",form="unformatted",access="sequential",status='old')
    !     read(01) (((f(alpha,i,j), alpha=0,8), i=1,nx), j=1,ny)
    !     read(01) (((g(alpha,i,j), alpha=0,4), i=1,nx), j=1,ny)
    !     read(01) ((u(i,j),i=1,nx), j=1,ny)
    !     read(01) ((v(i,j),i=1,nx), j=1,ny)
    !     reaD(01) ((T(i,j),i=1,nx), j=1,ny)
    !     close(01)
        
    !     dev = 0.0d0
    !     do j=1,ny
    !         do i=1,nx
    !             rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    !             dev = dev+g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    !         enddo
    !     enddo
    !     write(00,*) "RELOAD: Deviation in temperature: ", real(dev)
    !     if(dev.GT.1.0d0) then
    !         write(00,*) "Error: too large Deviation when reload data!"
    !         stop
    !     endif
    !     write(00,*) "Raw data is loaded from the file: backupFile-",trim(reloadFileName),".bin"
    else
        if (rank == 0) then
            write(00,*) "Error: initial field is not properly set"
        endif
        stop
    endif
    

#ifdef steadyFlow
    up = 0.0d0
    vp = 0.0d0
    Tp = 0.0d0
#endif
        
    f_post = 0.0d0
    g_post = 0.0d0
    
    binFileNum = 0
    pltFileNum = 0
    particleFileNum = 0
    dimensionlessTime = 0

    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 1
    statisticallyStationaryTime = 0
#ifdef unsteadyFlow
    if(loadInitField.EQ.1) statisticallyStationaryState = 1
#endif
    statisticallyConvergeState = 0
    
    return
end subroutine initial