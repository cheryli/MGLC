subroutine initial()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: un(0:18), unT(0:6)
    real(kind=8) :: us2
    real(kind=8) :: omega(0:18), omegaT(0:6)
    
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0
    
    if (rank == 0) then
        write(*,*) 'Mesh:',nx,ny,nz
        write(*,*) 'Rayleigh=',real(Rayleigh), ', Prandtl=',real(Prandtl), ', Mach=',real(Mach)
        write(*,*) "Ekman=",real(Ekman)
        write(*,*) "   "
        write(*,*) 'tauf=',real(tauf)
        write(*,*) "paraA=",real(paraA)
        write(*,*) "viscosity=",real(viscosity), ", diffusivity=",real(diffusivity)
        write(*,*) "omegaRatating=", real(omegaRatating)
        write(*,*) "itc_max=",itc_max
        write(*,*) "Ouput will begin at", int(timeUnit)
        write(*,*) "Output interval is", int(timeUnit)
        write(*,*) "Time unit: Sqrt(L0/(gBeta*DeltaT))=",dsqrt(dble(Ny)/gBeta)
        write(*,*) "Velocity unit: Sqrt(gBeta*L0*DeltaT)=",dsqrt(gBeta*dble(Ny))
        write(*,*) "    "
    
#ifdef RBconvection
    write(*,*) "    "
    write(*,*) "I am RBconvection"
    if(flowReversal.EQ.1) then
        write(*,*) "I am also flow reversal"
    write(*,*) "    "
    elseif(flowReversal.EQ.0) then
        write(*,*) "I am NOT flow reversal"
    write(*,*) "    "
    else
        write(*,*) "Error: Please check flowReversal setting!"
        stop
    endif
#endif

#ifdef benchmarkCavity
    write(*,*) "I am benchmarkCavity"
    write(*,*) "    "
#endif
    endif

    xp(0) = 0.0d0
    xp(total_nx+1) = dble(total_nx)
    do i = 1, total_nx
        xp(i) = dble(i)-0.5d0
    enddo
    yp(0) = 0.0d0
    yp(total_ny+1) = dble(total_ny)
    do j = 1, total_ny
        yp(j) = dble(j)-0.5d0
    enddo
    zp(0) = 0.0d0
    zp(total_nz+1) = dble(nz)
    do k = 1, total_nz
        zp(k) = dble(k)-0.5d0
    enddo
    
    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (T(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    allocate (up(nx,ny,nz))
    allocate (vp(nx,ny,nz))
    allocate (wp(nx,ny,nz))
    allocate (Tp(nx,ny,nz))
    
    allocate (f(0:18,nx,ny,nz))
    allocate (f_post(0:18,0:nx+1,0:ny+1,0:nz+1))
    allocate (g(0:6,nx,ny,nz))
    allocate (g_post(0:6,0:nx+1,0:ny+1,0:nz+1))
    
    allocate (Fx(nx,ny,nz))
    allocate (Fy(nx,ny,nz))
    allocate (Fz(nx,ny,nz))
    
    rho = 1.0d0
    
    omega(0) = 1.0d0/3.0d0
    do alpha=1,6
        omega(alpha) = 1.0d0/18.0d0
    enddo
    do alpha=7,18
        omega(alpha) = 1.0d0/36.0d0
    enddo
    
    omegaT(0) = (1.0d0-paraA)/7.0d0
    do alpha=1,6
        omegaT(alpha) = (paraA+6.0d0)/42.0d0
    enddo

    
    if(loadInitField.EQ.0) then 
        u = 0.0d0
        v = 0.0d0
        w = 0.0d0
        T = 0.0d0
        
    if (rank == 0) then
        write(*,*) "Initial field is set exactly"
#ifdef noslipWalls
    write(*,*) "Velocity B.C. for left/right walls are: ===No-slip wall==="
#endif
#ifdef TopBottomPlatesConstT
    write(*,*) "Temperature B.C. for bottom plate is:==Hot bottom plate==="
    write(*,*) "Temperature B.C. for top plate is:====Cold top plate==="
#endif
#ifdef LeftRightWallsConstT
    write(*,*) "Temperature B.C. for left wall is:===Hot left wall==="
    write(*,*) "Temperature B.C. for right wall is:==Cold right wall==="
#endif
#ifdef LeftRightWallsAdiabatic
    write(*,*) "Temperature B.C. for left/right walls are:===Adiabatic walls==="
#endif
#ifdef TopBottomPlatesAdiabatic
    write(*,*) "Temperature B.C. for top/bottom plates are:===Adiabatic plates==="
#endif

#ifdef BackFrontWallsAdiabatic
    write(*,*) "Temperature B.C. for back/front walls are:===Adiabatic walls==="
#endif    
    endif


#ifdef LeftRightWallsConstT
    ! left wall y=1
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                T(i,1,k) = Thot
            enddo
        enddo
    endif

    ! right wall y = ny
    if (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                T(i,ny,k) = Tcold
            enddo
        enddo
    endif

#endif


#ifdef TopBottomPlatesConstT
    ! bottom wall k = 1
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                T(i,j,1) = Thot
            enddo
        enddo
    endif


    ! top wall k = nz
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                T(i,j,nz) = Tcold
            enddo
        enddo
    endif
#endif

        !$omp parallel do default(none) shared(f,g,u,v,w,T,ex,ey,ez,omega,omegaT,rho) private(i,j,k,alpha,us2,un,unT)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                    do alpha=0,18
                        un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        f(alpha,i,j,k) = rho(i,j,k)*omega(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                    enddo
                    do alpha=0,6
                        unT(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                        g(alpha,i,j,k) = omegaT(alpha)*T(i,j,k)*(1.0d0+21.0d0/(6.0d0+paraA)*unT(alpha))
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

    ! elseif(loadInitField.EQ.1) then
    !     write(*,*) "Load initial field from previous simulation"
    !     write(*,*) "Start to read raw data >>>>>>>>>>>>"
    !     open(unit=01,file='./reload/backupFile-1000.bin',form="unformatted",access="sequential",status='old')
    !     read(01) (((u(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) ((((f(alpha,i,j,k),alpha=0,18),i=1,nx),j=1,ny),k=1,nz)
    !     read(01) ((((g(alpha,i,j,k),alpha=0,6),i=1,nx),j=1,ny),k=1,nz)
    !     close(01)
    !     write(*,*) "Read raw data!"
    else
        write(*,*) "Error: initial field is not properly set"
    endif
    
    up = 0.0d0
    vp = 0.0d0
    wp = 0.0d0
    Tp = 0.0d0
    
    f_post = 0.0d0
    g_post = 0.0d0
        
    fileNum = 0
    dimensionlessTime = 0
    
    NuVolAvg = 0.0d0
    ReVolAvg = 0.0d0
    NuVolAvg_mean = 0.0d0
    ReVolAvg_mean = 0.0d0
    statisticallyStationaryState = 0
    statisticallyConvergeState = 0

    return
end subroutine initial