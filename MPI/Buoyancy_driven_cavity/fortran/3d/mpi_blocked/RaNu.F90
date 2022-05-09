subroutine calNuRe()
    use commondata
    implicit none
    integer :: i, j, k 
    integer :: halfTime
    real(kind=8) :: NuVolAvg_firstHalf, NuVolAvg_lastHalf
    real(kind=8) :: ReVolAvg_firstHalf, ReVolAvg_lastHalf
    real(kind=8) :: NuVolAvg_temp
    real(kind=8) :: ReVolAvg_temp
    
    dimensionlessTime = dimensionlessTime+1
    !-------------------------------------------------------------------------------------------------------
    NuVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(w,T) private(i,j,k) reduction(+:NuVolAvg_temp)
    do k=1,nz
        do j=1,ny
            do i=1,nx       
                NuVolAvg_temp = NuVolAvg_temp+w(i,j,k)*T(i,j,k)
            enddo
        enddo
    enddo
    !$omp end parallel do
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(nx*ny*nz)*dble(nz)/diffusivity+1.0d0
    
    open(unit=01,file="Nu_Bulk.dat",status='unknown',position='append')
    write(01,*) itc/timeUnit, NuVolAvg(dimensionlessTime)
    close(01)
    
    NuVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)+NuVolAvg(i) 
    enddo
    NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------
    ReVolAvg_temp = 0.0d0
    !$omp parallel do default(none) shared(u,v,w) private(i,j,k) reduction(+:ReVolAvg_temp)
    do k=1,nz
        do j=1,ny
            do i=1,nx       
                ReVolAvg_temp = ReVolAvg_temp+(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))
            enddo
        enddo
    enddo
    !$omp end parallel do
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(nx*ny*nz))*dble(nz)/viscosity
    
    open(unit=02,file="./Re_Bulk.dat",status='unknown',position='append')
    write(02,*) itc/timeUnit, ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)

    ReVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)+ReVolAvg(i) 
    enddo
    ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------    

    if(statisticallyStationaryState.EQ.1) then
        if(MOD(dimensionlessTime-statisticallyStationaryTime,200).EQ.0) then
            
            halfTime = (statisticallyStationaryTime+dimensionlessTime)/2
            
            NuVolAvg_firstHalf = 0.0d0
            ReVolAvg_firstHalf = 0.0d0
            do i=statisticallyStationaryTime,halfTime
                NuVolAvg_firstHalf = NuVolAvg_firstHalf+NuVolAvg(i)
                ReVolAvg_firstHalf = ReVolAvg_firstHalf+ReVolAvg(i)
            enddo
            NuVolAvg_firstHalf = NuVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
            ReVolAvg_firstHalf = ReVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
            
            NuVolAvg_lastHalf = 0.0d0
            ReVolAvg_lastHalf = 0.0d0
            do i=halfTime,dimensionlessTime
                NuVolAvg_lastHalf = NuVolAvg_lastHalf+NuVolAvg(i)
                ReVolAvg_lastHalf = ReVolAvg_lastHalf+ReVolAvg(i)
            enddo
            NuVolAvg_lastHalf = NuVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
            ReVolAvg_lastHalf = ReVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
            
            if( (dabs((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf).LE.0.01d0).AND. &
                 (dabs((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf).LE.0.01d0) ) then

                if(flowReversal.EQ.0) then
                    statisticallyConvergeState = 1
                elseif(flowReversal.EQ.1) then
                    if(dimensionlessTime.GE.(statisticallyStationaryTime+flowReversalTime)) statisticallyConvergeState = 1
                endif
                
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|---------Statistically Converge Reached at Time =",dimensionlessTime,"----------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif
            
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "statisticallyStationaryTime=",statisticallyStationaryTime
            write(01,*) "halfTime =",halfTime
            write(01,*) "dimensionlessTime =",dimensionlessTime
            write(01,*) "NuVolAvg diff. (first and last half):",dabs(((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (first and last half):",dabs(((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "               "
            close(01)        
        endif
    endif

    if(statisticallyStationaryState.EQ.0) then
        if(MOD(dimensionlessTime,100).EQ.0) then
            if( (dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-100))/NuVolAvg_mean(dimensionlessTime)).LE.0.01d0).AND. &
                (dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-100))/ReVolAvg_mean(dimensionlessTime)).LE.0.01d0) ) then
                statisticallyStationaryState = 1
                statisticallyStationaryTime = dimensionlessTime
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|-----------Statistically Stationary Reached at Time =",statisticallyStationaryTime,"------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                write(01,*) "   "
                close(01)
            endif 
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =",dimensionlessTime
            write(01,*) "NuVolAvg diff. (100 timeunit interval):",dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-100)) &
                                                                            /NuVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (100 timeunit interval):",dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-100)) &
                                                                            /ReVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "               "
            close(01)       
        endif
    endif
    
    return
end subroutine calNuRe