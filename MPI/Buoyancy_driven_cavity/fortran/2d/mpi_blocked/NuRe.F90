subroutine calNuRe()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: halfTime
    real(kind=8) :: NuVolAvg_firstHalf, NuVolAvg_lastHalf
    real(kind=8) :: ReVolAvg_firstHalf, ReVolAvg_lastHalf
    real(kind=8) :: NuVolAvg_temp
    real(kind=8) :: ReVolAvg_temp
    real(kind=8) :: angularMomentum

    dimensionlessTime = dimensionlessTime+1
    !-------------------------------------------------------------------------------------------------------
    angularMomentum = 0.0d0
    do j=1,ny 
        do i=1,nx          
                angularMomentum = angularMomentum+(i-nxHalf)*v(i,j)-(j-nyHalf)*u(i,j)
        enddo
    enddo
    angularMomentum = angularMomentum/dble(fluidNumMax)
    
    open(unit=01,file="angularMomentum.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), angularMomentum
    close(01)
    
    NuVolAvg_temp = 0.0d0
    
    do j=1,ny
        do i=1,nx       
                NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)
        enddo
    enddo
    NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(fluidNumMax)*lengthUnit/diffusivity+1.0d0
    
    open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
    write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), NuVolAvg(dimensionlessTime)
    close(01)
    
    NuVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)+NuVolAvg(i) 
    enddo
    NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------
    ReVolAvg_temp = 0.0d0
   
    do j=1,ny
        do i=1,nx       
                ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
        enddo
    enddo
    ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(fluidNumMax))*lengthUnit/viscosity
    
    open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
    write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), ReVolAvg(dimensionlessTime)  !!for print purpose only
    close(02)

    ReVolAvg_mean(dimensionlessTime) = 0.0d0
    do i=1,dimensionlessTime
        ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)+ReVolAvg(i) 
    enddo
    ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
    !-------------------------------------------------------------------------------------------------------    

    if(statisticallyStationaryState.EQ.1) then
        if( ((dimensionlessTime-statisticallyStationaryTime).GE.minAvgPeriod).AND.(MOD(dimensionlessTime-statisticallyStationaryTime,stationaryCheckInterval).EQ.0) ) then
            
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
            
            if( (dabs((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf).LE.statConvError).AND. &
                 (dabs((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf).LE.statConvError) ) then

                if(flowReversal.EQ.0) then
                    statisticallyConvergeState = 1
                elseif(flowReversal.EQ.1) then
                    if(dimensionlessTime.GE.(statisticallyStationaryTime+flowReversalTimeMax)) statisticallyConvergeState = 1
                endif
                
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|---------Statistically Converge Reached at Time =", int(outputFrequency*dimensionlessTime),"----------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif

            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
            write(01,*) "statisticallyStationaryTime=", int(statisticallyStationaryTime*outputFrequency)
            write(01,*) "halfTime =", int(halfTime*outputFrequency)
            write(01,*) "NuVolAvg diff. (first and last half):",dabs(((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (first and last half):",dabs(((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf)*100.0d0),"%"
            write(01,*) "               "
            close(01)        
        endif
    endif

    if(statisticallyStationaryState.EQ.0) then
        if(MOD(dimensionlessTime,stationaryCheckInterval).EQ.0) then
            if( (dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/NuVolAvg_mean(dimensionlessTime)).LE.statStationError).AND. &
                (dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/ReVolAvg_mean(dimensionlessTime)).LE.statStationError) ) then
                statisticallyStationaryState = 1
                statisticallyStationaryTime = dimensionlessTime
                open(unit=01,file="statistically.txt",status='unknown',position="append")
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "|-----------Statistically Stationary Reached at Time =", int(statisticallyStationaryTime*outputFrequency),"------|"
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "calling backupData..."
                call backupData()
                write(01,*) "|-------------------------------------------------------------------------|"
                write(01,*) "   "
                close(01)
            endif 
            open(unit=01,file="statistically.txt",status='unknown',position="append")
            write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
            write(01,*) "NuVolAvg diff. (stationaryCheckInterval):",dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
                                                                            /NuVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "ReVolAvg diff. (stationaryCheckInterval):",dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
                                                                            /ReVolAvg_mean(dimensionlessTime)*100.0d0),"%"
            write(01,*) "               "
            close(01)       
        endif
    endif
    
    return
end subroutine calNuRe



! subroutine calNuRe()
!     use commondata
!     implicit none
!     integer(kind=4) :: i, j
!     integer(kind=4) :: halfTime
!     real(kind=8) :: NuVolAvg_firstHalf, NuVolAvg_lastHalf
!     real(kind=8) :: ReVolAvg_firstHalf, ReVolAvg_lastHalf
!     real(kind=8) :: NuVolAvg_temp
!     real(kind=8) :: ReVolAvg_temp
!     real(kind=8) :: angularMomentum

!     dimensionlessTime = dimensionlessTime+1
!     !-------------------------------------------------------------------------------------------------------
!     angularMomentum = 0.0d0
!     !$omp parallel do default(none) shared(u,v,obst) private(i,j) reduction(+:angularMomentum)
!     do j=1,ny 
!         do i=1,nx 
!             if(obst(i,j).EQ.0) then
!                 angularMomentum = angularMomentum+(i-nxHalf)*v(i,j)-(j-nyHalf)*u(i,j)
!             endif
!         enddo
!     enddo
!     !$omp end parallel do
!     angularMomentum = angularMomentum/dble(fluidNumMax)
    
!     open(unit=01,file="angularMomentum.dat",status='unknown',position='append')
!     write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), angularMomentum
!     close(01)
    
!     NuVolAvg_temp = 0.0d0
!     !$omp parallel do default(none) shared(v,T,obst) private(i,j) reduction(+:NuVolAvg_temp)
!     do j=1,ny
!         do i=1,nx       
!             if(obst(i,j).EQ.0) then
!                 NuVolAvg_temp = NuVolAvg_temp+v(i,j)*T(i,j)
!             endif
!         enddo
!     enddo
!     !$omp end parallel do
!     NuVolAvg(dimensionlessTime) = NuVolAvg_temp/dble(fluidNumMax)*lengthUnit/diffusivity+1.0d0
    
!     open(unit=01,file="Nu_VolAvg.dat",status='unknown',position='append')
!     write(01,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), NuVolAvg(dimensionlessTime)
!     close(01)
    
!     NuVolAvg_mean(dimensionlessTime) = 0.0d0
!     do i=1,dimensionlessTime
!         NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)+NuVolAvg(i) 
!     enddo
!     NuVolAvg_mean(dimensionlessTime) = NuVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
!     !-------------------------------------------------------------------------------------------------------
    
!     !-------------------------------------------------------------------------------------------------------
!     ReVolAvg_temp = 0.0d0
!     !$omp parallel do default(none) shared(u,v,obst) private(i,j) reduction(+:ReVolAvg_temp)
!     do j=1,ny
!         do i=1,nx       
!             if(obst(i,j).EQ.0) then
!                 ReVolAvg_temp = ReVolAvg_temp+(u(i,j)*u(i,j)+v(i,j)*v(i,j))
!             endif
!         enddo
!     enddo
!     !$omp end parallel do
!     ReVolAvg(dimensionlessTime) = dsqrt(ReVolAvg_temp/dble(fluidNumMax))*lengthUnit/viscosity
    
!     open(unit=02,file="Re_VolAvg.dat",status='unknown',position='append')
!     write(02,*) real(reloadDimensionlessTime+dimensionlessTime*outputFrequency), ReVolAvg(dimensionlessTime)  !!for print purpose only
!     close(02)

!     ReVolAvg_mean(dimensionlessTime) = 0.0d0
!     do i=1,dimensionlessTime
!         ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)+ReVolAvg(i) 
!     enddo
!     ReVolAvg_mean(dimensionlessTime) = ReVolAvg_mean(dimensionlessTime)/dble(dimensionlessTime)
!     !-------------------------------------------------------------------------------------------------------    

!     if(statisticallyStationaryState.EQ.1) then
!         if( ((dimensionlessTime-statisticallyStationaryTime).GE.minAvgPeriod).AND.(MOD(dimensionlessTime-statisticallyStationaryTime,stationaryCheckInterval).EQ.0) ) then
            
!             halfTime = (statisticallyStationaryTime+dimensionlessTime)/2
            
!             NuVolAvg_firstHalf = 0.0d0
!             ReVolAvg_firstHalf = 0.0d0
!             do i=statisticallyStationaryTime,halfTime
!                 NuVolAvg_firstHalf = NuVolAvg_firstHalf+NuVolAvg(i)
!                 ReVolAvg_firstHalf = ReVolAvg_firstHalf+ReVolAvg(i)
!             enddo
!             NuVolAvg_firstHalf = NuVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
!             ReVolAvg_firstHalf = ReVolAvg_firstHalf/dble(halfTime-statisticallyStationaryTime+1)
            
!             NuVolAvg_lastHalf = 0.0d0
!             ReVolAvg_lastHalf = 0.0d0
!             do i=halfTime,dimensionlessTime
!                 NuVolAvg_lastHalf = NuVolAvg_lastHalf+NuVolAvg(i)
!                 ReVolAvg_lastHalf = ReVolAvg_lastHalf+ReVolAvg(i)
!             enddo
!             NuVolAvg_lastHalf = NuVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
!             ReVolAvg_lastHalf = ReVolAvg_lastHalf/dble(dimensionlessTime-halfTime+1)
            
!             if( (dabs((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf).LE.statConvError).AND. &
!                  (dabs((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf).LE.statConvError) ) then

!                 if(flowReversal.EQ.0) then
!                     statisticallyConvergeState = 1
!                 elseif(flowReversal.EQ.1) then
!                     if(dimensionlessTime.GE.(statisticallyStationaryTime+flowReversalTimeMax)) statisticallyConvergeState = 1
!                 endif
                
!                 open(unit=01,file="statistically.txt",status='unknown',position="append")
!                 write(01,*) "|-------------------------------------------------------------------------|"
!                 write(01,*) "|---------Statistically Converge Reached at Time =", int(outputFrequency*dimensionlessTime),"----------|"
!                 write(01,*) "|-------------------------------------------------------------------------|"
!                 write(01,*) "   "
!                 close(01)
!             endif

!             open(unit=01,file="statistically.txt",status='unknown',position="append")
!             write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
!             write(01,*) "statisticallyStationaryTime=", int(statisticallyStationaryTime*outputFrequency)
!             write(01,*) "halfTime =", int(halfTime*outputFrequency)
!             write(01,*) "NuVolAvg diff. (first and last half):",dabs(((NuVolAvg_firstHalf-NuVolAvg_lastHalf)/NuVolAvg_lastHalf)*100.0d0),"%"
!             write(01,*) "ReVolAvg diff. (first and last half):",dabs(((ReVolAvg_firstHalf-ReVolAvg_lastHalf)/ReVolAvg_lastHalf)*100.0d0),"%"
!             write(01,*) "               "
!             close(01)        
!         endif
!     endif

!     if(statisticallyStationaryState.EQ.0) then
!         if(MOD(dimensionlessTime,stationaryCheckInterval).EQ.0) then
!             if( (dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/NuVolAvg_mean(dimensionlessTime)).LE.statStationError).AND. &
!                 (dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval))/ReVolAvg_mean(dimensionlessTime)).LE.statStationError) ) then
!                 statisticallyStationaryState = 1
!                 statisticallyStationaryTime = dimensionlessTime
!                 open(unit=01,file="statistically.txt",status='unknown',position="append")
!                 write(01,*) "|-------------------------------------------------------------------------|"
!                 write(01,*) "|-----------Statistically Stationary Reached at Time =", int(statisticallyStationaryTime*outputFrequency),"------|"
!                 write(01,*) "|-------------------------------------------------------------------------|"
!                 write(01,*) "calling backupData..."
!                 ! call backupData()
!                 write(01,*) "|-------------------------------------------------------------------------|"
!                 write(01,*) "   "
!                 close(01)
!             endif 
!             open(unit=01,file="statistically.txt",status='unknown',position="append")
!             write(01,*) "dimensionlessTime =", int(dimensionlessTime*outputFrequency)
!             write(01,*) "NuVolAvg diff. (stationaryCheckInterval):",dabs((NuVolAvg_mean(dimensionlessTime)-NuVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
!                                                                             /NuVolAvg_mean(dimensionlessTime)*100.0d0),"%"
!             write(01,*) "ReVolAvg diff. (stationaryCheckInterval):",dabs((ReVolAvg_mean(dimensionlessTime)-ReVolAvg_mean(dimensionlessTime-stationaryCheckInterval)) &
!                                                                             /ReVolAvg_mean(dimensionlessTime)*100.0d0),"%"
!             write(01,*) "               "
!             close(01)       
!         endif
!     endif
    
!     return
! end subroutine calNuRe
