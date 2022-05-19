subroutine calForce()
    use commondata
    implicit none
    integer :: alpha
    integer :: i, j
    integer :: ip, jp
    real(8) :: temp1, temp2
    real(8) :: q
    real(8) :: x0, y0
    integer :: myFlag
    integer :: cNum
    real(8) :: tempForceX, tempForceY, tempTorque
    integer :: cNum2
    real(8) :: dij
    real(8) :: Fxij(cNumMax), Fyij(cNumMax)
    real(8) :: dw
    real(8) :: Fwxij(cNumMax), Fwyij(cNumMax)
    real(8) :: forceScale
    character(len=100) :: filename

    do cNum=1,cNumMax
        wallTotalForceX(cNum) = 0.0d0
        wallTotalForceY(cNum) = 0.0d0
        totalTorque(cNum) = 0.0d0
    enddo
                            
    !$omp parallel do default(none) &
    !$omp shared(ex,ey,r,obst,xCenter,yCenter,rationalOmega,radius,f,f_post,Uc,Vc,omega) &
    !$omp private(i,j,alpha,ip,jp,myFlag,cNum,x0,y0,q,temp1,temp2,tempForceX,tempForceY,tempTorque) &
    !$omp reduction(+:wallTotalForceX,wallTotalForceY,totalTorque)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
            
                do alpha=1,8 !neighboring nodes, i.e., those are fluid nodes which is near solid wall
                    ip = i+ex(alpha)
                    jp = j+ey(alpha)
                        
                    if( obst(ip,jp).EQ.1 ) then

                        myFlag = 0
                    
                        do cNum=1,cNumMax
                            if( ((ip-xCenter(cNum))**2.0d0+(jp-yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i),dble(j),alpha,x0,y0,q) 

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
    
!----------------------------------------------------------------------------------------------------------------------------------
    
    do cNum=1,cNumMax
        
        !---------------particle-particle------------------
        forceScale = Pi*radius(cNum)**2.0d0*(rhoSolid-rhoAvg)*gravity/stiffParticle
        Fxij(cNum) = 0.0d0
        Fyij(cNum) = 0.0d0
        do cNum2=1,cNumMax
            if(cNum2.NE.cNum) then
                dij = dsqrt( (xCenter(cNum)-xCenter(cNum2))**2.0d0+(yCenter(cNum)-yCenter(cNum2))**2.0d0 )
                if(dij.GE.(radius(cNum)+radius(cNum2)+thresholdParticle)) then
                    Fxij(cNum) = Fxij(cNum)
                    Fyij(cNum) = Fyij(cNum)
                elseif( (dij.LT.(radius(cNum)+radius(cNum2)+thresholdParticle)).AND.(dij.GE.(radius(cNum)+radius(cNum2))) ) then
Fxij(cNum) = Fxij(cNum)+forceScale*((dij-radius(cNum)-radius(cNum2)-thresholdParticle)/thresholdParticle)**2.0d0 &
                                  *(xCenter(cNum)-xCenter(cNum2))/dij
Fyij(cNum) = Fyij(cNum)+forceScale*((dij-radius(cNum)-radius(cNum2)-thresholdParticle)/thresholdParticle)**2.0d0 &
                                  *(yCenter(cNum)-yCenter(cNum2))/dij
                    if(MOD(itc,1000).EQ.0) then
                        write(*,*) "Particle-particle interaction actived!"
                        write(*,*) "itc=",itc," ,cNum=",cNum
                    endif
                else
                    write(*,*) 'Particle-particle interpenetration!'
                    write(*,*) "itc=",itc," ,cNum=",cNum, cNum2
                    call output_Tecplot()
                    stop
                endif
            endif
        enddo
        !---------------particle-particle------------------

        !---------------particle-wall (spring force model)------------------
        Fwxij(cNum) = 0.0d0
        Fwyij(cNum) = 0.0d0
        forceScale = Pi*radius(cNum)**2.0d0*(rhoSolid-rho0)*gravity/stiffWall

        ! near the bottom wall
        dw = yCenter(cNum)-radius(cNum)-1.0d0
        if(dw .LT. 0) then
            write(*,*) "WARNING: particle penetrates the bottom wall!"
            write(*,*) "yCenter(cNum)=",yCenter(cNum)
            write(*,*) "Error: simulation stoped!"
            call output_Tecplot()
            stop
        elseif(dw .LT. thresholdWall) then
            Fwyij(cNum) = Fwyij(cNum) + forceScale*((dw-thresholdWall)/thresholdWall)**2.0d0
            if(MOD(itc,1000).EQ.0) then
                    write(*,*) "Repulsive force activated (near bottom wall)!!!!!! itc=",itc
            endif
        endif

        ! near the left wall
        dw = xCenter(cNum)-radius(cNum)-1.0d0
        if(dw .LT. 0) then
            write(*,*) "WARNING: particle penetrates the left wall!"
            write(*,*) "xCenter(cNum)=",xCenter(cNum)
            write(*,*) "Error: simulation stoped!"
            call output_Tecplot()
            stop
        elseif( dw .LT. thresholdWall ) then
            Fwxij(cNum) = Fwxij(cNum) + forceScale*((dw-thresholdWall)/thresholdWall)**2.0d0
            if(MOD(itc,1000).EQ.0) then
                    write(*,*) "Repulsive force activated (near left wall)!!!!!! itc=",itc
            endif
        endif

        ! near the right wall
        dw = dble(nx) - xCenter(cNum) - radius(cNum)
        if(dw .LT. 0) then  
            write(*,*) "WARNING: particle penetrates the right wall!"
            write(*,*) "xCenter(cNum)=",xCenter(cNum)
            write(*,*) "Error: simulation stoped!"
            call output_Tecplot()
            stop
        elseif( dw .LT. thresholdWall ) then
            Fwxij(cNum) = Fwxij(cNum) - forceScale*((dw-thresholdWall)/thresholdWall)**2.0d0
            if(MOD(itc,1000).EQ.0) then
                write(*,*) "Repulsive force activated (near right wall)!!!!!! itc=",itc
            endif
        endif
        !---------------particle-wall (spring force model)------------------

    enddo
!------------------------------------------------------------------------------------------------------------------------------

    do cNum=1,cNumMax
        wallTotalForceX(cNum) = wallTotalForceX(cNum)+Fxij(cNum)+Fwxij(cNum)
        wallTotalForceY(cNum) = wallTotalForceY(cNum)-(rhoSolid-rhoAvg)*Pi*radius0**2.0d0*gravity+Fyij(cNum)+Fwyij(cNum)
        totalTorque(cNum) = totalTorque(cNum)
        
        xCenterOld(cNum) = xCenter(cNum)
        yCenterOld(cNum) = yCenter(cNum)
        UcOld(cNum) = Uc(cNUm)
        VcOld(cNum) = Vc(cNum)
        rationalOmegaOld(cNum) = rationalOmega(cNum)
    enddo
    
    do cNum=1,cNumMax
        write(filename,*) cNum
        filename = adjustl(filename)

        open(unit=01,file='forceX-'//trim(filename)//'.dat',status='unknown',position='append')
        write(01,*) itc*t0, wallTotalForceX(cNum), Fxij(cNum)+Fwxij(cNum)
        close(01)

        open(unit=01,file='forceY-'//trim(filename)//'.dat',status='unknown',position='append')
        write(01,*) itc*t0, wallTotalForceY(cNum), Fyij(cNum)+Fwyij(cNum)
        close(01)
    enddo

    return
    end subroutine calForce