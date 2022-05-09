subroutine calForce()
    use mpi
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
    real(8) :: force_x(cNumMax), force_y(cNumMax), torque(cNumMax)
    real(8) :: dw
    real(8) :: Fwxij(cNumMax), Fwyij(cNumMax)
    real(8) :: forceScale

    ! do cNum=1,cNumMax
    !     wallTotalForceX(cNum) = 0.0d0
    !     wallTotalForceY(cNum) = 0.0d0
    !     totalTorque(cNum) = 0.0d0
    ! enddo

    do cNum=1,cNumMax
        force_x(cNumMax) = 0.0d0
        force_y(cNumMax) = 0.0d0
        torque(cNumMax) = 0.0d0
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
                            if( ((ip + i_start_global - xCenter(cNum))**2.0d0 &
                            +(jp + j_start_global - yCenter(cNum))**2.0d0) &
                            .LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i + i_start_global),dble(j + j_start_global),alpha,x0,y0,q) 
                                
                                temp1 = -( y0-yCenter(cNum) )*rationalOmega(cNum)
                                temp2 = ( x0-xCenter(cNum) )*rationalOmega(cNum)

    !---------No. 2: Wen JCP 2014----------------------------------------
    tempForceX = (ex(alpha)-Uc(cNum)-temp1)*f_post(alpha,i,j)-(ex(r(alpha))-Uc(cNum)-temp1)*f(r(alpha),i,j)
    tempForceY = (ey(alpha)-Vc(cNum)-temp2)*f_post(alpha,i,j)-(ey(r(alpha))-Vc(cNum)-temp2)*f(r(alpha),i,j)
    tempTorque = (x0-xCenter(cNum))*tempForceY-(y0-yCenter(cNum))*tempForceX 
                           
                                force_x(cNum) = force_x(cNum) + tempForceX 
                                force_y(cNum) = force_y(cNum) + tempForceY
                                torque(cNum) = torque(cNum) + tempTorque
                                ! wallTotalForceX(cNum) = wallTotalForceX(cNum)+tempForceX 
                                ! wallTotalForceY(cNum) = wallTotalForceY(cNum)+tempForceY 
                                ! totalTorque(cNum) = totalTorque(cNum)+tempTorque 
                                
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

    call MPI_ALLreduce(force_x, wallTotalForceX, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(force_y, wallTotalForceY, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(torque, totalTorque, cNumMax, MPI_REAL8, MPI_SUM, comm2d, rc)

 
 !---------------particle-wall (spring force model)------------------
    do cNum=1,cNumMax
        forceScale = Pi*radius(cNum)**2.0d0*(rhoSolid-rho0)*gravity/wallStiff

        dw = yCenter(cNum)-radius(cNum)-1.0d0

        if (rank == 0) then   
            if(MOD(itc,2000).EQ.0) then
                open(unit=01,file='dw.dat',status='unknown',position='append')
                write(01,*) itc, dw
                close(01)
            endif
        endif
            
        if(dw.GE.thresholdR) then
            Fwxij(cNum) = 0.0d0
            Fwyij(cNum) = 0.0d0
        elseif( (dw.LT.thresholdR).AND.(dw.GE.0.0d0) ) then
            Fwxij(cNum) = 0.0d0
            Fwyij(cNum) = forceScale*((dw-thresholdR)/thresholdR)**2.0d0
            if(MOD(itc,1000).EQ.0) then
                if (rank == 0) then   
                    write(*,*) "    "
                    write(*,*) "Repulsive force activated (near bottom wall)!!!!!! itc=",itc
                    write(*,*) "    "
                endif
            endif
        else
            if (rank == 0) then   
            write(*,*) "WARNING: particle penetrates the bottom wall!"
            write(*,*) "yCenter(cNum)=",yCenter(cNum)
            write(*,*) "Error: simulation stoped!"
            endif
            stop
        endif
    enddo
 !---------------particle-wall (spring force model)------------------
 
    do cNum=1,cNumMax        
        wallTotalForceX(cNum) = wallTotalForceX(cNum)+Fwxij(cNum)
        wallTotalForceY(cNum) = wallTotalForceY(cNum)-(rhoSolid-rhoAvg)*Pi*radius0**2.0d0*gravity+Fwyij(cNum)
        totalTorque(cNum) = totalTorque(cNum)
        
        xCenterOld(cNum) = xCenter(cNum)
        yCenterOld(cNum) = yCenter(cNum)
        UcOld(cNum) = Uc(cNUm)
        VcOld(cNum) = Vc(cNum)
        rationalOmegaOld(cNum) = rationalOmega(cNum)
    enddo

    return
end subroutine calForce