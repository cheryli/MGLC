subroutine bounceback_particle()
    use mpi
    use commondata
    implicit none
    integer :: i, j, alpha
    integer :: ip, jp
    real(8) :: q
    integer :: cNum
    integer :: myFlag
    integer :: fluidNum, total_fluidNum = 0
    real(8) :: temp1, temp2
    real(8) :: x0, y0
    real(8) :: total_rho = 0.0d0

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

    call MPI_Barrier(comm2d, rc)

    call MPI_ALLreduce(rhoAvg, total_rho, 1, MPI_REAL8, MPI_SUM, comm2d, rc)
    call MPI_ALLreduce(fluidNum, total_fluidNum, 1, MPI_INTEGER, MPI_SUM, comm2d, rc)

    rhoAvg = total_rho / dble(total_fluidNum)

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
                            if( ((ip + i_start_global-xCenter(cNum))**2.0d0+(jp + j_start_global - yCenter(cNum))**2.0d0).LE.radius(cNum)**2.0d0 ) then

                                myFlag = 1

                                call calQ(cNum,dble(i + i_start_global),dble(j + j_start_global),alpha,x0,y0,q) 

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
end subroutine bounceback_particle





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