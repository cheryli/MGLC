subroutine output()
    use mpi
    use commondata
    integer :: i, j
    integer :: p_rank, num(0:1) ,dx = 0, dy = 0, new_coords(0:1)
    real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :)
    real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :)

    if (rank2d > 0) then
        ! rank != 0 send data
        num(0) = nx
        num(1) = ny
        ! send to rank 0
        call MPI_Send(num, 2, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size
        call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
        call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
        call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
    else
        ! rank 0 collect data
        allocate(total_u(total_nx, total_ny))
        allocate(total_v(total_nx, total_ny))
        allocate(total_rho(total_nx, total_ny))

        ! determine the origin
        if (nx > total_nx / dims(0)) then ! --- 5 5 '5' 4 4 4
            dx = nx * coords(0)
        else                    ! --- 5 5 5 4 '4' 4
            dx = nx * coords(0) + mod(total_nx, dims(0))
        endif

        if (ny > total_ny / dims(1)) then ! --- 5 5 '5' 4 4 4
            dy = ny * coords(1)
        else                    ! --- 5 5 5 4 '4' 4
            dy = ny * coords(1) + mod(total_ny, dims(1))
        endif

        ! collect data from rank 0
        do j = 1, ny
            do i = 1, nx
                total_u(dx + i, dy + j) = u(i, j)
                total_v(dx + i, dy + j) = v(i, j)
                total_rho(dx + i, dy + j) = rho(i, j)
            enddo
        enddo

        ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) - 1

            call MPI_Cart_coords(comm2d, p_rank, 2, new_coords, rc)

            ! receive the block size
            call MPI_Recv(num, 2, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)

            ! creat buffer
            allocate(tmp_u(num(0), num(1)))
            allocate(tmp_v(num(0), num(1)))
            allocate(tmp_rho(num(0), num(1)))
            
            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)


            ! determine the origin
            if (num(0) > total_nx / dims(0)) then ! --- 5 5 '5' 4 4 4
                dx = num(0) * new_coords(0)
            else                    ! --- 5 5 5 4 '4' 4
                dx = num(0) * new_coords(0) + mod(total_nx, dims(0))
            endif

            if (num(1) > total_ny / dims(1)) then ! --- 5 5 '5' 4 4 4
                dy = num(1) * new_coords(1)
            else                    ! --- 5 5 5 4 '4' 4
                dy = num(1) * new_coords(1) + mod(total_ny, dims(1))
            endif

            ! assign data
            do j = 1, num(1)
                do i = 1, num(0)
                    total_u(dx + i, dy + j) = tmp_u(i, j)
                    total_v(dx + i, dy + j) = tmp_v(i, j)
                    total_rho(dx + i, dy + j) = tmp_rho(i, j)
                enddo
            enddo

            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_rho)
        enddo

        call output_ASCII(xp, yp, total_u, total_v, total_rho, total_nx, total_ny, itc)
        call output_Tecplot(xp, yp, total_u, total_v, total_rho, total_nx, total_ny, itc)
        call output_binary(total_u, total_v, total_rho, total_nx, total_ny, itc)
        call getVelocity(xp, yp, total_u, total_v, total_rho, total_nx, total_ny, U0, itc)

        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_rho)

    endif

end subroutine output


subroutine output_ASCII(xp, yp, u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
        write(02,*) 'TITLE="Lid Driven Cavity"'
        write(02,*) 'VARIABLES="X" "Y" "U" "V" "Pressure" '
        write(02,101) nx, ny
        do j=1,ny
            do i=1,nx
                write(02,100) xp(i), yp(j), u(i,j), v(i,j), rho(i,j)/3.0d0
            enddo
        enddo
100     format(1x,2(e11.4,' '),10(e13.6,' '))
101     format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII


!!!c--------------------------------
!!!c--------------------------------
subroutine output_Tecplot(xp, yp, u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5
    integer, parameter :: kmax=1
    character(len=40) :: zoneName
    
    write(B2,'(i9.9)') itc

    open(41,file='MRTcavity-'//B2//'.plt',form='binary')
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

    !c-- Number of variables in this data file (here 5 variables)
    write(41) 5

    !c-- Variable names.
    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='Pressure'
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

    !---------Specify Var Location. 0 = Do not specify, all data
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
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(u(i,j))
                write(41) real(v(i,j))
                write(41) real(rho(i,j)/3.0d0)
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
!!!c--------------------------------
!!!c-------------------------------- 


subroutine getVelocity(xp, yp, u, v, rho, nx, ny, U0, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), U0
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer :: i, j
    integer :: nxHalf, nyHalf
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)

    nxHalf = (nx-1)/2 + 1
    nyHalf = (ny-1)/2 + 1

    open(unit=02,file='./u-y'//trim(filename)//'.dat',status='unknown')
    do j=1,ny
        write(02,*) u(nxHalf,j)/U0, yp(j)/dble(nx)
    enddo
    close(02)
    
    open(unit=03,file='./x-v'//trim(filename)//'.dat',status='unknown')
    do i=1,nx
        write(03,*) xp(i)/dble(nx), v(i,nyHalf)/U0
    enddo
    close(03)

    return
end subroutine getVelocity


subroutine output_binary(u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer :: i, j
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)
    
    open(unit=01,file='MRTcavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) ((u(i,j),i=1,nx),j=1,ny)
    write(01) ((v(i,j),i=1,nx),j=1,ny)
    write(01) ((rho(i,j),i=1,nx),j=1,ny)
    close(01)

    return
end subroutine output_binary