subroutine output()
    use mpi
    use commondata
    implicit none
    integer :: i, j, k
    integer :: p_rank, num(0:2), new_coords(0:2), dx = 0, dy = 0, dz = 0
    real(8), allocatable :: total_u(:, :, :), total_v(:, :, :), total_w(:, :, :), total_rho(:, :, :)
    real(8), allocatable :: tmp_u(:, :, :), tmp_v(:, :, :), tmp_w(:, :, :), tmp_rho(:, :, :)
    
    ! use rank 0 to receive data and output the results

    if (rank3d > 0) then  !!! ----  rank != 0 send data
        ! collect the rank information
        num(0) = nx
        num(1) = ny
        num(2) = nz
        ! send to rank 0
        call MPI_Send(num, 3, MPI_INTEGER, 0, 0, comm3d, rc)    ! rank information
        call MPI_Send(u, nx*ny*nz, MPI_REAL8, 0, 1, comm3d, rc)
        call MPI_Send(v, nx*ny*nz, MPI_REAL8, 0, 2, comm3d, rc)
        call MPI_Send(w, nx*ny*nz, MPI_REAL8, 0, 3, comm3d, rc)
        call MPI_Send(rho, nx*ny*nz, MPI_REAL8, 0, 4, comm3d, rc)
    else    
        !!! ---- rank 0 collect data
        ! allocate array
        allocate(total_u(total_nx, total_ny, total_nz))
        allocate(total_v(total_nx, total_ny, total_nz))
        allocate(total_w(total_nx, total_ny, total_nz))
        allocate(total_rho(total_nx, total_ny, total_nz))

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

        if (nz > total_nz / dims(2)) then
            dz = nz * coords(2)
        else
            dz = nz * coords(2) + mod(total_nz, dims(2))
        endif

        ! collect data from rank 0
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    total_u(dx + i, dy + j, dz + k) = u(i, j, k)
                    total_v(dx + i, dy + j, dz + k) = v(i, j, k)
                    total_w(dx + i, dy + j, dz + k) = w(i, j, k)
                    total_rho(dx + i, dy + j, dz + k) = rho(i, j, k)
                enddo
            enddo
        enddo

        ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) * dims(2) - 1

            call MPI_Cart_coords(comm3d, p_rank, 3, new_coords, rc)

            ! receive the rank information
            call MPI_Recv(num, 3, MPI_INTEGER, p_rank, 0, comm3d, MPI_STATUS_IGNORE, rc)

            ! creat buffer
            allocate(tmp_u(num(0), num(1), num(2)))
            allocate(tmp_v(num(0), num(1), num(2)))
            allocate(tmp_w(num(0), num(1), num(2)))
            allocate(tmp_rho(num(0), num(1), num(2)))

            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 1, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 2, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_w, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 3, comm3d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1) * num(2), MPI_REAL8, p_rank, 4, comm3d, MPI_STATUS_IGNORE, rc)


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

            if (num(2) > total_nz / dims(2)) then
                dz = num(2) * new_coords(2)
            else
                dz = num(2) * new_coords(2) + mod(total_nz, dims(2))
            endif

            ! assign data
            do k = 1, num(2)
                do j = 1, num(1)
                    do i = 1, num(0)
                        total_u(i + dx, j + dy, k + dz) = tmp_u(i, j, k)
                        total_v(i + dx, j + dy, k + dz) = tmp_v(i, j, k)
                        total_w(i + dx, j + dy, k + dz) = tmp_w(i, j, k)
                        total_rho(i + dx, j + dy, k + dz) = tmp_rho(i, j, k)
                    enddo
                enddo
            enddo

            ! de-allocate buffer arrays
            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_w)
            deallocate(tmp_rho)
        enddo

        ! after collect total_* data, then output
        ! call output_ASCII(xp, yp, zp, total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        ! call output_binary(total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        call output_Tecplot(xp, yp, zp, total_u, total_v, total_w, total_rho, total_nx, total_ny, total_nz, itc)
        call getVelocity(xp, yp, zp, total_u, total_v, total_w, total_nx, total_ny, total_nz, U0, itc)

        ! de-allocate total arrays
        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_w)
        deallocate(total_rho)
    endif

end subroutine output








subroutine output_ASCII(xp, yp, zp, u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=77,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')
        write(77,*) 'TITLE="Lid Driven Cavity"'
        write(77,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "Pressure" '
        write(77,101) nx, ny, nz
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(77,100) xp(i), yp(j), zp(k), u(i,j,k), v(i,j,k), w(i,j,k), rho(i,j,k)/3.0d0
                enddo
            enddo
        enddo
100     format(1x,3(e11.4,' '),10(e13.6,' '))
101     format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'K=',1x,i5,1x,'F=POINT')
    close(77)

    return
end subroutine output_ASCII


!!!c--------------------------------
!!!c--------------------------------
subroutine output_Tecplot(xp, yp, zp, u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=9) :: B2
    REAL(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    
    write(B2,'(i9.9)') itc

    open(41,file='MRTcavity-'//B2//'.plt', access='stream', form='unformatted')
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
    write(41) 7

    !c-- Variable names.
    V1 = 'X'
    call dumpstring(V1)
    V2 = 'Y'
    call dumpstring(V2)
    V3 = 'Z'
    call dumpstring(V3)
    V4 = 'U'
    call dumpstring(V4)
    V5 = 'V'
    call dumpstring(V5)
    V6 = 'W'
    call dumpstring(V6)
    V7 = 'Pressure'
    call dumpstring(V7)

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
    write(41) nz

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
    write(41) 1
    write(41) 1

    !--------Has variable sharing 0 = no, 1 = yes.
    write(41) 0

    !----------Zone number to share connectivity list with (-1 = no
    ! sharing).
    write(41) -1

    !---------------------------------------------------------------------
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(41) real(xp(i))
                write(41) real(yp(j))
                write(41) real(zp(k))
                write(41) real(u(i,j,k))
                write(41) real(v(i,j,k))
                write(41) real(w(i,j,k))
                write(41) real(rho(i,j,k)/3.0d0)
            enddo
        enddo
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


subroutine getVelocity(xp, yp, zp, u, v, w, nx, ny, nz, U0, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: xp(0:nx+1), yp(0:ny+1), zp(0:nz+1)
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz)
    real(8), intent(in) :: U0
    integer :: i, j, k
    integer :: nxHalf, nyHalf, nzHalf
    character(len=9) :: B2

    write(B2,'(i9.9)') itc

    nxHalf = (nx - 1) / 2 + 1
    nyHalf = (ny - 1) / 2 + 1
    nzHalf = (nz - 1) / 2 + 1

    open(unit=02,file='./u-z_'//B2//'.dat',status='unknown')
    do k=1,nz
        write(02,*) u(nxHalf, nyHalf, k)/U0, zp(k)/dble(nz)
    enddo
    close(02)

    open(unit=03,file='./x-w_'//B2//'.dat',status='unknown')
    do i=1,nx
        write(03,*) xp(i)/dble(nx), w(i, nyHalf, nzHalf)/U0
    enddo
    close(03)

    return
end subroutine getVelocity


subroutine output_binary(u, v, w, rho, nx, ny, nz, itc)
    implicit none
    integer, intent(in) :: nx, ny, nz, itc
    real(8), intent(in) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), rho(nx, ny, nz)
    integer :: i, j, k
    character(len=100) :: filename
    
    write(filename,*) itc
    filename = adjustl(filename)
    
    open(unit=01,file='MRTcavity-'//trim(filename)//'.bin',form="unformatted",access="sequential")
    write(01) (((u(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((v(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    write(01) (((rho(i,j,k), i = 1, nx), j = 1, ny), k = 1, nz)
    close(01)

    return
end subroutine output_binary