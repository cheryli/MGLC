subroutine output()
    use mpi
    use commondata
    integer :: i, j
    integer :: p_rank, num(0:3) ,dx = 0, dy = 0
    real(8), allocatable :: total_u(:, :), total_v(:, :), total_rho(:, :)
    real(8), allocatable :: tmp_u(:, :), tmp_v(:, :), tmp_rho(:, :)
    ! real(8), allocatable :: stream(:, :), vorticity(:, :)
    ! integer :: total_obst(:, :), tmp_obst(:, :)

    if (rank2d > 0) then
        ! rank != 0 send data
        num(0) = nx
        num(1) = ny
        num(2) = i_start_global
        num(3) = j_start_global
        ! send to rank 0
        call MPI_Send(num, 4, MPI_INTEGER, 0, 0, comm2d, rc)    ! block size and origion
        call MPI_Send(u, nx*ny, MPI_REAL8, 0, 1, comm2d, rc)
        call MPI_Send(v, nx*ny, MPI_REAL8, 0, 2, comm2d, rc)
        call MPI_Send(rho, nx*ny, MPI_REAL8, 0, 3, comm2d, rc)
        ! call MPI_Send(T, nx*ny, MPI_REAL8, 0, 4, comm2d, rc)
    else
        ! rank 0 collect data
        allocate(total_u(total_nx, total_ny))
        allocate(total_v(total_nx, total_ny))
        allocate(total_rho(total_nx, total_ny))
        ! allocate(total_T(total_nx, total_ny))

        dx = i_start_global
        dy = j_start_global

        ! collect data from rank 0
        do j = 1, ny
            do i = 1, nx
                total_u(dx + i, dy + j) = u(i, j)
                total_v(dx + i, dy + j) = v(i, j)
                total_rho(dx + i, dy + j) = rho(i, j)
                ! total_T(dx + i, dy + j) = T(i, j)
            enddo
        enddo

         ! collect data from all other processors
        do p_rank = 1, dims(0) * dims(1) - 1
            ! receive the block size and origion
            call MPI_Recv(num, 4, MPI_INTEGER, p_rank, 0, comm2d, MPI_STATUS_IGNORE, rc)
            ! creat buffer
            allocate(tmp_u(num(0), num(1)))
            allocate(tmp_v(num(0), num(1)))
            allocate(tmp_rho(num(0), num(1)))
            ! allocate(tmp_T(num(0), num(1)))
            
            ! receive data
            call MPI_Recv(tmp_u, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 1, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_v, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 2, comm2d, MPI_STATUS_IGNORE, rc)
            call MPI_Recv(tmp_rho, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 3, comm2d, MPI_STATUS_IGNORE, rc)
            ! call MPI_Recv(tmp_T, num(0) * num(1), MPI_DOUBLE_PRECISION, p_rank, 4, comm2d, MPI_STATUS_IGNORE, rc)

            ! determine the origin
            dx = num(2)
            dy = num(3)

            ! assign data
            do j = 1, num(1)
                do i = 1, num(0)
                    total_u(dx + i, dy + j) = tmp_u(i, j)
                    total_v(dx + i, dy + j) = tmp_v(i, j)
                    total_rho(dx + i, dy + j) = tmp_rho(i, j)
                    ! total_T(dx + i, dy + j) = tmp_T(i, j)
                enddo
            enddo

            deallocate(tmp_u)
            deallocate(tmp_v)
            deallocate(tmp_rho)
            ! deallocate(tmp_T)
        enddo


        ! allocate(stream(total_nx, total_ny))
        ! allocate(vorticity(total_nx, total_ny))

        ! call compute_stream_vorticity(stream, vorticity, total_u, total_v, total_nx, total_ny)
        call output_Tecplot(X, Y, total_u, total_v, total_rho, total_nx, total_ny, itc)
        ! call output_binary(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        ! call output_ASCII(total_u, total_v, total_rho, total_T, total_nx, total_ny, itc)
        ! call out_Velocity_Nu(total_u, total_v, total_T, total_nx, total_ny, diffusivity, lengthUnit)

        deallocate(total_u)
        deallocate(total_v)
        deallocate(total_rho)
        ! deallocate(total_T)
        ! deallocate(stream)
        ! deallocate(vorticity)

    endif


end subroutine output




subroutine output_Tecplot(xp, yp, u, v, rho, nx, ny, itc)
    implicit none
    integer, intent(in) :: nx, ny, itc
    real(8), intent(in) :: xp(nx), yp(ny)
    real(8), intent(in) :: u(nx, ny), v(nx, ny), rho(nx, ny)
    integer(kind=4) :: i, j, k
    character(len=9) :: B2
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7,V8
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName
    

    write(B2,'(i9.9)') itc
    open(unit=41,file='movingCylinder-'//B2//'.plt', access='stream', form='unformatted')

    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0

    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) '#!TDV101'

    !c--Integer value of 1
    write(41) 1

    Title='MyFirst'
    call dumpstring(title)

    !c-- Number of variables in this data file

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
    V5='rho'
    call dumpstring(V5)
    ! V6='T'
    ! call dumpstring(V6)
    ! V7='stream'
    ! call dumpstring(V7)
    ! V8='vorticity'
    ! call dumpstring(V8)


    !c-----Zones-----------------------------

    !c--------Zone marker. Value = 299.0
    write(41) zoneMarker

    !--------Zone name.
    zoneName='ZONE 001'
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
    ! write(41) 1
    ! write(41) 1
    ! write(41) 1

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
                write(41) real(rho(i,j))
                ! write(41) real(T(i,j))
                ! write(41) real(stream(i,j))
                ! write(41) real(vorticity(i,j))
            end do
        end do
    enddo
    close(41)
    !---------------------------------------------------------------------

    return
end subroutine output_Tecplot


subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
end subroutine dumpstring