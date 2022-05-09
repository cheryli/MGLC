subroutine collision_with_message_exchange()
    use mpi
    use commondata
    implicit none
    integer :: req(30)

    ! caculate the boundary first
    call collision(1, nx, 1, 1)     ! bottom voundary
    call collision(1, nx, ny, ny)   ! top boundary
    call collision(1, 1, 1, ny)     ! left boundary
    call collision(nx, nx, 1, ny)   ! right boundary

    ! start meassage exchange
    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(f_post(2, 1, ny), 1, row_x, nbr_top, 2, comm2d, req(1), rc)
    call MPI_Irecv(f_post(2, 1, 0), 1, row_x, nbr_bottom, 2, comm2d, req(2), rc)

    call MPI_Isend(f_post(5, 1, ny), 1, row_x, nbr_top, 5, comm2d, req(3), rc)
    call MPI_Irecv(f_post(5, 1, 0), 1, row_x, nbr_bottom, 5, comm2d, req(4), rc)

    call MPI_Isend(f_post(6, 1, ny), 1, row_x, nbr_top, 6, comm2d, req(5), rc)
    call MPI_Irecv(f_post(6, 1, 0), 1, row_x, nbr_bottom, 6, comm2d, req(6), rc)

    ! message passing to bottom(j--)
    call MPI_Isend(f_post(4, 1, 1), 1, row_x, nbr_bottom, 4, comm2d, req(7), rc)
    call MPI_Irecv(f_post(4, 1, ny+1), 1, row_x, nbr_top, 4, comm2d, req(8), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, row_x, nbr_bottom, 7, comm2d, req(9), rc)
    call MPI_Irecv(f_post(7, 1, ny+1), 1, row_x, nbr_top, 7, comm2d, req(10), rc)

    call MPI_Isend(f_post(8, 1, 1), 1, row_x, nbr_bottom, 8, comm2d, req(11), rc)
    call MPI_Irecv(f_post(8, 1, ny+1), 1, row_x, nbr_top, 8, comm2d, req(12), rc)


    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(f_post(1, nx, 1), 1, column_y, nbr_right, 1, comm2d, req(13), rc)
    call MPI_Irecv(f_post(1, 0, 1), 1, column_y, nbr_left, 1, comm2d, req(14), rc)

    call MPI_Isend(f_post(5, nx, 1), 1, column_y, nbr_right, 5, comm2d, req(15), rc)
    call MPI_Irecv(f_post(5, 0, 1), 1, column_y, nbr_left, 5, comm2d, req(16), rc)

    call MPI_Isend(f_post(8, nx, 1), 1, column_y, nbr_right, 8, comm2d, req(17), rc)
    call MPI_Irecv(f_post(8, 0, 1), 1, column_y, nbr_left, 8, comm2d, req(18), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(f_post(3, 1, 1), 1, column_y, nbr_left, 3, comm2d, req(19), rc)
    call MPI_Irecv(f_post(3, nx+1, 1), 1, column_y, nbr_right, 3, comm2d, req(20), rc)

    call MPI_Isend(f_post(6, 1, 1), 1, column_y, nbr_left, 6, comm2d, req(21), rc)
    call MPI_Irecv(f_post(6, nx+1, 1), 1, column_y, nbr_right, 6, comm2d, req(22), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, column_y, nbr_left, 7, comm2d, req(23), rc)
    call MPI_Irecv(f_post(7, nx+1, 1), 1, column_y, nbr_right, 7, comm2d, req(24), rc)   


    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Isend(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, comm2d, req(25), rc)
    call MPI_Irecv(f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, comm2d, req(26), rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Isend(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, comm2d, req(27), rc)
    call MPI_Irecv(f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, comm2d, req(28), rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Isend(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, comm2d, req(29), rc)
    call MPI_Irecv(f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, comm2d, req(30), rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Isend(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, comm2d, req(31), rc)
    call MPI_Irecv(f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, comm2d, req(32), rc)


    ! then calculate the inner points
    call collision(2, nx-1, 2, ny-1)

    ! wait the message to finish
    call MPI_Waitall(32, req, MPI_STATUSES_IGNORE, rc)


end subroutine collision_with_message_exchange


subroutine collision(i_start, i_end, j_start, j_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer :: i, j
    integer :: alpha
    !---------------------------
    real(8) :: s(0:8)
    real(8) :: m(0:8)
    real(8) :: m_post(0:8)
    real(8) :: meq(0:8)
    !---------------------------

    do j = j_start, j_end
        do i = i_start, i_end
        
    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
    m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
    m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
    m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
    m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*( u(i,j)*u(i,j)-v(i,j)*v(i,j) )
            meq(8) = rho(i,j)*( u(i,j)*v(i,j) )

            s(0) = 0.0d0      !!s_{\rho}
            s(1) = Snu !!s_{e}
            s(2) = Snu !!s_{\epsilon}
            s(3) = 0.0d0      !!s_{j} 
            s(4) = Sq !!s_{q}
            s(5) = 0.0d0      !!s_{j}
            s(6) = Sq       !!s_{q}
            s(7) = Snu !!s_{\nu}
            s(8) = Snu       !!s_{\nu}
            
            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0 
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)*0.25d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0
        enddo
    enddo

    return
end subroutine collision