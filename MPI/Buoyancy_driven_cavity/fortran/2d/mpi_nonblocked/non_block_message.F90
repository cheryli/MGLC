subroutine f_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none

    ! calculate the boundary first
    call collision(1, nx, 1, 1)     ! bottom voundary
    call collision(1, nx, ny, ny)   ! top boundary
    call collision(1, 1, 1, ny)     ! left boundary
    call collision(nx, nx, 1, ny)   ! right boundary

    ! start meassage exchange
    ! message tag --- discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(f_post(2, 1, ny), 1, f_row_x, nbr_top, 2, comm2d, f_req(1), rc)
    call MPI_Irecv(f_post(2, 1, 0), 1, f_row_x, nbr_bottom, 2, comm2d, f_req(2), rc)

    call MPI_Isend(f_post(5, 1, ny), 1, f_row_x, nbr_top, 5, comm2d, f_req(3), rc)
    call MPI_Irecv(f_post(5, 1, 0), 1, f_row_x, nbr_bottom, 5, comm2d, f_req(4), rc)

    call MPI_Isend(f_post(6, 1, ny), 1, f_row_x, nbr_top, 6, comm2d, f_req(5), rc)
    call MPI_Irecv(f_post(6, 1, 0), 1, f_row_x, nbr_bottom, 6, comm2d, f_req(6), rc)

    ! message passing to bottom(j--)
    call MPI_Isend(f_post(4, 1, 1), 1, f_row_x, nbr_bottom, 4, comm2d, f_req(7), rc)
    call MPI_Irecv(f_post(4, 1, ny+1), 1, f_row_x, nbr_top, 4, comm2d, f_req(8), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, f_row_x, nbr_bottom, 7, comm2d, f_req(9), rc)
    call MPI_Irecv(f_post(7, 1, ny+1), 1, f_row_x, nbr_top, 7, comm2d, f_req(10), rc)

    call MPI_Isend(f_post(8, 1, 1), 1, f_row_x, nbr_bottom, 8, comm2d, f_req(11), rc)
    call MPI_Irecv(f_post(8, 1, ny+1), 1, f_row_x, nbr_top, 8, comm2d, f_req(12), rc)


    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(f_post(1, nx, 1), 1, f_column_y, nbr_right, 1, comm2d, f_req(13), rc)
    call MPI_Irecv(f_post(1, 0, 1), 1, f_column_y, nbr_left, 1, comm2d, f_req(14), rc)

    call MPI_Isend(f_post(5, nx, 1), 1, f_column_y, nbr_right, 5, comm2d, f_req(15), rc)
    call MPI_Irecv(f_post(5, 0, 1), 1, f_column_y, nbr_left, 5, comm2d, f_req(16), rc)

    call MPI_Isend(f_post(8, nx, 1), 1, f_column_y, nbr_right, 8, comm2d, f_req(17), rc)
    call MPI_Irecv(f_post(8, 0, 1), 1, f_column_y, nbr_left, 8, comm2d, f_req(18), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(f_post(3, 1, 1), 1, f_column_y, nbr_left, 3, comm2d, f_req(19), rc)
    call MPI_Irecv(f_post(3, nx+1, 1), 1, f_column_y, nbr_right, 3, comm2d, f_req(20), rc)

    call MPI_Isend(f_post(6, 1, 1), 1, f_column_y, nbr_left, 6, comm2d, f_req(21), rc)
    call MPI_Irecv(f_post(6, nx+1, 1), 1, f_column_y, nbr_right, 6, comm2d, f_req(22), rc)

    call MPI_Isend(f_post(7, 1, 1), 1, f_column_y, nbr_left, 7, comm2d, f_req(23), rc)
    call MPI_Irecv(f_post(7, nx+1, 1), 1, f_column_y, nbr_right, 7, comm2d, f_req(24), rc)   


    ! exchange message with corners --- message tag: discrete velocity
    ! exchange message along top-right (i++, j++)
    call MPI_Isend(f_post(5, nx, ny), 1, MPI_REAL8, cnr_top_right, 5, comm2d, f_req(25), rc)
    call MPI_Irecv(f_post(5, 0, 0), 1, MPI_REAL8, cnr_bottom_left, 5, comm2d, f_req(26), rc)

    ! exchange message along top-left (i--, j++)
    call MPI_Isend(f_post(6, 1, ny), 1, MPI_REAL8, cnr_top_left, 6, comm2d, f_req(27), rc)
    call MPI_Irecv(f_post(6, nx+1, 0), 1, MPI_REAL8, cnr_bottom_right, 6, comm2d, f_req(28), rc)

    ! exchange message along bottom-left (i--, j--)
    call MPI_Isend(f_post(7, 1, 1), 1, MPI_REAL8, cnr_bottom_left, 7, comm2d, f_req(29), rc)
    call MPI_Irecv(f_post(7, nx+1, ny+1), 1, MPI_REAL8, cnr_top_right, 7, comm2d, f_req(30), rc)

    ! exchange message along bottom-right (i++, j--)
    call MPI_Isend(f_post(8, nx, 1), 1, MPI_REAL8, cnr_bottom_right, 8, comm2d, f_req(31), rc)
    call MPI_Irecv(f_post(8, 0, ny+1), 1, MPI_REAL8, cnr_top_left, 8, comm2d, f_req(32), rc)


    ! then calculate the inner points
    call collision(2, nx-1, 2, ny-1)


end subroutine f_collision_with_message_exchange


subroutine g_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none

    ! calculate the boundary first
    call collisionT(1, nx, 1, 1)     ! bottom voundary
    call collisionT(1, nx, ny, ny)   ! top boundary
    call collisionT(1, 1, 1, ny)     ! left boundary
    call collisionT(nx, nx, 1, ny)   ! right boundary

    ! message tag:  0, 1, 2, discrete velocty

    ! ------------ exchange message along y ----------------
    ! message passing to top(j++)
    call MPI_Isend(g_post(2, 1, ny), 1, g_row_x, nbr_top, 2, comm2d, g_req(1), rc)
    call MPI_Irecv(g_post(2, 1, 0), 1, g_row_x, nbr_bottom, 2, comm2d, g_req(2), rc)


   ! message passing to bottom(j--)
    call MPI_Isend(g_post(4, 1, 1), 1, g_row_x, nbr_bottom, 4, comm2d, g_req(7), rc)
    call MPI_Irecv(g_post(4, 1, ny+1), 1, g_row_x, nbr_top, 4, comm2d, g_req(8), rc)

    ! ------------ exchange message along x ----------------
    ! message passing to right(i++)
    call MPI_Isend(g_post(1, nx, 1), 1, g_column_y, nbr_right, 1, comm2d, g_req(13), rc)
    call MPI_Irecv(g_post(1, 0, 1), 1, g_column_y, nbr_left, 1, comm2d, g_req(14), rc)
    
    ! message passing to left(i--)
    call MPI_Isend(g_post(3, 1, 1), 1, g_column_y, nbr_left, 3, comm2d, g_req(19), rc)
    call MPI_Irecv(g_post(3, nx+1, 1), 1, g_column_y, nbr_right, 3, comm2d, g_req(20), rc)
                
    ! then calculate the inner points
    call collision(2, nx-1, 2, ny-1)

end subroutine g_collision_with_message_exchange


subroutine collision(i_start, i_end, j_start, j_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:8), m_post(0:8), meq(0:8)
    real(kind=8) :: s(0:8)
    real(kind=8) :: fSource(0:8)


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

            Fx(i,j) = 0.0d0
            Fy(i,j) = rho(i,j)*gBeta*(T(i,j)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = (1.0d0-0.5d0*s(3))*Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = (1.0d0-0.5d0*s(5))*Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1.0d0-0.5d0*s(8))*(u(i,j)*Fy(i,j)+v(i,j)*Fx(i,j))

            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)
            enddo

    f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
    f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                    +m_post(7)/4.0d0
    f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                    -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)/4.0d0
    f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)/4.0d0
    f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)/4.0d0
    f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                    -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)/4.0d0

        enddo
    enddo
    
    return
end subroutine collision


subroutine collisionT(i_start, i_end, j_start, j_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    !------------------------
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: q(0:4)

    do j = j_start, j_end
        do i = i_start, i_end
            
    n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(1) = g(1,i,j)-g(3,i,j)
    n(2) = g(2,i,j)-g(4,i,j)
    n(3) = -4.0d0*g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
    n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)
        
            neq(0) = T(i,j)
            neq(1) = T(i,j)*u(i,j)
            neq(2) = T(i,j)*v(i,j)
            neq(3) = T(i,j)*paraA
            neq(4) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qnu
            q(4) = Qnu
        
            do alpha=0,4
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j) = 0.2d0*n_post(0)-0.2d0*n_post(3)
    g_post(1,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(2,i,j) = 0.2d0*n_post(0)+0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4)
    g_post(3,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(1)+0.05d0*n_post(3)+0.25d0*n_post(4)
    g_post(4,i,j) = 0.2d0*n_post(0)-0.5d0*n_post(2)+0.05d0*n_post(3)-0.25d0*n_post(4) 
        
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collisionT
