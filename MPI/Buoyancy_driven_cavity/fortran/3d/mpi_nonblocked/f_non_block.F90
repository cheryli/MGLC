subroutine f_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none
    integer :: tag(5), idx


    !!! caculate the boundary first 
    call collision(1, 1, 1, ny, 1, nz)
    call collision(nx, nx, 1, ny, 1, nz)
    call collision(1, nx, 1, 1, 1, nz)
    call collision(1, nx, ny, ny, 1, nz)
    call collision(1, nx, 1, ny, 1, 1)
    call collision(1, nx, 1, ny, nz, nz)

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)s
    tag = (/ 1, 7, 9, 11, 13 /)
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), nx, 1, 1), 1, f_surface_x, nbr_surface(1), tag(idx), comm3d, f_req(0+idx), rc)
        call MPI_Irecv(f_post(tag(idx), 0, 1, 1), 1, f_surface_x, nbr_surface(2), tag(idx), comm3d, f_req(5+idx), rc)
    enddo

    ! message passing to (i--)
    tag = (/ 2, 8, 10 ,12, 14 /) 
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), 1, 1, 1), 1, f_surface_x, nbr_surface(2), tag(idx), comm3d, f_req(10+idx), rc)
        call MPI_Irecv(f_post(tag(idx), nx+1, 1, 1), 1, f_surface_x, nbr_surface(1), tag(idx), comm3d, f_req(15+idx), rc)
    enddo

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    tag = (/ 3, 7, 8, 15, 17 /)
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), 1, ny, 1), 1, f_surface_y, nbr_surface(3), tag(idx), comm3d, f_req(20+idx), rc)
        call MPI_Irecv(f_post(tag(idx), 1, 0, 1), 1, f_surface_y, nbr_surface(4), tag(idx), comm3d, f_req(25+idx), rc)
    enddo

    ! message passing to (j--)
    tag = (/ 4, 9, 10, 16, 18 /)
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), 1, 1, 1), 1, f_surface_y, nbr_surface(4), tag(idx), comm3d, f_req(30+idx), rc)
        call MPI_Irecv(f_post(tag(idx), 1, ny+1, 1), 1, f_surface_y, nbr_surface(3), tag(idx), comm3d, f_req(35+idx), rc)
    enddo

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    tag = (/ 5, 11, 12, 15, 16 /)
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), 1, 1, nz), 1, f_surface_z, nbr_surface(5), tag(idx), comm3d, f_req(40+idx), rc)
        call MPI_Irecv(f_post(tag(idx), 1, 1, 0), 1, f_surface_z, nbr_surface(6), tag(idx), comm3d, f_req(45+idx), rc)
    enddo

    ! message passing to (k--)
    tag = (/ 6, 13, 14, 17, 18 /)
    do idx = 1, 5
        call MPI_Isend(f_post(tag(idx), 1, 1, 1), 1, f_surface_z, nbr_surface(6), tag(idx), comm3d, f_req(50+idx), rc)
        call MPI_Irecv(f_post(tag(idx), 1, 1, nz+1), 1, f_surface_z, nbr_surface(5), tag(idx), comm3d, f_req(55+idx), rc)
    enddo

    ! /* ----------------------------- line data ----------------------------- */
    ! ------------ exchange message perpendicular to z ----------------
    ! message passing to (i++, j++) (7)
    call MPI_Isend(f_post(7, nx, ny, 1), 1, f_line_z, nbr_line(7), 7, comm3d, f_req(61), rc)
    call MPI_Irecv(f_post(7, 0, 0, 1), 1, f_line_z, nbr_line(10), 7, comm3d, f_req(62), rc)

    ! message passing to (i--, j--) (10)
    call MPI_Isend(f_post(10, 1, 1, 1), 1, f_line_z, nbr_line(10), 10, comm3d, f_req(63), rc)
    call MPI_Irecv(f_post(10, nx+1, ny+1, 1), 1, f_line_z, nbr_line(7), 10, comm3d, f_req(64), rc)

    ! message passing to (i++, j--) (9)
    call MPI_Isend(f_post(9, nx, 1, 1), 1, f_line_z, nbr_line(9), 9, comm3d, f_req(65), rc)
    call MPI_Irecv(f_post(9, 0, ny+1, 1), 1, f_line_z, nbr_line(8), 9, comm3d, f_req(66), rc)

    ! message passing to (i--, j++) (8)
    call MPI_Isend(f_post(8, 1, ny, 1), 1, f_line_z, nbr_line(8), 8, comm3d, f_req(67), rc)
    call MPI_Irecv(f_post(8, nx+1, 0, 1), 1, f_line_z, nbr_line(9), 8, comm3d, f_req(68), rc)

    ! ------------ exchange message perpendicular to y ----------------
    ! message passing to (i++, k++) (11)
    call MPI_Isend(f_post(11, nx, 1, nz), 1, f_line_y, nbr_line(11), 11, comm3d, f_req(69), rc)
    call MPI_Irecv(f_post(11, 0, 1, 0), 1, f_line_y, nbr_line(14), 11, comm3d, f_req(70), rc)

    ! message passing to (i--, k--) (14)
    call MPI_Isend(f_post(14, 1, 1, 1), 1, f_line_y, nbr_line(14), 14, comm3d, f_req(71), rc)
    call MPI_Irecv(f_post(14, nx+1, 1, nz+1), 1, f_line_y, nbr_line(11), 14, comm3d, f_req(72), rc)

    ! message passing to (i++, k--) (13)
    call MPI_Isend(f_post(13, nx, 1, 1), 1, f_line_y, nbr_line(13), 13, comm3d, f_req(73), rc)
    call MPI_Irecv(f_post(13, 0, 1, nz+1), 1, f_line_y, nbr_line(12), 13, comm3d, f_req(74), rc)

    ! message passing to (i--, k++) (12)
    call MPI_Isend(f_post(12, 1, 1, nz), 1, f_line_y, nbr_line(12), 12, comm3d, f_req(75), rc)
    call MPI_Irecv(f_post(12, nx+1, 1, 0), 1, f_line_y, nbr_line(13), 12, comm3d, f_req(76), rc)

    ! ------------ exchange message perpendicular to x ----------------
    ! message passing to (j++, k++) (15)\
    call MPI_Isend(f_post(15, 1, ny, nz), 1, f_line_x, nbr_line(15), 15, comm3d, f_req(77), rc)
    call MPI_Irecv(f_post(15, 1, 0, 0), 1, f_line_x, nbr_line(18), 15, comm3d, f_req(78), rc)

    ! message passing to (j--, k--) (18)
    call MPI_Isend(f_post(18, 1, 1, 1), 1, f_line_x, nbr_line(18), 18, comm3d, f_req(79), rc)
    call MPI_Irecv(f_post(18, 1, ny+1, nz+1), 1, f_line_x, nbr_line(15), 18, comm3d, f_req(80), rc)

    ! message passing to (j++, k--) (17)
    call MPI_Isend(f_post(17, 1, ny, 1), 1, f_line_x, nbr_line(17), 17, comm3d, f_req(81), rc)
    call MPI_Irecv(f_post(17, 1, 0, nz+1), 1, f_line_x, nbr_line(16), 17, comm3d, f_req(82), rc)

    ! message passing to (j--, k++) (16)
    call MPI_Isend(f_post(16, 1, 1, nz), 1, f_line_x, nbr_line(16), 16, comm3d, f_req(83), rc)
    call MPI_Irecv(f_post(16, 1, ny+1, 0), 1, f_line_x, nbr_line(17), 16, comm3d, f_req(84), rc)
    
    ! don't need to echange point data

    !!! then compute the inner points
    call collision(2, nx-1, 2, ny-1, 2, nz-1)
    
    ! wait the message to finish
    call MPI_Waitall(84, f_req, MPI_STATUSES_IGNORE, rc)

end subroutine f_collision_with_message_exchange




subroutine collision(i_start, i_end, j_start, j_end, k_start, k_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,w,Fx,Fy,Fz,T) private(i,j,k,alpha,s,m,m_post,meq,fSource) 
    do k = k_start, k_end
        do j = j_start, j_end
            do i = i_start, i_end
    !--------------------------------------------------------------------------------------------------------------------
    !---m0    
    m(0) =f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m1
    m(1) = -30.0d0*f(0,i,j,k)-11.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
    +8.0d0*( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m2
    m(2) = 12.0d0*f(0,i,j,k)-4.0d0*( f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) ) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !---m3
    m(3) = f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m4
    m(4) = -4.0d0*(f(1,i,j,k)-f(2,i,j,k))+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k) &
                    +f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)
    !---m5
    m(5) = f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m6
    m(6) = -4.0d0*(f(3,i,j,k)-f(4,i,j,k))+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k) &
                +f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m7
    m(7) = f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m8
    m(8) = -4.0d0*(f(5,i,j,k)-f(6,i,j,k))+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k) &
            +f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)
    !---m9
    m(9) = 2.0d0*(f(1,i,j,k)+f(2,i,j,k))-f(3,i,j,k)-f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m10
     m(10) = -4.0d0*(f(1,i,j,k)+f(2,i,j,k))+2.0d0*(f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k)) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)-2.0d0*( f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k) )
    !---m11
    m(11) = f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k)+f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m12
     m(12) = -2.0d0*(f(3,i,j,k)+f(4,i,j,k)-f(5,i,j,k)-f(6,i,j,k))+( f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k) )-( f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k) )
    !---m13
    m(13) = f(7,i,j,k)-f(8,i,j,k)-f(9,i,j,k)+f(10,i,j,k)
    !---m14
    m(14) = f(15,i,j,k)-f(16,i,j,k)-f(17,i,j,k)+f(18,i,j,k)
    !---m15
    m(15) = f(11,i,j,k)-f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m16
    m(16) = f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)-f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)+f(14,i,j,k)
    !---m17
    m(17) = -f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)
    !---m18
    m(18) = f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)-f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
    !--------------------------------------------------------------------------------------------------------------------

    meq(0) = rho(i,j,k)
    meq(1) = -11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(2) = 3.0d0*rho(i,j,k)-11.0d0/2.0d0*rho(i,j,k)*( u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k) )
    meq(3) = rho(i,j,k)*u(i,j,k)
    meq(4) = -2.0d0/3.0d0*meq(3)
    meq(5) = rho(i,j,k)*v(i,j,k)
    meq(6) = -2.0d0/3.0d0*meq(5)
    meq(7) = rho(i,j,k)*w(i,j,k)
    meq(8) = -2.0d0/3.0d0*meq(7)
    meq(9) = rho(i,j,k)*(2.0d0*u(i,j,k)*u(i,j,k)-v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(10) = -0.5d0*meq(9)
    meq(11) = rho(i,j,k)*(v(i,j,k)*v(i,j,k)-w(i,j,k)*w(i,j,k))
    meq(12) = -0.5d0*meq(11)
    meq(13) = rho(i,j,k)*(u(i,j,k)*v(i,j,k))
    meq(14) = rho(i,j,k)*(v(i,j,k)*w(i,j,k))
    meq(15) = rho(i,j,k)*(w(i,j,k)*u(i,j,k))
    meq(16) = 0.0d0
    meq(17) = 0.0d0
    meq(18) = 0.0d0

            s(0) = 0.0d0 
            s(1) = Snu  !!!s_{e}
            s(2) = Snu   !!! s_{epsilon}
            s(3) = 0.0d0 
            s(4) = Sq   !!! s_{q}
            s(5) = 0.0d0 
            s(6) = Sq   !!! s_{q}
            s(7) = 0.0d0 
            s(8) = Sq   !!! s_{q}
            s(9) = Snu !!! s_{nu}
            s(10) = Snu   !!! s_{pi}
            s(11) = Snu   !!! s_{nu}
            s(12) = Snu !!! s_{pi}
            s(13) = Snu !!! s_{nu}
            s(14) = Snu   !!! s_{nu}
            s(15) = Snu   !!! s_{nu}
            s(16) = Sq   !!! s_{m}
            s(17) = Sq   !!! s_{m}
            s(18) = Sq   !!! s_{m}

            Fx(i,j,k) = -2.0d0*rho(i,j,k)*v(i,j,k)*omegaRatating
            Fy(i,j,k) = 2.0d0*rho(i,j,k)*u(i,j,k)*omegaRatating
            Fz(i,j,k) = rho(i,j,k)*gBeta*(T(i,j,k)-Tref)

            fSource(0) = 0.0d0
            fSource(1) = 38.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(2) = -11.0d0*(u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k))
            fSource(3) = Fx(i,j,k)
            fSource(4) = -2.0d0/3.0d0*Fx(i,j,k)
            fSource(5) = Fy(i,j,k)
            fSource(6) = -2.0d0/3.0d0*Fy(i,j,k)
            fSource(7) = Fz(i,j,k)
            fSource(8) = -2.0d0/3.0d0*Fz(i,j,k)
            fSource(9) = 4.0d0*u(i,j,k)*Fx(i,j,k)-2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(10) = -2.0d0*u(i,j,k)*Fx(i,j,k)+v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k) 
            fSource(11) = 2.0d0*v(i,j,k)*Fy(i,j,k)-2.0d0*w(i,j,k)*Fz(i,j,k)
            fSource(12) = -v(i,j,k)*Fy(i,j,k)+w(i,j,k)*Fz(i,j,k)
            fSource(13) = u(i,j,k)*Fy(i,j,k)+v(i,j,k)*Fx(i,j,k)
            fSource(14) = v(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fy(i,j,k)
            fSource(15) = u(i,j,k)*Fz(i,j,k)+w(i,j,k)*Fx(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+(1.0d0-0.5d0*s(alpha))*fSource(alpha)
            enddo

    f_post(0,i,j,k) = m_post(0)/19.0d0-5.0d0/399.0d0*m_post(1)+m_post(2)/21.0d0

    f_post(1,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(3)-m_post(4) )*0.1d0+( m_post(9)-m_post(10) )/18.0d0

    f_post(3,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(4,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(5)-m_post(6) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
+( m_post(11)-m_post(12) )/12.0d0

    f_post(5,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
+( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0-11.0d0/2394.0d0*m_post(1)-m_post(2)/63.0d0 &
-( m_post(7)-m_post(8) )*0.1d0-( m_post(9)-m_post(10) )/36.0d0 &
-( m_post(11)-m_post(12) )/12.0d0

!---
    f_post(7,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0+( m_post(16)-m_post(17) )*0.125d0

    f_post(8,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0-( m_post(16)+m_post(17) )*0.125d0

    f_post(9,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(5)-m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
-m_post(13)*0.25d0+( m_post(16)+m_post(17) )*0.125d0

    f_post(10,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(5)+m_post(6) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0+m_post(11)/12.0d0+m_post(12)/24.0d0 &
+m_post(13)*0.25d0-( m_post(16)-m_post(17) )*0.125d0

!---
    f_post(11,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)-0.1250d0*( m_post(16)-m_post(18) )

    f_post(12,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)+0.125d0*( m_post(16)+m_post(18) )

    f_post(13,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+0.025d0*( 4.0d0*m_post(3)+m_post(4)-4.0d0*m_post(7)-m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
-0.25d0*m_post(15)-0.125d0*( m_post(16)+m_post(18) )

    f_post(14,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-0.025d0*( 4.0d0*m_post(3)+m_post(4)+4.0d0*m_post(7)+m_post(8) ) &
+m_post(9)/36.0d0+m_post(10)/72.0d0-m_post(11)/12.0d0-m_post(12)/24.0d0 &
+0.25d0*m_post(15)+0.125d0*( m_post(16)-m_post(18) )

!---
    f_post(15,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)+0.125d0*( m_post(17)-m_post(18) )

    f_post(16,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)-0.125d0*( m_post(17)+m_post(18) )

    f_post(17,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
+( 4.0d0*m_post(5)+m_post(6)-4.0d0*m_post(7)-m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
-0.25d0*m_post(14)+0.125d0*( m_post(17)+m_post(18) )

    f_post(18,i,j,k) = m_post(0)/19.0d0+4.0d0/1197.0d0*m_post(1)+m_post(2)/252.0d0 &
-( 4.0d0*m_post(5)+m_post(6)+4.0d0*m_post(7)+m_post(8) )*0.025d0 &
-( m_post(9)+m_post(10)*0.5d0 )/18.0d0 &
+0.25d0*m_post(14)-0.125d0*( m_post(17)-m_post(18) )

            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collision