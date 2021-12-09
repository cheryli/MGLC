subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    !---------------------------
    real(8) :: s(0:18)
    real(8) :: m(0:18)
    real(8) :: m_post(0:18)
    real(8) :: meq(0:18)
    !---------------------------
    ! real(8) :: us2, un, feq
    !---------------------------

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
    ! loop body, re-indented for readability 
        
    m(0) = f(0,i,j,k) &
        + f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(1) = -30.0d0 * f(0,i,j,k) &
        - 11.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + 8.0d0 * (f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k)) &
        + 8.0d0 * (f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(2) = 12.0d0 * f(0,i,j,k) &
        - 4.0d0 * (f(1,i,j,k) + f(2,i,j,k) + f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
    m(3) = f(1,i,j,k) - f(2,i,j,k) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)
    m(4) = -4.0d0 * (f(1,i,j,k) - f(2,i,j,k)) &
        + f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) - f(12,i,j,k) &
        + f(13,i,j,k) - f(14,i,j,k)       
    m(5) = f(3,i,j,k) - f(4,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(6) = -4.0d0 * (f(3,i,j,k) - f(4,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) - f(9,i,j,k) - f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(7) = f(5,i,j,k) - f(6,i,j,k) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(8) = -4.0d0 * (f(5,i,j,k) - f(6,i,j,k)) &
        + f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) + f(15,i,j,k) + f(16,i,j,k) - f(17,i,j,k) - f(18,i,j,k)
    m(9) = 2.0d0 * (f(1,i,j,k) + f(2,i,j,k)) - f(3,i,j,k) - f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(10) = -4.0d0 * (f(1,i,j,k) + f(2,i,j,k)) + 2.0d0 * (f(3,i,j,k) + f(4,i,j,k) + f(5,i,j,k) + f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) + f(11,i,j,k) + f(12,i,j,k) &
        + f(13,i,j,k) + f(14,i,j,k) - 2.0 * (f(15,i,j,k) + f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k))
    m(11) = f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(12) = -2.0d0 * (f(3,i,j,k) + f(4,i,j,k) - f(5,i,j,k) - f(6,i,j,k)) &
        + f(7,i,j,k) + f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k)
    m(13) = f(7,i,j,k) - f(8,i,j,k) - f(9,i,j,k) + f(10,i,j,k)
    m(14) = f(15,i,j,k) - f(16,i,j,k) - f(17,i,j,k) + f(18,i,j,k)
    m(15) = f(11,i,j,k) - f(12,i,j,k) - f(13,i,j,k) + f(14,i,j,k)
    m(16) = f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) - f(10,i,j,k) - f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) + f(14,i,j,k)
    m(17) = -f(7,i,j,k) - f(8,i,j,k) + f(9,i,j,k) + f(10,i,j,k) &
        + f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) - f(18,i,j,k)
    m(18) = f(11,i,j,k) + f(12,i,j,k) &
        - f(13,i,j,k) - f(14,i,j,k) - f(15,i,j,k) - f(16,i,j,k) + f(17,i,j,k) + f(18,i,j,k)
        
        
        meq(0) = rho(i,j,k) 
        meq(1) = rho(i,j,k) * (-11.0d0 + 19.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(2) = rho(i,j,k) * (3.0d0 - 11.0d0/2.0d0 * (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)))
        meq(3) = rho(i,j,k) * u(i,j,k)
        meq(4) = -2.0d0/3.0d0 * rho(i,j,k) * u(i,j,k)
        meq(5) = rho(i,j,k) * v(i,j,k)
        meq(6) = -2.0d0/3.0d0 * rho(i,j,k) * v(i,j,k)
        meq(7) = rho(i,j,k) * w(i,j,k)
        meq(8) = -2.0d0/3.0d0 * rho(i,j,k) * w(i,j,k)
        meq(9) = rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(10) = -1.0d0/2.0d0 * rho(i,j,k) * (2.0d0 * u(i,j,k)*u(i,j,k) - v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(11) = rho(i,j,k) * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(12) = -1.0d0/2.0d0 * (v(i,j,k)*v(i,j,k) - w(i,j,k)*w(i,j,k))
        meq(13) = rho(i,j,k) * u(i,j,k) * v(i,j,k)
        meq(14) = rho(i,j,k) * v(i,j,k) * w(i,j,k)
        meq(15) = rho(i,j,k) * u(i,j,k) * w(i,j,k)
        meq(16) = 0.0d0
        meq(17) = 0.0d0
        meq(18) = 0.0d0


        s(0) = 0.0d0    !! s_{\rho}
        s(1) = Snu      !! s_{e}
        s(2) = Snu      !! s_{epsilon}
        s(3) = 0.0d0    !! s_{j}
        s(4) = Sq       !! s_{q}
        s(5) = 0.0d0    !! s_{j}
        s(6) = Sq       !! s_{q}
        s(7) = 0.0d0    !! s_{j}
        s(8) = Sq       !! s_{q}
        s(9) = Snu      !! s_{\nu}
        s(10) = Snu     !! s_{\pi}
        s(11) = Snu     !! s_{\nu}
        s(12) = Snu     !! s_{\pi}
        s(13) = Snu     !! s_{nu}
        s(14) = Snu     !! s_{nu}
        s(15) = Snu     !! s_{nu}
        s(16) = Sq      !! s_{m}
        s(17) = Sq      !! s_{m}
        s(18) = Sq      !! s_{m}
        
        do alpha=0,18
            m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))
        enddo
 
    f_post(0,i,j,k) = m_post(0)/19.0d0 - 5.0d0/399.0d0*m_post(1) + m_post(2)/21.0d0

    ! ------------
    f_post(1,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(3)/10.0d0 &
        - m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0

    f_post(2,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(3)/10.0d0 &
        + m_post(4)/10.0d0 + m_post(9)/18.0d0 - m_post(10)/18.0d0
    
    ! ------------
    f_post(3,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(5)/10.0d0 &
        - m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(4,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(5)/10.0d0 &
        + m_post(6)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 + m_post(11)/12.0d0 - m_post(12)/12.0d0
    
    f_post(5,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 + m_post(7)/10.0d0 &
        - m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    f_post(6,i,j,k) = m_post(0)/19.0d0 - 11.0d0/2394.0d0*m_post(1) - m_post(2)/63.0d0 - m_post(7)/10.0d0 &
        + m_post(8)/10.0d0 - m_post(9)/36.0d0 + m_post(10)/36.0d0 - m_post(11)/12.0d0 + m_post(12)/12.0d0

    ! ------------
    f_post(7,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 + m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(8,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(5)/10.0d0 + m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 - m_post(16)/8.0d0 - m_post(17)/8.0d0

    f_post(9,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 - m_post(13)/4.0d0 + m_post(16)/8.0d0 + m_post(17)/8.0d0

    f_post(10,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(5)/10.0d0 - m_post(6)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        + m_post(11)/12.0d0 + m_post(12)/24.0d0 + m_post(13)/4.0d0 - m_post(16)/8.0d0 + m_post(17)/8.0d0

    ! -----------
    f_post(11,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 - m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(12,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 + m_post(16)/8.0d0 + m_post(18)/8.0d0

    f_post(13,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(3)/10.0d0 &
        + m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 - m_post(15)/4.0d0 - m_post(16)/8.0d0 - m_post(18)/8.0d0
    
    f_post(14,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(3)/10.0d0 &
        - m_post(4)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 + m_post(9)/36.0d0 + m_post(10)/72.0d0 &
        - m_post(11)/12.0d0 - m_post(12)/24.0d0 + m_post(15)/4.0d0 + m_post(16)/8.0d0 - m_post(18)/8.0d0

    ! -----------
    f_post(15,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 + m_post(17)/8.0d0 - m_post(18)/8.0d0
    
    f_post(16,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 + m_post(7)/10.0d0 + m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 - m_post(17)/8.0d0 - m_post(18)/8.0d0

    f_post(17,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 + m_post(5)/10.0d0 &
        + m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        - m_post(14)/4.0d0 + m_post(17)/8.0d0 + m_post(18)/8.0d0

    f_post(18,i,j,k) = m_post(0)/19.0d0 + 4.0d0/1197.0d0*m_post(1) + m_post(2)/252.0d0 - m_post(5)/10.0d0 &
        - m_post(6)/40.0d0 - m_post(7)/10.0d0 - m_post(8)/40.0d0 - m_post(9)/18.0d0 - m_post(10)/36.0d0 &
        + m_post(14)/4.0d0 - m_post(17)/8.0d0 + m_post(18)/8.0d0


                ! us2 = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
                ! do alpha = 0, 18
                !     un = u(i,j,k)*ex(alpha) + v(i,j,k)*ey(alpha) + w(i,j,k)*ez(alpha)
                !     feq = rho(i,j,k) * omega(alpha) * (1.0d0 + 3.0d0*un + 4.5d0*un*un - 1.5d0*us2)
                    
                !     f_post(alpha,i,j,k) = f(alpha,i,j,k) - Snu*(f(alpha,i,j,k) - feq)
                ! enddo
    
            enddo
        enddo
    enddo

    return
end subroutine collision
