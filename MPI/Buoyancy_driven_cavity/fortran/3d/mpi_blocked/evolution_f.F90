subroutine collision()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,w,Fx,Fy,Fz,T) private(i,j,k,alpha,s,m,m_post,meq,fSource) 
    do k=1,nz
        do j=1,ny
            do i=1,nx
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


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    integer :: alpha
    
    !$omp parallel do default(none) shared(f,f_post,ex,ey,ez) private(i,j,k,ip,jp,kp,alpha)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do alpha=0,18
                    ip = i-ex(alpha)
                    jp = j-ey(alpha)
                    kp = k-ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j, k

#ifdef noslipWalls
    !Back plane (i = 1)
    if (coords(0) == 0) then
        do k=1,nz
            do j=1,ny
                f(1,1,j,k) = f_post(2,1,j,k)
                f(7,1,j,k) = f_post(10,1,j,k)
                f(9,1,j,k) = f_post(8,1,j,k)
                f(11,1,j,k) = f_post(14,1,j,k)
                f(13,1,j,k) = f_post(12,1,j,k)
            enddo
        enddo
    endif

    ! !Front plane (i=nx)
    if (coords(0) == dims(0) - 1) then
        do k=1,nz
            do j=1,ny
                f(2,nx,j,k) = f_post(1,nx,j,k)
                f(8,nx,j,k) = f_post(9,nx,j,k)
                f(10,nx,j,k) = f_post(7,nx,j,k)
                f(12,nx,j,k) = f_post(13,nx,j,k)
                f(14,nx,j,k) = f_post(11,nx,j,k)
            enddo
        enddo
    endif

    ! Left plane (j=1)
    if (coords(1) == 0) then
        do k=1,nz
            do i=1,nx
                f(3,i,1,k) = f_post(4,i,1,k)
                f(7,i,1,k) = f_post(10,i,1,k)
                f(8,i,1,k) = f_post(9,i,1,k)
                f(15,i,1,k) = f_post(18,i,1,k)
                f(17,i,1,k) = f_post(16,i,1,k)
            enddo
        enddo
    endif

    ! Right plane (j=ny)
    if  (coords(1) == dims(1) - 1) then
        do k=1,nz
            do i=1,nx
                !Right plane (j=ny)
                f(4,i,ny,k) = f_post(3,i,ny,k)
                f(9,i,ny,k) = f_post(8,i,ny,k)
                f(10,i,ny,k) = f_post(7,i,ny,k)
                f(16,i,ny,k) = f_post(17,i,ny,k)
                f(18,i,ny,k) = f_post(15,i,ny,k)
            enddo
        enddo
    endif
#endif

    ! Bottom side (k=1)
    if (coords(2) == 0) then
        do j=1,ny
            do i=1,nx
                f(5,i,j,1) = f_post(6,i,j,1)
                f(11,i,j,1) = f_post(14,i,j,1)
                f(12,i,j,1) = f_post(13,i,j,1)
                f(15,i,j,1) = f_post(18,i,j,1)
                f(16,i,j,1) = f_post(17,i,j,1)
            enddo
        enddo
    endif


    !Top side (k=nz)
    if (coords(2) == dims(2) - 1) then
        do j=1,ny
            do i=1,nx
                f(6,i,j,nz) = f_post(5,i,j,nz)
                f(13,i,j,nz) = f_post(12,i,j,nz)
                f(14,i,j,nz) = f_post(11,i,j,nz)
                f(17,i,j,nz) = f_post(16,i,j,nz)
                f(18,i,j,nz) = f_post(15,i,j,nz)
            enddo
        enddo
    endif

    return
end subroutine bounceback
    

subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k

    !$omp parallel do default(none) shared(f,rho,u,v,w,Fx,Fy,Fz) private(i,j,k)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k) = f(0,i,j,k)+f(1,i,j,k)+f(2,i,j,k)+f(3,i,j,k)+f(4,i,j,k)+f(5,i,j,k)+f(6,i,j,k) &
            +f(7,i,j,k)+f(8,i,j,k)+f(9,i,j,k)+f(10,i,j,k)+f(11,i,j,k)+f(12,i,j,k)+f(13,i,j,k)+f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)+f(17,i,j,k)+f(18,i,j,k)
            
                u(i,j,k) = ( f(1,i,j,k)-f(2,i,j,k)+f(7,i,j,k)-f(8,i,j,k)+f(9,i,j,k)-f(10,i,j,k)+f(11,i,j,k)-f(12,i,j,k)+f(13,i,j,k)-f(14,i,j,k)+0.5d0*Fx(i,j,k) )/rho(i,j,k)
            
                v(i,j,k) = ( f(3,i,j,k)-f(4,i,j,k)+f(7,i,j,k)+f(8,i,j,k)-f(9,i,j,k)-f(10,i,j,k)+f(15,i,j,k)-f(16,i,j,k)+f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fy(i,j,k) )/rho(i,j,k)
            
                w(i,j,k) = ( f(5,i,j,k)-f(6,i,j,k)+f(11,i,j,k)+f(12,i,j,k)-f(13,i,j,k)-f(14,i,j,k)+f(15,i,j,k)+f(16,i,j,k)-f(17,i,j,k)-f(18,i,j,k)+0.5d0*Fz(i,j,k) )/rho(i,j,k)
            
            enddo
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macro