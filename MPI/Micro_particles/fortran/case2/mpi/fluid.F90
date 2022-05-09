subroutine collision()
    use commondata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: m(0:8), m_post(0:8), meq(0:8)
    real(8) :: s(0:8)

    !$omp parallel do default(none) shared(f,f_post,rho,u,v,obst) private(i,j,alpha,s,m,m_post,meq) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then

    m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
    m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*f(5,i,j)+2.0d0*f(6,i,j)+2.0d0*f(7,i,j)+2.0d0*f(8,i,j)
    m(2) = 4.0d0*f(0,i,j)-2.0d0*f(1,i,j)-2.0d0*f(2,i,j)-2.0d0*f(3,i,j)-2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
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
            meq(4) = -meq(3)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -meq(5)
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

            endif
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine collision


subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, alpha
    integer :: ip, jp
    real(8) :: q
    integer :: cNum
    integer :: myFlag
    integer :: fluidNum
    real(8) :: temp1, temp2
    real(8) :: x0, y0
    
    !$omp parallel do default(none) shared(ex,ey,f,f_post,obst) private(i,j,alpha,ip,jp)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                do alpha=0,8
                    ip = i+ex(alpha)
                    jp = j+ey(alpha)

                    f(alpha,ip,jp) = f_post(alpha,i,j)
                
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine streaming


subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j

    do j=1,ny
        !Left side
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !Right side
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo

    do i=1,nx
        !Top side
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)
        f(8,i,ny) = f_post(6,i,ny)

        !Bottom side
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)
        f(6,i,1) = f_post(8,i,1)
    enddo

    return
end subroutine bounceback


subroutine macro()
    use commondata
    implicit none
    integer :: i, j

    !$omp parallel do default(none) shared(f,rho,u,v,obst) private(i,j) 
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0)  then
            
                rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
                u(i,j) = ( f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j) )/rho(i,j)
                v(i,j) = ( f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j) )/rho(i,j)
            
            endif
        enddo
    enddo
    !$omp end parallel do

    return
end subroutine macro


subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2

    error1 = 0.0d0
    error2 = 0.0d0
    
    !$omp parallel do default(none) shared(u,up,v,vp,obst) private(i,j) reduction(+:error1,error2)
    do j=1,ny
        do i=1,nx
            if(obst(i,j).EQ.0) then
                error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
                error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)
                
                up(i,j) = u(i,j)
                vp(i,j) = v(i,j) 
            endif
        enddo
    enddo
    !$omp end parallel do
    
    errorU = sqrt(error1)/sqrt(error2)

    write(*,*) itc,' ',errorU

    return
    end subroutine check