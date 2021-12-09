subroutine macro()
    use commondata
    implicit none
    integer :: i, j, k, alpha

    rho = 0.0d0
    u = 0.0d0
    v = 0.0d0
    w = 0.0d0
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha = 0, 18
                    rho(i,j,k) = rho(i,j,k) + f(alpha,i,j,k)
                    u(i,j,k) = u(i,j,k) + f(alpha,i,j,k) * ex(alpha)
                    v(i,j,k) = v(i,j,k) + f(alpha,i,j,k) * ey(alpha)
                    w(i,j,k) = w(i,j,k) + f(alpha,i,j,k) * ez(alpha)
                enddo

                u(i,j,k) = u(i,j,k)/rho(i,j,k)
                v(i,j,k) = v(i,j,k)/rho(i,j,k)
                w(i,j,k) = w(i,j,k)/rho(i,j,k)
            enddo
        enddo
    enddo

    return
end subroutine macro