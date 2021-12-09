subroutine streaming()
    use commondata
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp   
    integer :: alpha
    
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                do alpha=0,18
                    ip = i - ex(alpha)
                    jp = j - ey(alpha)
                    kp = k - ez(alpha)
                    
                    f(alpha,i,j,k) = f_post(alpha,ip,jp,kp)
                enddo
            enddo
        enddo
    enddo
    
    return
end subroutine streaming