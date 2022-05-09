subroutine bounceback()
    use commondata
    implicit none
    integer :: i, j


    !Left side (i=1)
    if (coords(0) == 0) then
        do j=1,ny 
            f(1,1,j) = f_post(3,1,j)
            f(5,1,j) = f_post(7,1,j)
            f(8,1,j) = f_post(6,1,j)
        enddo
    endif

    if (coords(0) == dims(0)-1) then
        !Right side (i=nx)
        do j=1,ny
            f(3,nx,j) = f_post(1,nx,j)
            f(6,nx,j) = f_post(8,nx,j)
            f(7,nx,j) = f_post(5,nx,j)
        enddo
    endif


    if (coords(1) == 0) then
        do i=1,nx
            !Bottom side (j=1)
            f(2,i,1) = f_post(4,i,1)
            f(5,i,1) = f_post(7,i,1)
            f(6,i,1) = f_post(8,i,1)
        enddo
    endif

    if (coords(1) == dims(1) - 1) then
        do i=1,nx
            !Top side (j=ny)
            f(4,i,ny) = f_post(2,i,ny)
            f(7,i,ny) = f_post(5,i,ny)-rho(i,ny)*(Uwall)/6.0d0
            f(8,i,ny) = f_post(6,i,ny)-rho(i,ny)*(-Uwall)/6.0d0
        enddo
    endif
    
    return
end subroutine bounceback