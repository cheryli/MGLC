subroutine allocate_all()
    use commondata

    allocate(X(total_nx))
    allocate(Y(total_ny))
    
    allocate(u(nx,ny))
    allocate(v(nx,ny))
    allocate(rho(nx,ny))
    allocate(up(nx,ny))
    allocate(vp(nx,ny))

    allocate(obst(0:nx+1, 0:ny+1))
    allocate(obstNew(0:nx+1, 0:ny+1))

    allocate(f(0:8, -2:nx+3, -2:ny+3))
    allocate(f_post(0:8, -1:nx+2, -1:ny+2))

end subroutine allocate_all    



subroutine free_all()
    use commondata
    implicit none

    deallocate(X)
    deallocate(Y)

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)

    deallocate(f)
    deallocate(f_post)

    deallocate(obst)
    deallocate(obstNew)

end subroutine free_all