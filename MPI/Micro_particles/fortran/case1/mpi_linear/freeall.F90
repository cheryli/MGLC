subroutine freeall()
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

end subroutine freeall