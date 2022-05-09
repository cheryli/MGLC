subroutine g_collision_with_message_exchange()
    use mpi
    use commondata
    implicit none

    !!! caculate the boundary first 
    call collisionT(1, 1, 1, ny, 1, nz)
    call collisionT(nx, nx, 1, ny, 1, nz)
    call collisionT(1, nx, 1, 1, 1, nz)
    call collisionT(1, nx, ny, ny, 1, nz)
    call collisionT(1, nx, 1, ny, 1, 1)
    call collisionT(1, nx, 1, ny, nz, nz)

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)s
    call MPI_Isend(g_post(1, nx, 1, 1), 1, g_surface_x, nbr_surface(1), 1, comm3d, g_req(1), rc)
    call MPI_Irecv(g_post(1, 0, 1, 1), 1, g_surface_x, nbr_surface(2), 1, comm3d, g_req(2), rc)

    ! message passing to (i--)
    call MPI_Isend(g_post(2, 1, 1, 1), 1, g_surface_x, nbr_surface(2), 2, comm3d, g_req(3), rc)
    call MPI_Irecv(g_post(2, nx+1, 1, 1), 1, g_surface_x, nbr_surface(1), 2, comm3d, g_req(4), rc)


    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    call MPI_Isend(g_post(3, 1, ny, 1), 1, g_surface_y, nbr_surface(3), 3, comm3d, g_req(5), rc)
    call MPI_Irecv(g_post(3, 1, 0, 1), 1, g_surface_y, nbr_surface(4), 3, comm3d, g_req(6), rc)

    ! message passing to (j--)
    call MPI_Isend(g_post(4, 1, 1, 1), 1, g_surface_y, nbr_surface(4), 4, comm3d, g_req(7), rc)
    call MPI_Irecv(g_post(4, 1, ny+1, 1), 1, g_surface_y, nbr_surface(3), 4, comm3d, g_req(8), rc)

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    call MPI_Isend(g_post(5, 1, 1, nz), 1, g_surface_z, nbr_surface(5), 5, comm3d, g_req(9), rc)
    call MPI_Irecv(g_post(5, 1, 1, 0), 1, g_surface_z, nbr_surface(6), 5, comm3d, g_req(10), rc)

    ! message passing to (k--)
    call MPI_Isend(g_post(6, 1, 1, 1), 1, g_surface_z, nbr_surface(6), 6, comm3d, g_req(11), rc)
    call MPI_Irecv(g_post(6, 1, 1, nz+1), 1, g_surface_z, nbr_surface(5), 6, comm3d, g_req(12), rc)


    ! don't need to echange line data and point data

    !!! then compute the inner points
    call collisionT(2, nx-1, 2, ny-1, 2, nz-1)
    
    ! wait the message to finish
    call MPI_Waitall(12, g_req, MPI_STATUSES_IGNORE, rc)

end subroutine g_collision_with_message_exchange




subroutine g_message_passing_sendrecv()
    use mpi
    use commondata
    implicit none

    ! message tag --- discrete velocty

    ! /* ----------------------------- surface data ----------------------------- */
    ! ------------ exchange message along x ----------------
    ! message passing to (i++)
    call MPI_Sendrecv(g_post(1, nx, 1, 1), 1, g_surface_x, nbr_surface(1), 1, &
                g_post(1, 0, 1, 1), 1, g_surface_x, nbr_surface(2), 1, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! message passing to (i--)
    call MPI_Sendrecv(g_post(2, 1, 1, 1), 1, g_surface_x, nbr_surface(2), 2, &
                g_post(2, nx+1, 1, 1), 1, g_surface_x, nbr_surface(1), 2, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along y ----------------
    ! message passing to (j++)
    call MPI_Sendrecv(g_post(3, 1, ny, 1), 1, g_surface_y, nbr_surface(3), 3, &
                g_post(3, 1, 0, 1), 1, g_surface_y, nbr_surface(4), 3, &
                comm3d, MPI_STATUS_IGNORE, rc)


    ! message passing to (j--)
    call MPI_Sendrecv(g_post(4, 1, 1, 1), 1, g_surface_y, nbr_surface(4), 4, &
                g_post(4, 1, ny+1, 1), 1, g_surface_y, nbr_surface(3), 4, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! ------------ exchange message along z ----------------
    ! message passing to (k++)
    call MPI_Sendrecv(g_post(5, 1, 1, nz), 1, g_surface_z, nbr_surface(5), 5, &
                g_post(5, 1, 1, 0), 1, g_surface_z, nbr_surface(6), 5, &
                comm3d, MPI_STATUS_IGNORE, rc)


    ! message passing to (k--)
    call MPI_Sendrecv(g_post(6, 1, 1, 1), 1, g_surface_z, nbr_surface(6), 6, &
                g_post(6, 1, 1, nz+1), 1, g_surface_z, nbr_surface(5), 6, &
                comm3d, MPI_STATUS_IGNORE, rc)

    ! don't need to echange line data and point data


end subroutine g_message_passing_sendrecv



subroutine collisionT(i_start, i_end, j_start, j_end, k_start, k_end)
    use commondata
    implicit none
    integer, intent(in) :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: i, j, k
    integer :: alpha
    !------------------------
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: q(0:6)

    !$omp parallel do default(none) shared(g,g_post,u,v,w,T) private(i,j,k,alpha,n,neq,q,n_post) 
    do k = k_start, k_end
        do j = j_start, j_end
            do i = i_start, i_end
            
    n(0) = g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(1) = g(1,i,j,k)-g(2,i,j,k)
    n(2) = g(3,i,j,k)-g(4,i,j,k)
    n(3) = g(5,i,j,k)-g(6,i,j,k)
    n(4) = -6.0d0*g(0,i,j,k)+g(1,i,j,k)+g(2,i,j,k)+g(3,i,j,k)+g(4,i,j,k)+g(5,i,j,k)+g(6,i,j,k)
    n(5) = 2.0d0*g(1,i,j,k)+2.0d0*g(2,i,j,k)-g(3,i,j,k)-g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
    n(6) = g(3,i,j,k)+g(4,i,j,k)-g(5,i,j,k)-g(6,i,j,k)
        
            neq(0) = T(i,j,k)
            neq(1) = T(i,j,k)*u(i,j,k)
            neq(2) = T(i,j,k)*v(i,j,k)
            neq(3) = T(i,j,k)*w(i,j,k)
            neq(4) = T(i,j,k)*paraA
            neq(5) = 0.0d0
            neq(6) = 0.0d0
        
            q(0) = 0.0d0
            q(1) = Qd
            q(2) = Qd
            q(3) = Qd
            q(4) = Qnu
            q(5) = Qnu
            q(6) = Qnu
        
            do alpha=0,6
                n_post(alpha) = n(alpha)-q(alpha)*(n(alpha)-neq(alpha))
            enddo
        
    g_post(0,i,j,k) = n_post(0)/7.0d0-n_post(4)/7.0d0
    g_post(1,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(2,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(1)+n_post(4)/42.0d0+n_post(5)/6.0d0
    g_post(3,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(4,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(2)+n_post(4)/42.0d0-n_post(5)/12.0d0+0.25d0*n_post(6)
    g_post(5,i,j,k) = n_post(0)/7.0d0+0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
    g_post(6,i,j,k) = n_post(0)/7.0d0-0.5d0*n_post(3)+n_post(4)/42.0d0-n_post(5)/12.0d0-0.25d0*n_post(6)
        
            enddo
        enddo
    enddo
    !$omp end parallel do
    
    return
end subroutine collisionT