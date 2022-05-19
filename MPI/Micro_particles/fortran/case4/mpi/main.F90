program main
    use mpi
    use omp_lib     
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads
    real(8) :: start_time, end_time

    call mpi_starts()

    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------
#ifdef _OPENMP
    call OMP_set_num_threads(24)
    write(*,*) "Start OpenMP......"
    myMaxThreads = OMP_get_max_threads()
    write(*,*) "Running threads =",myMaxThreads
#endif
    !--------------------OpenMP--------------------OpenMP--------------------OpenMP----------

    call allocate_all()

    call initial()
    call output()
    
    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    start_time = MPI_Wtime()

    do while( (errorU.GT.eps).AND.(itc.LT.itc_max) )

        call collision()

        ! call message_passing_sendrecv()
        call send_all_fp()

        call streaming()

        call bounceback()

        ! call message_particle_f_post()

        call bounceback_particle()

        call macro()

        call calForce()

        call send_all_f()    
        ! call message_particle_f()
        
        itc = itc+1

        if(MOD(itc,500).EQ.0) then
            call check()
            call output()
        endif
        !~ if(MOD(itc,5000).EQ.0) then
            !~ call output_Tecplot()
        !~ endif

        call updateCenter()

    enddo

    call MPI_Barrier(MPI_COMM_WORLD, rc)
    end_time = MPI_Wtime()
    
    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif

    if (rank == 0) then
        write(*,*) "Time (CPU) = ",finish-start, "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", finish2-start2, "s"
#endif
        write(*,*) "Time (MPI) = ", real(end_time - start_time), "s"
    endif

    itc = itc+1
    call output()
    call free_all()
    
    call MPI_Finalize(rc)

    stop
end program main