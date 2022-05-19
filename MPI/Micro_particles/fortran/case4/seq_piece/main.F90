program main
    use omp_lib     
    use commondata
    implicit none
    real(8) :: start, finish
    real(8) :: start2, finish2
    integer :: myMaxThreads
    
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
    call output_Tecplot()
    
    call CPU_TIME(start)
#ifdef _OPENMP
    start2 = OMP_get_wtime()
#endif

    do while( (errorU.GT.eps).AND.(itc.LT.itc_max) )

        call collision()

        call streaming()

        call macro()

        call calForce()
        
        itc = itc+1

        if(MOD(itc,500).EQ.0) then
            call check()
        endif
        
        if(MOD(itc,500).EQ.0) then
            call output_Tecplot()
        endif

        call updateCenter()

    enddo

    call CPU_TIME(finish)
#ifdef _OPENMP
    finish2 = OMP_get_wtime()
#endif

    write(*,*) "Time (CPU) = ",finish-start, "s"
#ifdef _OPENMP
    write(*,*) "Time (OMP) = ", finish2-start2, "s"
#endif

    itc = itc+1
    call output_Tecplot()

    call free_all()

    stop
end program main