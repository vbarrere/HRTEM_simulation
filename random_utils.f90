    module random_utils

    implicit none

    contains

    subroutine stable_seed

        use constants, only: base_seed
        use variable, only: xyz_files, snapshot_index, augmentation_index, seed

        character(len=255)  ::  filename
        integer             ::  i
        integer(kind=8)     ::  h

        filename = trim(xyz_files(snapshot_index))
        i = len_trim(filename)
        do while(i .gt. 0 .and. filename(i:i) .ne. '/')
            i = i - 1
        enddo
        h = int(base_seed, kind=8)
        do i = i + 1, len_trim(filename)
            h = modulo(33_8*h + int(iachar(filename(i:i)), kind=8), 2147483000_8)
        enddo
        h = modulo(33_8*h + int(augmentation_index, kind=8), 2147483000_8)
        seed = int(modulo(h, 2147483000_8))
        if (seed .le. 0) seed = seed + 2147483000
    
    endsubroutine

    
    subroutine seed_random

        use variable, only: seed

        integer                 ::  i, n_seed, clock_count
        integer, allocatable    ::  seed_values(:)
        integer(kind=8)         ::  state

        call random_seed(size=n_seed)
        allocate(seed_values(n_seed))

        if (seed .ge. 0) then
            state = int(seed, kind=8)
            do i= 1, n_seed
                state = modulo(1103515245_8 * state + 12345_8 + 1009_8 * int(i, kind=8), 2147483000_8)
                seed_values(i) = int(state)
            enddo
        else
            call system_clock(count=clock_count)
            state = int(clock_count, kind=8)
            do i= 1, n_seed
                state = modulo(1103515245_8 * state + 12345_8 + 1009_8 * int(i, kind=8), 2147483000_8)
                seed_values(i) = int(state)
            enddo
        endif
        call random_seed(put=seed_values)
        deallocate(seed_values)
    endsubroutine


    double precision function random_uniform(xmin, xmax)

        double precision, intent(in)    ::  xmin, xmax
        double precision                ::  u

        call random_number(u)
        random_uniform = xmin + (xmax - xmin) * u
    
    endfunction

    subroutine sample_aberration(index, xmin, xmax)
        
        use constants, only: aberr_periodicity, pi
        use variable, only: aberr_re, aberr_im

        integer                         ::  n
        double precision, intent(in)    ::  xmin, xmax
        integer, intent(in)             ::  index
        double precision                ::  magnitude, phi

        n = aberr_periodicity(index)
        magnitude = random_uniform(xmin, xmax)
        if (n .eq. 0) then
            aberr_re(index) = magnitude
            aberr_im(index) = 0.0d0
            return
        endif
        phi = random_uniform(0.0d0, 2.0d0 * pi)
        aberr_re(index) = magnitude * cos(dble(n) * phi)
        aberr_im(index) = magnitude * sin(dble(n) * phi)
        
    endsubroutine


endmodule