module utils_io

    implicit none
    
    character(len=10)   :: id_sim_tab(48000)
    integer             :: n_atoms_tab(48000), n_steps_tab(48000)
    double precision    :: composition_tab(48000), initial_temperature_tab(48000)

    contains


    subroutine read_file_list

        use constants, only: file_list
        use variable, only: max_files, xyz_files

        integer :: ierr, i_file, n_files

        n_files = 0
        open(10, file=file_list, status='old', action='read')
        do i_file = 1, size(xyz_files)
            read(10, '(A)', iostat=ierr) xyz_files(i_file)
            if (ierr /= 0) exit
            n_files = n_files + 1
        end do
        close(10)
        if (n_files .eq. 0) stop 'error: no xyz files found'
        if (max_files .gt. 0) then
            max_files = min(max_files, n_files) 
        else
            max_files = n_files
        end if

    endsubroutine


    subroutine load_data

        use variable, only: xyz_files, data_file

        integer :: i_file, ierr

        open(10, file=data_file, status='old', action='read')
        do i_file = 1, size(xyz_files)
            read(10, *, iostat=ierr) id_sim_tab(i_file), n_atoms_tab(i_file), composition_tab(i_file), n_steps_tab(i_file), &
                        initial_temperature_tab(i_file)
            if (ierr /= 0) exit
        enddo
        close(10)

    endsubroutine


    subroutine read_xyz
        
        use variable, only: xyz_files, snapshot_index, pos, species, epot, box
        use descriptor, only: n_atoms

        integer             :: i_atom, index0, index1
        double precision    :: lattice(9)
        character(len=255)  :: line

        open(10, file=xyz_files(snapshot_index), status='old', action='read')
        read(10, *) n_atoms
        read(10, '(A)') line
        do i_atom = 1, n_atoms
            read(10, *) species(i_atom)(1:2), pos(1, i_atom), pos(2, i_atom), pos(3, i_atom), epot(i_atom)
        end do
        close(10)
        index0 = index(line, 'Lattice="') + len('Lattice="')
        index1 = index(line(index0:), '"')
        index1 = index1 + index0 - 2
        read(line(index0:index1), *) lattice
        box = (/ lattice(1), lattice(5), lattice(9) /)
        pos(1, :) = modulo(pos(1, :), box(1))
        pos(2, :) = modulo(pos(2, :), box(2))
        pos(3, :) = modulo(pos(3, :), box(3))
        
    endsubroutine


    subroutine read_data
        
        use variable, only: xyz_files, snapshot_index, found
        use descriptor, only: id_sim, n_atoms, composition, n_steps, initial_temperature

        integer             :: i_file, idx
        character(len=255)  :: id_from_file, filename, basename

        filename = xyz_files(snapshot_index)
        idx = index(filename, '/', back=.true.)
        if (idx .gt. 0) then
            basename = filename(idx+1:)
        else
            basename = filename
        end if
        id_from_file = basename(1:index(basename, '.xyz')-1)
        found = .false.
        do i_file = 1, size(xyz_files)
            if (trim(id_sim_tab(i_file)) .eq. id_from_file) then
                if (n_atoms_tab(i_file) .ne. n_atoms) then
                    write(*,*) "error: n_atoms mismatch for snapshot ", snapshot_index, ": n_atoms in xyz file = ", n_atoms, &
                        ", n_atoms in data file = ", n_atoms_tab(i_file), " (id_sim = ", trim(id_sim_tab(i_file)), ")"
                else
                    found = .true.
                    id_sim = id_sim_tab(i_file)
                    n_atoms = n_atoms_tab(i_file)
                    composition = composition_tab(i_file)
                    n_steps = n_steps_tab(i_file)
                    initial_temperature = initial_temperature_tab(i_file)
                end if
                exit
            end if
        enddo

    endsubroutine


    subroutine read_input

        use variable, only: fs, edge, sc_mrad, vib1, vib2, vibdir, oapr, dose_e_per_a2, readout_noise_e, &
                            doptc, dopsc, dovib, aberr_re, aberr_im
        use random_utils, only: random_uniform, sample_aberration

        doptc = 1
        fs = 0.7d0
        dopsc = 1
        dovib = 2
        edge = 0.0d0

        sc_mrad = random_uniform(0.05d0, 0.25d0)
        vib1 = random_uniform(0.005d0, 0.04d0)
        vib2 = random_uniform(0.005d0, 0.04d0)
        vibdir = random_uniform(0.0d0, 180.0d0)
        oapr = random_uniform(20.0d0, 30.0d0)
        dose_e_per_a2 = random_uniform(3000.0d0, 15000.0d0)
        readout_noise_e = random_uniform(0.0d0, 1.0d0)

        aberr_re = 0.0d0
        aberr_im = 0.0d0
        
        !call sample_aberration(2, -2.0d0, -1.0d0) ! Defocus
        call sample_aberration(2, -20.0d0, -10.0d0) ! Defocus
        call sample_aberration(3, 0.0d0, 6.0d0) ! A1 2-fold astigmatism
        call sample_aberration(4, 0.0d0, 50.0d0) ! B2 Axial coma
        call sample_aberration(5, 0.0d0, 50.0d0) ! A2 3-fold astigmatism
        call sample_aberration(6, -15000.0d0, 0.0d0) ! C3 Spherical aberration (Cs)
        call sample_aberration(7, 0.0d0, 700.0d0) ! S3 Star aberration
        call sample_aberration(8, 0.0d0, 700.0d0) ! A3 4-fold astigmatism
        call sample_aberration(9, 0.0d0, 1500.0d0) ! B4 5th-order term
        call sample_aberration(10, 0.0d0, 1500.0d0) ! D4 5th-order term
        call sample_aberration(11, 0.0d0, 1500.0d0) ! A4 5th-order term

    endsubroutine


    subroutine save_data

        use descriptor
        use variable, only: size, nx, ny

        integer             :: i_file, ierr, i_px
        character(len=10)   :: rank_suffix
        integer             :: pixel_row(nx*ny)

        open(10, file="data.dat", status='replace')
        open(11, file="images.dat", status='replace')
        write(10, '(A)') 'id_sim, n_atoms, n_steps, initial_temperature, epot_total, composition, gyration_radius, '&
                    'nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in, d_com, coreshell_index'
        do i_file = 0, size-1
            write(rank_suffix, '(I0)') i_file
            open(12, file='descriptors_rank_' // trim(adjustl(rank_suffix)) // '.tmp', status='old')
            open(13, file='images_rank_' // trim(adjustl(rank_suffix)) // '.tmp', status='old')
            do 
                read(12, *, iostat=ierr) id_sim_bis, n_atoms, n_steps, initial_temperature, epot_total, composition, & 
                    gyration_radius, nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in, d_com, coreshell_index
                if (ierr.ne.0) exit         
                write(10, *) id_sim_bis, n_atoms, n_steps, initial_temperature, epot_total, composition, gyration_radius, &
                    nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in, d_com, coreshell_index
                
                    read(13, *, iostat=ierr) id_sim_bis, pixel_row
                if (ierr.ne.0) exit
                write(11, '(A)', advance='no') id_sim_bis
                do i_px = 1, nx*ny
                    write(11, '(1X,I0)', advance='no') pixel_row(i_px)
                enddo
                write(11, *)
            enddo
            close(12)
            close(13)
        enddo    
        close(11)
        close(10)

    endsubroutine


    subroutine save_data_row
    
        use constants, only: image_unit, descriptor_unit
        use descriptor
        use variable, only: nx, ny, box, image

        integer             :: i_px, j_px, qimage(nx, ny)
        double precision    :: scaled, image_min, image_max

        image_min = minval(image)
        image_max = maxval(image)
        if (image_max .le. image_min) then
            qimage = 0
            return
        endif
        write(image_unit, '(A)', advance='no') id_sim_bis
        do j_px = 1, ny
            do i_px = 1, nx
                scaled = -128.0d0 + 255.0d0 * (image(i_px, j_px) - image_min) / (image_max - image_min)
                qimage(i_px, j_px) = max(-128, min(127, nint(scaled)))
                write(image_unit, '(1X,I0)', advance='no') qimage(i_px, j_px)
            enddo
        enddo
        write(image_unit, *)
        write(descriptor_unit, *) id_sim_bis, n_atoms, n_steps, initial_temperature, epot_total, composition, gyration_radius, &
            nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in, d_com, box(1)
    
        endsubroutine


endmodule