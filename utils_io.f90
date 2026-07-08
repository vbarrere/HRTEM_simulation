module utils_io

    implicit none

    contains

    subroutine read_file_list

        use constants,  only: file_list
        use variable,   only: max_files, xyz_files

        integer             :: ierr, i_file, n_files

        n_files = 0
        open(10, file=file_list, status='old', action='read')
        do
            read(10, '(A)', iostat=ierr)
            if(ierr.ne.0) exit
            n_files = n_files + 1
        end do
        close(10)
        allocate(xyz_files(n_files))
        open(10, file=file_list, status='old', action='read')
        do i_file = 1, n_files
            read(10, '(A)') xyz_files(i_file)
        end do
        close(10)
        if (n_files .eq. 0) stop 'error: no xyz files found'
        if (max_files .gt. 0) then
            max_files = min(max_files, n_files) 
        else
            max_files = n_files
        end if

    endsubroutine


    subroutine read_xyz
        
        use variable,   only: xyz_files, snapshot_index, pos, species, epot, box
        use descriptor, only: n_atoms

        integer                     :: i_atom, index0, index1
        double precision            :: lattice(9)
        character(len=255)          :: line

        open(10, file=xyz_files(snapshot_index), status='old', action='read')
        read(10, *) n_atoms
        read(10, '(A)') line
        !allocate(species(n_atoms), pos(3, n_atoms), epot(n_atoms))
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
        
        use variable,   only: xyz_files, snapshot_index
        use descriptor, only: id_sim, n_atoms, composition, n_steps, initial_temperature

        integer                     :: slash_index, ext_pos, ierr, n_atoms_tmp, n_steps_tmp
        double precision            :: composition_tmp, initial_temperature_tmp
        logical                     :: found
        character(len=10)           :: id_tmp
        character(len=255)          :: filename, basename, id_from_file

        filename = trim(xyz_files(snapshot_index))
        slash_index = index(filename, '/', back=.true.)
        basename = filename(slash_index+1:)
        ext_pos = index(basename, '.xyz')
        if (ext_pos .gt. 0) then
            id_from_file = basename(1:ext_pos-1)
        else
            id_from_file = basename
        end if

        found = .false.
        open(10, file='../data2.dat', status='old', action='read')
        do
            read(10, *, iostat=ierr) id_tmp, n_atoms_tmp, composition_tmp, n_steps_tmp, initial_temperature_tmp
            if (ierr.ne.0) exit
            if (trim(id_tmp) .eq. trim(id_from_file)) then
                !write(*,*) "id_from_file: ", trim(id_from_file), " id_tmp: ", trim(id_tmp)
                found = .true.
                id_sim = id_tmp
                if (n_atoms_tmp .ne. n_atoms) then
                    write(*,*) "error: n_atoms mismatch for snapshot ", snapshot_index, ": n_atoms in xyz file = ", n_atoms, &
                        ", n_atoms in data file = ", n_atoms_tmp
                    stop
                end if
                !n_atoms = n_atoms_tmp
                composition = composition_tmp
                n_steps = n_steps_tmp
                initial_temperature = initial_temperature_tmp
                exit
            end if
        enddo
        close(10)
        if (.not. found) then
            write(*, *) "error: no matching data found for snapshot ", snapshot_index
            stop
        end if
        !write(*,*) "Data for snapshot ", snapshot_index, ": id_sim=", trim(id_sim), ", n_atoms=", n_atoms, &
        !    ", composition=", composition, ", n_steps=", n_steps, ", initial_temperature=", initial_temperature
    endsubroutine


    subroutine save_xyz
        
        use variable,   only: snapshot_index, augmentation_index, box, pos_cluster, epot, species
        use descriptor, only: n_atoms

        integer                                 :: i_atom
        character(len=256)                      :: filename, snap_str, aug_str

        write(snap_str,*) snapshot_index
        write(aug_str,*) augmentation_index
        filename = '../xyz_process/snapshot_' // trim(adjustl(snap_str)) // '_' // trim(adjustl(aug_str)) // '.xyz'
        open(10, file=filename, action='write', status='replace')
        write(10, *) n_atoms
        write(10, '(A, 9ES22.15, A)') 'Lattice="', box(1), 0.0d0, 0.0d0, 0.0d0, box(2), 0.0d0, 0.0d0, 0.0d0, box(3), &
            '" Properties=species:S:1:pos:R:3:epot:R:1'
        do i_atom = 1, n_atoms
            write(10, *) species(i_atom), pos_cluster(1, i_atom), pos_cluster(2, i_atom), &
                    pos_cluster(3, i_atom), epot(i_atom)
        enddo
        close(10)

    endsubroutine

    subroutine read_input

        use variable,  only:    fs, edge, sc_mrad, vib1, vib2, vibdir, oapr, dose_e_per_a2, readout_noise_e, &
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
        
        call sample_aberration(2, -2.0d0, -1.0d0) ! Defocus
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

endmodule