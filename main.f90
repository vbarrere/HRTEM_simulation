program main

    use mpi
    use constants, only: file_list, image_unit, descriptor_unit
    use variable, only: max_files, snapshot_index, species,augmentation_index, accepted, atom_typ1, atom_typ2, n_accepted, placed, &
                        box, nx, ny, nz, ht
    use utils_io
    use nano_process
    use random_utils
    use slc, only: run_slc
    use msa
    use wavimg, only: run_wavimg
    use image

    implicit none

    integer             :: ierr, rank, size, first_file, last_file, accepted_local, accepted_total, i_file
    character(len=255)  :: rank_suffix, rank_images, rank_descriptors
    character(len=255)  :: xyz_dir, images_data, descriptors_data, env_var

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, size, ierr)

    call get_environment_variable('xyz_dir', xyz_dir)
    call get_environment_variable('images_data', images_data)
    call get_environment_variable('descriptors_data', descriptors_data)
    call get_environment_variable('max_files', env_var)
    read(env_var, *) max_files
    call get_environment_variable('atom_typ1', atom_typ1)
    call get_environment_variable('atom_typ2', atom_typ2)
    call get_environment_variable('n_px', env_var)
    read(env_var, *) nx
    read(env_var, *) ny
    call get_environment_variable('nz', env_var)
    read(env_var, *) nz
    call get_environment_variable('ht', env_var)
    read(env_var, *) ht

    if(rank.eq.0) call execute_command_line('find ' // xyz_dir // & 
            ' -maxdepth 1 -name "*.xyz" | sort -V > ' // file_list)
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call read_file_list
    first_file = rank * max_files / size + 1
    last_file  = (rank + 1) * max_files / size
    write(rank_suffix, *) rank
    rank_images = 'images_rank_' // trim(adjustl(rank_suffix)) // '.tmp'
    rank_descriptors = 'descriptors_rank_' // trim(adjustl(rank_suffix)) // '.tmp'

    open(image_unit, file=rank_images, action='write', status='replace')
    open(descriptor_unit, file=rank_descriptors, action='write', status='replace')

    accepted_local = 0
    do snapshot_index = first_file, last_file
        n_accepted = 0
        call read_xyz
        call read_data
        call extract_largest_cluster
        if (count(species .eq. atom_typ1) .eq. 0 .or. count(species .ne. atom_typ1) .eq. 0) cycle
        call compute_descriptors
        do augmentation_index = 1, 10
            call stable_seed
            accepted = .false.
            call seed_random
            call compute_rotation
            call place_in_box
            if (.not. placed) cycle
            !call save_xyz
            call read_input
            box = box * 0.1d0
            call prepare_hrtem_particle
            
            call run_slc
            call run_msa
            call run_wavimg

            call quantize_image_int8
            call write_image_row
            call write_descriptor_row
            box = box * 10.0d0
            accepted = .true.
            if (accepted) n_accepted = n_accepted + 1
        enddo
        accepted_local = accepted_local + n_accepted
    enddo

    close(image_unit)
    close(descriptor_unit)

    call mpi_reduce(accepted_local, accepted_total, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)

    if(rank.eq.0) then
        open(image_unit, file=images_data, action='write', status='replace')
        open(descriptor_unit, file=descriptors_data, action='write', status='replace')
        call write_descriptor_header(descriptor_unit)
        do i_file = 0, size - 1
            write(rank_suffix, *) i_file
            call append_file_to_unit('images_rank_' // trim(adjustl(rank_suffix)) // '.tmp', image_unit)
            call append_file_to_unit('descriptors_rank_' // trim(adjustl(rank_suffix)) // '.tmp', descriptor_unit)
        enddo
        close(image_unit)
        close(descriptor_unit)
        write(*, '(A,1X,I0,1X,A,1X,I0)') 'prep: ', accepted_total, 'accepted images from', max_files, 'snapshots'
        write(*, '(A,1X,A)') 'prep: wrote images to', images_data
        write(*, '(A,1X,A)') 'prep: wrote descriptors to', descriptors_data
    endif

    call mpi_barrier(mpi_comm_world, ierr)
    call execute_command_line('rm -f ' // trim(rank_images) // ' ' // trim(rank_descriptors), exitstat=ierr)
    if(rank.eq.0) call execute_command_line('rm -f ' // file_list, exitstat=ierr)


    call mpi_finalize(ierr)

endprogram