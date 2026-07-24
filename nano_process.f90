module nano_process
    
    implicit none
    
    contains

    subroutine extract_largest_cluster
        
        use descriptor, only: n_clusters, n_atoms
        use constants, only: cluster_cutoff
        use variable, only: pos, box, species, epot

        integer             ::  i_atom, j_atom, stack_top, current_atom, largest_cluster, i_cluster, largest_size
        integer             ::  cluster_id(n_atoms), cluster_size(n_atoms), stack(n_atoms)
        double precision    ::  cutoff2, delta(3), pos_tmp(3, n_atoms), epot_tmp(n_atoms)
        logical             ::  visited(n_atoms)
        character(len=2)    ::  species_tmp(n_atoms)

        cutoff2 = cluster_cutoff * cluster_cutoff
        cluster_id = 0
        cluster_size = 0
        visited = .false.
        pos_tmp = 0.0d0
        n_clusters = 0

        do i_atom = 1, n_atoms
            if (visited(i_atom)) cycle
            n_clusters = n_clusters + 1
            stack_top = 1
            stack(1) = i_atom
            visited(i_atom) = .true.
            cluster_id(i_atom) = n_clusters
            pos_tmp(:, i_atom) = pos(:, i_atom)

            do while(stack_top .gt. 0)
                current_atom = stack(stack_top)
                stack_top = stack_top - 1
                cluster_size(n_clusters) = cluster_size(n_clusters) + 1
            
                do j_atom = 1, n_atoms
                    if (visited(j_atom)) cycle
                    delta = pos(:,j_atom) - pos(:,current_atom)
                    delta = delta - box * nint(delta / box)
                    if(dot_product(delta, delta) .gt. cutoff2) cycle
                    stack_top = stack_top + 1
                    stack(stack_top) = j_atom
                    visited(j_atom) = .true.
                    cluster_id(j_atom) = n_clusters
                    pos_tmp(:, j_atom) = pos(:, j_atom)
                enddo
            enddo
        enddo

        largest_cluster = 1
        largest_size = cluster_size(1)
        do i_cluster = 2, n_clusters
            if (cluster_size(i_cluster) .gt. largest_size) then
                largest_cluster = i_cluster
                largest_size = cluster_size(i_cluster)
            endif
        enddo
        epot_tmp = epot(1:n_atoms)
        species_tmp = species(1:n_atoms)
        i_cluster = 0
        pos = 0.0d0
        species = '  '
        epot = 0.0d0
        do i_atom = 1, n_atoms
            if (cluster_id(i_atom) .ne. largest_cluster) cycle
            i_cluster = i_cluster + 1
            species(i_cluster) = species_tmp(i_atom)
            pos(:, i_cluster) = pos_tmp(:, i_atom)
            epot(i_cluster) = epot_tmp(i_atom)
        enddo
        n_atoms = largest_size
        
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) * nint((pos(1, 1:n_atoms) - pos(1, 1) )/ box(1))
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) * nint((pos(2, 1:n_atoms) - pos(2, 1) )/ box(2))
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) * nint((pos(3, 1:n_atoms) - pos(3, 1) )/ box(3))
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - sum(pos(1, 1:n_atoms)) / n_atoms
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - sum(pos(2, 1:n_atoms)) / n_atoms
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - sum(pos(3, 1:n_atoms)) / n_atoms
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) * nint(pos(1, 1:n_atoms) / box(1))
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) * nint(pos(2, 1:n_atoms) / box(2))
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) * nint(pos(3, 1:n_atoms) / box(3))
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - sum(pos(1, 1:n_atoms)) / n_atoms
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - sum(pos(2, 1:n_atoms)) / n_atoms
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - sum(pos(3, 1:n_atoms)) / n_atoms
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) / 2.0d0
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) / 2.0d0
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) / 2.0d0

    endsubroutine
    

    subroutine compute_descriptors
        
        use, intrinsic :: ieee_arithmetic
        use descriptor
        use variable, only: species, pos, atom_typ1
        
        integer             ::  i_atom
        double precision    ::  dist2, radius_limit, term1, term2, term3, radius_limit2
        
        nat1 = 0
        nat2 = 0
        nat1_out = 0
        nat1_in = 0
        nat2_out = 0
        nat2_in = 0

        mass_center = 0.0d0
        mass_center_nat1 = 0.0d0
        mass_center_nat2 = 0.0d0

        do i_atom = 1, n_atoms
            mass_center = mass_center + pos(:, i_atom)
            if (species(i_atom) .eq. atom_typ1) then
                nat1 = nat1 + 1
                mass_center_nat1 = mass_center_nat1 + pos(:, i_atom)
            else
                nat2 = nat2 + 1
                mass_center_nat2 = mass_center_nat2 + pos(:, i_atom)
            end if
        enddo
        mass_center = mass_center / dble(n_atoms)
        mass_center_nat1 = mass_center_nat1 / dble(nat1)
        mass_center_nat2 = mass_center_nat2 / dble(nat2)

        gyration_radius = 0.0d0
        do i_atom = 1, n_atoms
            gyration_radius = gyration_radius + dot_product(pos(:, i_atom) - mass_center, &
                    pos(:, i_atom) - mass_center)
        enddo
        gyration_radius = sqrt(gyration_radius / dble(n_atoms))
        radius_limit = 0.8d0 * gyration_radius
        radius_limit2 = radius_limit * radius_limit
        do i_atom = 1, n_atoms
            dist2 = dot_product(pos(:, i_atom) - mass_center, pos(:, i_atom) - mass_center)
            if(dist2 .ge. radius_limit2) then
                if(species(i_atom) .eq. atom_typ1) then
                    nat1_out = nat1_out + 1
                else
                    nat2_out = nat2_out + 1
                endif
            else
                if(species(i_atom) .eq. atom_typ1) then
                    nat1_in = nat1_in + 1
                else
                    nat2_in = nat2_in + 1
                endif
            endif
        enddo
        d_com = sqrt(dot_product(mass_center_nat1 - mass_center_nat2, mass_center_nat1 - mass_center_nat2))
        if(gyration_radius .gt. 0.0d0) then
            term1 = dble(nat1_out) / dble(nat1_out + nat2_out) - composition
            term2 = dble(nat1_in) / dble(nat1_in + nat2_in) - composition
            term3 = d_com / (2.0d0 * gyration_radius)
            coreshell_index = 2.0d0 * abs(term1) + 2.0d0 * abs(term2) - term3
        else
            coreshell_index = ieee_value(coreshell_index, ieee_quiet_nan)
        endif

    endsubroutine


    subroutine compute_rotation
    
        use descriptor, only: n_atoms, id_sim, id_sim_bis, mass_center
        use variable, only: pos, pos_cluster, augmentation_index
        use constants, only: pi

        integer             ::  i_atom
        double precision    ::  u, position_scale, phi, theta, rot_matrix(3, 3)
        character(len=10)   ::  str_augmentation_index

        do i_atom = 1, n_atoms
            pos_cluster(:, i_atom) = pos(:, i_atom) - mass_center
        enddo
        call random_number(u)
        position_scale = 0.95d0 + 0.1d0 * u
        pos_cluster = pos_cluster * position_scale
        call random_number(phi)
        call random_number(theta)
        phi = 2.0d0 * pi * phi
        theta = pi * theta

        rot_matrix(1, 1) = cos(phi) * cos(theta)
        rot_matrix(1, 2) = -sin(phi)
        rot_matrix(1, 3) = cos(phi) * sin(theta)
        rot_matrix(2, 1) = sin(phi) * cos(theta)
        rot_matrix(2, 2) = cos(phi)
        rot_matrix(2, 3) = sin(phi) * sin(theta)
        rot_matrix(3, 1) = -sin(theta)
        rot_matrix(3, 2) = 0.0d0
        rot_matrix(3, 3) = cos(theta)

        do i_atom = 1, n_atoms
            pos_cluster(:, i_atom) = matmul(rot_matrix, pos_cluster(:, i_atom))
        enddo
        write(str_augmentation_index, '(I0)') augmentation_index
        id_sim_bis = trim(id_sim) // "_" // trim(adjustl(str_augmentation_index))

    endsubroutine


    subroutine place_in_box
    
        use constants, only: lateral_border_margin, box_hrtem
        use variable, only: placed, shift, pos_cluster
        use descriptor, only: n_atoms

        integer             ::  i_axis
        double precision    ::  min_coord(3), max_coord(3), span(3), free_lateral(2), u(2), box_size(3)

        placed = .false.
        shift = 0.0d0
        box_size = box_hrtem * 10.0d0
        do i_axis = 1, 3
            min_coord(i_axis) = minval(pos_cluster(i_axis, 1:n_atoms))
            max_coord(i_axis) = maxval(pos_cluster(i_axis, 1:n_atoms))
        enddo
        span = max_coord - min_coord
        if (any(span .ge. box_size - 1.0d-8)) return
        free_lateral = box_size(1:2) - span(1:2) - 2.0d0 * lateral_border_margin
        if (any(free_lateral .lt. 0.0d0)) return

        call random_number(u)
        shift(1:2) = -min_coord(1:2) + lateral_border_margin + u * free_lateral
        shift(3) = -min_coord(3) + 0.5d0 * (box_size(3) - span(3))
        pos_cluster(:, 1:n_atoms) = pos_cluster(:, 1:n_atoms) + spread(shift, dim=2, ncopies=n_atoms)
        placed = .true.
    
    endsubroutine
    

    subroutine prepare_hrtem_particle
        
        use constants, only: dbf_ag, dbf_co, box_hrtem
        use descriptor, only: n_atoms
        use variable, only: species, atom_typ1, pos_cluster, atomic_number, biso

        integer ::  i_atom

        do i_atom = 1, n_atoms
            if (species(i_atom) .eq. atom_typ1) then
                atomic_number(i_atom) = 47
                biso(i_atom) = dbf_ag
            else
                atomic_number(i_atom) = 27
                biso(i_atom) = dbf_co
            endif
            pos_cluster(:, i_atom) = pos_cluster(:, i_atom) / (box_hrtem*10.0d0)
        enddo

    endsubroutine
    
end module