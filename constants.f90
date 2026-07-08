module constants

    implicit none

    character(len=255), parameter           :: file_list = "xyz_file_list.tmp"
    integer, parameter                      :: image_unit = 11, descriptor_unit = 12
    double precision, parameter             :: interatomic_distance = 2.889d0
    double precision, parameter             :: cluster_cutoff = 1.5d0 * interatomic_distance
    double precision, parameter             :: pi = acos(-1.0d0)
    double precision, parameter             :: lateral_border_margin = 2.0d0
    double precision, parameter             :: dbf_ag = 0.019d0, dbf_co = 0.015d0
    double precision, parameter             :: hc = 1.2398419843320026d0 ! Planck's constant times speed of light
    double precision, parameter             :: e0 = 510.99895069d0, sigma0 = 2.0886573497d0
    double precision, parameter             :: r8pi2 = 1.0d0/(8.0d0*pi*pi)
    double precision, parameter             :: euler = 0.5772156649015328606065120d0
    double precision, parameter             :: v0 = 0.03809982119d0
    integer, parameter                      :: base_seed = 1357911
    integer, parameter                      :: aberr_periodicity(24) = (/ &
                                                & 1, 0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, &
                                                & 2, 4, 6, 1, 3, 5, 7, 0, 2, 4, 6, 8 /)
    integer, parameter                      :: nx_max = 96, ny_max = 96, nz_max = 15
    integer, parameter                      :: n_atoms_max = 5000
    integer, parameter                      :: n_types_max = 5
    double precision, parameter             :: box_hrtem(3) = (/ 12.0d0, 12.0d0, 12.0d0 /)
    
endmodule


module variable

    use constants, only: n_atoms_max, n_types_max, nx_max, ny_max, nz_max

    implicit none

    integer                                 :: max_files, snapshot_index, augmentation_index, seed
    integer                                 :: n_accepted, n_seed
    double precision                        :: box(3), shift(3)
    character(len=2)                        :: atom_typ1, atom_typ2
    logical                                 :: accepted, placed
    double precision                        :: pos(3, n_atoms_max), epot(n_atoms_max), pos_cluster(3, n_atoms_max)
    character(len=2)                        :: species(n_atoms_max)
    character(len=255), allocatable         :: xyz_files(:)
    integer                                 :: nx, ny, nz, doptc, dopsc, dovib
    double precision                        :: ht, fs, edge, sc_mrad, vib1, vib2, vibdir, oapr, dose_e_per_a2
    double precision                        :: readout_noise_e, aberr_re(24), aberr_im(24)
    character(len=255)                      :: img_file
    integer                                 :: atomic_number(n_atoms_max), atomic_number_type(n_types_max)
    double precision                        :: biso(n_atoms_max), biso_type(n_types_max), lambda, sigma_lambda, dx, dy, dz, &
                                                volume_slc, g2, gx(nx_max), gy(ny_max), f_re, f_im, gthr2
    integer                                 :: type_index(n_atoms_max), n_types, i_type, index_slice, slice_count(nz_max)
    double complex                          :: ftab(n_types_max, nx_max, ny_max), trans(nx_max, ny_max, nz_max), &
                                                uhat(nx_max, ny_max), trans_fft(nx_max, ny_max), wave(nx_max, ny_max), &
                                                wave_fft(nx_max, ny_max), image_c(nx_max, ny_max)

    double precision                        :: a(2), b(6), image_min, image_max
    double precision                        :: image(nx_max, ny_max)
    double precision                        :: dchix, dchiy, wx, wy
    integer                                 :: qimage(nx_max, ny_max)

endmodule





module descriptor

    implicit none

    character(len=10)                       :: id_sim               ! Simulation ID
    integer                                 :: n_atoms              ! Number of atoms
    double precision                        :: composition          ! Composition
    integer                                 :: n_steps              ! Number of steps
    double precision                        :: initial_temperature  ! Initial temperature
    double precision                        :: epot_total           ! Total potential energy
    double precision                        :: gyration_radius      ! Gyration radius
    integer                                 :: nat1                 ! Number of atoms of type 1
    integer                                 :: nat2                 ! Number of atoms of type 2
    integer                                 :: nat1_out             ! Number of atoms of type 1 outside the virtual sphere
    integer                                 :: nat2_out             ! Number of atoms of type 2 outside the virtual sphere
    integer                                 :: nat1_in              ! Number of atoms of type 1 inside the virtual sphere
    integer                                 :: nat2_in              ! Number of atoms of type 2 inside the virtual sphere
    integer                                 :: n_clusters           ! Number of clusters
    double precision                        :: d_com                ! Distance between centers of mass of the two clusters
    double precision                        :: coreshell_index      ! Core-shell index

endmodule