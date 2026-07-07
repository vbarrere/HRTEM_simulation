module image

    implicit none

    contains

    subroutine quantize_image_int8

        use variable, only: image, image_min, image_max, qimage, nx, ny

        integer :: i_px, j_px
        double precision :: scaled

        image_min = minval(image)
        image_max = maxval(image)
        if (image_max .le. image_min) then
            qimage = 0
            return
        endif
        do j_px = 1, ny
            do i_px = 1, nx
                scaled = -128.0d0 + 255.0d0 * (image(i_px, j_px) - image_min) / (image_max - image_min)
                qimage(i_px, j_px) = max(-128, min(127, nint(scaled)))
            enddo
        enddo

    endsubroutine

    subroutine write_image_row

        use variable, only: qimage, nx, ny

        integer             :: i_px, j_px

        do j_px = 1, ny
            do i_px = 1, nx
                if(i_px.eq.1 .and. j_px.eq.1) then
                    write(11, '(I0)', advance='no') qimage(i_px, j_px)
                else
                    write(11, '(1X,I0)', advance='no') qimage(i_px, j_px)
                endif
            enddo
        enddo
        write(11, *)

    endsubroutine

    
    subroutine write_descriptor_row
        
        use descriptor
        

        write(12,*) id_sim, n_atoms, n_steps, initial_temperature, epot_total, composition, gyration_radius, nat1, nat2, nat1_out, &
                    nat2_out, nat1_in, nat2_in, d_com

        !write(12, '(I0,1X,ES16.8E3,1X,ES16.8E3,1X,I0,1X,I0,1X,A,1X,I0)', advance='no') &
        !    n_atoms, composition, coreshell_index, snapshot_index, augmentation_index, &
        !    trim(input_xyz), seed
        !write(12, '(1X,I0,1X,I0,1X,I0,1X,ES16.8E3,1X,ES16.8E3,1X,I0,1X,I0,1X,I0,1X,I0)', advance='no') &
        !    nat1, nat2, n_clusters, gyration_radius, d_com, &
        !    descriptors%n_shell, descriptors%n_ag_shell, descriptors%n_core, descriptors%n_ag_core
        !write(12, '(1X,3(ES16.8E3,1X),3(I0,1X),12(ES16.8E3,1X),3(I0,1X))', advance='no') &
        !    position_scale, image_min, image_max, nx, ny, nz, ht, fs, sc_mrad, vib1, vib2, vibdir, oapr, edge, &
        !    dose_e_per_a2, readout_noise_e, dbf_ag, dbf_co, doptc, dopsc, dovib
        !write(12, '(9(ES16.8E3,1X),3(ES16.8E3,1X))', advance='no') rot, shift
       ! write(12, '(24(ES16.8E3,1X),24(ES16.8E3,1X))') aberr_re, aberr_im
    
    endsubroutine

    subroutine append_file_to_unit(input_file, output_unit)

        integer, intent(in)          :: output_unit
        integer                      :: error
        character(len=*), intent(in) :: input_file
        character(len=131072)        :: line

        open(22, file=input_file, action='read', status='old', iostat=error)
        do
            read(22, '(A)', iostat=error) line
            if(error.ne.0) exit
            write(output_unit, '(A)') trim(line)
        enddo
        close(22)
        
    
    endsubroutine append_file_to_unit


    subroutine write_descriptor_header(unit_no)

        integer, intent(in) :: unit_no

        write(unit_no, '(A)', advance='no') 'nat xag coreshell_index snapshot_index augmentation_index xyz_file seed'
        write(unit_no, '(A)', advance='no') ' n_ag n_co n_clusters rg_ang d_com_ang n_shell n_ag_shell n_core n_ag_core'
        write(unit_no, '(A)', advance='no') ' position_scale image_min image_max nx ny nz ht fs sc_mrad vib1 vib2 vibdir oapr edge'
        write(unit_no, '(A)', advance='no') ' dose_e_per_a2 readout_noise_e dbf_ag dbf_co doptc dopsc dovib'
        write(unit_no, '(A)', advance='no') ' rot11 rot21 rot31 rot12 rot22 rot32 rot13 rot23 rot33'
        write(unit_no, '(A)', advance='no') ' shift_x_ang shift_y_ang shift_z_ang'
        write(unit_no, '(A)', advance='no') ' aberr_re01 aberr_re02 aberr_re03 aberr_re04 aberr_re05 aberr_re06'
        write(unit_no, '(A)', advance='no') ' aberr_re07 aberr_re08 aberr_re09 aberr_re10 aberr_re11 aberr_re12'
        write(unit_no, '(A)', advance='no') ' aberr_re13 aberr_re14 aberr_re15 aberr_re16 aberr_re17 aberr_re18'
        write(unit_no, '(A)', advance='no') ' aberr_re19 aberr_re20 aberr_re21 aberr_re22 aberr_re23 aberr_re24'
        write(unit_no, '(A)', advance='no') ' aberr_im01 aberr_im02 aberr_im03 aberr_im04 aberr_im05 aberr_im06'
        write(unit_no, '(A)', advance='no') ' aberr_im07 aberr_im08 aberr_im09 aberr_im10 aberr_im11 aberr_im12'
        write(unit_no, '(A)', advance='no') ' aberr_im13 aberr_im14 aberr_im15 aberr_im16 aberr_im17 aberr_im18'
        write(unit_no, '(A)', advance='no') ' aberr_im19 aberr_im20 aberr_im21 aberr_im22 aberr_im23 aberr_im24'
        write(unit_no, *)

    endsubroutine write_descriptor_header

endmodule