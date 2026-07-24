module msa

    implicit none

    integer ::   index_slice

    contains

    subroutine run_msa

        use constants, only: pi
        use variable, only: wave, dx, dy, nz, trans, nx, ny, lambda, dz, wave_fft, gmax
        use fft, only: fft2, fft_index

        integer             ::  i_px, j_px, mx, my
        double precision    ::  aperture2, chi, gx, gy, g2

        wave = dcmplx(1.0d0, 0.0d0)
        !gmax= min(0.5d0/dx, 0.5d0/dy)
        aperture2 = (2.0d0*gmax/3.0d0)**2
        do index_slice = 1, nz
            wave(1:nx, 1:ny) = wave(1:nx, 1:ny) * trans(1:nx, 1:ny, index_slice)
            call fft2(wave(1:nx, 1:ny), wave_fft(1:nx, 1:ny), -1)
            do j_px = 1, ny
                my = fft_index(j_px, ny)
                gy = my / (dy * ny)
                do i_px = 1, nx
                    mx = fft_index(i_px, nx)
                    gx = mx / (dx * nx)
                    g2 = gx**2 + gy**2
                    if (g2 .gt. aperture2) then
                        wave_fft(i_px, j_px) = dcmplx(0.0d0, 0.0d0)
                    else
                        chi = pi * lambda * g2 * dz
                        wave_fft(i_px, j_px) = wave_fft(i_px, j_px) * dcmplx(cos(chi), -sin(chi))
                    endif
                enddo
            enddo
            call fft2(wave_fft(1:nx, 1:ny), wave(1:nx, 1:ny), 1)
            wave(1:nx, 1:ny) = wave(1:nx, 1:ny) / (nx * ny)
        enddo

    endsubroutine

endmodule