module fft

    use, intrinsic :: iso_c_binding

    implicit none

    interface
        
        function fftw_plan_dft_2d(n0, n1, in, out, sign, flags) bind(C, name='fftw_plan_dft_2d')
        
            use, intrinsic :: iso_c_binding
            integer(c_int), value             :: n0, n1, sign
            integer(c_int), value             :: flags
            complex(c_double_complex)         :: in(*), out(*)
            type(c_ptr)                       :: fftw_plan_dft_2d
        
        endfunction fftw_plan_dft_2d

        subroutine fftw_execute_dft(plan, in, out) bind(C, name='fftw_execute_dft')
        
            use, intrinsic :: iso_c_binding
            type(c_ptr), value                :: plan
            complex(c_double_complex)         :: in(*), out(*)
        
        endsubroutine fftw_execute_dft
    
    endinterface

    integer(c_int), parameter               :: fftw_estimate = 64, fftw_unaligned = 2
    type(c_ptr)                             :: fft_plan(2) = c_null_ptr

    contains

    subroutine fft2(fft_in, fft_out, fft_sign)

        use variable, only: nx, ny


        integer, intent(in)             :: fft_sign
        integer                         :: iplan
        double complex, intent(in)      :: fft_in(nx, ny)
        double complex, intent(out)     :: fft_out(nx, ny)

        iplan = (fft_sign + 3) / 2
        if(.not.c_associated(fft_plan(iplan))) then
            fft_plan(iplan) = fftw_plan_dft_2d(ny, nx, fft_in, fft_out, int(fft_sign, c_int), fftw_estimate + fftw_unaligned)
        endif

        call fftw_execute_dft(fft_plan(iplan), fft_in, fft_out)
    
    endsubroutine

    integer function fft_index(index, n)

        integer, intent(in)     :: index, n
        integer                 :: n2    

        n2 = n / 2
        fft_index = modulo(index - 1 + n2, n) - n2

    endfunction

endmodule