      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         character trans
         integer m, n, kl, ku, lda, incx, incy
         double precision alpah, beta, yc(9), xc(9), a(3, 9), y(9), x(9)
         data((a(i, j), i=1, 3), j=1, 9)/-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0,&
         &-1.0, 2.0, -1.0/
         data(yc(i), i=1, 9)/1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8, 1.0_8/
         data(x(i), i=1, 9)/1.0_8, 1.0_8, 2.0_8, 2.0_8, 3.0_8, 3.0_8, 4.0_8, 4.0_8, 5.0_8/
         call json%initialize()
         call json%create_array(root, '')

         y = yc
         call addcase(root, '0', 'n', 6, 9, 1, 1, 1.5_8, a, 3, x, 1, 2.5_8, y, 1)

         y = yc
         call addcase(root, '1', 'n', 6, 0, 1, 1, 1.5_8, a, 3, x, 1, 2.5_8, yc, 1)

         y = yc
         call addcase(root, '2', 'n', 6, 9, 1, 1, 0._8, a, 3, x, 1, 1._8, y, 1)

         y = yc
         call addcase(root, '3', 'n', 3, 9, 1, 1, 1.5_8, a, 3, x, 1, 2.5_8, y, 2)

         y = yc
         call addcase(root, '4', 'n', 3, 9, 1, 1, 1.5_8, a, 3, x, 1, 0._8, y, 2)

         y = yc
         call addcase(root, '5', 'n', 6, 9, 1, 1, 1.5_8, a, 3, x, 1, 0._8, y, 1)

         y = yc
         call addcase(root, '6', 't', 6, 9, 1, 1, 1.5_8, a, 3, x, 1, 1._8, y, 1)

         y = yc
         call addcase(root, '7', 't', 6, 9, 1, 1, 0._8, a, 3, x, 1, 2.5_8, y, 1)

         y = yc
         x(9) = 3.0_8
         y(1) = 9
         y(2) = 0
         y(3) = 0
         y(4) = 0
         y(5) = 0
         y(6) = 0
         y(7) = 0
         y(8) = 0
         y(9) = -9
         call addcase(root, '8', 'n', 6, 9, 1, 1, 1._8, a, 3, x, -1, 0._8, y, -1)

         call print(root)
      contains
         subroutine addcase(root, ncase, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dgbmv
            character(len=*) ncase, trans
            integer m, n, kl, ku, lda, incx, incy
            double precision alpha, beta, y(9), x(9), a(3, 9)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'trans', trans)
            call json%add(case, 'm', m)
            call json%add(case, 'n', n)
            call json%add(case, 'kl', kl)
            call json%add(case, 'ku', ku)
            call json%add(case, 'alpha', alpha)
            call json%create_array(array, 'a')
            do j = 1, 9
               do i = 1, 3
                  call json%add(array, '', a(i, j))
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'lda', lda)
            call json%add(case, 'x', x)
            call json%add(case, 'incx', incx)
            call json%add(case, 'beta', beta)
            call json%add(case, 'y', y)
            call json%add(case, 'incy', incy)
            call dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            call json%add(case, 'expect', y)
            nullify (case)
         end subroutine addcase

         subroutine print(root)
            use, intrinsic :: iso_fortran_env, only: real64
            use json_module, CK => json_CK, CK => json_CK
            implicit none
            type(json_core) :: json
            type(json_value), pointer :: root
            logical :: status_ok
            character(kind=CK, len=:), allocatable :: error_msg
            call json%print(root, './tests/fixtures/level2/gbmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
