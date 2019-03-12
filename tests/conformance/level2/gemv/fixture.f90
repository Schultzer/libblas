      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision yc(9), y(9), x(9), a(6, 9)
         data((a(i, j), i=1, 6), j=1, 9)/-0.054877470000,&
        & 0.250141320000, 0.618243290000, -0.172623500000,&
        & -2.223900270000, -1.26361440000,&
        &  0.358728900000, -0.011045480000, -0.940649160000,&
        & -0.115825320000, -0.814968710000,&
        &  0.24226350000,&
        &  -1.425098390000, 0.365941120000, 0.248412650000,&
        &  0.065288180000, 0.019156390000,&
        &  0.25733840000, -0.649010080000,&
        &  -0.119168760000, 0.664135700000, 1.100969100000,&
        &  0.143771480000, -0.11775360000,&
        &  -0.912068370000, -1.437586240000, -0.797089530000,&
        &  1.254083110000, 0.772142190000, -0.21951560000,&
        &  -0.424810280000, -0.418980100000, 0.996986860000,&
        &  -0.275778030000, 1.256018820000, 0.64667440000,&
        &  1.299312300000, -0.873262110000, 0.008370960000,&
        &  -0.880871720000,&
        &  0.596259020000, 0.11971760000,&
        &  -0.282173880000, 1.455988400000, 0.229019590000, 0.99654393000,&
        &  0.78185918000, -0.7767766000,&
        &  -0.61598991000, 0.04658030000, -1.130385780000, 0.576718780000,&
        &  -1.28074943000, 1.6254473000/
         data(x(i), i=1, 9)/1, 1, 2, 2, 3, 3, 4, 4, 5/
         data(yc(i), i=1, 9)/1, 1, 1, 1, 1, 1, 4, 4, 5/
         call json%initialize()
         call json%create_array(root, '')

         y = yc
         call addcase(root, '0', 'n', 6, 9, 1.5_8, a, 6, x, 1, 2.5_8, y, 1)

         y = yc
         call addcase(root, '1', 'n', 6, 9, 1.5_8, a, 6, x, -1, 0._8, y, -1)

         y = yc
         call addcase(root, '2', 'n', 6, 9, 1.5_8, a, 6, x, -1, 0._8, y, 1)

         y = yc
         call addcase(root, '3', 't', 6, 9, 1.5_8, a, 6, x, -1, 1._8, y, 1)

         y = yc
         call addcase(root, '4', 't', 6, 9, 0._8, a, 6, x, -1, 1._8, y, 1)

         y = yc
         call addcase(root, '5', 't', 6, 9, 0._8, a, 6, x, -1, 1.5_8, y, 1)

         y = yc
         call addcase(root, '6', 't', 3, 4, 1._8, a, 6, x, -2, 0._8, y, -2)
         call print(root)
      contains
         subroutine addcase(root, ncase, trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dgemv
            character(len=*) ncase, trans
            integer m, n, lda, incx, incy
            double precision alpha, beta, y(9), x(9), a(lda, 9)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'trans', trans)
            call json%add(case, 'm', m)
            call json%add(case, 'n', n)
            call json%add(case, 'alpha', alpha)
            call json%create_array(array, 'a')
            do j = 1, 9
               do i = 1, 6
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
            call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
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
            call json%print(root, './tests/fixtures/level2/gemv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end

