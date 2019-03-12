      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision alpha, beta, yc(6), y(6), x(6), upper(6, 6), lower(6, 6)
         data((upper(i, j), i=1, 6), j=1, 6)/0, 0, 0, 0, 0, 1.2629542848807933,&
        &0, 0, 0, 0, -0.9285670347135381, -0.2947204467905602,&
        &0, 0, 0, -1.1476570092363514, -0.28946157368822334,&
        &-0.29921511789731614,&
        &0, 0, 0.43568329935571865,&
        &-1.237538421929958,&
        &-0.22426788527830935, 0.37739564598170106,&
        &0, -0.057106774383808755, 0.5036079722337261, 1.085769362145687,&
        &-0.6909538396968303, -1.2845993538721883,&
        &-0.23570655643950122, -0.5428882550102544, -0.4333103174567822,&
        &-0.6494716467962331, 0.726750747385451, 1.1519117540872/
         data(yc(i), i=1, 6)/&
        &-0.17262350264585732,&
        &-2.2239002740099374,&
        &-1.263614384970583,&
        &0.3587288959713519,&
        &-0.011045478465663564,&
        &-0.9406491626186084/
         data(x(i), i=1, 6)/&
        &-0.17262350264585732,&
        &-2.2239002740099374,&
        &-1.263614384970583,&
        &0.3587288959713519,&
        &-0.011045478465663564,&
        &-0.9406491626186084/
         data((lower(i, j), i=1, 6), j=1, 6)/0.42224244493991137,&
        &-0.7566161579452455, -0.5090229720808566,&
        &-0.7133912411518395, -0.5207411698065698, -0.8821312454529107,&
        &0.7525384253822267, 0.5578293548896909, 0.5946176517754793,&
        &-0.08945109276100993, -0.17983183590695262, 0,&
        &0.3094478561542928, -0.29360545612871647, -0.459479708224535,&
        &0.9853681223466992, 0, 0,&
        &-0.0437639313749969, 0.8481489396654069, 0.19752193428575993,&
        &0, 0, 0,&
        &-0.703576878644526, -0.9738448490388691, 0,&
        &0, 0, 0,&
        &-0.008812844287604094, 0, 0,&
        &0, 0, 0/
         call json%initialize()
         call json%create_array(root, '')

         y = yc
         call addcase(root, '0', 'u', 6, 5, 0.75_8, upper, 6, x, 1, 0.25_8, y, 1)

         y = yc
         call addcase(root, '1', 'l', 6, 5, 0.75_8, lower, 6, x, 1, 0.25_8, y, 1)

         y = yc
         call addcase(root, '2', 'l', 6, 5, 0.75_8, lower, 6, x, -1, 0.25_8, y, -1)

         y = yc
         call addcase(root, '3', 'l', 3, 5, 0.75_8, lower, 6, x, 2, 0.25_8, y, -2)

         y = yc
         call addcase(root, '3', 'l', 6, 5, 0._8, lower, 6, x, 1, 1._8, y, 1)

         y = yc
         call addcase(root, '5', 'l', 6, 5, 0.25_8, lower, 6, x, 1, 1._8, y, 1)

         y = yc
         call addcase(root, '6', 'l', 6, 5, 0._8, lower, 6, x, 1, 0._8, y, 1)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dsbmv
            character(len=*) ncase, uplo
            integer m, n, lda, incx, incy
            double precision alpha, beta, y(6), x(6), a(6, 6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'n', n)
            call json%add(case, 'k', k)
            call json%add(case, 'alpha', alpha)
            call json%create_array(array, 'a')
            do j = 1, 6
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
            call dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
            call json%add(case, 'expect', y)
         end subroutine addcase

         subroutine print(root)
            use, intrinsic :: iso_fortran_env, only: real64
            use json_module, CK => json_CK, CK => json_CK
            implicit none
            type(json_core) :: json
            type(json_value), pointer :: root
            logical :: status_ok
            character(kind=CK, len=:), allocatable :: error_msg
            call json%print(root, './tests/fixtures/level2/sbmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end

