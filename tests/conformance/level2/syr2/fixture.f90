      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision x(6), y(6), a(6, 6), ac(6, 6)
         data((ac(i, j), i=1, 6), j=1, 6)/1.053750863028617, 0, 0, 0, 0, 0,&
         &0.197684262345795, -1.068692711254786, 0, 0, 0, 0,&
         &0.2626454586627623, -1.2329011995712644,&
         &-0.00372353379218051,&
         &0, 0, 0,&
         &-0.9740025611125269,&
         &0.6893726977654734,&
         &-0.955839103276798,&
         &-1.2317070584140966,&
         &0,&
         &0,&
         &-0.9106806824932887,&
         &0.7412763052602079,&
         &0.06851153327714439,&
         &-0.3237507545879617,&
         &-1.0865030469936974,&
         &0,&
         &-0.767790184730859,&
         &-1.1197200611269833,&
         &-0.4481742366033955,&
         &0.47173637445323024,&
         &-1.180490682884277,&
         &1.4702569970829857/
         data(x(i), i=1, 6)/&
         &-0.08252376201716412,&
         &0.6060734308621007,&
         &-0.8874201453170976,&
         &0.10542139019376515,&
         &0.3528744733184766,&
         &0.5503933584550523/
         data(y(i), i=1, 6)/&
         &-7.7995958197861910,&
         &-0.90987740363925695,&
         &-1.3999755270779133,&
         &14.987896066159010,&
         &6.2202929537743330,&
         &11.557443037629128/
         call json%initialize()
         call json%create_array(root, '')

         a = ac
         call addcase(root, '0', 'u', 6, 1._8, x, 1, y, 1, a, 6)

         a = ac
         call addcase(root, '1', 'l', 6, 1._8, x, -1, y, -1, a, 6)

         x(2) = 0
         x(5) = 0
         y(1) = 0
         y(2) = 0
         y(4) = 0
         a = ac
         call addcase(root, '2', 'l', 6, 1._8, x, -1, y, -1, a, 6)

         a = ac
         call addcase(root, '3', 'u', 6, 1._8, x, -1, y, -1, a, 6)

         a = ac
         call addcase(root, '4', 'u', 0, 1._8, x, -1, y, -1, a, 6)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, n, alpha, x, incx, y, incy, a, lda)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dsyr2
            character(len=*) ncase, uplo
            integer n, lda, incx, incy
            double precision alpha, y(6), x(6), a(6, 6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'n', n)
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
            call json%add(case, 'y', y)
            call json%add(case, 'incy', incy)
            call dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
            call json%create_array(array, 'expect')
            do j = 1, 6
               do i = 1, 6
                  call json%add(array, '', a(i, j))
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
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
            call json%print(root, './tests/fixtures/level2/syr2.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
