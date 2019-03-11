      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision x(6), xc(6), a(6, 6), upper(6, 6), lower(6, 6)
         data((upper(i, j), i=1, 6), j=1, 6)/&
         &0,&
         &0,&
         &0,&
         &0,&
         &0,&
         &1.053750863028617,&
         &0,&
         &0,&
         &0,&
         &0,&
         &0.197684262345795,&
         &-1.068692711254786,&
         &0,&
         &0,&
         &0,&
         &0.2626454586627623,&
         &-1.2329011995712644,&
         &-0.00372353379218051,&
         &0,&
         &0,&
         &-0.9740025611125269,&
         &0.6893726977654734,&
         &-0.955839103276798,&
         &-1.2317070584140966,&
         &0,&
         &-0.9106806824932887,&
         &0.7412763052602079,&
         &0.06851153327714439,&
         &-0.3237507545879617,&
         &-1.0865030469936974,&
         &-0.767790184730859,&
         &-1.1197200611269833,&
         &-0.4481742366033955,&
         &0.47173637445323024,&
         &-1.180490682884277,&
         &1.4702569970829857/
         data((lower(i, j), i=1, 6), j=1, 6)/&
          &1.2629542848807933,&
          &-0.3262333607056494,&
          &1.3297992629225006,&
          &1.2724293214294047,&
          &0.4146414344564082,&
          &-1.5399500419037095,&
          &-0.2947204467905602,&
          &-0.005767172747536955,&
          &2.404653388857951,&
          &0.7635934611404596,&
          &-0.7990092489893682,&
          &0,&
          &-0.29921511789731614,&
          &-0.411510832795067,&
          &0.2522234481561323,&
          &-0.8919211272845686,&
          &0,&
          &0,&
          &0.37739564598170106,&
          &0.1333363608148414,&
          &0.8041895097449078,&
          &0,&
          &0,&
          &0,&
          &-1.2845993538721883,&
          &0.04672617218835198,&
          &0,&
          &0,&
          &0,&
          &0,&
          &1.1519117540872,&
          &0,&
          &0,&
          &0,&
          &0,&
          &0/
         data(xc(i), i=1, 6)/&
         &-0.08252376201716412,&
         &0.6060734308621007,&
         &0,&
         &0.10542139019376515,&
         &0.3528744733184766,&
         &0.5503933584550523/
         call json%initialize()
         call json%create_array(root, '')

         x = xc
         call addcase(root, '0', 'u', 'n', 'n', 6, 5, upper, 6, x, 1)

         x = xc
         call addcase(root, '1', 'u', 'n', 'u', 6, 5, upper, 6, x, 1)

         x = xc
         call addcase(root, '2', 'l', 'n', 'n', 6, 5, lower, 6, x, 1)

         x = xc
         call addcase(root, '3', 'l', 'n', 'u', 6, 5, lower, 6, x, 1)

         x = xc
         call addcase(root, '4', 'u', 't', 'n', 6, 5, upper, 6, x, -1)

         x = xc
         call addcase(root, '5', 'u', 't', 'u', 6, 5, upper, 6, x, -1)

         x = xc
         call addcase(root, '6', 'l', 't', 'n', 6, 5, lower, 6, x, -1)

         x = xc
         call addcase(root, '7', 'l', 't', 'u', 6, 5, lower, 6, x, -1)

         x = xc
         x(3) = 0
         call addcase(root, '8', 'l', 't', 'u', 0, 5, lower, 6, x, -1)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, trans, diag, n, k, a, lda, x, incx)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dtbmv
            character(len=*) ncase, uplo, trans, diag
            integer n, k, lda, incx
            double precision x(6), a(6, 6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'trans', trans)
            call json%add(case, 'diag', diag)
            call json%add(case, 'n', n)
            call json%add(case, 'k', k)
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
            call dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
            call json%add(case, 'expect', x)
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
            call json%print(root, './tests/fixtures/level2/tbmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
