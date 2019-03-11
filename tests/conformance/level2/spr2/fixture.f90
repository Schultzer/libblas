      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision alpha, x(8), xc(8), y(8), yc(8), xcase1(8), ycase1(8)
         ! apsize = n * (k + 1) - k * (k + 1) / 2
         double precision a(36), ac(36)
         data(ac(i), i=1, 36)/2.024761390077735,&
         &-1.0457176517867426,&
         &-0.8962112639605789,&
         &-0.060634778120552485,&
         &-0.5013783177950893,&
         &0.9260627253330198,&
         &-1.4577072101223822,&
         &0.09505622698864298,&
         &0.8476649635960255,&
         &-1.6243645297605886,&
         &1.5761581812187344,&
         &-1.475547635260523,&
         &-0.14460820731054116,&
         &-1.0750101908181724,&
         &0.40654273194492646,&
         &-0.14727079003896576,&
         &1.5415930688269568,&
         &-0.9818556688038707,&
         &0.4965781726616992,&
         &1.6969478807230851,&
         &-0.26073630856812263,&
         &0.5013218277237272,&
         &-1.0135396704947914,&
         &1.6147522354681292,&
         &0.005641984852495504,&
         &-2.9048990603455724,&
         &-1.1071648189687495,&
         &1.5475669326182715,&
         &-0.10150344763172572,&
         &0.042650249796697896,&
         &-1.5967180142971973,&
         &0.490967372597059,&
         &0.421603365384753,&
         &1.8739038985953016,&
         &1.0345143239443348,&
         &0.08181031035401386/
         data(yc(i), i=1, 8)/&
         &0.7021167106675735,&
         &2.5071111484833684,&
         &-1.890027143624024,&
         &-0.5898127901911715,&
         &-1.7145022968458246,&
         &-0.4209978978166964,&
         &0.310141376504687,&
         &1.7025705860144955/
         data(xc(i), i=1, 8)/&
         &-0.08252376201716412,&
         &0.6060734308621007,&
         &-0.8874201453170976,&
         &0.10542139019376515,&
         &0.3528744733184766,&
         &0.5503933584550523,&
         &-1.1343309685168443,&
         &1.4623515387464268/
         data(xcase1(i), i=1, 8)/1, 2, 3, 4, 5, 6, 7, 8/
         data(ycase1(i), i=1, 8)/9, 12, 13, 14, 15, 16, 17, 18/
         call json%initialize()
         call json%create_array(root, '')

         x = xc
         y = yc
         a = ac
         call addcase(root, '0', 'u', 8, 1.0_8, x, 1, y, 1, a)

         x = xcase1
         y = ycase1
         x(3) = 0
         x(5) = 0
         y(5) = 0
         a = ac
         call addcase(root, '1', 'l', 8, 0.25_8, x, -1, y, -1, a)

         a = ac
         call addcase(root, '2', 'u', 8, 0.25_8, x, -1, y, -1, a)

         a = ac
         call addcase(root, '3', 'u', 8, 0._8, x, -1, y, -1, a)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, n, alpha, x, incx, y, incy, ap)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dspr2
            character(len=*) ncase, uplo
            integer n, incx, incy
            double precision alpha, y(8), x(8), ap(36)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'n', n)
            call json%add(case, 'alpha', alpha)
            call json%add(case, 'x', x)
            call json%add(case, 'incx', incx)
            call json%add(case, 'y', y)
            call json%add(case, 'incy', incy)
            call json%add(case, 'ap', ap)
            call dspr2(uplo, n, alpha, x, incx, y, incy, ap)
            call json%add(case, 'expect', ap)
         end subroutine addcase

         subroutine print(root)
            use, intrinsic :: iso_fortran_env, only: real64
            use json_module, CK => json_CK, CK => json_CK
            implicit none
            type(json_core) :: json
            type(json_value), pointer :: root
            logical :: status_ok
            character(kind=CK, len=:), allocatable :: error_msg
            call json%print(root, './tests/fixtures/level2/spr2.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
