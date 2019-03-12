      PROGRAM fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         complex*16 ap(21), v6(6), y(6), x(6)
         data(ap(i), i=1, 21)/&
         &(1.2629542848807933, 0.9921603654457979),&
         &(-0.9285670347135381, -0.4527839725531578),&
         &(-0.2947204467905602, -0.8320432961178319),&
         &(-1.1476570092363514, 0.8320471285723897),&
         &(-0.28946157368822334, -0.22732869142475534),&
         &(-0.29921511789731614, 0.2661373616721048),&
         &(0.43568329935571865, -0.054877473711578625),&
         &(-1.237538421929958, 0.2501413228541527),&
         &(-0.22426788527830935, 0.6182432935662469),&
         &(0.37739564598170106, -0.17262350264585732),&
         &(-0.057106774383808755, 0.3587288959713519),&
         &(0.5036079722337261, -0.011045478465663564),&
         &(1.085769362145687, -0.9406491626186084),&
         &(-0.6909538396968303, -0.11582532215695436),&
         &(-1.2845993538721883, -0.8149687088699175),&
         &(-0.23570655643950122, -1.4250983947324998),&
         &(-0.5428882550102544, 0.36594112304921983),&
         &(-0.4333103174567822, 0.2484126488725964),&
         &(-0.6494716467962331, 0.06528818167162072),&
         &(0.726750747385451, 0.01915639166027384),&
         &(1.1519117540872, 0.2573383771555333)/
         data(v6(i), i=1, 6)/&
         &(-0.6490100777088978, 0.7721421858045301),&
         &(-0.11916876241803812, -0.21951562675343952),&
         &(0.6641356998941105, -0.4248102833772871),&
         &(1.100969102194087, -0.418980099421959),&
         &(0.14377148075806995, 0.9969868609091059),&
         &(-0.11775359816595128, -0.27577802908802723)/

         call json%initialize()
         call json%create_array(root, '')

         x = v6
         x(4) = 0
         x(2) = 0

         y = v6
         y(4) = 0
         call addcase(root, '0', 'u', 6, (0.2_8, -0.8_8), ap, x, 1, (0.3_8, -0.7_8), y, 1)

         y = v6
         y(4) = 0
         call addcase(root, '1', 'l', 6, (0.2_8, -0.8_8), ap, x, -1, (0.3_8, -0.7_8), y, -1)

         y = v6
         y(4) = 0
         call addcase(root, '2', 'l', 6, (0._8, 0._8), ap, x, -1, (1._8, 0._8), y, -1)

         y = v6
         y(4) = 0
         call addcase(root, '3', 'l', 6, (0._8, 0._8), ap, x, -1, (0._8, 0._8), y, -1)

         y = v6
         y(4) = 0
         call addcase(root, '4', 'l', 6, (1._8, 0._8), ap, x, -1, (1._8, 0._8), y, -1)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, n, alpha, ap, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external zhpmv
            character(len=*) ncase, uplo
            integer n, incx, incy
            complex*16 alpha, beta, ap(21), x(6), y(6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'n', n)
            call json%create_array(data, 'alpha')
            call json%add(data, '', dreal(alpha))
            call json%add(data, '', dimag(alpha))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'ap')
            do i = 1, 21
               call json%create_array(data, '')
               call json%add(data, '', dreal(ap(i)))
               call json%add(data, '', dimag(ap(i)))
               call json%add(array, data)
               nullify (data)
            enddo
            call json%add(case, array)
            nullify (array)
            call json%create_array(array, 'x')
            do i = 1, 6
               call json%create_array(data, '')
               call json%add(data, '', dreal(x(i)))
               call json%add(data, '', dimag(x(i)))
               call json%add(array, data)
               nullify (data)
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'incx', incx)
            call json%create_array(data, 'beta')
            call json%add(data, '', dreal(beta))
            call json%add(data, '', dimag(beta))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'y')
            do i = 1, 6
               call json%create_array(data, '')
               call json%add(data, '', dreal(y(i)))
               call json%add(data, '', dimag(y(i)))
               call json%add(array, data)
               nullify (data)
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'incy', incy)
            call zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            call json%create_array(array, 'expect')
            do i = 1, 6
               call json%create_array(data, '')
               call json%add(data, '', dreal(y(i)))
               call json%add(data, '', dimag(y(i)))
               call json%add(array, data)
               nullify (data)
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
            call json%print(root, './tests/fixtures/level2/complex/hpmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
