      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         complex*16 lower(6, 6), upper(6, 6), v6(6), x(6)
         data((upper(i, j), i=1, 6), j=1, 6)/&
         &(1.2629542848807933, -0.42951310949188126),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(-0.9285670347135381, -0.8320432961178319),&
         &(-0.2947204467905602, -1.166570547084707),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(-1.1476570092363514, -0.22732869142475534),&
         &(-0.28946157368822334, 0.2661373616721048),&
         &(-0.29921511789731614, -0.3767027185836281),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0.43568329935571865, 0.2501413228541527),&
         &(-1.237538421929958, 0.6182432935662469),&
         &(-0.22426788527830935, -0.17262350264585732),&
         &(0.37739564598170106, -2.2239002740099374),&
         &(0, 0),&
         &(0, 0),&
         &(-0.057106774383808755, -0.011045478465663564),&
         &(0.5036079722337261, -0.9406491626186084),&
         &(1.085769362145687, -0.11582532215695436),&
         &(-0.6909538396968303, -0.8149687088699175),&
         &(-1.2845993538721883, 0.24226348085968588),&
         &(0, 0),&
         &(-0.23570655643950122, 0.36594112304921983),&
         &(-0.5428882550102544, 0.2484126488725964),&
         &(-0.4333103174567822, 0.06528818167162072),&
         &(-0.6494716467962331, 0.01915639166027384),&
         &(0.726750747385451, 0.2573383771555333),&
         &(1.1519117540872, 0)/
         data((lower(i, j), i=1, 6), j=1, 6)/&
         &(1.2629542848807933, -0.42951310949188126),&
         &(-0.3262333607056494, 1.2383041008533804),&
         &(1.3297992629225006, -0.2793462818542693),&
         &(1.2724293214294047, 1.7579030898107073),&
         &(0.4146414344564082, 0.5607460908880562),&
         &(-1.5399500419037095, -0.4527839725531578),&
         &(0, 0),&
         &(-0.2947204467905602, -1.166570547084707),&
         &(-0.005767172747536955, -1.0655905803882961),&
         &(2.404653388857951, -1.563782051071005),&
         &(0.7635934611404596, 1.1565369971501793),&
         &(-0.7990092489893682, 0.8320471285723897),&
         &(0, 0),&
         &(0, 0),&
         &(-0.29921511789731614, -0.3767027185836281),&
         &(-0.411510832795067, 2.4413646288945894),&
         &(0.2522234481561323, -0.7953391172553718),&
         &(-0.8919211272845686, -0.054877473711578625),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0.37739564598170106, -2.2239002740099374),&
         &(0.1333363608148414, -1.263614384970583),&
         &(0.8041895097449078, 0.3587288959713519),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(-1.2845993538721883, 0.24226348085968588),&
         &(0.04672617218835198, -1.4250983947324998),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(1.1519117540872, 0)/
         data(v6(i), i=1, 6)/&
         &(-0.6490100777088978, 0.7721421858045301),&
         &(-0.11916876241803812, -0.21951562675343952),&
         &(0.6641356998941105, -0.4248102833772871),&
         &(1.100969102194087, -0.418980099421959),&
         &(0.14377148075806995, 0.9969868609091059),&
         &(-0.11775359816595128, -0.27577802908802723)/

         !FIXME change the input matrix, this is related the the bug in the slice function [https://github.com/R-js/blasjs/issues/3].

         call json%initialize()
         call json%create_array(root, '')

         x = v6
         x(6) = 0
         call addcase(root, '0', 'u', 'n', 'n', 6, lower, 6, x, 1)

         x = v6
         x(6) = 0
         call addcase(root, '1', 'u', 'n', 'u', 6, upper, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '2', 'l', 'n', 'n', 6, lower, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '3', 'l', 'n', 'u', 6, lower, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '4', 'u', 't', 'n', 6, upper, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '5', 'u', 't', 'u', 6, upper, 6, x, 1)

         x = v6
         x(6) = 0
         call addcase(root, '6', 'l', 't', 'n', 6, lower, 6, x, 1)

         x = v6
         x(6) = 0
         call addcase(root, '7', 'l', 't', 'u', 6, lower, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '8', 'u', 'c', 'n', 6, upper, 6, x, 1)

         x = v6
         x(1) = 0
         call addcase(root, '9', 'u', 'c', 'u', 6, upper, 6, x, 1)

         x = v6
         x(6) = 0
         call addcase(root, '10', 'l', 'c', 'n', 6, lower, 6, x, 1)

         x = v6
         x(6) = 0
         call addcase(root, '11', 'l', 'c', 'u', 6, lower, 6, x, -1)

         x = v6
         call addcase(root, '12', 'l', 'c', 'u', 0, lower, 6, x, -1)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, trans, diag, n, a, lda, x, incx)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external ztrsv
            character(len=*) ncase, uplo, trans, diag
            integer n, incx, lda
            complex*16 a(6, 6), x(6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'trans', trans)
            call json%add(case, 'diag', diag)
            call json%add(case, 'n', n)
            call json%create_array(array, 'a')
            do j = 1, 6
               do i = 1, 6
                  call json%create_array(data, '')
                  call json%add(data, '', dreal(a(i, j)))
                  call json%add(data, '', dimag(a(i, j)))
                  call json%add(array, data)
                  nullify (data)
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'lda', lda)
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
            nullify (array)
            call ztrsv(uplo, trans, diag, n, a, lda, x, incx)
            call json%create_array(array, 'expect')
            do i = 1, 6
               call json%create_array(data, '')
               call json%add(data, '', dreal(x(i)))
               call json%add(data, '', dimag(x(i)))
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
            call json%print(root, './tests/fixtures/level2/complex/trsv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
