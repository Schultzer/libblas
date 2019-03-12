      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         complex*16 upper(4, 6), lower(4, 6), x(6), y(6), v6(6)
         data((upper(i, j), i=1, 4), j=1, 6)/&
         &(0, 0),&
         &(0, 0),&
         &(0, 0),&
         &(1.2629542350769043, 0.9921603654457979),&
         &(0, 0),&
         &(0, 0),&
         &(-0.9285670518875122, -0.4527839725531578),&
         &(-0.29472044110298157, -0.8320432961178319),&
         &(0, 0),&
         &(-1.147657036781311, 0.8320471285723897),&
         &(-0.28946158289909363, -0.22732869142475534),&
         &(-0.2992151081562042, 0.2661373616721048),&
         &(0.43568331003189087, -0.054877473711578625),&
         &(-1.237538456916809, 0.2501413228541527),&
         &(-0.2242678850889206, 0.6182432935662469),&
         &(0.3773956596851349, -0.17262350264585732),&
         &(0.503607988357544, -0.011045478465663564),&
         &(1.0857694149017334, -0.9406491626186084),&
         &(-0.6909538507461548, -0.11582532215695436),&
         &(-1.2845993041992188, -0.8149687088699175),&
         &(-0.433310329914093, 0.2484126488725964),&
         &(-0.649471640586853, 0.06528818167162072),&
         &(0.7267507314682007, 0.01915639166027384),&
         &(1.151911735534668, 0.2573383771555333)/
         data((lower(i, j), i=1, 4), j=1, 6)/&
         &(1.2629542848807933, 0.9921603654457979),&
         &(-0.3262333607056494, -0.42951310949188126),&
         &(1.3297992629225006, 1.2383041008533804),&
         &(1.2724293214294047, -0.2793462818542693),&
         &(-0.2947204467905602, -0.8320432961178319),&
         &(-0.005767172747536955, -1.166570547084707),&
         &(2.404653388857951, -1.0655905803882961),&
         &(0.7635934611404596, -1.563782051071005),&
         &(-0.29921511789731614, 0.2661373616721048),&
         &(-0.411510832795067, -0.3767027185836281),&
         &(0.2522234481561323, 2.4413646288945894),&
         &(-0.8919211272845686, -0.7953391172553718),&
         &(0.37739564598170106, -0.17262350264585732),&
         &(0.1333363608148414, -2.2239002740099374),&
         &(0.8041895097449078, -1.263614384970583),&
         &(0, 0),&
         &(-1.2845993538721883, -0.8149687088699175),&
         &(0.04672617218835198, 0.24226348085968588),&
         &(0, 0),&
         &(0, 0),&
         &(1.1519117540872, 0.2573383771555333),&
         &(0, 0),&
         &(0, 0),&
         &(0, 0)/
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
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '0', 'u', 6, 3, (0.2_8, 0.8_8), upper, 4, x, 1, (0.3_8, -0.7_8), y, 1)

         x = v6
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '1', 'u', 6, 3, (0.2_8, 0.8_8), upper, 4, x, -1, (0._8, 0._8), y, -1)

         x = v6
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '2', 'u', 6, 3, (0._8, 0._8), upper, 4, x, -1, (0._8, 0._8), y, -1)

         x = v6
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '3', 'u', 6, 3, (0._8, 0._8), upper, 4, x, -1, (1._8, 0._8), y, -1)

         x = v6
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '4', 'u', 6, 3, (0.2_8, 0.8_8), upper, 4, x, -1, (1._8, 0._8), y, -1)

         x = v6
         y = x
         y(4) = 0
         y(2) = 0
         call addcase(root, '5', 'l', 6, 3, (0.2_8, 0.8_8), lower, 4, x, 1, (1._8, 0._8), y, 1)

         call print(root)
      contains
         subroutine addcase(root, ncase, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external zhbmv
            character(len=*) ncase, uplo
            integer n, k, lda, incx, incy
            complex*16 alpha, beta, a(4, 6), x(6), y(6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'n', n)
            call json%add(case, 'k', k)
            call json%create_array(data, 'alpha')
            call json%add(data, '', dreal(alpha))
            call json%add(data, '', dimag(alpha))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'a')
            do j = 1, 6
               do i = 1, 4
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
            call zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
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
            call json%print(root, './tests/fixtures/level2/complex/hbmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
