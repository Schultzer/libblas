
      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         complex*16 a(6, 8), x(8), v(6), y(6), y2(8)
         double precision are(6, 8), aim(6, 8)
         data((are(i, j), i=1, 6), j=1, 8)/&
         &1.2629542848807933,&
         &-0.3262333607056494,&
         &1.3297992629225006,&
         &1.2724293214294047,&
         &0,&
         &0,&
         &0.4146414344564082,&
         &-1.5399500419037095,&
         &-0.9285670347135381,&
         &-0.2947204467905602,&
         &-0.005767172747536955,&
         &0,&
         &2.404653388857951,&
         &0.7635934611404596,&
         &-0.7990092489893682,&
         &-1.1476570092363514,&
         &-0.28946157368822334,&
         &-0.29921511789731614,&
         &0,&
         &-0.411510832795067,&
         &0.2522234481561323,&
         &-0.8919211272845686,&
         &0.43568329935571865,&
         &-1.237538421929958,&
         &0,&
         &0,&
         &-0.22426788527830935,&
         &0.37739564598170106,&
         &0.1333363608148414,&
         &0.8041895097449078,&
         &0,&
         &0,&
         &0,&
         &-0.057106774383808755,&
         &0.5036079722337261,&
         &1.085769362145687,&
         &0,&
         &0,&
         &0,&
         &0,&
         &-0.6909538396968303,&
         &-1.2845993538721883,&
         &0,&
         &0,&
         &0,&
         &0,&
         &0,&
         &0.04672617218835198/
         data((aim(i, j), i=1, 6), j=1, 8)/&
         &0.9921603654457979,&
         &-0.42951310949188126,&
         &1.2383041008533804,&
         &-0.2793462818542693,&
         &0,&
         &0,&
         &1.7579030898107073,&
         &0.5607460908880562,&
         &-0.4527839725531578,&
         &-0.8320432961178319,&
         &-1.166570547084707,&
         &0,&
         &-1.0655905803882961,&
         &-1.563782051071005,&
         &1.1565369971501793,&
         &0.8320471285723897,&
         &-0.22732869142475534,&
         &0.2661373616721048,&
         &0,&
         &-0.3767027185836281,&
         &2.4413646288945894,&
         &-0.7953391172553718,&
         &-0.054877473711578625,&
         &0.2501413228541527,&
         &0,&
         &0,&
         &0.6182432935662469,&
         &-0.17262350264585732,&
         &-2.2239002740099374,&
         &-1.263614384970583,&
         &0,&
         &0,&
         &0,&
         &0.3587288959713519,&
         &-0.011045478465663564,&
         &-0.9406491626186084,&
         &0,&
         &0,&
         &0,&
         &0,&
         &-0.11582532215695436,&
         &-0.8149687088699175,&
         &0,&
         &0,&
         &0,&
         &0,&
         &0,&
         &0.24226348085968588/
         data(x(i), i=1, 8)/&
         &(-0.6490100777088978, 0.7721421858045301),&
         &(-0.11916876241803812, -0.21951562675343952),&
         &(0.6641356998941105, -0.4248102833772871),&
         &(1.100969102194087, -0.418980099421959),&
         &(0.14377148075806995, 0.9969868609091059),&
         &(-0.11775359816595128, -0.27577802908802723),&
         &(-0.9120683669483379, 1.2560188173061),&
         &(-1.4375862408299789, 0.6466743904953449)/
         data(v(i), i=1, 6)/&
         &(-0.6490100777088978, 0.7721421858045301),&
         &(-0.11916876241803812, -0.21951562675343952),&
         &(0.6641356998941105, -0.4248102833772871),&
         &(1.100969102194087, -0.418980099421959),&
         &(0.14377148075806995, 0.9969868609091059),&
         &(-0.11775359816595128, -0.27577802908802723)/
         call json%initialize()
         call json%create_array(root, '')

         call copy(are, aim, a, 8, 6)

         y = v
         call addcase(root, '0', 'n', 6, 8, 3, 2, (0._8, 0._8), a, 6, x, 1, (2.5_8, 0.5_8), y, 1)

         y = v
         call addcase(root, '1', 'n', 6, 8, 3, 2, (0.2_8, 0.8_8), a, 6, x, 1, (2.5_8, 0.5_8), y, 1)

         y2 = x
         call addcase(root, '2', 't', 6, 8, 3, 2, (0.2_8, 0.8_8), a, 6, v, 1, (0._8, 0._8), y2, 1)

         y2 = x
         call addcase(root, '3', 'c', 6, 8, 3, 2, (0.2_8, 0.8_8), a, 6, v, -1, (1._8, 0._8), y2, -1)

         y2 = x
         call addcase(root, '4', 'c', 6, 8, 3, 2, (0._8, 0._8), a, 6, v, -1, (1._8, 0._8), y2, -1)

         call print(root)
      contains
         subroutine addcase(root, ncase, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external zgbmv, lsame
            logical lsame
            character(len=*) ncase, trans
            integer m, n, kl, ku, lda, incx, incy, y_end, x_end
            complex*16 alpha, beta, a(6, 8), x(*), y(*)
            if (lsame(trans, 't') .or. lsame(trans, 'c')) then
               y_end = 8
               x_end = 6
            else
               y_end = 6
               x_end = 8
            end if
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'trans', trans)
            call json%add(case, 'm', m)
            call json%add(case, 'n', n)
            call json%add(case, 'kl', kl)
            call json%add(case, 'ku', ku)
            call json%create_array(data, 'alpha')
            call json%add(data, '', dreal(alpha))
            call json%add(data, '', dimag(alpha))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'a')
            do j = 1, 8
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
            do i = 1, int(x_end)
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
            do i = 1, int(y_end)
               call json%create_array(data, '')
               call json%add(data, '', dreal(y(i)))
               call json%add(data, '', dimag(y(i)))
               call json%add(array, data)
               nullify (data)
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'incy', incy)
            call zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            call json%create_array(array, 'expect')
            do i = 1, int(y_end)
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
            call json%print(root, './tests/fixtures/level2/complex/gbmv.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print

         subroutine copy(re, im, a, n, lda)
            integer n, lda, cur
            complex*16 a(lda, n)
            double precision re(lda*n), im(lda*n)
            do j = 1, n
               do i = 1, lda
                  cur = (j - 1)*lda + i
                  a(i, j) = complex(re(cur), im(cur))
               enddo
            enddo
         end subroutine copy

      end
