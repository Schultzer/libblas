      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         complex*16 m8x8(8, 8), c(8, 8)
         data((m8x8(i, j), i=1, 8), j=1, 8)/&
         &(1.2629542848807933, 0.9921603654457979),&
         &(-0.3262333607056494, -0.42951310949188126),&
         &(1.3297992629225006, 1.2383041008533804),&
         &(1.2724293214294047, -0.2793462818542693),&
         &(0.4146414344564082, 1.7579030898107073),&
         &(-1.5399500419037095, 0.5607460908880562),&
         &(-0.9285670347135381, -0.4527839725531578),&
         &(-0.2947204467905602, -0.8320432961178319),&
         &(-0.005767172747536955, -1.166570547084707),&
         &(2.404653388857951, -1.0655905803882961),&
         &(0.7635934611404596, -1.563782051071005),&
         &(-0.7990092489893682, 1.1565369971501793),&
         &(-1.1476570092363514, 0.8320471285723897),&
         &(-0.28946157368822334, -0.22732869142475534),&
         &(-0.29921511789731614, 0.2661373616721048),&
         &(-0.411510832795067, -0.3767027185836281),&
         &(0.2522234481561323, 2.4413646288945894),&
         &(-0.8919211272845686, -0.7953391172553718),&
         &(0.43568329935571865, -0.054877473711578625),&
         &(-1.237538421929958, 0.2501413228541527),&
         &(-0.22426788527830935, 0.6182432935662469),&
         &(0.37739564598170106, -0.17262350264585732),&
         &(0.1333363608148414, -2.2239002740099374),&
         &(0.8041895097449078, -1.263614384970583),&
         &(-0.057106774383808755, 0.3587288959713519),&
         &(0.5036079722337261, -0.011045478465663564),&
         &(1.085769362145687, -0.9406491626186084),&
         &(-0.6909538396968303, -0.11582532215695436),&
         &(-1.2845993538721883, -0.8149687088699175),&
         &(0.04672617218835198, 0.24226348085968588),&
         &(-0.23570655643950122, -1.4250983947324998),&
         &(-0.5428882550102544, 0.36594112304921983),&
         &(-0.4333103174567822, 0.2484126488725964),&
         &(-0.6494716467962331, 0.06528818167162072),&
         &(0.726750747385451, 0.01915639166027384),&
         &(1.1519117540872, 0.2573383771555333),&
         &(0.9921603654457979, 1.2629542848807933),&
         &(-0.42951310949188126, -0.3262333607056494),&
         &(1.2383041008533804, 1.3297992629225006),&
         &(-0.2793462818542693, 1.2724293214294047),&
         &(1.7579030898107073, 0.4146414344564082),&
         &(0.5607460908880562, -1.5399500419037095),&
         &(-0.4527839725531578, -0.9285670347135381),&
         &(-0.8320432961178319, -0.2947204467905602),&
         &(-1.166570547084707, -0.005767172747536955),&
         &(-1.0655905803882961, 2.404653388857951),&
         &(-1.563782051071005, 0.7635934611404596),&
         &(1.1565369971501793, -0.7990092489893682),&
         &(0.8320471285723897, -1.1476570092363514),&
         &(-0.22732869142475534, -0.28946157368822334),&
         &(0.2661373616721048, -0.29921511789731614),&
         &(-0.3767027185836281, -0.411510832795067),&
         &(2.4413646288945894, 0.2522234481561323),&
         &(-0.7953391172553718, -0.8919211272845686),&
         &(-0.054877473711578625, 0.43568329935571865),&
         &(0.2501413228541527, -1.237538421929958),&
        &(0.6182432935662469, -0.22426788527830935),&
        &(-0.17262350264585732, 0.37739564598170106),&
        &(-2.2239002740099374, 0.1333363608148414),&
        &(-1.263614384970583, 0.8041895097449078),&
        &(0.3587288959713519, -0.057106774383808755),&
        &(-0.011045478465663564, 0.5036079722337261),&
        &(-0.9406491626186084, 1.085769362145687),&
        &(-0.11582532215695436, -0.6909538396968303)/
         call json%initialize()
         call json%create_array(root, '')

         c = m8x8
         call addcase(root, '0', 'n', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '1', 'n', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '2', 'n', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (1._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '3', 'n', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '4', 'n', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '5', 'n', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (1._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '6', 't', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '7', 't', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '8', 't', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '9', 't', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '10', 'n', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (1._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '11', 'n', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '12', 'n', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '13', 'c', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '14', 'c', 'n', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '15', 't', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '16', 't', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '17', 'c', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '18', 'c', 'c', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '19', 'c', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0.2_8, 0.8_8), c, 8)

         c = m8x8
         call addcase(root, '20', 'c', 't', 6, 8, 4, (0.3_8, -0.7_8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '21', 'c', 't', 6, 8, 4, (0._8, 0._8), m8x8, 8, m8x8, 8, (1._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '22', 'c', 't', 6, 8, 4, (0._8, 0._8), m8x8, 8, m8x8, 8, (0.2_8, 0.2_8), c, 8)

         c = m8x8
         call addcase(root, '23', 'c', 't', 6, 8, 4, (0._8, 0._8), m8x8, 8, m8x8, 8, (0._8, 0._8), c, 8)

         c = m8x8
         call addcase(root, '24', 'c', 't', 6, 8, 4, (0._8, 0._8), m8x8, 8, m8x8, 8, (0.2_8, 0._8), c, 8)
         call print(root)
      contains
         subroutine addcase(root, ncase, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external zgemm
            character(len=*) ncase, transa, transb
            integer n, m, lda, ldb, ldc
            complex*16 alpha, beta, a(8, 8), b(8, 8), c(8, 8)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'transa', transa)
            call json%add(case, 'transb', transb)
            call json%add(case, 'm', m)
            call json%add(case, 'n', n)
            call json%add(case, 'k', k)
            call json%create_array(data, 'alpha')
            call json%add(data, '', dreal(alpha))
            call json%add(data, '', dimag(alpha))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'a')
            do j = 1, 8
               do i = 1, 8
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
            call json%create_array(array, 'b')
            do j = 1, 8
               do i = 1, 8
                  call json%create_array(data, '')
                  call json%add(data, '', dreal(b(i, j)))
                  call json%add(data, '', dimag(b(i, j)))
                  call json%add(array, data)
                  nullify (data)
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'ldb', lda)
            call json%create_array(data, 'beta')
            call json%add(data, '', dreal(beta))
            call json%add(data, '', dimag(beta))
            call json%add(case, data)
            nullify (data)
            call json%create_array(array, 'c')
            do j = 1, 8
               do i = 1, 8
                  call json%create_array(data, '')
                  call json%add(data, '', dreal(c(i, j)))
                  call json%add(data, '', dimag(c(i, j)))
                  call json%add(array, data)
                  nullify (data)
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'ldc', ldc)
            call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            call json%create_array(array, 'expect')
            do j = 1, 8
               do i = 1, 8
                  call json%create_array(data, '')
                  call json%add(data, '', dreal(c(i, j)))
                  call json%add(data, '', dimag(c(i, j)))
                  call json%add(array, data)
                  nullify (data)
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
            call json%print(root, './tests/fixtures/level3/complex/gemm.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
