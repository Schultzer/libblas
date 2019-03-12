        program fixture
           use, intrinsic :: iso_fortran_env, only: wp => real64
           use json_module
           type(json_core) :: json
           type(json_value), pointer :: root
           complex*16 m6x6(6, 6), a(6, 6), b(6, 6)
           data((m6x6(i, j), i=1, 6), j=1, 6)/&
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
           &(1.1519117540872, 0.2573383771555333)/

           call json%initialize()
           call json%create_array(root, '')

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           call addcase(root, '0', 'l', 'u', 'n', 'n', 0, 6, (0.2_8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           call addcase(root, '1', 'l', 'u', 'n', 'n', 4, 0, (0.2_8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           call addcase(root, '2', 'l', 'u', 'n', 'n', 4, 6, (0._8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           call addcase(root, '3', 'l', 'u', 'n', 'n', 4, 6, (0.2_8, 0.6_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '4', 'l', 'u', 'n', 'u', 4, 6, (1._8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '5', 'l', 'l', 'n', 'n', 4, 6, (1._8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '6', 'l', 'l', 'n', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '7', 'l', 'u', 't', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '8', 'l', 'u', 'c', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '9', 'l', 'l', 'c', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '10', 'l', 'l', 'c', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           b(1, 1) = 0
           call addcase(root, '11', 'l', 'l', 't', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 2) = 0
           call addcase(root, '12', 'r', 'u', 'n', 'n', 4, 6, (1._8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 2) = 0
           call addcase(root, '13', 'r', 'u', 'n', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(5, 4) = 0
           call addcase(root, '14', 'r', 'l', 'n', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(5, 4) = 0
           call addcase(root, '15', 'r', 'l', 'n', 'u', 4, 6, (1.0_8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 6) = 0
           call addcase(root, '16', 'r', 'u', 't', 'n', 4, 6, (1.0_8, 0._8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 6) = 0
           call addcase(root, '17', 'r', 'u', 'c', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 6) = 0
           call addcase(root, '18', 'r', 'u', 'c', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(1, 6) = 0
           call addcase(root, '19', 'r', 'l', 'c', 'u', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           a(2, 1) = 0
           call addcase(root, '20', 'r', 'l', 'c', 'n', 4, 6, (0.5_8, 0.5_8), a, 6, b, 6)

           call copy(m6x6, a, 6, 6)
           call copy(m6x6, b, 6, 6)
           call addcase(root, '21', 'r', 'l', 't', 'n', 4, 6, (1._8, 0._8), a, 6, b, 6)
           call print(root)
        contains
           subroutine addcase(root, ncase, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
              type(json_core) :: json
              type(json_value), pointer :: root, case, array, data
              external ztrsm
              character side, uplo, transa, diag
              character(len=*) ncase
              integer n, m, lda, ldb
              complex*16 alpha, a(lda, n), b(ldb, n)
              call json%create_object(case, '')
              call json%add(root, case)
              call json%add(case, 'side', side)
              call json%add(case, 'uplo', uplo)
              call json%add(case, 'transa', transa)
              call json%add(case, 'diag', diag)
              call json%add(case, 'm', m)
              call json%add(case, 'n', n)
              call json%create_array(data, 'alpha')
              call json%add(data, '', dreal(alpha))
              call json%add(data, '', dimag(alpha))
              call json%add(case, data)
              nullify (data)
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
              call json%create_array(array, 'b')
              do j = 1, 6
                 do i = 1, 6
                    call json%create_array(data, '')
                    call json%add(data, '', dreal(b(i, j)))
                    call json%add(data, '', dimag(b(i, j)))
                    call json%add(array, data)
                    nullify (data)
                 enddo
              enddo
              call json%add(case, array)
              nullify (array)
              call json%add(case, 'ldb', ldb)
              call ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
              call json%create_array(array, 'expect')
              do j = 1, 6
                 do i = 1, 6
                    call json%create_array(data, '')
                    call json%add(data, '', dreal(b(i, j)))
                    call json%add(data, '', dimag(b(i, j)))
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
              call json%print(root, './tests/fixtures/level3/complex/trsm.json')
              call json%destroy(root)
              if (json%failed()) then
                 call json%check_for_errors(status_ok, error_msg)
                 write (*, *) 'Error: '//error_msg
                 call json%clear_exceptions()
                 call json%destroy(root)
              end if
           end subroutine print
           subroutine copy(x, a, n, lda)
              integer n, lda
              complex*16 x(lda, n), a(lda, n)
              do j = 1, n
                 do i = 1, lda
                    a(i, j) = x(i, j)
                 enddo
              enddo
           end subroutine copy
        end
