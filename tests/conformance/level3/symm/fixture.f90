      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision m6x6(6, 6), c(6, 6)
         data((m6x6(i, j), i=1, 6), j=1, 6)/&
         &1.2629542848807933,&
         &-0.3262333607056494,&
         &1.3297992629225006,&
         &1.2724293214294047,&
         &0.4146414344564082,&
         &-1.5399500419037095,&
         &-0.9285670347135381,&
         &-0.2947204467905602,&
         &-0.005767172747536955,&
         &2.404653388857951,&
         &0.7635934611404596,&
         &-0.7990092489893682,&
         &-1.1476570092363514,&
         &-0.28946157368822334,&
         &-0.29921511789731614,&
         &-0.411510832795067,&
         &0.2522234481561323,&
         &-0.8919211272845686,&
         &0.43568329935571865,&
         &-1.237538421929958,&
         &-0.22426788527830935,&
         &0.37739564598170106,&
         &0.1333363608148414,&
         &0.8041895097449078,&
         &-0.057106774383808755,&
         &0.5036079722337261,&
         &1.085769362145687,&
         &-0.6909538396968303,&
         &-1.2845993538721883,&
         &0.04672617218835198,&
         &-0.23570655643950122,&
         &-0.5428882550102544,&
         &-0.4333103174567822,&
         &-0.6494716467962331,&
         &0.726750747385451,&
         &1.1519117540872/
         call json%initialize()
         call json%create_array(root, '')

         c = m6x6
         call addcase(root, '0', 'l', 'u', 6, 0, 0.3_8, m6x6, 6, m6x6, 6, 0.2_8, c, 6)

         c = m6x6
         call addcase(root, '1', 'l', 'u', 0, 0, 0.3_8, m6x6, 6, m6x6, 6, 0.2_8, c, 6)

         c = m6x6
         call addcase(root, '3', 'l', 'u', 6, 4, 0._8, m6x6, 6, m6x6, 6, 1._8, c, 6)

         c = m6x6
         call addcase(root, '4', 'l', 'u', 6, 4, 0._8, m6x6, 6, m6x6, 6, 0._8, c, 6)

         c = m6x6
         call addcase(root, '5', 'l', 'u', 4, 6, 0._8, m6x6, 6, m6x6, 6, 0.2_8, c, 6)

         c = m6x6
         call addcase(root, '6', 'l', 'u', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 0.2_8, c, 6)

         c = m6x6
         call addcase(root, '7', 'l', 'u', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 0._8, c, 6)

         c = m6x6
         call addcase(root, '8', 'l', 'l', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 0._8, c, 6)

         c = m6x6
         call addcase(root, '9', 'l', 'l', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 0.2_8, c, 6)

         c = m6x6
         call addcase(root, '10', 'r', 'u', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 0._8, c, 6)

         c = m6x6
         call addcase(root, '11', 'r', 'u', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 1._8, c, 6)

         c = m6x6
         call addcase(root, '12', 'r', 'l', 4, 6, 0.3_8, m6x6, 6, m6x6, 6, 1._8, c, 6)

         call print(root)
      contains
         subroutine addcase(root, ncase, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array
            external dsymm
            character side, uplo
            character(len=*) ncase
            integer n, m, lda, ldb, ldc
            double precision alpha, beta, a(6, 6), b(6, 6), c(6, 6)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'side', side)
            call json%add(case, 'uplo', uplo)
            call json%add(case, 'm', m)
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
            call json%create_array(array, 'b')
            do j = 1, 6
               do i = 1, 6
                  call json%add(array, '', b(i, j))
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'ldb', ldb)
            call json%add(case, 'beta', beta)
            call json%create_array(array, 'c')
            do j = 1, 6
               do i = 1, 6
                  call json%add(array, '', c(i, j))
               enddo
            enddo
            call json%add(case, array)
            nullify (array)
            call json%add(case, 'ldc', ldc)
            call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            call json%create_array(array, 'expect')
            do j = 1, 6
               do i = 1, 6
                  call json%add(array, '', c(i, j))
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
            call json%print(root, './tests/fixtures/level3/symm.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
