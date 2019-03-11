      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision dd1, dd2, dx1, dy1, dparam(5), dcparam(5)
         data(dcparam(i), i=1, 5)/.0, .0, .0, .0, .0/

         call json%initialize()
         call json%create_array(root, '')

         dparam = dcparam
         dd1 = -4
         dd2 = 2
         dx1 = 3
         dy1 = 9
         call addcase(root, '0', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 4
         dd2 = 0
         dx1 = 3
         dy1 = 9
         call addcase(root, '1', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 1
         dd2 = 2
         dx1 = 3
         dy1 = 1
         call addcase(root, '2', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 2
         dd2 = -1
         dx1 = 3
         dy1 = 8
         call addcase(root, '3', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 2
         dd2 = 1
         dx1 = 3
         dy1 = 8
         call addcase(root, '4', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 5.960464477539063e-8
         dd2 = 2.9802322387695312e-8
         dx1 = 3
         dy1 = 8
         call addcase(root, '5', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 5.960464477539063e-8
         dd2 = 2.9802322387695312e-8
         dx1 = 3
         dy1 = 2
         call addcase(root, '6', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 1.1920928955078125e-7
         dd2 = 5.960464477539063e-8
         dx1 = 3
         dy1 = 2
         call addcase(root, '7', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 16777216
         dd2 = 33554432
         dx1 = 3
         dy1 = 2
         call addcase(root, '8', dd1, dd2, dx1, dy1, dparam)

         dparam = dcparam
         dd1 = 33554432
         dd2 = 16777216
         dx1 = 3
         dy1 = 2
         call addcase(root, '9', dd1, dd2, dx1, dy1, dparam)
         call print(root)
      contains
         subroutine addcase(root, ncase, dd1, dd2, dx1, dy1, dparam)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external drotmg
            character(len=*) ncase
            integer n, incx, incy
            double precision dd1, dd2, dx1, dy1, dparam(5)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'dd1', dd1)
            call json%add(case, 'dd2', dd2)
            call json%add(case, 'dx1', dx1)
            call json%add(case, 'dy1', dy1)
            call json%add(case, 'dparam', dparam)
            call drotmg(dd1, dd2, dx1, dy1, dparam)
            call json%create_object(data, 'expect')
            call json%add(data, 'dd1', dd1)
            call json%add(data, 'dd2', dd2)
            call json%add(data, 'dx1', dx1)
            call json%add(data, 'dy1', dy1)
            call json%add(data, 'dparam', dparam)
            call json%add(case, data)
            nullify (data)
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
            call json%print(root, './tests/fixtures/level1/rotmg.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end

