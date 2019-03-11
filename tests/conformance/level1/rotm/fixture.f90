      program fixture
         use, intrinsic :: iso_fortran_env, only: wp => real64
         use json_module
         type(json_core) :: json
         type(json_value), pointer :: root
         double precision dx(4), dy(4), dparam(5), nan
         call json%initialize()
         call json%create_array(root, '')

         !FIXME fortran json produce incorrect json, in json NaN == "null", serde-json does not deserialize "null" into NaN [https://github.com/serde-rs/json/issues/202

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! nan = 0
         ! nan = 0/nan
         ! dparam(1) = -2
         ! dparam(2) = nan
         ! dparam(3) = nan
         ! dparam(4) = nan
         ! dparam(5) = nan
         ! call addcase(root, '0', 4, dx, 1, dy, 1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = -2
         ! dparam(2) = nan
         ! dparam(3) = nan
         ! dparam(4) = nan
         ! dparam(5) = nan
         ! call addcase(root, '1', 0, dx, 1, dy, 1, dparam)

         dx(1) = 1
         dx(2) = 2
         dx(3) = 3
         dx(4) = 4
         dy(1) = 5
         dy(2) = 6
         dy(3) = 7
         dy(4) = 8
         dparam(1) = -1
         dparam(2) = 2
         dparam(3) = 3
         dparam(4) = 4
         dparam(5) = 5
         call addcase(root, '2', 4, dx, 1, dy, 1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = 0
         ! dparam(2) = nan
         ! dparam(3) = 3
         ! dparam(4) = 4
         ! dparam(5) = nan
         ! call addcase(root, '3', 4, dx, 1, dy, 1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = 1
         ! dparam(2) = 2
         ! dparam(3) = nan
         ! dparam(4) = nan
         ! dparam(5) = 3
         ! call addcase(root, '4', 4, dx, 1, dy, 1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = 1
         ! dparam(2) = 2
         ! dparam(3) = nan
         ! dparam(4) = nan
         ! dparam(5) = 3
         ! call addcase(root, '5', 4, dx, 1, dy, -1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = 1
         ! dparam(2) = 2
         ! dparam(3) = nan
         ! dparam(4) = nan
         ! dparam(5) = 3
         ! call addcase(root, '6', 4, dx, 1, dy, -1, dparam)

         ! dx(1) = 1
         ! dx(2) = 2
         ! dx(3) = 3
         ! dx(4) = 4
         ! dy(1) = 5
         ! dy(2) = 6
         ! dy(3) = 7
         ! dy(4) = 8
         ! dparam(1) = 0
         ! dparam(2) = nan
         ! dparam(3) = 4
         ! dparam(4) = 5
         ! dparam(5) = nan
         ! call addcase(root, '7', 2, dx, 2, dy, -2, dparam)

         dx(1) = 1
         dx(2) = 2
         dx(3) = 3
         dx(4) = 4
         dy(1) = 5
         dy(2) = 6
         dy(3) = 7
         dy(4) = 8
         dparam(1) = -1
         dparam(2) = 1
         dparam(3) = 2
         dparam(4) = 3
         dparam(5) = 4
         call addcase(root, '8', 2, dx, 2, dy, -2, dparam)

         call print(root)
      contains
         subroutine addcase(root, ncase, n, dx, incx, dy, incy, dparam)
            type(json_core) :: json
            type(json_value), pointer :: root, case, array, data
            external drotm
            character(len=*) ncase
            integer n, incx, incy
            double precision dx(4), dy(4), dparam(5)
            call json%create_object(case, '')
            call json%add(root, case)
            call json%add(case, 'n', n)
            call json%add(case, 'dx', dx)
            call json%add(case, 'incx', incx)
            call json%add(case, 'dy', dy)
            call json%add(case, 'incy', incy)
            call json%add(case, 'dparam', dparam)
            call drotm(n, dx, incx, dy, incy, dparam)
            call json%create_object(data, 'expect')
            call json%add(data, 'dx', dx)
            call json%add(data, 'dy', dy)
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
            call json%print(root, './tests/fixtures/level1/rotm.json')
            call json%destroy(root)
            if (json%failed()) then
               call json%check_for_errors(status_ok, error_msg)
               write (*, *) 'Error: '//error_msg
               call json%clear_exceptions()
               call json%destroy(root)
            end if
         end subroutine print
      end
