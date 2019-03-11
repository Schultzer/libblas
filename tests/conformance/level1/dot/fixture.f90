program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision x(10), y(10)
   call json%initialize()
   call json%create_array(root, '')

   x = (/1, 2, 3, 4, 0, 0, 0, 0, 0, 0/)
   y = (/5, 6, 7, 8, 0, 0, 0, 0, 0, 0/)
   call addcase(root, '0', 4, x, 1, y, 1)

   x = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 2/)
   y = (/5, 6, 7, 8, 5, 4, 3, 3, -1, -9/)
   call addcase(root, '1', 10, x, 1, y, 1)

   x = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 2/)
   y = (/5, 6, 7, 8, 5, 4, 3, 3, -1, -9/)
   call addcase(root, '2', 0, x, 1, y, 1)

   x = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 2/)
   y = (/5, 6, 7, 8, 5, 4, 3, 3, -1, -9/)
   call addcase(root, '3', 5, x, 2, y, 2)

   x = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 2/)
   y = (/5, 6, 7, 8, 5, 4, 3, 3, -1, -9/)
   call addcase(root, '4', 5, x, -2, y, -2)

   x = (/1, 2, 3, 4, 5, 6, 7, 0, 0, 0/)
   y = (/5, 6, 7, 8, 5, 4, 3, 0, 0, 0/)
   call addcase(root, '5', 7, x, 1, y, 1)

   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx, y, incy)
      type(json_core) :: json
      type(json_value), pointer :: root, case
      external ddot
      character(len=*) ncase
      integer n, incx, incy
      double precision x(10), y(10), dot, ddot
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'x', x)
      call json%add(case, 'incx', incx)
      call json%add(case, 'y', y)
      call json%add(case, 'incy', incy)
      dot = ddot(n, x, incx, y, incy)
      call json%add(case, 'expect', dot)
      nullify (case)
   end subroutine addcase

   subroutine print(root)
      use, intrinsic :: iso_fortran_env, only: real64
      use json_module, CK => json_CK
      implicit none

      type(json_core) :: json
      type(json_value), pointer :: root
      logical :: status_ok
      character(kind=CK, len=:), allocatable :: error_msg
      call json%print(root, './tests/fixtures/level1/dot.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
