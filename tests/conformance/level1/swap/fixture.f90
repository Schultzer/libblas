program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision x(4), y(4)
   call json%initialize()
   call json%create_array(root, '')

   x = (/1, 2, 3, 4/)
   y = (/1, 2, 3, 4/)
   call addcase(root, '0', 4, x, 1, y, -1)

   x = (/1, 2, 3, 4/)
   y = (/1, 2, 3, 4/)
   call addcase(root, '1', 4, x, -1, y, 1)

   x = (/1, 2, 3, 4/)
   y = (/7, 8, 9, 5/)
   call addcase(root, '2', 4, x, 1, y, 1)

   x = (/1, 2, 3, 0/)
   y = (/7, 8, 9, 0/)
   call addcase(root, '3', 3, x, 1, y, 1)

   x = (/1, 2, 0, 0/)
   y = (/7, 8, 0, 0/)
   call addcase(root, '4', 2, x, 1, y, 1)

   x = (/1, 2, 0, 0/)
   y = (/7, 8, 0, 0/)
   call addcase(root, '5', 0, x, 1, y, 1)

   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx, y, incy)
      type(json_core) :: json
      type(json_value), pointer :: root, case, data
      external dswap
      character(len=*) ncase
      integer n, incx
      double precision x(4), y(4)
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'x', x)
      call json%add(case, 'incx', incx)
      call json%add(case, 'y', y)
      call json%add(case, 'incy', incy)
      call dswap(n, x, incx, y, incy)
      call json%create_object(data, 'expect')
      call json%add(data, 'x', x)
      call json%add(data, 'y', y)
      call json%add(case, data)
      nullify (data)
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
      call json%print(root, './tests/fixtures/level1/swap.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
