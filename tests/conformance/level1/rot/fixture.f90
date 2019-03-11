program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision x(2), y(2), pi
   call json%initialize()
   call json%create_array(root, '')
   pi = 4.d0*datan(1.d0)

   x = (/1, 0/)
   y = (/0, 1/)
   call addcase(root, '0', 2, x, 1, y, 1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))

   x = (/1, 0/)
   y = (/0, 1/)
   call addcase(root, '0', 2, x, -1, y, -1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))

   x = (/1, 0/)
   y = (/0, 1/)
   call addcase(root, '0', 0, x, -1, y, -1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))
   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx, y, incy, c, s)
      type(json_core) :: json
      type(json_value), pointer :: root, case, data
      external drot
      character(len=*) ncase
      integer n, incx, incy
      double precision x(2), y(2), c, s
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'x', x)
      call json%add(case, 'incx', incx)
      call json%add(case, 'y', y)
      call json%add(case, 'incy', incy)
      call json%add(case, 'c', c)
      call json%add(case, 's', s)
      call drot(n, x, incx, y, incy, c, s)
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
      call json%print(root, './tests/fixtures/level1/rot.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
