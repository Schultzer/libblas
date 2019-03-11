program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(4), y(4), x2(1), y2(1)
   call json%initialize()
   call json%create_array(root, '')

   x = (/(1, 5), (2, 6), (3, 7), (4, 8)/)
   y = (/(5, 9), (6, 10), (7, 11), (8, 12)/)
   call addcase(root, '0', 4, x, 1, y, 1)

   x = (/(1, 5), (2, 6), (3, 7), (4, 8)/)
   y = (/(5, 9), (6, 10), (7, 11), (8, 12)/)
   call addcase(root, '1', 4, x, -1, y, -1)

   x = (/(1, 5), (2, 6), (3, 7), (4, 8)/)
   y = (/(5, 9), (6, 10), (7, 11), (8, 12)/)
   call addcase(root, '2', 0, x, -1, y, -1)

   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx, y, incy)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external zdotu
      character(len=*) ncase
      integer n, incx, incy
      complex*16 x(4), y(4), dot, zdotu
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(array, 'x')
      do i = 1, 4
         call json%create_array(data, '')
         call json%add(data, '', dreal(x(i)))
         call json%add(data, '', dimag(x(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incx', incx)
      call json%create_array(array, 'y')
      do i = 1, 4
         call json%create_array(data, '')
         call json%add(data, '', dreal(y(i)))
         call json%add(data, '', dimag(y(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incy', incy)
      dot = zdotu(n, x, incx, y, incy)
      call json%create_array(data, 'expect')
      call json%add(data, '', dreal(dot))
      call json%add(data, '', dimag(dot))
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
      call json%print(root, './tests/fixtures/level1/complex/dotu.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
