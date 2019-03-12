program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(2), y(2)
   double precision pi
   call json%initialize()
   call json%create_array(root, '')
   pi = 4.d0*datan(1.d0)

   x = (/(1, 0), (0, 0)/)
   y = (/(0, 1), (0, 0)/)
   call addcase(root, '0', 2, x, 1, y, 1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))

   x = (/(1, 0), (0, 0)/)
   y = (/(0, 1), (0, 0)/)
   call addcase(root, '0', 2, x, -1, y, -1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))

   x = (/(1, 0), (0, 0)/)
   y = (/(0, 1), (0, 0)/)
   call addcase(root, '0', 0, x, -1, y, -1, cos(pi*1.d0/6.d0), sin(pi*1.d0/6.d0))
   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx, y, incy, c, s)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data, d
      external zdrot
      character(len=*) ncase
      integer n, incx, incy
      complex*16 x(2), y(2)
      double precision c, s
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(array, 'x')
      do i = 1, 2
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
      call json%add(case, 'c', c)
      call json%add(case, 's', s)
      do i = 1, 2
         call json%create_array(data, '')
         call json%add(data, '', dreal(y(i)))
         call json%add(data, '', dimag(y(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incy', incy)
      call zdrot(n, x, incx, y, incy, c, s)
      call json%create_object(data, 'expect')
      call json%create_array(array, 'x')
      do i = 1, 2
         call json%create_array(d, '')
         call json%add(d, '', dreal(x(i)))
         call json%add(d, '', dimag(x(i)))
         call json%add(array, d)
         nullify (d)
      enddo
      call json%add(data, array)
      nullify (array)
      call json%create_array(array, 'y')
      do i = 1, 2
         call json%create_array(d, '')
         call json%add(d, '', dreal(y(i)))
         call json%add(d, '', dimag(y(i)))
         call json%add(array, d)
         nullify (d)
      enddo
      call json%add(data, array)
      nullify (array)
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
      call json%print(root, './tests/fixtures/level1/complex/rot.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
