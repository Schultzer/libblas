program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision a, b, c, s
   call json%initialize()
   call json%create_array(root, '')

   a = 4
   b = 2
   c = 0
   s = 0
   call addcase(root, '0', a, b, c, s)

   a = 0
   b = 0
   c = 0
   s = 0
   call addcase(root, '1', a, b, c, s)

   a = -4
   b = -2
   c = 0
   s = 0
   call addcase(root, '2', a, b, c, s)
   call print(root)
contains
   subroutine addcase(root, ncase, a, b, c, s)
      type(json_core) :: json
      type(json_value), pointer :: root, case, data
      external drotg
      character(len=*) ncase
      double precision a, b, c, s
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'a', a)
      call json%add(case, 'b', b)
      call json%add(case, 'c', c)
      call json%add(case, 's', s)
      call drotg(a, b, c, s)
      call json%create_object(data, 'expect')
      call json%add(data, 'a', a)
      call json%add(data, 'b', b)
      call json%add(data, 'c', c)
      call json%add(data, 's', s)
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
      call json%print(root, './tests/fixtures/level1/rotg.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
