program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 a(8, 8)
   data((a(i, j), i=1, 8), j=1, 8)/&
   &(1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),&
   &(0, 0), (1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),&
   &(0, 0), (0, 0), (1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0),&
   &(0, 0), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0), (0, 0), (0, 0),&
   &(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0), (0, 0),&
   &(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0),&
   &(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 0),&
   &(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0)/
   call json%initialize()
   call json%create_array(root, '')

   call addcase(root, 'pack', 8, a, 8)

   call print(root)
contains
   subroutine addcase(root, ncase, n, a, lda)
      use, intrinsic :: iso_fortran_env, only: wp => real64
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      character(len=*) ncase
      integer n, lda
      complex*16 a(lda, n), c(8)
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(array, 'mat')
      do j = 1, 8
         do i = 1, 8
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
      c = pack(a, a /= 0)
      call json%create_array(array, 'expect')
      do i = 1, lda
         call json%create_array(data, '')
         call json%add(data, '', dreal(c(i)))
         call json%add(data, '', dimag(c(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
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
      call json%print(root, './tests/fixtures/matrix/pack_test.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
