program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 m6x6(6, 6), a(6, 6), b(6, 6), c(6, 6)
   data((m6x6(i, j), i=1, 6), j=1, 6)/&
   &(1.2629542848807933, 0.9921603654457979),&
   &(-0.3262333607056494, -0.42951310949188126),&
   &(1.3297992629225006, 1.2383041008533804),&
   &(1.2724293214294047, -0.2793462818542693),&
   &(0.4146414344564082, 1.7579030898107073),&
   &(-1.5399500419037095, 0.5607460908880562),&
   &(-0.9285670347135381, -0.4527839725531578),&
   &(-0.2947204467905602, -0.8320432961178319),&
   &(-0.005767172747536955, -1.166570547084707),&
   &(2.404653388857951, -1.0655905803882961),&
   &(0.7635934611404596, -1.563782051071005),&
   &(-0.7990092489893682, 1.1565369971501793),&
   &(-1.1476570092363514, 0.8320471285723897),&
   &(-0.28946157368822334, -0.22732869142475534),&
   &(-0.29921511789731614, 0.2661373616721048),&
   &(-0.411510832795067, -0.3767027185836281),&
   &(0.2522234481561323, 2.4413646288945894),&
   &(-0.8919211272845686, -0.7953391172553718),&
   &(0.43568329935571865, -0.054877473711578625),&
   &(-1.237538421929958, 0.2501413228541527),&
   &(-0.22426788527830935, 0.6182432935662469),&
   &(0.37739564598170106, -0.17262350264585732),&
   &(0.1333363608148414, -2.2239002740099374),&
   &(0.8041895097449078, -1.263614384970583),&
   &(-0.057106774383808755, 0.3587288959713519),&
   &(0.5036079722337261, -0.011045478465663564),&
   &(1.085769362145687, -0.9406491626186084),&
   &(-0.6909538396968303, -0.11582532215695436),&
   &(-1.2845993538721883, -0.8149687088699175),&
   &(0.04672617218835198, 0.24226348085968588),&
   &(-0.23570655643950122, -1.4250983947324998),&
   &(-0.5428882550102544, 0.36594112304921983),&
   &(-0.4333103174567822, 0.2484126488725964),&
   &(-0.6494716467962331, 0.06528818167162072),&
   &(0.726750747385451, 0.01915639166027384),&
   &(1.1519117540872, 0.2573383771555333)/
   call json%initialize()
   call json%create_array(root, '')

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '0', 'u', 'n', 4, 6, 0._8, a, 6, 1._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '1', 'u', 'n', 4, 6, 0._8, a, 6, 0._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '2', 'l', 'n', 4, 6, 0._8, a, 6, 0._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '3', 'l', 'n', 4, 6, 0._8, a, 6, 0.5_8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '4', 'u', 'n', 4, 6, 0.2_8, a, 6, 0.5_8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '5', 'u', 'n', 4, 6, 0.2_8, a, 6, 0._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '6', 'u', 'n', 4, 6, 0.2_8, a, 6, 1._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '7', 'l', 'n', 4, 6, 0.2_8, a, 6, 1._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '8', 'u', 'c', 4, 6, 0.2_8, a, 6, 0._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '9', 'u', 'c', 4, 6, 0.2_8, a, 6, 0.2_8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '10', 'l', 'c', 4, 6, 0.2_8, a, 6, 0._8, c, 6)

   call copy(m6x6, a, 6, 6)
   call copy(m6x6, c, 6, 6)
   a(1, 1) = (0._8, 0._8)
   call addcase(root, '11', 'l', 'c', 4, 6, 0.2_8, a, 6, 0.2_8, c, 6)

   call print(root)
contains
   subroutine addcase(root, ncase, uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external zherk
      character uplo, trans
      character(len=*) ncase
      integer n, m, lda, ldc
      double precision alpha, beta
      complex*16 a(lda, n), c(ldc, n)
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'uplo', uplo)
      call json%add(case, 'trans', trans)
      call json%add(case, 'n', n)
      call json%add(case, 'k', k)
      call json%add(case, 'alpha', alpha)
      call json%create_array(array, 'a')
      do j = 1, 6
         do i = 1, 6
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
      call json%add(case, 'beta', beta)
      call json%create_array(array, 'c')
      do j = 1, 6
         do i = 1, 6
            call json%create_array(data, '')
            call json%add(data, '', dreal(c(i, j)))
            call json%add(data, '', dimag(c(i, j)))
            call json%add(array, data)
            nullify (data)
         enddo
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'ldc', ldc)
      call zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
      call json%create_array(array, 'expect')
      do j = 1, 6
         do i = 1, 6
            call json%create_array(data, '')
            call json%add(data, '', dreal(c(i, j)))
            call json%add(data, '', dimag(c(i, j)))
            call json%add(array, data)
            nullify (data)
         enddo
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
      call json%print(root, './tests/fixtures/level3/complex/herk.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
   subroutine copy(x, a, n, lda)
      integer n, lda
      complex*16 x(lda, n), a(lda, n)
      do j = 1, n
         do i = 1, lda
            a(i, j) = x(i, j)
         enddo
      enddo
   end subroutine copy
end
