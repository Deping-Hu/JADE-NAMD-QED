!...Kronecker delta function
function Kdelta(i,j) result(retval)
implicit none
integer,intent(in)::i,j
real*8::retval

if(i==j)then
   retval=1.d0
else
   retval=0.d0
endif
end function Kdelta

subroutine diag_real(A,n,e,v)
  implicit none
  !Arguments(input):
  real(8),intent(in)::   A(n,n)          !Symmetric matrix
  integer,intent(in)::   n               !Eigenvalue problem size
  !Arguments(output):
  real(8),intent(out)::   e(n)   !Eigenvalues
  real(8),intent(out)::   v(n,n) !Corresponding eigenvectors

  !local variables
  integer,allocatable,dimension(:)   :: IWORK
  character(1)                       :: JOBZ, RANGE, UPLO
  real(8)                            :: VL, VU, ABSTOL
  real(8),allocatable,dimension(:)   :: WORK
  integer                            :: LDA, IL, IU, M, LDZ, NB, LWORK,INFO
  integer,allocatable,dimension(:)   :: IFAIL
  integer                            :: ILAENV
  real(8),allocatable:: Atemp(:,:)   !=A
  !write(*,*)"Hello diag_real"
  !write(*,*)"A=",A
  allocate(Atemp(n,n))
  Atemp=A

  JOBZ='V'       !'V' means compute both eigenvalues and eigenvectors
  RANGE='I'      !'I' means only selected eigenvalues/eigenvectors will be found
  UPLO='U'       !'U' means the upper triangle of A is provided.
  LDA=n          !Leading dimension of A is N
  VL=0.0D0       !VL can be set to anything because it is not used when
  RANGE='I'
  VU=0.0D0       !VU can be set to anything because it is not used when
  RANGE='I'
  IL=1           !The smallest eigenvalue number
  IU=n        !The largest eigenvalue number.
  ABSTOL=2*tiny(ABSTOL) !Set tolerance that yields highest possible accuracy in the
  !calculations; equivalent to  ABSTOL=2*DLAMCH('S')
  LDZ=n          !Leading dimension of v is N
  NB=ILAENV(1,'DSYTRD','VIU',n,n,n,n) !Determine the optimal block size
  LWORK=(NB+3)*N !Set the size of array WORK
  !write(*,*)"NB,LWORK=",NB,LWORK
  LWORK=8*n   !for mkl
  allocate(WORK(LWORK)) !Allocate array WORK, which is a real(8) work array used by DSYGVX
  allocate(IWORK(5*N))  !Allocate array IWORK, which is an integer work array used by DSYGVX
  allocate(IFAIL(N)) !Allocate IFAIL, an integer flag array used by DSYGVX

  call DSYEVX(JOBZ, RANGE, UPLO, n, Atemp, LDA, VL, VU, IL, IU, &
&   ABSTOL, M, e, v, LDZ, WORK, LWORK, IWORK, &
&   IFAIL, INFO )

  if (INFO/=0) then
    write(6,*) 'Error in diag4: subroutine DSYEVX failed with INFO=',INFO
    stop
  endif

  deallocate(IFAIL)
  deallocate(IWORK)
  deallocate(WORK)
  deallocate(Atemp)
  !write(*,*)"Bye diag_real"
end subroutine diag_real

!Similarity Transformation of matrix
subroutine Simil_Trans_real(mode,n,U,A)
implicit none
character(len=*),intent(in)::mode
integer,intent(in)::n
real*8,intent(in)::U(n,n)
real*8,intent(inout)::A(n,n)
real*8,allocatable::Ut(:,:)
allocate(Ut(n,n))

!write(*,*)"Warning: If U has the phase factor, after rotation"
!write(*,*)"The matrix A also has the phase factor."

Ut=transpose(U)
if(trim(adjustl(mode))=="UtAU")then
  A=matmul(Ut,A)
  A=matmul(A,U)
else if (trim(adjustl(mode))=="UAUt")then
  A=matmul(U,A)
  A=matmul(A,Ut)
else if (trim(adjustl(mode))=="UAU")then
  A=matmul(U,A)
  A=matmul(A,U)
else if (trim(adjustl(mode))=="UtAUt")then
  A=matmul(Ut,A)
  A=matmul(A,Ut)
else 
  stop "mode we do not know."
endif

deallocate(Ut)
end subroutine Simil_Trans_real
!----------------------------------------------------------------------
subroutine Simil_Trans_complex(mode,n,U,A)
implicit none
integer,intent(in)::n
complex(kind=8),intent(in)::U(n,n)
complex(kind=8),intent(inout)::A(n,n)
character(len=*),intent(in)::mode
complex(kind=8),allocatable::Ut(:,:)
allocate(Ut(n,n))

!write(*,*)"Warning: If U has the phase factor, after rotation"
!write(*,*)"The matrix A also has the phase factor."

Ut=transpose(conjg(U))  !U^t=conjugate transpose of U
if(trim(adjustl(mode))=="UtAU")then
  A=matmul(Ut,A)
  A=matmul(A,U)
else if (trim(adjustl(mode))=="UAUt")then
  A=matmul(U,A)
  A=matmul(A,Ut)
else if (trim(adjustl(mode))=="UAU")then
  A=matmul(U,A)
  A=matmul(A,U)
else if (trim(adjustl(mode))=="UtAUt")then
  A=matmul(Ut,A)
  A=matmul(A,Ut)
else
  stop "mode we do not know."
endif

deallocate(Ut)
end subroutine Simil_Trans_complex
