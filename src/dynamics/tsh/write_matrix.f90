!
!...write matrix into a file
module write_matrix_2file
implicit none

interface write_matrix
  module procedure write_real_v1,write_real_v2,write_v_mat,write_1real,write_2real,&
                   write_real_matrix
end interface

interface print_matrix
  module procedure print_real_v,print_real_mat,print_complex_mat
end interface

CONTAINS

!...write a real 1d array
subroutine write_real_v1(nrow,x,filename)
implicit none
integer,intent(in)::nrow
real*8,intent(in)::x(nrow)
character(len=*),intent(in)::filename
logical::od
integer::i

inquire(unit=99,opened=od)
if(od)then
  write(*,*)"unit=99 has been used."
  stop
else
  open(unit=99,file=filename)
  do i=1,nrow
     write(99,'(f16.8)')x(i)
  enddo
  close(99)
endif

end subroutine write_real_v1
!
!...write two real 1d arraies
subroutine write_real_v2(nrow,x,y,filename)
implicit none
integer,intent(in)::nrow
real*8,intent(in)::x(nrow),y(nrow)
character(len=*),intent(in)::filename
logical::od
integer::i

inquire(unit=99,opened=od)
if(od)then
  write(*,*)"unit=99 has been used."
  stop
else
  open(unit=99,file=filename)
  do i=1,nrow
     write(99,'(f16.8,1x,f16.8)')x(i),y(i)
  enddo
  close(99)
endif

end subroutine write_real_v2
!
!write 1d arraies with a square matrix A(ns,ns,nrow)
subroutine write_v_mat(nrow,x,ns,A,filename)
implicit none
integer,intent(in)::nrow,ns
real*8,intent(in)::x(nrow),A(ns,ns,nrow)
character(len=*),intent(in)::filename
logical::od
integer::i,j,k

inquire(unit=99,opened=od)
if(od)then
  write(*,*)"unit=99 has been used."
  stop
else
  open(unit=99,file=filename)
  do i=1,nrow
     write(99,'(1000(f16.8,1x))')x(i),((A(j,k,i),k=1,ns),j=1,ns)
  enddo
  close(99)
endif

end subroutine write_v_mat
!
!...write 1 real numbers to a file
subroutine write_1real(y,filename)
implicit none
real*8,intent(in)::y
character(len=*),intent(in)::filename
logical::od

inquire(unit=99,opened=od)
if(od)then
  write(*,*)"unit=99 has been used."
  stop
else
  open(unit=99,file=filename)
  write(99,'(e23.12)')y
  close(99)
endif

end subroutine write_1real
!
!...write 2 real numbers to a file
subroutine write_2real(x,y,filename)
implicit none
real*8,intent(in)::x,y
character(len=*),intent(in)::filename
logical::od

inquire(unit=99,opened=od)
if(od)then
  write(*,*)"unit=99 has been used."
  stop
else
  open(unit=99,file=filename)
  write(99,'(e23.12,1x,e23.12)')x,y
  close(99)
endif

end subroutine write_2real
!...write a matrix into a file
subroutine write_real_matrix(nrow,ncol,Arr,filename)
integer,intent(in)::nrow,ncol
character(len=*),intent(in)::filename
real*8,intent(in)::Arr(nrow,ncol)
integer i,io
logical::alive
character(len=80)::cha,fmt1

call i2char(ncol,cha)
!write(*,*)"cha=",cha

!fmt1="(1x,"//trim(adjustl(cha))//"(e23.13,1x))"
fmt1="(1x,"//trim(adjustl(cha))//"(f16.8,1x))"
write(*,*)"fmt1=",fmt1
!write(*,fmt1)arr(1,:)

open(unit=1,file=filename,iostat=io)
if (io/=0) then
    write(*,*) 'Could not open file "',trim(filename),'"!'
    stop 1
endif

do i=1,nrow
   write(1,fmt=fmt1)Arr(i,:)
enddo
close(1)

end subroutine write_real_matrix
!============================================================================
!
!...print a real 1d array on the screen
subroutine print_real_v(vname,n,v)
implicit none
character(len=*),intent(in)::vname
integer,intent(in)::n
real*8,intent(in)::v(n)
integer::i

write(*,*)trim(adjustl(vname)),"="
do i=1,n
   write(*,*)v(i)
enddo

end subroutine print_real_v

!print a real matrix on the screen
subroutine print_real_mat(vname,n1,n2,v)
implicit none
character(len=*),intent(in)::vname
integer,intent(in)::n1,n2
real*8,intent(in)::v(n1,n2)
integer::i

write(*,*)trim(adjustl(vname)),"="
do i=1,n1
   write(*,'(1000(f16.8,1x))')v(i,:)
enddo

end subroutine print_real_mat
!print a real matrix on the screen
subroutine print_complex_mat(vname,n1,n2,v)
implicit none
character(len=*),intent(in)::vname
integer,intent(in)::n1,n2
complex*16,intent(in)::v(n1,n2)
integer::i,j

write(*,*)"v=",v
write(*,*)trim(adjustl(vname)),"="
do i=1,n1
   write(*,'(1000(e23.12e3,1x,e23.12e3))')v(i,:)
enddo

end subroutine print_complex_mat
!============================================================================



end module write_matrix_2file
