!...Some old subroutines used to deal with strings
!-----------------------------------------------------------------------
      logical function match(str1,str2,m)
      implicit none
      character*(*) str1,str2
      integer m,i

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      i=len(str1)+1
      if(str1.eq.'#')then
         if (str1.eq.str2(1:len(str1))) match = .TRUE. !old
      else
         if(str1.eq.str2(1:len(str1)))then
           if(str2(i:i).eq.' '.or.ichar(str2(i:i)).eq.9) then
             match=.true.
           endif
         endif
      endif
      return
      end function match
!-----------------------------------------------------------------------
      integer function intread(str,m)
      implicit none
      character*(*) str
      integer m

      integer n1,n2,nmax

      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.&
               (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' ' &
               .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) intread
      m = n2

      return
      end function intread
! -----------------------------------------------------------------------
! return the real*8 that starts at loc m in str
! skip initial spaces and tabs
! scans to end of number marked by space or tab
! also return m = loc of char after end of number in str
      double precision function realread(str,m)
      implicit none
      character*(*) str
      integer m

      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.&
           (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' ' &
           .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      read (str(n1:n2-1),*) realread
      m = n2

      return
      end function realread
! -----------------------------------------------------------------------
! return the substring that starts at loc m in str
! skip initial spaces and tabs
! substring is any chars up to next space or tab
! also return m = loc of char after end of substr in str
      subroutine strread(str,m,substr)
      implicit none
      character*(*) str,substr
      integer m

      integer n1,n2,nmax
      
      n1 = m
      nmax = len(str)
      do while (n1.le.nmax.and.&
           (str(n1:n1).eq.' '.or.ichar(str(n1:n1)).eq.9))
        n1 = n1 + 1
      enddo
      n2 = n1
      do while (n2.le.nmax.and.str(n2:n2).ne.' '&
           .and.ichar(str(n2:n2)).ne.9)
        n2 = n2 + 1
      enddo
      substr = str(n1:n2-1)
      m = n2

      return
      end subroutine strread
! -----------------------------------------------------------------------
! returns the actual length of str
! backtracks from end of string skipping over spaces
      integer function length(str)
      implicit none
      character*(*) str

      integer n

      n = len(str)
      do while (n.gt.0.and.str(n:n).eq.' ')
        n = n - 1
      enddo
      length = n

      return
      end function length
! -----------------------------------------------------------------------
!...convert integer i to a character cha
subroutine i2char(i,cha)
implicit none
integer,intent(in)::i
character(len=*),intent(out)::cha
integer::len_cha,len_i
character(len=20)::cha2

len_cha=len(cha)
write(*,*)"length(cha)=",len_cha
len_i=int(dlog10(dble(i)))+1
write(*,*)"length(i)=",len_i

if(len_i>len_cha)then
  write(*,*)"length of",trim(cha)," is too short for the integer."
  stop
else
  write(cha2,'(i10)')i
endif
!write(*,*)"cha2=",cha2
cha=trim(adjustl(cha2))

end subroutine i2char
! -----------------------------------------------------------------------
!...convert a character cha to an integer i
subroutine char2i(cha,i)
implicit none
character(len=*),intent(in)::cha
integer,intent(out)::i

write(*,*)"# The default value of maximum integer is 2147483647."
write(*,*)"# You can not ask computer to convert a string to an & 
integer number larger than it."
if(len_trim(adjustl(cha))>10)then
  write(*,*)"length of ",trim(cha)," is too long."
  stop
else
  read(cha,*)i
endif

end subroutine char2i
! -----------------------------------------------------------------------
!...combine two characters: ch3=ch1+ch2
subroutine add2cha(ch1,ch2,ch3)
implicit none
character(len=*),intent(in)::ch1,ch2
character(len=*),intent(out)::ch3
ch3=trim(adjustl(ch1))//trim(adjustl(ch2))
end subroutine add2cha
! -----------------------------------------------------------------------
