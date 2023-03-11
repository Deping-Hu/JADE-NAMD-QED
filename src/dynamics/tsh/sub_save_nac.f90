      subroutine  sub_save_nac (n_atom, nstate, index_state, &
                                nac_x, &
                                nac_y, &
                                nac_z, &
                                atom_label, &
                                it, time, &
                                file_save_nac )

      implicit none
      include 'param.def'

      integer, intent(in) :: n_atom, it, nstate, index_state
      integer, intent(in) :: file_save_nac
      double precision, intent(in) :: time
      double precision, intent(in), dimension(nstate,nstate,n_atom) ::   &
                                                nac_x, &
                                                nac_y, &
                                                nac_z
                                                
      character*2, intent(in), dimension(n_atom)  ::  atom_label

      integer ::  i,j,k
     
      do i=1, nstate
      do j=1, nstate

         if (i .lt. j) then
            write (file_save_nac, *) n_atom

            write (file_save_nac, 9998)  "AU", "Step:", it, "State:", i, "State:", j, &
                                           "Time:", time

            do  k=1, n_atom  
                write (file_save_nac, 9999) atom_label(k), &
                                             nac_x(i,j,k), &
                                             nac_y(i,j,k), &
                                             nac_z(i,j,k)
            enddo

         end if
      enddo
      enddo


9999   format(a, 1x, 3(f20.10, 1x))
9998   format(a, 1x, a,  1x, i10, 1x, a, 1x, i10, 1x, a, 1x, i10, 1x, a, f20.10)

       return
 
       end 

