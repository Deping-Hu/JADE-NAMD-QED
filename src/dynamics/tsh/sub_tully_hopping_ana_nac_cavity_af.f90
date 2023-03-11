    subroutine sub_tully_hopping_ana_nac_cavity_af(n_atom, &
                                                   mass, &
                                                   old_vel_x, &
                                                   old_vel_y, &
                                                   old_vel_z, &
                                                   vel_x, &
                                                   vel_y, &
                                                   vel_z, &
                                                   n_state, &
                                                   n_state_ele, &
                                                   old_Hcav, &
                                                   Hcav, &
                                                   pes_all, &
                                                   pes_ref, &
                                                   NACv_cav_x, &
                                                   NACv_cav_y, &
                                                   NACv_cav_z, &
                                                   old_NACv_cav_x, &
                                                   old_NACv_cav_y, &
                                                   old_NACv_cav_z, &
                                                   rho, &
                                                   U_ss, &
                                                   old_U_ss, &
                                                   time, &
                                                   it, &
                                                   dtime, &
                                                   ntime_ele, &
                                                   index_state, &
                                                   cor_dec, &
                                                   file_save_ele, &
                                                   file_md_out, &
                                                   label_hop_prob)

       use write_matrix_2file
       implicit none

!
!     ----- Parameters -----
!
       include 'param.def'
!
!     ----- Argument -----

       integer, intent(in) :: n_atom, n_state, ntime_ele, it, n_state_ele, &

                              file_save_ele, file_md_out, label_hop_prob

       integer, intent(inout) :: index_state

       double precision, intent(in) :: dtime, pes_ref, time, cor_dec

       double precision, intent(in), dimension(n_atom) :: mass, vel_x, vel_y, vel_z, &

                                                          old_vel_x, old_vel_y, old_vel_z

       double precision, intent(in), dimension(n_state) ::  pes_all

       double precision, intent(in), dimension(n_state, n_state) :: Hcav, old_Hcav, old_U_ss, U_ss

       double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_x, NACv_cav_y, NACv_cav_z, &

                                                              old_NACv_cav_x, old_NACv_cav_y, old_NACv_cav_z

       complex(kind=8), intent(inout), dimension(n_state, n_state) :: rho

!
!  Local

       integer :: i, j, k, l, it_ele

       double precision :: delt, num_random, rho_decress_all

       double precision, dimension(n_state) :: pro_hop, pro_sum, pes_eff

       double precision, dimension(n_atom) :: vel_x_eff, vel_y_eff, vel_z_eff

       double precision, dimension(n_state, n_state) :: U_ss_eff, Hcav_eff, U_ss_trans

       double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_x_eff, NACv_cav_y_eff, NACv_cav_z_eff

       complex(kind=8), dimension(n_state, n_state) :: old_rho, rho_af, old_rho_af, tmp1, tmp2, tmp3, tmp4, &

                                                       h_vv, h_vv_conjg, k1, k2, k3, k4, rho_k_tmp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       delt = dtime/ntime_ele

       pro_hop = 0.d0

       U_ss_trans = transpose(old_U_ss)

       rho_af = matmul(old_U_ss,rho) 
       rho_af = matmul(rho_af,U_ss_trans) 

       old_rho_af = rho_af

       old_rho = rho

       do it_ele = 1, ntime_ele

!!!!! Begin density matrix evolution 

          Hcav_eff = old_Hcav + it_ele*(Hcav - old_Hcav)/ntime_ele

          NACv_cav_x_eff = old_NACv_cav_x + &
                        it_ele*(NACv_cav_x - old_NACv_cav_x)/ntime_ele

          NACv_cav_y_eff = old_NACv_cav_y + &
                        it_ele*(NACv_cav_y - old_NACv_cav_y)/ntime_ele

          NACv_cav_z_eff = old_NACv_cav_z + &
                        it_ele*(NACv_cav_z - old_NACv_cav_z)/ntime_ele

          vel_x_eff = old_vel_x + &
                      it_ele*(vel_x - old_vel_x)/ntime_ele

          vel_y_eff = old_vel_y + &
                      it_ele*(vel_y - old_vel_y)/ntime_ele

          vel_z_eff = old_vel_z + &
                      it_ele*(vel_z - old_vel_z)/ntime_ele

          do i = 1, n_state
          do j = 1, n_state
             if (i .eq. j) then
                Hcav_eff(i, j) = Hcav_eff(i, j) - pes_ref
             end if
          end do
          end do

          h_vv(:, :) = cmplx(0.d0, 0.d0)

          do i = 1, n_state
          do j = 1, n_state
             do k = 1, n_atom
                h_vv(i, j) = h_vv(i, j) &
                             + NACv_cav_x_eff(i, j, k)*vel_x_eff(k) &
                             + NACv_cav_y_eff(i, j, k)*vel_y_eff(k) &
                             + NACv_cav_z_eff(i, j, k)*vel_z_eff(k)
             end do
          end do
          end do

!!!!!!! 4th Runge-Kutta to propagate the density matrix

!!!!!!!!!!!!! k1
          rho_k_tmp = rho_af
          tmp1 = MATMUL(Hcav_eff, rho_k_tmp)
          tmp2 = MATMUL(rho_k_tmp, Hcav_eff)

          tmp3 = MATMUL(h_vv, rho_k_tmp)
          tmp4 = MATMUL(rho_k_tmp, h_vv)

          k1 = -c1*(tmp1-tmp2)-(tmp3-tmp4)

          rho_k_tmp = rho_af + k1*delt/2.d0

!!!!!!!!!!!!! k2
          tmp1 = MATMUL(Hcav_eff, rho_k_tmp)
          tmp2 = MATMUL(rho_k_tmp, Hcav_eff)

          tmp3 = MATMUL(h_vv, rho_k_tmp)
          tmp4 = MATMUL(rho_k_tmp, h_vv)

          k2 = -c1*(tmp1-tmp2)-(tmp3-tmp4)

          rho_k_tmp = rho_af + k2*delt/2.d0

!!!!!!!!!!!!! k3
          tmp1 = MATMUL(Hcav_eff, rho_k_tmp)
          tmp2 = MATMUL(rho_k_tmp, Hcav_eff)

          tmp3 = MATMUL(h_vv, rho_k_tmp)
          tmp4 = MATMUL(rho_k_tmp, h_vv)

          k3 = -c1*(tmp1-tmp2)-(tmp3-tmp4)

          rho_k_tmp = rho_af + k3*delt

!!!!!!!!!!!!! k4
          tmp1 = MATMUL(Hcav_eff, rho_k_tmp)
          tmp2 = MATMUL(rho_k_tmp, Hcav_eff)

          tmp3 = MATMUL(h_vv, rho_k_tmp)
          tmp4 = MATMUL(rho_k_tmp, h_vv)

          k4 = -c1*(tmp1-tmp2)-(tmp3-tmp4)

!!!!! Sum all
          rho_af = rho_af +(1/6.d0)*delt*(k1+2*k2+2*k3+k4)
!!!!! End RK


       end do

       U_ss_trans = transpose(U_ss)

       rho = matmul(U_ss_trans,rho_af) 
       rho = matmul(rho,U_ss) 

!!!! Global Flux Surface Hopping (GFSH), calculate hop probability at each nuclear step
       if (label_hop_prob .eq. 3) then

          if (real(rho(index_state, index_state)) .lt. real(old_rho(index_state, index_state))) then

             rho_decress_all = 0.d0

             do i = 1, n_state
                if (real(rho(i, i)) .lt. real(old_rho(i, i))) then
                   rho_decress_all = rho_decress_all + (real(old_rho(i, i)) - real(rho(i, i)))
                end if
             end do

             do i = 1, n_state

                pro_hop(i) = (real(rho(i, i)) - real(old_rho(i, i))) &
                             /real(old_rho(index_state, index_state)) &
                             *(real(old_rho(index_state, index_state)) - real(rho(index_state, index_state))) &
                             /rho_decress_all

             end do

          end if

       endif

!   Decoherence corrections after electron propagation
!       if (cor_dec .gt. 0) then
!          if (real(rho(index_state,index_state)) .gt. 0.d0) then
!             call sub_decoherent_corr(n_atom, &
!                                      mass, &
!                                      vel_x, &
!                                      vel_y, &
!                                      vel_z, &
!                                      n_state, &
!                                      pes_all, &
!                                      rho, &
!                                      dtime, &
!                                      index_state, &
!                                      cor_dec, &
!                                      file_md_out)
!
!
!          endif
!       endif

       do i = 1, n_state
          if (pro_hop(i) .lt. 0.d0) then
             pro_hop(i) = 0.d0
          end if
       end do

       pro_hop(index_state) = 1
       do i = 1, n_state
          if (index_state .ne. i) then
             pro_hop(index_state) = pro_hop(index_state) - pro_hop(i)
          end if
       end do

       pro_sum = 0.d0
       do i = 1, n_state
          do j = 1, i
             pro_sum(i) = pro_sum(i) + pro_hop(j)
          end do
       end do

       write (file_save_ele, *) "------------------------------------"

       do i = 1, n_state
       do j = 1, n_state
          write (file_save_ele, 9999) it, time*TOFS, &
             "rho", i, j, &
             REAL(rho(i, j)), '+i', AIMAG(rho(i, j))
       end do
       end do

       write (file_save_ele, 9996) it, time*TOFS, &
          "The current state", index_state
       write (file_save_ele, 9997) it, time*TOFS, &
          "Hopping probability", pro_hop(:)
       write (file_save_ele, 9997) it, time*TOFS, &
          "Area for hopping", pro_sum(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!       generate random number
!       the init_random_seed may be set before the dynamic cycle.
       call get_random_number(num_random)

!       Surface-hopping check

       if ((0 .le. num_random) &
           .and. &
           (num_random .lt. pro_sum(1)) &
           ) then
          index_state = 1
          goto 7000
       else
          do i = 2, n_state
             if ((pro_sum(i - 1) .le. num_random) &
                 .and. &
                 (num_random .lt. pro_sum(i)) &
                 ) then
                index_state = i
                goto 7000
             end if
          end do
       end if

7000   write (file_save_ele, 9997) it, time*TOFS, &
          "Random number", num_random
       write (file_save_ele, 9996) it, time*TOFS, &
          "The new state", index_state
       write (file_save_ele, *) "------------------------------------"

9996   format(i10, 1x, f15.8, 3x, a, 1x, i3)
9997   format(i10, 1x, f15.8, 3x, a, 1x, 10(f15.8, 1x))
9999   format(i10, 1x, f15.8, 3x, a, 1x, i3, 1x, i3, 1x, f10.7, a, f10.7)

    end
