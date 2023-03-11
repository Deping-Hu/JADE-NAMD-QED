      subroutine sub_sh_ana_nac_cavity(n_atom, &
                                       i_state, &
                                       n_state, n_state_ele, n_photon_state, md_state, &
                                       ntime, dtime, ntime_ele, &
                                       n_sav_stat, n_sav_traj, &
                                       label_ZN, qm_method, seed_random, &
                                       cor_dec, &
                                       label_nac_phase, &
                                       label_reject_hops, &
                                       hop_e, &
                                       label_hop_prob, &
                                       wc, gc, & 
                                       photon_loss_rate, &
                                       af_cut, &
                                       kkr_mu, &
                                       efield_x, &
                                       efield_y, &
                                       efield_z, &
                                       Hamitonian_Compoment, &
                                       label_permanent_dipole, &
                                       label_initial_basis)

         implicit none
         include 'param.def'

         integer, intent(in) :: label_hop_prob
         integer, intent(in) :: n_atom, i_state, n_state, n_state_ele, n_photon_state
         integer, intent(in) :: ntime, ntime_ele
         integer, intent(in) :: n_sav_stat, n_sav_traj
         integer, intent(in) :: qm_method, label_ZN, kkr_mu
         integer, intent(inout) :: seed_random
         double precision, intent(in) :: cor_dec, photon_loss_rate, af_cut
         integer, intent(in) :: label_nac_phase
         integer, intent(in) :: label_reject_hops
         integer, intent(in) :: efield_x, efield_y, efield_z
         double precision, intent(in) ::  hop_e, wc, gc
         character(len=20), intent(in) :: Hamitonian_Compoment, label_permanent_dipole, label_initial_basis

         integer, intent(inout), dimension(n_state) :: md_state
         double precision, intent(inout) :: dtime

         double precision, dimension(n_state, n_state) :: U_ss, U_ss_trans, old_U_ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  - Define local variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         double precision, dimension(n_atom) :: &
            coor_x, &
            coor_y, &
            coor_z

         double precision, dimension(n_atom) :: &
            vel_x, &
            vel_y, &
            vel_z

         character*2, dimension(n_atom)  ::  atom_label

         integer :: file_save_energy, file_save_traj, &
                    file_save_grad, &
                    file_save_nac, &
                    file_save_vel, file_md_out, &
                    file_save_state, file_save_all, &
                    file_save_ele, file_save_pe, &
                    file_save_dipole, label_ml_pes, file_save_mu, &
                    file_save_hcav

         integer :: i, j, it

         double precision, dimension(n_atom) ::   mass, charge
         integer :: input_numer

         double precision, dimension(n_atom) ::   gradient_x, &
            gradient_y, &
            gradient_z

         double precision, dimension(n_state, n_atom) :: &
            gra_all_x, &
            gra_all_y, &
            gra_all_z

         double precision, dimension(n_state, n_atom) :: &
            old_gra_all_x, &
            old_gra_all_y, &
            old_gra_all_z

         double precision, dimension(n_state, n_state, n_atom) :: &
            nac_x, &
            nac_y, &
            nac_z

         double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_x, NACv_cav_y, NACv_cav_z

         double precision, dimension(n_state,n_state,n_atom) :: old_NACv_cav_x, old_NACv_cav_y, old_NACv_cav_z

         double precision, dimension(n_state, n_state, n_atom) :: &
            old_nac_x, &
            old_nac_y, &
            old_nac_z

         double precision, dimension(n_atom) :: &
            old_vel_x, &
            old_vel_y, &
            old_vel_z

         double precision :: time, pes_ref, pes_current

         double precision, dimension(n_state) ::  pes_all, &
            old_pes_all, &
            pop_tmp, &
            pop_tmp_1

         double precision, dimension(n_state_ele) ::  pes_all_ele

         complex(kind=8), dimension(n_state, n_state) :: rho

         double precision,  dimension(n_state, n_state) ::   Hcav, old_Hcav

         integer :: index_state, old_index_state, new_index_state

!       integer :: i_state, j_state, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
         double precision :: tmp_mass
         integer          :: tmp_charge
         character*2      :: tmp_label

         double precision :: num_random_initial

         integer::   nuclear_step_start, nuclear_step_finish
         integer::   count_rate, count_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! random number initialized
         call init_random_seed(seed_random)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   -       Find all masses !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         file_md_out = 23
         open (unit=file_md_out, file="dynamics.out")

         write (file_md_out, *) "--------------------------------------"
         write (file_md_out, *) "Run the surface hopping dynamics!"
         write (file_md_out, *) "--------------------------------------"
         write (file_md_out, *)
         write (file_md_out, *)
         write (file_md_out, *) "Initialization"
         write (file_md_out, *) "Parameters for the dynamics calculations"
         write (file_md_out, *) " Total time:", dtime*ntime
         write (file_md_out, *) " Time step:", dtime
         write (file_md_out, *) " The number of time steps:", ntime
         write (file_md_out, *) "The frequency to save the energy:", n_sav_stat
         write (file_md_out, *) "The frequency to save the geometry:", n_sav_traj
         write (file_md_out, *)
         write (file_md_out, *) "The quantum calculations"
         write (file_md_out, *) "Which methods? ", qm_method
         write (file_md_out, *) "(1: MNDO, 2: CASSCF-MOLPRO, 11: TDDFT-Turbomole,...)"
         write (file_md_out, *)
         write (file_md_out, *) "which state to start?", i_state
         write (file_md_out, *) "How many states are involved?", n_state
         write (file_md_out, *) "How many electronic states are involved?", n_state_ele
         write (file_md_out, 9998) "Which states are involved?", md_state(:)
         write (file_md_out, *) " Time step for electronic motions", ntime_ele
         if (0 .ge. cor_dec) then
            write (file_md_out, *) "No decoherent correction"
         else
            write (file_md_out, *) "Decoherent correction factor", cor_dec
         end if
         if (label_nac_phase .eq. 0) then
            write (file_md_out, *) "Do not follow the phase of NAC"
         end if
         if (label_nac_phase .eq. 1) then
            write (file_md_out, *) "Follow the phase of NAC"
         end if
         if ((label_nac_phase .ne. 0) &
             .and. &
             (label_nac_phase .ne. 1) &
             ) then
            write (file_md_out, *) "How to treat the phase of NAC?"
            stop
         end if

         if (label_reject_hops .eq. 0) then
            write (file_md_out, *) &
               "Do not reject the hops with large energy gap"
         else
            write (file_md_out, *) "Reject the hops when energy gap is larger than", hop_e, "(eV)"
         end if

9998     format(a, 10(i3, ";"))

         write (file_md_out, *) "--------------------------------------"

         file_save_energy = 24
         file_save_traj = 25
         file_save_vel = 26
         file_save_state = 27
         file_save_all = 28
         file_save_ele = 29
         file_save_pe = 30
         file_save_dipole = 31
         file_save_grad = 32
         file_save_mu = 33
         file_save_hcav = 34
         file_save_nac = 37

         open (unit=file_save_energy, file="energy_time.out")
         open (unit=file_save_traj, file="traj_time.out")
         open (unit=file_save_grad, file="grad_time.out")
         open (unit=file_save_nac, file="nac_time.out")
         open (unit=file_save_vel, file="vel_time.out")
         open (unit=file_save_state, file="current_state.out")
         open (unit=file_save_all, file="hop_all_time.out")
         open (unit=file_save_ele, file="ele_time.out")
         open (unit=file_save_pe, file="pe_time.out")
         open (unit=file_save_dipole, file="di_time.out")
         open (unit=file_save_mu, file="mu_time.out")
         open (unit=file_save_hcav, file="hcav_time.out")


!      Obttain the gradient and pes at time zero
         write (file_md_out, *) "Calculate the gradient"

         label_ml_pes = 0

         it = 0
         time = 0

         write (file_md_out, *) "Time steps:", it
         write (*, *) "Time steps:", it

!      Read initial coordinates

         call sub_read_current_geom(n_atom, &
                                    coor_x, coor_y, coor_z, &
                                    vel_x, vel_y, vel_z, &
                                    atom_label, &
                                    n_state, &
                                    n_state_ele, &
                                    rho, &
                                    index_state, &
                                    it, time)

         do i = 1, n_atom
            tmp_label = atom_label(i)
            call sub_get_mass(tmp_label, tmp_mass, tmp_charge)
            mass(i) = tmp_mass
            charge(i) = tmp_charge
         end do

         call sub_many_ana_cavity_pes(n_atom, &
                                      n_state, n_state_ele, md_state, &
                                      atom_label, &
                                      coor_x, coor_y, coor_z, &
                                      label_ZN, &
                                      label_ml_pes, &
                                      qm_method, &
                                      index_state, &
                                      gra_all_x, &
                                      gra_all_y, &
                                      gra_all_z, &
                                      nac_x, &
                                      nac_y, &
                                      nac_z, &
                                      NACv_cav_x, &
                                      NACv_cav_y, &
                                      NACv_cav_z, &
                                      Hcav, &
                                      pes_all, &
                                      pes_all_ele, &
                                      it, &
                                      time, &
                                      wc, &
                                      gc, &
                                      af_cut, &
                                      kkr_mu, &
                                      U_ss, &
                                      efield_x, &
                                      efield_y, &
                                      efield_z, &
                                      Hamitonian_Compoment, &
                                      label_permanent_dipole, &
                                      file_md_out, &
                                      file_save_dipole, &
                                      file_save_mu, &
                                      file_save_hcav)

         if (trim(adjustl(label_initial_basis)) == "adibatic_fock") then

            U_ss_trans = transpose(U_ss)

            rho = matmul(U_ss_trans,rho)
            rho = matmul(rho, U_ss)

            do i = 1, n_state
               pop_tmp(i) = rho(i, i)
            end do

            pop_tmp_1(1) = pop_tmp(1)

            do i = 2, n_state
               pop_tmp_1(i) = pop_tmp_1(i - 1) + pop_tmp(i)
            end do

            call get_random_number(num_random_initial)
            call get_random_number(num_random_initial)
            call get_random_number(num_random_initial)
            call get_random_number(num_random_initial)
            call get_random_number(num_random_initial)

            if (num_random_initial .lt. pop_tmp_1(1)) then
               index_state = 1
            end if

            do i = 2, n_state
               if ((num_random_initial .gt. pop_tmp_1(i - 1)) .and. (num_random_initial .lt. pop_tmp_1(i))) then
                  index_state = i
               end if
            end do

         end if

         call sub_save_restart(n_atom, &
                               coor_x, coor_y, coor_z, &
                               vel_x, vel_y, vel_z, &
                               atom_label, &
                               n_state, &
                               rho, &
                               index_state, &
                               it, time)

!        Save the energy
         pes_current = pes_all(index_state)
         pes_ref = pes_all_ele(1)

         write (file_save_energy, *) &
                              "# Step,Time, &
                                 Pot_ene, &
                                 Kin_ene, &
                                 Tot_ene, &
                                 lin_mom_x, &
                                 lin_mom_y, &
                                 lin_mom_z"

         call sub_save_energy(n_atom, &
                              pes_current, &
                              vel_x, vel_y, vel_z, &
                              mass, &
                              it, time, &
                              file_save_energy)

         call sub_save_pe(it, time, &
                          n_state, &
                          index_state, &
                          pes_all, &
                          file_save_pe)

         call sub_save_traj(n_atom, &
                            coor_x, coor_y, coor_z, &
                            vel_x, vel_y, vel_z, &
                            atom_label, &
                            it, time, &
                            file_save_traj, file_save_vel)

         gradient_x(:) = gra_all_x(index_state, :)
         gradient_y(:) = gra_all_y(index_state, :)
         gradient_z(:) = gra_all_z(index_state, :)

         call sub_save_grad(n_atom, n_state, index_state, &
                            gra_all_x, &
                            gra_all_y, &
                            gra_all_z, &
                            atom_label, &
                            it, time, &
                            file_save_grad)

         call sub_save_nac(n_atom, n_state, index_state, &
                            nac_x, &
                            nac_y, &
                            nac_z, &
                            atom_label, &
                            it, time, &
                            file_save_nac)
         write (file_save_state, *) &
                               "# Step, Time, &
                                  Current state, &
                                  Current potential, &
                                  Velocity norm, &
                                  current gradient"

         call sub_save_hopping_ana_nac_cavity(n_atom, n_state, md_state, &
                                              it, time, &
                                              index_state, &
                                              pes_all, &
                                              gra_all_x, &
                                              gra_all_y, &
                                              gra_all_z, &
                                              nac_x, &
                                              nac_y, &
                                              nac_z, &
                                              vel_x, vel_y, vel_z, &
                                              rho, &
                                              U_ss, &
                                              file_save_state, &
                                              file_save_all &
                                              )

!      Begin the propagation
         do it = 1, ntime

            call system_clock(nuclear_step_start,count_rate,count_max)
            nuclear_step_start = nuclear_step_start*1.0/count_rate

            time = it*dtime

            write (file_md_out, *) "Time steps:", it
            write (*, *) "Time steps:", it

            old_gra_all_x = gra_all_x
            old_gra_all_y = gra_all_y
            old_gra_all_z = gra_all_z

            old_nac_x = nac_x
            old_nac_y = nac_y
            old_nac_z = nac_z

            old_NACv_cav_x = NACv_cav_x

            old_NACv_cav_y = NACv_cav_y

            old_NACv_cav_z = NACv_cav_z

            old_Hcav = Hcav

            old_U_ss = U_ss
            old_pes_all = pes_all

            old_vel_x = vel_x
            old_vel_y = vel_y
            old_vel_z = vel_z

            gradient_x(:) = gra_all_x(index_state, :)
            gradient_y(:) = gra_all_y(index_state, :)
            gradient_z(:) = gra_all_z(index_state, :)

            call sub_newton_vel(n_atom, &
                                vel_x, vel_y, vel_z, &
                                gradient_x, &
                                gradient_y, &
                                gradient_z, &
                                mass, &
                                dtime)

!          Propagate the velocity full time step dt
            call sub_newton_coor(n_atom, &
                                 coor_x, coor_y, coor_z, &
                                 vel_x, vel_y, vel_z, &
                                 dtime)

!          Get all gradients

            call sub_many_ana_cavity_pes(n_atom, &
                                         n_state, n_state_ele, md_state, &
                                         atom_label, &
                                         coor_x, coor_y, coor_z, &
                                         label_ZN, &
                                         label_ml_pes, &
                                         qm_method, &
                                         index_state, &
                                         gra_all_x, &
                                         gra_all_y, &
                                         gra_all_z, &
                                         nac_x, &
                                         nac_y, &
                                         nac_z, &
                                         NACv_cav_x, &
                                         NACv_cav_y, &
                                         NACv_cav_z, &
                                         Hcav, &
                                         pes_all, &
                                         pes_all_ele, &
                                         it, &
                                         time, &
                                         wc, &
                                         gc, &
                                         af_cut, &
                                         kkr_mu, &
                                         U_ss, &
                                         efield_x, &
                                         efield_y, &
                                         efield_z, &
                                         Hamitonian_Compoment, &
                                         label_permanent_dipole, &
                                         file_md_out, &
                                         file_save_dipole, &
                                         file_save_mu, &
                                         file_save_hcav)

            pes_current = pes_all(index_state)

!          Propagate the velocity half time step dt/2
            gradient_x(:) = gra_all_x(index_state, :)
            gradient_y(:) = gra_all_y(index_state, :)
            gradient_z(:) = gra_all_z(index_state, :)

            call sub_newton_vel(n_atom, &
                                vel_x, vel_y, vel_z, &
                                gradient_x, &
                                gradient_y, &
                                gradient_z, &
                                mass, &
                                dtime)

            if (label_nac_phase .eq. 1) then
               call sub_nac_phase(n_atom, &
                                  n_state, &
                                  old_nac_x, &
                                  old_nac_y, &
                                  old_nac_z, &
                                  nac_x, &
                                  nac_y, &
                                  nac_z)
            end if

            old_index_state = index_state

            write(*,*) "Propagate in the AF state"

            call sub_tully_hopping_ana_nac_cavity_af(n_atom, &
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


            new_index_state = index_state

!         Reject all hops when the energy gap is large than hop_e if necessary

            if (label_reject_hops .ne. 0) then
               if (old_index_state .ne. new_index_state) then
                  if (abs(pes_all(old_index_state) - pes_all(new_index_state)) &
                      .ge. &
                      (hop_e/TOEV) &
                      ) then
                     new_index_state = old_index_state
                     index_state = old_index_state
                  end if
               end if
            end if

            if (old_index_state .ne. new_index_state) then
               write (file_md_out, *) "Possible hop happens at Step", &
                  it
               write (file_md_out, *) md_state(old_index_state), &
                  "-->", &
                  md_state(new_index_state)

               call sub_velocity_scaling(n_atom, &
                                         n_state, &
                                         mass, &
                                         pes_all, &
                                         nac_x, &
                                         nac_y, &
                                         nac_z, &
                                         old_index_state, &
                                         new_index_state, &
                                         vel_x, vel_y, vel_z, &
                                         file_md_out)
               index_state = new_index_state
            end if

            if (old_index_state .ne. new_index_state) then

               write (file_md_out, *) "Real hop due to TDSE happens at Step", it

               write (file_md_out, *) md_state(old_index_state), &
                  "-->", &
                  md_state(new_index_state)
               write (file_md_out, *) "Time(fs), energy gap (au)", &
                  time*TOFS, &
                  pes_all(old_index_state) - pes_all(new_index_state)
            end if

            old_index_state = index_state

            call sub_tully_hopping_ana_nac_cavity_lb(n_atom, &
                                                     n_state, &
                                                     n_state_ele, &
                                                     rho, &
                                                     U_ss, &
                                                     time, &
                                                     it, &
                                                     dtime, &
                                                     ntime_ele, &
                                                     index_state, &
                                                     photon_loss_rate, &
                                                     file_save_ele, &
                                                     file_md_out, &
                                                     label_hop_prob)

            new_index_state = index_state

            if (old_index_state .ne. new_index_state) then

               write (file_md_out, *) "Real hop due to Linblad happens at Step", it

               write (file_md_out, *) md_state(old_index_state), &
                  "-->", &
                  md_state(new_index_state)
               write (file_md_out, *) "Time(fs), energy gap (au)", &
                  time*TOFS, &
                  pes_all(old_index_state) - pes_all(new_index_state)
            end if

!   Decoherence corrections after electron propagation
            if (cor_dec .gt. 0) then
               if (real(rho(index_state,index_state)) .gt. 0.d0) then
                  call sub_decoherent_corr(n_atom, &
                                          mass, &
                                          vel_x, &
                                          vel_y, &
                                          vel_z, &
                                          n_state, &
                                          pes_all, &
                                          rho, &
                                          dtime, &
                                          index_state, &
                                          cor_dec, &
                                          file_md_out)


              endif
           endif

!!!! I forget why I do the electronic calulation again after hopping

!           if (old_index_state .ne. new_index_state) then
!
!               call sub_many_ana_cavity_pes(n_atom, &
!                                            n_state, n_state_ele, md_state, &
!                                            atom_label, &
!                                            coor_x, coor_y, coor_z, &
!                                            label_ZN, &
!                                            label_ml_pes, &
!                                            qm_method, &
!                                            index_state, &
!                                            gra_all_x, &
!                                            gra_all_y, &
!                                            gra_all_z, &
!                                            nac_x, &
!                                            nac_y, &
!                                            nac_z, &
!                                            NACv_cav_x, &
!                                            NACv_cav_y, &
!                                            NACv_cav_z, &
!                                            Hcav, &
!                                            pes_all, &
!                                            pes_all_ele, &
!                                            it, &
!                                            time, &
!                                            wc, &
!                                            gc, &
!                                            af_cut, &
!                                            kkr_mu, &
!                                            U_ss, &
!                                            efield_x, &
!                                            efield_y, &
!                                            efield_z, &
!                                            Hamitonian_Compoment, &
!                                            label_permanent_dipole, &
!                                            file_md_out, &
!                                            file_save_dipole, &
!                                            file_save_mu, &
!                                            file_save_hcav)
!           end if

!       save restart file
           call sub_save_restart(n_atom, &
                                  coor_x, coor_y, coor_z, &
                                  vel_x, vel_y, vel_z, &
                                  atom_label, &
                                  n_state, &
                                  rho, &
                                  index_state, &
                                  it, time)

           call sub_save_hopping_ana_nac_cavity(n_atom, n_state, md_state, &
                                                 it, time, &
                                                 index_state, &
                                                 pes_all, &
                                                 gra_all_x, &
                                                 gra_all_y, &
                                                 gra_all_z, &
                                                 nac_x, &
                                                 nac_y, &
                                                 nac_z, &
                                                 vel_x, vel_y, vel_z, &
                                                 rho, &
                                                 U_ss, &
                                                 file_save_state, &
                                                 file_save_all &
                                                 )

            call sub_save_pe(it, time, &
                             n_state, &
                             index_state, &
                             pes_all, &
                             file_save_pe)

!!!       Save the current geometry and velocity
            call sub_write_current_geom(n_atom, &
                                        coor_x, coor_y, coor_z, &
                                        vel_x, vel_y, vel_z, &
                                        atom_label, &
                                        n_state, &
                                        n_state_ele, &
                                        rho, &
                                        index_state, &
                                        it, time)

!!!!!!    Save the energy
            if (mod(it, n_sav_stat) .eq. 0) then
               pes_current = pes_all(index_state)
               call sub_save_energy(n_atom, &
                                    pes_current, &
                                    vel_x, vel_y, vel_z, &
                                    mass, &
                                    it, time, &
                                    file_save_energy)

            end if

!!!!!     Save the trajectory
            if (mod(it, n_sav_traj) .eq. 0) then
               call sub_save_traj(n_atom, &
                                  coor_x, coor_y, coor_z, &
                                  vel_x, vel_y, vel_z, &
                                  atom_label, &
                                  it, time, &
                                  file_save_traj, file_save_vel)
            end if

!!!!!     Save the gradient
            if (mod(it, n_sav_traj) .eq. 0) then
               call sub_save_grad(n_atom, n_state, index_state, &
                                  gra_all_x, &
                                  gra_all_y, &
                                  gra_all_z, &
                                  atom_label, &
                                  it, time, &
                                  file_save_grad)

               call sub_save_nac(n_atom, n_state, index_state, &
                                  nac_x, &
                                  nac_y, &
                                  nac_z, &
                                  atom_label, &
                                  it, time, &
                                  file_save_nac)
            end if

!          write (28,9999) it, coor_x(:), gradient_x(:), vel_x(:)

         call system_clock(nuclear_step_finish,count_rate,count_max)
         nuclear_step_finish = nuclear_step_finish*1.0/count_rate

         write(*,*) "total time for one nulcear step", nuclear_step_finish - nuclear_step_start

         end do

         close (file_md_out)
         close (file_save_energy)
         close (file_save_traj)
         close (file_save_vel)
         close (file_save_state)
         close (file_save_all)
         close (file_save_ele)
         close (file_save_pe)
         close (file_save_grad)
         close (file_save_nac)
         close (file_save_dipole)
         close (file_save_mu)
         close (file_save_hcav)

9999     format(i10, 10(f20.10, 1x))
9997     format(a)

         return

      end
