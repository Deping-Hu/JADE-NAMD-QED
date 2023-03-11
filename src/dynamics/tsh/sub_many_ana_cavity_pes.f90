      subroutine    sub_many_ana_cavity_pes (n_atom, &
                                n_state, n_state_ele, md_state, &
                                atom_label, &
                                coor_x, coor_y, coor_z, &
                                label_no_nac, &
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
                                file_save_hcav )
 

      use write_matrix_2file

      implicit none
      include 'param.def'

      integer, intent(in) :: n_atom, qm_method, n_state, n_state_ele, label_no_nac, label_ml_pes, kkr_mu
      double precision, intent(in) :: time, wc, gc, af_cut
      integer, intent(in) :: efield_x, efield_y, efield_z
      character(len=20), intent(in) :: Hamitonian_Compoment, label_permanent_dipole
      integer, intent(in) :: file_md_out, file_save_dipole, file_save_mu, file_save_hcav, &
                             index_state, it
      integer, intent(in), dimension(n_state) :: md_state(n_state)
      character*2, intent(in), dimension(n_atom)  ::  atom_label 
      double precision, intent(in), dimension(n_atom) ::   &
                                                      coor_x, &
                                                      coor_y, &
                                                      coor_z


      double precision, intent(inout), dimension(n_state) :: pes_all
      double precision, intent(inout), dimension(n_state, n_atom) ::   &
                                                       gra_all_x, &
                                                       gra_all_y, &
                                                       gra_all_z
 
      double precision, intent(inout), dimension(n_state, n_state, n_atom) ::   &
                                                       nac_x, &
                                                       nac_y, &
                                                       nac_z

      double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_x
      double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_y
      double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_z

      integer  ::  i, j, k, i_status, system, i_state, j_state
      integer  ::  file_coor
      character*200 :: string

!      real ::   cpu_time_tmp1, cpu_time_tmp2, cpu_time_tmp3, cpu_time_tmp4
      integer::   cpu_time_tmp1, cpu_time_tmp2, cpu_time_tmp3, cpu_time_tmp4
      integer::   count_rate, count_max

!!!!! cavity

      double precision,  dimension(n_state, n_state) ::   Hcav

      double precision,  dimension(n_state, n_state, n_atom) :: gHcav_x
      double precision,  dimension(n_state, n_state, n_atom) :: gHcav_y
      double precision,  dimension(n_state, n_state, n_atom) :: gHcav_z

      double precision,  dimension(n_state, n_state) :: U_ss

!!!!! without cavity

      integer  ::  np

      double precision, dimension(n_state_ele):: pes_all_ele

      double precision, allocatable, dimension(:,:):: &
                                      gra_all_ele_x, &
                                      gra_all_ele_y, &
                                      gra_all_ele_z

      double precision, allocatable, dimension(:,:,:):: &
                                      nac_ele_x, &
                                      nac_ele_y, &
                                      nac_ele_z

      double precision, allocatable, dimension(:,:):: &
                                      mu_x, &
                                      mu_y, &
                                      mu_z, &
                                      mu_all, &
                                      mu_learned_x, &
                                      mu_learned_y, &
                                      mu_learned_z, &
                                      mu_learned_all

      double precision, allocatable, dimension(:,:,:):: &
                                      g_mu_xx, &
                                      g_mu_xy, &
                                      g_mu_xz, &
                                      g_mu_yx, &
                                      g_mu_yy, &
                                      g_mu_yz, &
                                      g_mu_zx, &
                                      g_mu_zy, &
                                      g_mu_zz, &
                                      g_mu_all_x, &
                                      g_mu_all_y, &
                                      g_mu_all_z
                                     

      double precision, allocatable, dimension(:,:) :: self_int_all

      double precision, allocatable, dimension(:,:,:) :: g_self_int_x
      double precision, allocatable, dimension(:,:,:) :: g_self_int_y
      double precision, allocatable, dimension(:,:,:) :: g_self_int_z

      np = 2

      allocate(gra_all_ele_x(n_state_ele,n_atom))
      allocate(gra_all_ele_y(n_state_ele,n_atom))
      allocate(gra_all_ele_z(n_state_ele,n_atom))

      allocate(nac_ele_x(n_state_ele,n_state_ele,n_atom))
      allocate(nac_ele_y(n_state_ele,n_state_ele,n_atom))
      allocate(nac_ele_z(n_state_ele,n_state_ele,n_atom))

      allocate(mu_x(n_state_ele,n_state_ele))
      allocate(mu_y(n_state_ele,n_state_ele))
      allocate(mu_z(n_state_ele,n_state_ele))
      allocate(mu_all(n_state_ele,n_state_ele))

      allocate(mu_learned_x(n_state_ele,n_state_ele))
      allocate(mu_learned_y(n_state_ele,n_state_ele))
      allocate(mu_learned_z(n_state_ele,n_state_ele))
      allocate(mu_learned_all(n_state_ele,n_state_ele))

      allocate(self_int_all(n_state_ele,n_state_ele))

      allocate(g_mu_xx(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_xy(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_xz(n_atom,n_state_ele,n_state_ele))

      allocate(g_mu_yx(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_yy(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_yz(n_atom,n_state_ele,n_state_ele))

      allocate(g_mu_zx(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_zy(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_zz(n_atom,n_state_ele,n_state_ele))

      allocate(g_mu_all_x(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_all_y(n_atom,n_state_ele,n_state_ele))
      allocate(g_mu_all_z(n_atom,n_state_ele,n_state_ele))

      allocate(g_self_int_x(n_atom,n_state_ele,n_state_ele))
      allocate(g_self_int_y(n_atom,n_state_ele,n_state_ele))
      allocate(g_self_int_z(n_atom,n_state_ele,n_state_ele))

      file_coor=22
      open(unit=file_coor, file="coor_temp.xyz")
      write (file_coor, *)  n_atom
      write (file_coor, *)
      do i=1,n_atom
      write (file_coor, 9999)  atom_label(i),  &
                       coor_x(i),      &
                       coor_y(i),      &
                       coor_z(i)
      enddo
      close(file_coor)

9999  format (a2,1x, 3(f24.14,1x))

      gra_all_ele_x = 0.d0
      gra_all_ele_y = 0.d0
      gra_all_ele_z = 0.d0
      nac_ele_x     = 0.d0
      nac_ele_y     = 0.d0
      nac_ele_z     = 0.d0

!     Write the interface for QC calculations
        !write(*,*) "qm_method:", qm_method
        open(unit=102, file="qm_interface")
        write (102,*)  "The interface between Fortran and Python"

        if ( qm_method .ne. 3 ) then 
           write (file_md_out, *) "Quan-Chem Theory:", qm_method
           write (102,*)  "Quan-Chem package: ", qm_method, it
        else if ( qm_method .eq. 3 ) then

           if (label_ml_pes .eq. 1) then
              write (file_md_out, *) "kkr regression", qm_method
              write (102,*)  "Quan-Chem package: ", qm_method, it
           else if (label_ml_pes .eq. 0) then
              write (file_md_out, *) "Quan-Chem Theory:", qm_method
              write (102,*)  "Quan-Chem package: ", qm_method, it
           endif

        endif    

        write (102,*)  "-----------------------------------"
        write (102,*)  "Current Geometry"
        write (102, *)  n_atom
        write (102, *)  "UNIT(au)"
        do i=1,n_atom
           write (102, 9999)  atom_label(i),  &
                      coor_x(i),      &
                      coor_y(i),      &
                      coor_z(i)
        enddo
        write (102,*)  "-----------------------------------"
        write (102,*)  "Number of atom:", n_atom
        write (102,*)  "Number of electronic states involved in the dynamics:", n_state_ele
        write (102,*)  "Current state:", index_state

        close(102)
  
        ! start electronic structure calculation.
        i_status = 0
        

        call system_clock(cpu_time_tmp1,count_rate,count_max)
        cpu_time_tmp1 = cpu_time_tmp1*1.0/count_rate

        if ( qm_method .ne. 3 ) then 
           i_status=system ("quantum.py")
        else if ( qm_method .eq. 3 ) then
           if ( label_ml_pes .eq. 1) then
              i_status=system ("jade_kkr.py")

              if (i_status  .ne. 0) then
                  write (file_md_out, *) "kkr2qm"
                  i_status=system ("quantum.py")
              endif  

           else if ( label_ml_pes .eq. 0) then
              i_status=system ("quantum.py")
           endif

        endif


        if (i_status  .ne. 0) then
            write (*,*) "The interface errors!", i_status
            write (*,*) "Check QM calculations at t=0!"     
            stop   
        endif  

        call system_clock(cpu_time_tmp2,count_rate,count_max)
        cpu_time_tmp2 = cpu_time_tmp2*1.0/count_rate

        write(*,*) "total quantum calculation time", cpu_time_tmp2 - cpu_time_tmp1

        i_status=system ("jade_kkr.py")
        if (i_status  .ne. 0) then
            write (*,*) "The machine learning dipole errors!", i_status
            stop   
        endif  

        call system_clock(cpu_time_tmp3,count_rate,count_max)
        cpu_time_tmp3 = cpu_time_tmp3*1.0/count_rate

        write(*,*) "total machine learning time", cpu_time_tmp3 - cpu_time_tmp2


        i_status=system ("rotation_mu_gradient_back.py")
        if (i_status  .ne. 0) then
            write (*,*) "The rotation of mu and their gradients errors!", i_status
            stop   
        endif  

        call system_clock(cpu_time_tmp4,count_rate,count_max)
        cpu_time_tmp4 = cpu_time_tmp4*1.0/count_rate

        write(*,*) "total rotation time", cpu_time_tmp4 - cpu_time_tmp3


        ! direct read qm results dat & dump.
        
        open(unit=101, file="qm_results.dat")
        read (101,*)
        read (101,*)  
        do i=1,n_atom
           read (101, *)
        enddo
        read (101,*)


        do j=1,n_state_ele
           read (101, *) pes_all_ele(j)
        enddo

        read (101,*)
        do j=1,n_state_ele
           read (101,*)
           do i=1,n_atom
              read (101, *)    &
                             gra_all_ele_x(j,i),      &
                             gra_all_ele_y(j,i),      &
                             gra_all_ele_z(j,i)
           enddo
        enddo


        if (label_no_nac .eq. 0) then
           read (101,*)
           do j=1, n_state_ele
              do k=1, n_state_ele
                 read (101,*)
                 do i=1,n_atom
                    read (101, *)   &
                                nac_ele_x(j,k,i),      &
                                nac_ele_y(j,k,i),      &
                                nac_ele_z(j,k,i)
                 enddo
             enddo
           enddo

        endif

        read (101,*)
        do j=1, n_state_ele
           do k=1, n_state_ele
              read (101,*)
              read (101, *)   &
                             mu_x(j,k),      &
                             mu_y(j,k),      &
                             mu_z(j,k)
          enddo
        enddo

        close(101)

       do i=1, n_state_ele
          do j=1, n_state_ele
             mu_all(i,j) =  mu_x(i,j) * efield_x + mu_y(i,j) * efield_y + mu_z(i,j) * efield_z
         enddo
       enddo

        write (file_save_dipole, *) "----------------------------------"
        write (file_save_dipole, *) "----------------------------------"
        write (file_save_dipole, *)  "STEP",  it 
        write (file_save_dipole, *) "----------------------------------"
        open(unit=103, file="qm_other.dat")
        do while (.TRUE.)
           read (103, 8001, end=8998) string
           write (file_save_dipole, 8001) string
        enddo
8998    close(103)
        write (file_save_dipole, *) "----------------------------------"

8001    format(a100)


       open(unit=101, file="ml_dipole_rotation_back.dat")
       read (101,*)
       read (101,*)  
       do i=1,n_atom
          read (101, *)
       enddo
       read (101,*)


       do i=1,n_state_ele
          do j=1,n_state_ele
             if(i .le. j) then
                read (101, *) mu_learned_x(i,j)
                read (101, *) mu_learned_y(i,j)
                read (101, *) mu_learned_z(i,j)
             else if(i .gt. j) then
                mu_learned_x(i,j) = mu_learned_x(j,i)
                mu_learned_y(i,j) = mu_learned_y(j,i)
                mu_learned_z(i,j) = mu_learned_z(j,i)
             endif
          enddo
       enddo

       do i=1, n_state_ele
          do j=1, n_state_ele
             mu_learned_all(i,j) =  mu_learned_x(i,j) * efield_x + mu_learned_y(i,j) * efield_y + mu_learned_z(i,j) * efield_z
         enddo
       enddo

! without permanent dipole
       if(trim(adjustl(label_permanent_dipole))=="False")then
          do i=1,n_state_ele
             mu_all(i,i) = 0.d0
             mu_learned_all(i,i) = 0.d0
          enddo
       endif
       
!!!!!!!!


       if(kkr_mu==1)then
          self_int_all = matmul(mu_learned_all,mu_learned_all)
       else if(kkr_mu==0)then
          self_int_all = matmul(mu_all,mu_all)
       endif

!       mu_learned_all=0
!       self_int_all(1,1) = 0
!       self_int_all(2,2) = 0

!       write (file_save_mu, 9998) it, time*TOFS, &
!              ((mu_x(j,k),j=1,n_state_ele),k=1,n_state_ele), &
!              ((mu_y(j,k),j=1,n_state_ele),k=1,n_state_ele), &
!              ((mu_z(j,k),j=1,n_state_ele),k=1,n_state_ele)

       write (file_save_mu, 9998) it, time*TOFS, &
              ((mu_all(j,k),j=1,n_state_ele),k=1,n_state_ele), &
              ((mu_learned_all(j,k),j=1,n_state_ele),k=1,n_state_ele), &
              (((mu_learned_all(j,k)-mu_all(j,k))/mu_all(j,k),j=1,n_state_ele),k=1,n_state_ele)
   
9998    format(i10, 1x, 100(f20.10, 1x))

       read (101,*)
       do i=1,n_state_ele
         do j=1,n_state_ele

          if(i .le. j) then
             read (101,*)
             do k=1,n_atom
                read (101, *)    &
                               g_mu_xx(k,i,j),      &
                               g_mu_xy(k,i,j),      &
                               g_mu_xz(k,i,j)
             enddo

             read (101,*)
             do k=1,n_atom
                read (101, *)    &
                               g_mu_yx(k,i,j),      &
                               g_mu_yy(k,i,j),      &
                               g_mu_yz(k,i,j)
             enddo

             read (101,*)
             do k=1,n_atom
                read (101, *)    &
                               g_mu_zx(k,i,j),      &
                               g_mu_zy(k,i,j),      &
                               g_mu_zz(k,i,j)
             enddo

          else if(i .gt. j) then
             do k=1,n_atom

                g_mu_xx(k,i,j) = g_mu_xx(k,j,i)
                g_mu_xy(k,i,j) = g_mu_xy(k,j,i)
                g_mu_xz(k,i,j) = g_mu_xz(k,j,i)

                g_mu_yx(k,i,j) = g_mu_yx(k,j,i)
                g_mu_yy(k,i,j) = g_mu_yy(k,j,i)
                g_mu_yz(k,i,j) = g_mu_yz(k,j,i)

                g_mu_zx(k,i,j) = g_mu_zx(k,j,i)
                g_mu_zy(k,i,j) = g_mu_zy(k,j,i)
                g_mu_zz(k,i,j) = g_mu_zz(k,j,i)

             enddo

          endif

        enddo
       enddo

       close(101)

       do i=1,n_state_ele
         do j=1,n_state_ele
            do k=1,n_atom
               g_mu_all_x(k,i,j) = g_mu_xx(k,i,j)*efield_x + g_mu_yx(k,i,j)*efield_y + g_mu_zx(k,i,j)*efield_z
               g_mu_all_y(k,i,j) = g_mu_xy(k,i,j)*efield_x + g_mu_yy(k,i,j)*efield_y + g_mu_zy(k,i,j)*efield_z
               g_mu_all_z(k,i,j) = g_mu_xz(k,i,j)*efield_x + g_mu_yz(k,i,j)*efield_y + g_mu_zz(k,i,j)*efield_z
            enddo
          enddo
       enddo


! without permanent dipole
       if(trim(adjustl(label_permanent_dipole))=="False")then
          do i=1,n_state_ele
             do k=1,n_atom
                g_mu_all_x(k,i,i) = 0.d0
                g_mu_all_y(k,i,i) = 0.d0
                g_mu_all_z(k,i,i) = 0.d0
             enddo
          enddo
       endif
!!!!!!!!

       do i=1,n_atom
          g_self_int_x(i,:,:) = matmul(mu_learned_all,g_mu_all_x(i,:,:)) + matmul(g_mu_all_x(i,:,:),mu_learned_all) 
          g_self_int_y(i,:,:) = matmul(mu_learned_all,g_mu_all_y(i,:,:)) + matmul(g_mu_all_y(i,:,:),mu_learned_all) 
          g_self_int_z(i,:,:) = matmul(mu_learned_all,g_mu_all_z(i,:,:)) + matmul(g_mu_all_z(i,:,:),mu_learned_all) 
       enddo

       if(kkr_mu==1)then
          call construct_Hcav(Hcav,n_state,n_atom,n_state_ele,pes_all_ele,mu_learned_all,gc,wc,Hamitonian_Compoment,self_int_all)
       else if(kkr_mu==0)then
          call construct_Hcav(Hcav,n_state,n_atom,n_state_ele,pes_all_ele,mu_all,gc,wc,Hamitonian_Compoment,self_int_all)
       endif

       write (file_save_hcav, 9997) it, time*TOFS, &
               ((Hcav(i,j),i=1,n_state),j=1,n_state)
       call print_matrix("Hcav",n_state,n_state,Hcav)

9997    format(i10, 1x, 100(f20.10, 1x))

       call calc_gHcav (gHcav_x, &
                        gHcav_y, &
                        gHcav_z, &
                        pes_all_ele, &
                        Hcav, &
                        n_state, &
                        n_atom, &
                        n_state_ele, &
                        gra_all_ele_x, &
                        gra_all_ele_y, &
                        gra_all_ele_z, &
                        nac_ele_x, &
                        nac_ele_y, &
                        nac_ele_z, &
                        NACv_cav_x, &
                        NACv_cav_y, &
                        NACv_cav_z, &
                        g_mu_all_x, &
                        g_mu_all_y, &
                        g_mu_all_z, &
                        gc,   &
                        wc,   &
                        af_cut, &
                        Hamitonian_Compoment, &
                        g_self_int_x, &
                        g_self_int_y, &
                        g_self_int_z)


       call diag_real(Hcav,n_state,pes_all,U_ss)

       call calc_pola_NACv(n_state, &
                           n_atom, &
                           pes_all, &
                           gHcav_x,&
                           gHcav_y,&
                           gHcav_z,&
                           gra_all_x,&
                           gra_all_y,&
                           gra_all_z,&
                           nac_x, &
                           nac_y, &
                           nac_z, &
                           U_ss)

       deallocate(gra_all_ele_x)
       deallocate(gra_all_ele_y)
       deallocate(gra_all_ele_z)
       deallocate(nac_ele_x)
       deallocate(nac_ele_y)
       deallocate(nac_ele_z)
       deallocate(mu_x)
       deallocate(mu_y)
       deallocate(mu_z)
       deallocate(mu_learned_x)
       deallocate(mu_learned_y)
       deallocate(mu_learned_z)
       deallocate(mu_all)
       deallocate(mu_learned_all)

       deallocate(g_mu_xx)
       deallocate(g_mu_xy)
       deallocate(g_mu_xz)

       deallocate(g_mu_yx)
       deallocate(g_mu_yy)
       deallocate(g_mu_yz)

       deallocate(g_mu_zx)
       deallocate(g_mu_zy)
       deallocate(g_mu_zz)

       deallocate(g_mu_all_x)
       deallocate(g_mu_all_y)
       deallocate(g_mu_all_z)

       deallocate(self_int_all)
       deallocate(g_self_int_x)
       deallocate(g_self_int_y)
       deallocate(g_self_int_z)

       return
 
       end 
