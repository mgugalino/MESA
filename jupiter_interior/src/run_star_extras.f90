! ***********************************************************************
!  src/run_star_extras.f90 mugalino
!  Written for ASTR 635 Exoplanets, UMD College Park 
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         ! added a new prescription for heating in the planetary interior
         s% other_energy => other_interior_heating 
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 4
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: rcb_index
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'Prcb_ledoux'
         rcb_index = minloc(s% gradL - s% gradr, dim=1, mask=(s% gradL - s% gradr > 0.0))
         vals(1) = s%pgas(rcb_index)

         names(2) = 'Rrcb_ledoux'
         vals(2) = s%r(rcb_index)

         names(3) = 'Prcb_schw'
         rcb_index = minloc(s% grada - s% gradr, dim=1, mask=(s% grada - s% gradr > 0.0))
         vals(3) = s%pgas(rcb_index)

         names(4) = 'Rrcb_schw'
         vals(4) = s%r(rcb_index)
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step

      ! adds another source of heating in the interior
      ! follows prescription for depth dependent heating
      subroutine other_interior_heating(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: cell_idx, min_idx
         real(dp) :: Lirr, P_dep, sigma_dep, r_dep, gamma

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! Irradiation luminosity derived from dayside flux, dayside flux * (pi*r**2)
         Lirr = s%irradiation_flux * pi * s%r(1)**2 ! irradiation_flux is the dayside flux (erg cm-2 s-1)
         P_dep = s%x_ctrl(1) ! heating pressure in dyne cm-2
         gamma = s%x_ctrl(2) ! gamma = Gamma / Lstar

         min_idx = minloc(s%pgas, dim=1, mask=(s%pgas >= P_dep)) ! index of pressure closest to Pdep
         sigma_dep = s%scale_height(min_idx) / 2.0 ! scale height at roughly P = Pdep
         r_dep = s%r(min_idx)

         do cell_idx=1, s% nz
            s% extra_heat(cell_idx) = gamma * Lirr * (1./ sqrt(2.0 * pi * sigma_dep**2.0)) * (1.0 / (s%rho(cell_idx) * 4.0 * pi * s%r(cell_idx)**2.0)) &
            * exp(-1.0 * (s%r(cell_idx) - r_dep)**2.0 / (2 * pi * sigma_dep**2))
         end do

      end subroutine other_interior_heating
      
      end module run_star_extras
      
