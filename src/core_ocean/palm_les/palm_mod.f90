!> @file palm.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Large-Eddy Simulation (LES) model for the convective boundary layer,
!> optimized for use on parallel machines (implementation realized using the
!> Message Passing Interface (MPI)). The model can also be run on vector machines
!> (less well optimized) and workstations. Versions for the different types of
!> machines are controlled via cpp-directives.
!> Model runs are only feasible using the ksh-script mrun.
!>
!> @todo create routine last_actions instead of calling lsm_last_actions etc.
!> @todo move chem_init call to init_3d_model or to check_parameters
!------------------------------------------------------------------------------!

 subroutine palm(T_mpas,S_mpas,U_mpas,V_mpas,lt_mpas, &
             lat_mpas,nVertLevels,wtflux,wtflux_solar, wsflux,uwflux, &
             vwflux,fac,dep1,dep2,dzLES,nzLES,        &
             endTime, dtDisturb, tIncrementLES,sIncrementLES,             &
             uIncrementLES,vIncrementLES,tempLES,    &
             salinityLES, uLESout, vLESout, dtLS, first, zLES, &
             disturbMax, disturbAmp, disturbBot, disturbTop, disturbNblocks, &
             botDepth, timeAv)


    USE arrays_3d

    USE control_parameters

    USE configure_3D_MODEL

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices

!    USE netcdf_data_input_mod,                                                 &
!        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
!               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE fft_xy,                                                                 &
        ONLY: fft_finalize

    USE kinds

    USE ppr_1d

    USE pegrid

    USE random_generator_parallel, ONLY: deallocate_random_generator

    use surface_mod

    use tridia_solver, ONLY: tridia_deallocate

    use turbulence_closure_mod, ONLY: tcm_deallocate_arrays

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

    use statistics

    #if defined( __cudaProfiler )
        USE cudafor
    #endif

    IMPLICIT NONE

!
! -- Variables from MPAS
   integer(iwp) :: nVertLevels, il, jl, jloc, kl, knt, nzLES, iz, disturbNblocks
   integer(iwp) :: nzMPAS, zmMPASspot, zeMPASspot
   Real(wp),intent(in)                             :: dtLS
   logical,intent(in)                              :: first
   Real(wp),dimension(nVertLevels),intent(inout)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),dimension(nVertLevels),intent(inout)   :: tIncrementLES, sIncrementLES, &
                                                        uIncrementLES, vIncrementLES
   Real(wp),dimension(nzLES),intent(out)           :: tempLES, salinityLES, zLES
   Real(wp),dimension(nzLES),intent(out)           :: uLESout, vLESout
   Real(wp),dimension(nVertLevels),intent(in)      :: lt_mpas
   Real(wp),allocatable,dimension(:)   :: T_mpas2, S_mpas2, U_mpas2, V_mpas2
   Real(wp),allocatable,dimension(:)   :: Tles, Sles, Ules, Vles, zmid, zedge
   real(wp),allocatable,dimension(:)   :: zeLES, wtLES, wsLES, wuLES, wvLES
   Real(wp),allocatable,dimension(:)   :: zeLESInv
   Real(wp) :: wtflux, wsflux, uwflux, vwflux, dzLES, z_fac, z_frst, z_cntr
   real(wp) :: z_fac1, z_fac2, z_facn, tol, test, lat_mpas, fac, dep1, dep2
   real(wp) :: dtDisturb, endTime, thickDiff, disturbMax, disturbAmp
   real(wp) :: disturbBot, disturbTop, botDepth, timeAv
   real(wp) :: wtflux_solar, sumValT, sumValS, sumValU, sumValV, thickVal
   !
!-- Local variables
   CHARACTER(LEN=9)  ::  time_to_string  !<
   CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
   INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
   INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI
   Real(wp) :: coeff1, coeff2

!-- arrays and parameters for PPR remapping

   integer, parameter :: nvar = 4
   integer, parameter :: ndof = 1
   real(wp) :: fLES(ndof, nvar, nzLES)
   real(wp),allocatable :: fMPAS(:,:,:)
   type(rmap_work) :: work
   type(rmap_opts) :: opts
   type(rcon_ends) :: bc_l(nvar)
   type(rcon_ends) :: bc_r(nvar)

!-- this specifies options for the method, here is quartic interp
   opts%edge_meth = p5e_method
   opts%cell_meth = pqm_method
   opts%cell_lims = null_limit

   bc_l(:)%bcopt = bcon_loose
   bc_r(:)%bcopt = bcon_loose

   call init_control_parameters

   dt_disturb = dtDisturb
   end_time = endTime
   ideal_solar_division = fac
   ideal_solar_efolding1 = dep1
   ideal_solar_efolding2 = dep2
   wb_solar = wtflux_solar
   nz = nzLES
   disturb_nblocks = disturbNblocks
   dt_ls = dtLS
   dt_avg = timeAv

   disturbance_level_b = disturbBot
   disturbance_level_t = disturbTop
   disturbance_amplitude = disturbAmp
   disturbance_energy_limit = disturbMax

   allocate(zmid(nVertLevels),zedge(nVertLevels+1))
   allocate(T_mpas2(nVertLevels),S_mpas2(nVertLevels),U_mpas2(nVertLevels))
   allocate(V_mpas2(nVertLevels))

!   lt_mpas(:) = 50.0_wp
   zmid(1) = -0.5_wp*lt_mpas(1)
   zedge(1) = 0

   do il=2,nVertLevels
      zmid(il) = zmid(il-1) - 0.5*(lt_mpas(il-1) + lt_mpas(il))
      zedge(il) = zedge(il-1) - lt_mpas(il-1)
   enddo

   zedge(nvertLevels+1) = zedge(nVertLevels) - lt_mpas(nVertLevels)

   do il=1,nVertLevels
     if(zmid(il) < botDepth) then
       zmMPASspot = il
       nzMPAS = il
       exit
     endif
   enddo

   do il=1,nVertLevels
     if(zedge(il) < botDepth) then
       zeMPASspot = il
       exit
     endif
   enddo

   botDepth = zedge(zeMPASspot)

   allocate(fMPAS(ndof, nvar, nzMPAS))

!   U_mpas(:) = 0.0_wp
!   V_mpas(:) = 0.0_wp
!   S_mpas(:) = 34.0_wp!
!   do i=1,nVertLevels
!     T_mpas(i) = 293.15 + 0.005*zmid(i)
!   enddo

   ! will need to interpolate profiles to an LES grid, or maybe assume the same???
   ! at end assign fluxes tot mpas variables and end routine.

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because
!-- it will be defined in init_pegrid but is used before in cpu_log.
    CALL MPI_INIT( ierr )

    comm_palm = MPI_COMM_WORLD
    comm2d = comm_palm
!
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!
#endif

!TODO add check for right / acceptable range.
    top_momentumflux_u = uwflux
    top_momentumflux_v = vwflux
    top_heatflux =  wtflux
    top_salinityflux = -wsflux
    latitude = lat_mpas * 180.0 / pi

    !TODO ooverride the LES setting from a namelist
!
!-- Initialize measuring of the CPU-time remaining to the run
    CALL local_tremain_ini
!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )

    !
!-- Read control parameters from NAMELIST files and read environment-variables
    CALL parin

!-- Determine processor topology and local array indices
    CALL init_pegrid
    allocate(zu(nzb:nzt+1),zeLES(nzb-1:nzt+1),Tles(0:nzLES+1),Sles(0:nzLES+1))
    allocate(zw(nzb:nzt+1),Ules(0:nzLES+1),Vles(0:nzLES+1))
    allocate(zeLESinv(nzb-1:nzt+1))

    nzt = nzLES
    ! construct a stretched stretched grid
    z_cntr = botDepth
    z_frst = -dzLES
    z_fac1 = z_cntr / z_frst
    z_fac2 = 1.0_wp / REAL(nzt,kind=wp)
    z_fac = 1.10_wp
    tol = 1.0E-10_wp
    test = 10.00_wp
    knt = 0

    do while (test > tol)
      knt = knt + 1
      z_facn = (z_fac1*(z_fac - 1.0_wp) + 1.0_wp)**z_fac2
      test = abs(1.0 - z_facn / z_fac)
      if(knt .gt. 500) THEN
        print *, 'cannot find stretching factor,'
        print *, 'z_fac = ',z_fac, 'z_facn = ',z_facn, 'knt = ',knt
        stop
      ENDIF
      z_fac = z_facn
    enddo

    zeLES(nzt+1) = dzLES
    zeLES(nzt) = 0.0_wp
    zeLES(nzt-1) = -dzLES
    iz = 2
    do il = nzt-2,nzb,-1
      zeLES(il) = zeLES(nzt-1)*(z_fac**(real(iz,kind=wp)) - 1.0_wp) / (z_fac - 1.0_wp)
      iz = iz + 1
    enddo

!    zeLES(nzt+1) = 1.0_wp
!    zeLES(nzt) = 0.0_wp
!    do il = nzt-1,nzb,-1
!       zeLES(il) = zeLES(il+1) - 1.0_wp
!    enddo

    zeLES(nzb-1) = max(z_cntr,zeLES(nzb) - (zeLES(nzb+1) - zeLES(nzb)))

    do il = nzt,nzb,-1
      zu(il) = 0.5*(zeLES(il) + zeLES(il-1))
    enddo
    zu(nzt+1) = dzLES 
    zLES(nzb+1:nzt) = zu(nzb+1:nzt)

     call work%init(nzLES+1,nvar,opts) 
   !
!-- Generate grid parameters, initialize generic topography and further process
!-- topography information if required
    CALL init_grid
!-- Check control parameters and deduce further quantities
    CALL check_parameters
! interpolate mpas data to les and send to init variable
!
    CALL allocate_3d_arrays

!-- Initialize all necessary variables
    CALL init_3d_model

    il = 1
    jl=nzt
!    do while(zu(jl) > zmid(2))
!       tLSforcing(jl) = T_mpas(1) + 273.15_wp
!       uLSforcing(jl) = U_mpas(1)
!       vLSforcing(jl) = V_mpas(1)
!       sLSforcing(jl) = S_mpas(1)
!       jl = jl - 1
!    enddo

    fMPAS(1,1,:) = T_mpas(1:nzMPAS)
    fMPAS(1,2,:) = S_mpas(1:nzMPAS)
    fMPAS(1,3,:) = U_mpas(1:nzMPAS)
    fMPAS(1,4,:) = V_mpas(1:nzMPAS)

    jl=1
    do il = nzt,nzb-1,-1
      zeLESinv(jl) = zeLES(il)
      jl = jl + 1
    enddo

    call rmap1d(nzMPAS+1,nzLES+1,nvar,ndof,abs(zedge(1:nzMPAS+1)),abs(zeLESinv(1:nzLES+1)), &
                fMPAS, fLES, bc_l, bc_r, work, opts)

    jl = 1
    do il = nzt,nzb+1,-1
      tLSforcing(il) = fLES(1,1,jl) + 273.15
      sLSforcing(il) = fLES(1,2,jl)
      uLSforcing(il) = fLES(1,3,jl)
      vLSforcing(il) = fLES(1,4,jl)
      jl = jl + 1
    enddo

!    jloc = jl
!    il = 1
!    do while (il < nVertLevels-1)
!      do jl=nzt,nzb+1,-1
!         if(zu(jl) < zmid(il+1)) then
!           il = il+1
!           il = min(il,nVertLevels-1)
!         endif
!
!         coeff2 = (T_mpas(il) - T_mpas(il+1)) / (zmid(il) - zmid(il+1))
!         coeff1 = T_mpas(il+1) - coeff2*zmid(il+1)
!         tLSforcing(jl) = coeff2*zu(jl) + coeff1 + 273.15_wp
!         
!         coeff2 = (S_mpas(il) - S_mpas(il+1)) / (zmid(il) - zmid(il+1))
!         coeff1 = S_mpas(il+1) - coeff2*zmid(il+1)
!         sLSforcing(jl) = coeff2*zu(jl) + coeff1
!
!         coeff2 = (U_mpas(il) - U_mpas(il+1)) / (zmid(il) - zmid(il+1))
!         coeff1 = U_mpas(il+1) - coeff2*zmid(il+1)
!         uLSforcing(jl) = coeff2*zu(jl) + coeff1
!
!         coeff2 = (V_mpas(il) - V_mpas(il+1)) / (zmid(il) - zmid(il+1))
!         coeff1 = V_mpas(il+1) - coeff2*zmid(il+1)
!         vLSforcing(jl) = coeff2*zu(jl) + coeff1

!      enddo
!    enddo

if(first) then
       do jl=nzt,nzb+1,-1
          pt(jl,:,:) = tLSforcing(jl)
          sa(jl,:,:) = sLSforcing(jl)
          u(jl,:,:) = uLSforcing(jl)
          v(jl,:,:) = vLSforcing(jl)
          tempLES(jl) = tLSforcing(jl)
          salinityLES(jl) = sLSforcing(jl)
          uLESout(jl) = uLSforcing(jl)
          vLESout(jl) = vLSforcing(jl)
          uProfileInit(jl) = uLSforcing(jl)
          vProfileInit(jl) = vLSforcing(jl)
          tProfileInit(jl) = tLSforcing(jl)
          sProfileInit(jl) = sLSforcing(jl)
       enddo
else
        do jl = nzt,nzb+1,-1
          pt(jl,:,:) = tempLES(jl) + 273.15_wp
          sa(jl,:,:) = salinityLES(jl)
          u(jl,:,:) = uLESout(jl)
          v(jl,:,:) = vLESout(jl)
          uProfileInit(jl) = uLSforcing(jl)
          vProfileInit(jl) = vLSforcing(jl)
          tProfileInit(jl) = tLSforcing(jl)
          sProfileInit(jl) = sLSforcing(jl)
          enddo

  endif

    pt(nzb,:,:) = pt(nzb+1,:,:)
    sa(nzb,:,:) = sa(nzb+1,:,:)
    u(nzb,:,:) = u(nzb+1,:,:)
    v(nzb,:,:) = v(nzb+1,:,:)

    pt(nzt+1,:,:) = pt(nzt,:,:)
    sa(nzt+1,:,:) = sa(nzt,:,:)
    u(nzt+1,:,:) = u(nzt,:,:)
    v(nzt+1,:,:) = v(nzt,:,:)

    !
!-- Output of program header
!    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )

!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )

!    IF ( do3d_at_begin )  THEN
!       CALL data_output_3d( 0 )
!    ENDIF

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStart()
#endif
!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStop()
#endif

!
!-- If required, repeat output of header including the required CPU-time
!    IF ( myid == 0 )  CALL header
!
!-- If required, final  user-defined actions, and
!-- last actions on the open files and close files. Unit 14 was opened
!-- in wrd_local but it is closed here, to allow writing on this
!-- unit in routine user_last_actions.

    CALL cpu_log( log_point(4), 'last actions', 'start' )

    CALL close_file( 0 )

    CALL cpu_log( log_point(4), 'last actions', 'stop' )

!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
!    CALL cpu_statistics

    ! need tto interpolate back to mpas for fluxes, include sgs terms?
!    Tles = hom(:,1,4,statistic_regions) 
!    Sles = hom(:,1,23,statistic_regions)
!    Ules = hom(:,1,1,statistic_regions)
!    Vles = hom(:,1,2,statistic_regions)

    if(average_count_meanpr /= 0) then
 
       meanFields_avg(:,1) = meanFields_avg(:,1) / average_count_meanpr
       meanFields_avg(:,2) = meanFields_avg(:,2) / average_count_meanpr
       meanFields_avg(:,3) = meanFields_avg(:,3) / average_count_meanpr
       meanFields_avg(:,4) = meanFields_avg(:,4) / average_count_meanpr

    endif

    Tles = meanFields_avg(:,3)
    Sles = meanFields_avg(:,4)
    Ules = meanFields_avg(:,1)
    Vles = meanFields_avg(:,2)
 ! need to integrate over layers in mpas to get increments

 if(minval(tempLES) < 100.0_wp) tempLES(:) = tempLES(:) + 273.15_wp
    tProfileInit(1:) = tempLES(:)
    sProfileInit(1:) = salinityLES(:)
    uProfileInit(1:) = uLESout(:)
    vProfileInit(1:) = vLESout(:)

!    jl=1
!    do il=nzt,nzb+1,-1
!      fLES(1,1,jl) = Tles(il) !tProfileInit(il)
!      fLES(1,2,jl) = Sles(il) !sProfileInit(il)
!      fLES(1,3,jl) = Ules(il) !uProfileInit(il)
!      fLES(1,4,jl) = Vles(il) !vProfileInit(il)
!      jl = jl+1
!    enddo

!    call rmap1d(nzLES+1,nzMPAS+1,nvar,ndof,abs(zeLESinv(1:nzLES+1)),abs(zedge(1:nzMPAS+1)), &
!                fLES,fMPAS,bc_l,bc_r,work,opts)

!    print *, ' '
!    print *, fMPAS(1,1,:) - 273.15 - T_mpas(:nzMPAS)
!    print *, ' '
!    print *, ' '

    jl=1
    do il=nzt,nzb+1,-1
      fLES(1,1,jl) = (Tles(il) - tProfileInit(il)) / dtLS
      fLES(1,2,jl) = (Sles(il) - sProfileInit(il)) / dtLS
      fLES(1,3,jl) = (Ules(il) - uProfileInit(il)) / dtLS
      fLES(1,4,jl) = (Vles(il) - vProfileInit(il)) / dtLS
      jl = jl+1
    enddo

!    print *, fLES(1,1,:)
!    print *, ' '
    call rmap1d(nzLES+1,nzMPAS+1,nvar,ndof,abs(zeLESinv(1:nzLES+1)),abs(zedge(1:nzMPAS+1)), &
                fLES,fMPAS,bc_l,bc_r,work,opts)

    tIncrementLES(:) = 0.0_wp
    sIncrementLES(:) = 0.0_wp
    uIncrementLES(:) = 0.0_wp
    vIncrementLES(:) = 0.0_wp
    do jl=1,nzMPAS
      tIncrementLES(jl) = fMPAS(1,1,jl)
      sIncrementLES(jl) = fMPAS(2,1,jl)
      uIncrementLES(jl) = fMPAS(3,1,jl)
      vIncrementLES(jl) = fMPAS(4,1,jl)
    enddo

!    print *, fMPAS(1,1,:) 
!    stop
 
 !   il = nzt
 !   thickDiff = 0.0_wp
 !   !find stopping spot
 !   do knt = 1,nVertLevels
 !      if(zedge(knt) < botDepth) then
 !              exit
 !      endif
 !   enddo

 !   do jl=1,knt-3
 !     
 !      sumValT = thickDiff*(tProfileInit(il) - Tles(il))
 !      sumValS = thickDiff*(sProfileInit(il) - Sles(il))
 !      sumValU = thickDiff*(uProfileInit(il) - uLES(il))
 !      sumValV = thickDiff*(vProfileInit(il) - vLES(il))
!
!       thickVal = 0.0_wp
!       do while (zw(il-1) >= zedge(jl+1)-1e-6)
!          sumValT = sumValT - (zw(il) - zw(il-1))*(tProfileInit(il) - Tles(il))
!          sumValS = sumValS - (zw(il) - zw(il-1))*(sProfileInit(il) - Sles(il))
!          sumValU = sumValU - (zw(il) - zw(il-1))*(uProfileInit(il) - uLES(il))
!          sumValV = sumValV - (zw(il) - zw(il-1))*(vProfileInit(il) - vLES(il))
!          thickVal = thickVal + (zw(il) - zw(il-1))
!          il = max(il - 1,1)
!       enddo
!       if (thickVal < lt_mpas(jl) -1e-6) then
!          sumValT = sumValT - (lt_mpas(jl) - thickVal)*(tProfileInit(il) - Tles(il))
!          sumValS = sumValS - (lt_mpas(jl) - thickVal)*(sProfileInit(il) - Sles(il))
!          sumValU = sumValU - (lt_mpas(jl) - thickVal)*(uProfileInit(il) - uLES(il))
!          sumValV = sumValV - (lt_mpas(jl) - thickVal)*(vProfileInit(il) - vLES(il))
!          thickDiff = (zw(il) - zw(il-1)) - (lt_mpas(jl) - thickVal) 
!        else
!          thickDiff = 0.0_wp
!        endif
!
!       tIncrementLES(jl) = sumValT / (dtLS*lt_mpas(jl))
!       sIncrementLES(jl) = sumValS / (dtLS*lt_mpas(jl))
!       uIncrementLES(jl) = sumValU / (dtLS*lt_mpas(jl))
!       vIncrementLES(jl) = sumValV / (dtLS*lt_mpas(jl))
!
!    enddo

!    print *, 'il = ',il
!    print *, uProfileInit(nzt-20:nzt) - Ules(nzt-20:nzt)
!    print *, '****************************'
!    print *, '*******************************'
!    print *, uIncrementLES(:20)*dtLS
!    stop
    tempLES = Tles(1:nzLES) - 273.15_wp
    salinityLES = Sles(1:nzLES)
    uLESout = Ules(1:nzLES)
    vLESout = Vles(1:nzLES)
    DEALLOCATE( pt_init, q_init, s_init, ref_state, sa_init, ug,         &
                       u_init, v_init, vg, hom, hom_sum, meanFields_avg )

   deallocate(hor_index_bounds)

    deallocate(zu,zeLES,Tles,Sles)
    deallocate(hyp, Ules,Vles)
    deallocate(ddzu, ddzw, dd2zu, dzu, dzw, zw, ddzu_pres, nzb_s_inner,  &
               nzb_s_outer, nzb_u_inner, nzb_u_outer, nzb_v_inner,       &
               nzb_v_outer, nzb_w_inner, nzb_w_outer, nzb_diff_s_inner,  &
               nzb_diff_s_outer, wall_flags_0, advc_flags_1, advc_flags_2)

    deallocate(u_stk, v_stk, u_stk_zw, v_stk_zw )

    call deallocate_bc
    call deallocate_3d_variables
    call tcm_deallocate_arrays
    if (random_generator == 'random-parallel') call deallocate_random_generator
    call tridia_deallocate

    close(18)

    CALL fft_finalize
#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

END subroutine palm

subroutine init_control_parameters
    USE arrays_3d

    USE control_parameters

    USE statistics, only: flow_statistics_called
    USE kinds


    openfile = file_status(.FALSE.,.FALSE.)

    rayleigh_damping_factor = -1.0_wp
    rayleigh_damping_height = -1.0_wp
    timestep_count = 0
        poisfft_initialized = .FALSE.
        init_fft = .FALSE.
        psolver = 'poisfft'
        momentum_advec = 'ws-scheme'
        loop_optimization = 'vector'
        bc_e_b = 'neumann'
        bc_lr = 'cyclic'
        bc_ns = 'cyclic'
        bc_p_b = 'neumann'
        bc_p_t = 'neumann'
        bc_pt_b = 'neumann'
        bc_pt_t = 'neumann'
        bc_sa_t = 'neumann'
        bc_sa_b = 'neumann'
        bc_uv_b = 'neumann'
        bc_uv_t = 'neumann'
        ibc_uv_b = 1
        coupling_mode = 'uncoupled'
        fft_method = 'temperton-algorithm'
        topography = 'flat'
        initializing_actions = 'set_constant_profiles'
        random_generator = 'numerical-recipes'
        random_generator = 'random-parallel'
        reference_state = 'initial_profile'
        data_output = ' '
        data_output_user = ' '
        doav = ' '
        data_output_masks = ' ' 
        data_output_pr = ' '
        domask = ' '
        do2d = ' '
        do3d = ' '

        do3d_no(0:1) = 0
        meanFields_avg = 0.0_wp

        abort_mode = 1
        time_avg = 0.0_wp
        average_count_pr = 0
        average_count_meanpr = 0
        average_count_3d = 0
        current_timestep_number = 0
        coupling_topology = 0
        dist_range = 0
        doav_n = 0
        dopr_n = 0
        dopr_time_count = 0
        dopts_time_count = 0
        dots_time_count = 0
        dp_level_ind_b = 0 
        dvrp_filecount = 0
        ensemble_member_nr = 0

        iran = -1234567
        length = 0
        io_group = 0
        io_blocks = 1
        masks = 0
        maximum_parallel_io_streams = -1
        mgcycles = 0
        mg_cycles = 4
        mg_switch_to_pe0_level = -1
        ngsrb = 2
        nr_timesteps_this_run = 0
        nsor = 20
        nsor_ini = 100
        normalizing_region = 0
        num_leg = 0
        num_var_fl_user = 0
        nz_do3 = -9999
        y_shift = 0
        mask_size(max_masks,3) = -1
        mask_size_l(max_masks,3) = -1
        mask_start_l(max_masks,3) = -1
        pt_vertical_gradient_level_ind(10) = -9999
        sa_vertical_gradient_level_ind(10) = -9999
        stokes_drift_method = -9999

        dz(10) = -1.0_wp 
        dzconst = 2.5_wp
        dt_disturb = 20.0_wp 
        dt_do3d = 9999999.9_wp
        dt_3d = 0.01_wp

        simulated_time = 0.0_wp
        flow_statistics_called = .FALSE.
        disturbance_created = .FALSE.
        time_disturb = 0.0_wp
        time_dopr = 0.0_wp
        time_dopr_av = 0.0_wp
        time_dots = 0.0_wp
        time_do2d_xy = 0.0_wp
        time_do2d_xz = 0.0_wp
        time_do2d_yz = 0.0_wp
        time_do3d = 0.0_wp
        time_do_av = 0.0_wp
        time_run_control = 0.0_wp

end subroutine init_control_parameters

subroutine deallocate_memory
        
        use pegrid

        deallocate(hor_index_bounds)

end subroutine deallocate_memory