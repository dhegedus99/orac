subroutine get_cm_td(cyear, cdoy, cmonth, modis_surf_path, &
     modis_brdf_path, occci_path, imager_flags, imager_geolocation, &
     imager_angles, channel_info, ecmwf, assume_full_path, include_full_brdf, &
     use_occci, use_swansea_climatology, swan_g, verbose, mcd43_maxqa, &
     surface, source_atts, netcdf_info)

   use channel_structures_m
   use cox_munk_m
   use cox_munk_constants_m
   use ecmwf_m, only : ecmwf_t
   use fill_grid_m
   use imager_structures_m
   use interpol_m
   use mcd43c_m
   use ocean_colour_m
   use system_utils_m
   use preproc_constants_m
   use preproc_structures_m
   use ross_thick_li_sparse_r_m
   use source_attributes_m
   use surface_structures_m
   use netcdf
   use orac_ncdf_m
   use global_attributes_m
   use netcdf_output_m

   implicit none

   ! Input variables
   character(len=*),           intent(in)    :: cyear
   character(len=*),           intent(in)    :: cdoy
   character(len=*),           intent(in)    :: cmonth
   character(len=*),           intent(in)    :: modis_surf_path
   character(len=*),           intent(in)    :: modis_brdf_path
   character(len=*),           intent(in)    :: occci_path
   type(imager_flags_t),       intent(in)    :: imager_flags
   type(imager_geolocation_t), intent(in)    :: imager_geolocation
   type(imager_angles_t),      intent(in)    :: imager_angles
   type(channel_info_t),       intent(in)    :: channel_info
   type(ecmwf_t),              intent(in)    :: ecmwf
   logical,                    intent(in)    :: assume_full_path
   logical,                    intent(in)    :: include_full_brdf
   logical,                    intent(in)    :: use_occci
   logical,                    intent(in)    :: use_swansea_climatology
   real,                       intent(in)    :: swan_g
   logical,                    intent(in)    :: verbose
   integer,                    intent(in)    :: mcd43_maxqa
   type(surface_t),            intent(inout) :: surface
   type(source_attributes_t),  intent(inout) :: source_atts
   type(netcdf_output_info_t), intent(inout) :: netcdf_info

   ! Local variables

   ! General
   integer                           :: i, j,k, ii, jj, kk, i_oc, i_view
   logical                           :: flag
   integer                           :: nsea, nland
   integer                           :: seacount
   integer                           :: lndcount
   logical                           :: read_qa
   logical,            allocatable   :: mask(:,:)
   integer,            allocatable   :: bands(:)
   integer,            allocatable   :: band_to_sw_index(:)
   real,               allocatable   :: solza(:), satza(:)
   real,               allocatable   :: solaz(:), relaz(:)
   type(interpol_t),   allocatable   :: interp(:)

   ! Band definitions
   integer, parameter :: n_modbands = 7
   integer, target    :: modbands(n_modbands) = [1, 2, 3, 4, 5, 6, 7]
   integer, parameter :: n_swanbands = 4
   integer, target    :: swanbands(n_swanbands) = [4, 1, 2, 6]
   integer, parameter :: n_coxbands = 9
   integer, parameter :: coxbands(n_coxbands) = [1, 2, 3, 4, 5, 6, 7, 8, 9]

   ! Land surface reflectance
   character(len=path_length)        :: modis_surf_path_file
   character(len=path_length)        :: modis_brdf_path_file
   type(mcd43c1_t)                   :: mcdc1
   type(mcd43c3_t)                   :: mcdc3
   integer                           :: n_inbands
   integer,            pointer       :: in_bands(:)
   integer                           :: n_bands
   real,               allocatable   :: tmp_data(:,:)
   integer(kind=byte), allocatable   :: fg_mask(:,:)
   real,               allocatable   :: wsalnd(:,:)
   real,               allocatable   :: wgtlnd(:,:,:)
   real,               allocatable   :: rholnd(:,:,:), tmprho(:,:,:)
   real                              :: tmp_val, tmp_p(2)
   integer                           :: stat

   ! Sea surface reflectance
   real,               allocatable   :: u10sea(:), v10sea(:)
   real,               allocatable   :: latsea(:), lonsea(:)
   real,               allocatable   :: refsea(:,:)
   real,               allocatable   :: rhosea(:,:,:)
   type(cox_munk_shared_geo_wind_t)  :: cox_munk_shared_geo_wind
   type(ocean_colour_t), allocatable :: ocean_colour(:,:)
#ifdef __INTEL_COMPILER
   type(ocean_colour_t), allocatable :: ocean_colour2(:,:)
#endif
   integer(kind=lint), allocatable   :: dummy_chan_vec1d(:)
   integer,            allocatable   :: logical_column_as_int(:)

      if (verbose) write(*,*) &
         'get_surface_reflectance(): Beginning SEA SURFACE REFLECTANCE ', &
         'calculation'

      nsea=10000
      ! Allocate and populate the local arrays required for sea pixels
      allocate(solza(nsea))
      allocate(satza(nsea))
      allocate(solaz(nsea))
      allocate(relaz(nsea))
      allocate(u10sea(nsea))
      allocate(v10sea(nsea))
      allocate(refsea(n_coxbands,nsea))
      allocate(ocean_colour(n_coxbands,nsea))
      if (include_full_brdf) then
         allocate(rhosea(n_coxbands,nsea,4))
      end if
      if (use_occci) then
         allocate(latsea(nsea))
         allocate(lonsea(nsea))
      end if



      !call generate_uniform(nsea, 0., 20., satza)
      !call generate_uniform(nsea, 0., 20., solza)
      !call generate_uniform(nsea, 0., 20., solaz)
      !call generate_uniform(nsea, 0., 20., relaz)
      print*, minval(satza), maxval(satza)
      print*, minval(solza), maxval(solza)
      print*, minval(solaz), maxval(solaz)
      print*, minval(relaz), maxval(relaz)
      do i=1, nsea
          call correlated_angles(solza(i), satza(i), solaz(i), relaz(i))
          !print*, solza
      end do

      print*, minval(satza), maxval(satza)
      print*, minval(solza), maxval(solza)
      print*, minval(solaz), maxval(solaz)
      print*, minval(relaz), maxval(relaz)
       
      call generate_gaussian(0.0, 5., nsea, u10sea)
      call generate_gaussian(0.0, 5., nsea, v10sea)
      print*, minval(u10sea), maxval(u10sea)
      print*, minval(v10sea), maxval(v10sea)

      print*, solza
      call generate_correlated_lognormal(nsea, 0.0536, 0.1787, 0.0029, 0.0162, 0.57, ocean_colour(1,:)%totabs, ocean_colour(1,:)%totbsc) 
      call generate_correlated_lognormal(nsea, 0.0594, 0.0689, 0.0022, 0.0152, 0.67, ocean_colour(2,:)%totabs, ocean_colour(2,:)%totbsc) 
      call generate_correlated_lognormal(nsea, 0.5756, 0.4199, 0.0018, 0.0156, 0.02, ocean_colour(3,:)%totabs, ocean_colour(3,:)%totbsc)  
      print*, minval(ocean_colour(1,:)%totabs), maxval(ocean_colour(1,:)%totabs), minval(ocean_colour(1,:)%totbsc), maxval(ocean_colour(1,:)%totbsc)  
      print*, minval(ocean_colour(2,:)%totabs), maxval(ocean_colour(2,:)%totabs), minval(ocean_colour(2,:)%totbsc), maxval(ocean_colour(2,:)%totbsc)  
      print*, minval(ocean_colour(3,:)%totabs), maxval(ocean_colour(3,:)%totabs), minval(ocean_colour(3,:)%totbsc), maxval(ocean_colour(3,:)%totbsc)
     do i = 4, n_coxbands
         ocean_colour(i,:)%totabs = totalabs(i)
         ocean_colour(i,:)%totbsc = totalbsc(i) 
     end do
     do i = 1, 3
         ocean_colour(i,:)%totabs = totalabs(i) + baseabs(i)
         ocean_colour(i,:)%totbsc = totalbsc(i) + basebsc(i) 
     end do
         ocean_colour(:,1)%have_data = .true.

         do j = 1, nsea
            call cox_munk3_calc_shared_geo_wind(solza(j), satza(j), solaz(j), &
                 relaz(j), u10sea(j), v10sea(j), cox_munk_shared_geo_wind)

            do i = 1, n_coxbands
               call cox_munk3(coxbands(i), cox_munk_shared_geo_wind, &
                    ocean_colour(i,j), refsea(i, j))
            end do
         end do

         if (verbose) then
            do i = 1, n_coxbands
               write(*,*) 'Sea refl: sw_index, wvl, min, max = ', i, &
                    minval(refsea(i,:)), maxval(refsea(i,:))
            end do
         end if

         if (include_full_brdf) then
            allocate(tmprho(n_coxbands,nsea,4))
#ifdef __INTEL_COMPILER
            allocate(ocean_colour2(n_coxbands,nsea))
            
            do i = 1, n_coxbands
               ocean_colour2(i,:) = ocean_colour(coxbands(i),:)
            end do
            call cox_munk_rho_0v_0d_dv_and_dd(coxbands, solza(:), satza(:), &
                 solaz(:), relaz(:), ocean_colour2(:,:), &
                 u10sea, v10sea, sreal_fill_value, tmprho(:,:,1), &
                 tmprho(:,:,2), tmprho(:,:,3), tmprho(:,:,4), verbose)
            deallocate(ocean_colour2)
#else
            call cox_munk_rho_0v_0d_dv_and_dd(coxbands, solza(:), satza(:), &
                 solaz(:), relaz(:), ocean_colour(coxbands,:), &
                 u10sea, v10sea, sreal_fill_value, tmprho(:,:,1), &
                 tmprho(:,:,2), tmprho(:,:,3), tmprho(:,:,4), verbose)
#endif

            do i = 1, n_coxbands
               rhosea(coxbands(i),:,:) = tmprho(i,:,:)
            end do
            deallocate(tmprho)
         end if


   print*, netcdf_info%ncid_cox

   if (n_coxbands .ne. 0) then
      allocate(dummy_chan_vec1d(n_coxbands))
      dummy_chan_vec1d=0_lint
      ii = 1
      do i = 1, channel_info%nchannels_total
         if (channel_info%channel_sw_flag(i) .eq. 1) then
            dummy_chan_vec1d(ii)=i
            ii = ii+1
         end if
      end do

      !call ncdf_write_array( &
      !     netcdf_info%ncid_cox, &
      !     'alb_abs_ch_numbers', &
      !     netcdf_info%vid_cox_alb_abs_ch_numbers, &
      !     dummy_chan_vec1d, &
      !     1, 1, n_coxbands)
      deallocate(dummy_chan_vec1d)
   end if

   print*, size(refsea, dim=1), size(refsea, dim=2)
   !print*, ocean_colour%totabs
   if (n_coxbands .ne. 0) then
      call ncdf_write_array( &
           netcdf_info%ncid_cox, &
           'alb_data', &
           netcdf_info%vid_cox_alb_data, &
           refsea, &
           1, 1, n_coxbands, &
           1, 1, nsea)

      print*, 'here'
      call ncdf_write_array( &
           netcdf_info%ncid_cox, &
           'oc_totabs', &
           netcdf_info%vid_oc_ta, &
           ocean_colour%totabs, &
           1, 1, n_coxbands , &
           1, 1, nsea)
      print*, 'here'
      call ncdf_write_array( &
           netcdf_info%ncid_cox, &
           'oc_totbsc', &
           netcdf_info%vid_oc_tb, &
           ocean_colour%totbsc, &
           1, 1, n_coxbands , &
           1, 1, nsea)

      if (include_full_brdf) then
         call ncdf_write_array( &
              netcdf_info%ncid_cox, &
              'rho_0v', &
              netcdf_info%vid_cox_rho_0v_data, &
              rhosea(:,:,1), &
              1, 1, n_coxbands, &
              1, 1, nsea)

         call ncdf_write_array( &
              netcdf_info%ncid_cox, &
              'rho_0d', &
              netcdf_info%vid_cox_rho_0d_data, &
              rhosea(:,:,2), &
              1, 1, n_coxbands, &
              1, 1, nsea)

         call ncdf_write_array( &
              netcdf_info%ncid_cox, &
              'rho_dv', &
              netcdf_info%vid_cox_rho_dv_data, &
              rhosea(:,:,3), &
              1, 1, n_coxbands, &
              1, 1, nsea)

         call ncdf_write_array( &
              netcdf_info%ncid_cox, &
              'rho_dd', &
              netcdf_info%vid_cox_rho_dd_data, &
              rhosea(:,:,4), &
              1, 1, n_coxbands, &
              1, 1, nsea)
      end if
   end if

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'u10_data', &
        netcdf_info%vid_u10_data, &
        u10sea, &
        1, 1, nsea)

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'v10_data', &
        netcdf_info%vid_v10_data, &
        v10sea, &
        1, 1, nsea)

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'solzen_data', &
        netcdf_info%vid_cox_solzen, &
        solza, &
        1, 1, nsea)

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'satzen_data', &
        netcdf_info%vid_cox_satzen, &
        satza, &
        1, 1, nsea)

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'solaz_data', &
        netcdf_info%vid_cox_solaz, &
        solaz, &
        1, 1, nsea)

   call ncdf_write_array( &
        netcdf_info%ncid_cox, &
        'relazi_data', &
        netcdf_info%vid_cox_relazi, &
        relaz, &
        1, 1, nsea)

   call ncdf_close(netcdf_info%ncid_cox, 'netcdf_create_config(): ".cox.nc"')


      ! Tidy up cox_munk input arrays
      deallocate(solza)
      deallocate(satza)
      deallocate(solaz)
      deallocate(relaz)
      deallocate(u10sea)
      deallocate(v10sea)




end subroutine get_cm_td


subroutine generate_lognormal(mean, stddev, n, x)
    implicit none
    integer,  intent(in) :: n
    real,  intent(in) :: mean, stddev
    real :: log_mean, log_stddev
    real,  intent(out) :: x(n)
    integer :: i

    ! Define the parameters for the log-normal distribution
    log_mean = log(mean**2 / sqrt(stddev**2 + mean**2))
    log_stddev = sqrt(log(1.0 + (stddev**2 / mean**2)))

    ! Initialize random number generator
    call random_seed()

    ! Generate log-normal distribution data
    do i = 1, n
        call random_number(x(i))
        ! Generate a standard normal variable with Box-Muller method
        x(i) = sqrt(-2.0 * log(x(i))) * cos(2.0 * 3.141592653589793238 * x(i))
        ! Transform to log-normal
        x(i) = exp(log_mean + log_stddev * x(i))
    end do

end subroutine generate_lognormal

subroutine generate_gaussian(mean, stddev, n, x)
    implicit none
    integer,  intent(in) :: n
    real,  intent(in) :: mean, stddev
    real,  intent(out) :: x(n)
    real :: u1, u2, z0, z1
    integer :: i

    ! Initialize random number generator
    call random_seed()

    ! Generate Gaussian distribution data using Box-Muller transform
    do i = 1, n, 2
        ! Generate two uniform random numbers u1 and u2 in the range (0, 1)
        call random_number(u1)
        call random_number(u2)

        ! Apply Box-Muller transform to get two independent standard normal variables z0 and z1
        z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793238 * u2)
        z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * 3.141592653589793238 * u2)

        ! Scale and shift to get Gaussian variables with specified mean and stddev
        x(i) = mean + stddev * z0
        if (i + 1 <= n) x(i+1) = mean + stddev * z1
    end do


end subroutine generate_gaussian

subroutine generate_uniform(n, min_val, max_val, x)
    implicit none
    integer,  intent(in) :: n
    real, intent(in) :: min_val, max_val
    real,  intent(out) :: x(n)
    integer :: i

    ! Initialize random number generator
    call random_seed()

    ! Generate uniform distribution data in the range [min_val, max_val]
    do i = 1, n
        call random_number(x(i))              ! Generate a uniform random number in the range [0, 1)
        x(i) = min_val + x(i) * (max_val - min_val)  ! Scale to [min_val, max_val]
    end do

end subroutine generate_uniform

subroutine generate_correlated_lognormal(n,mean1, stddev1, mean2, stddev2, corr, x1, x2)
    implicit none
    integer,  intent(in) :: n
    real,  intent(in) ::  mean1, stddev1, mean2, stddev2, corr
    real,  intent(out) :: x1(n), x2(n)
    real :: log_mean1, log_stddev1, log_mean2, log_stddev2
    real :: u1, u2, z1, z2 
    integer :: i

    ! Calculate the parameters for the underlying normal distributions
    log_mean1 = log(mean1**2 / sqrt(stddev1**2 + mean1**2))
    log_stddev1 = sqrt(log(1.0 + (stddev1**2 / mean1**2)))
    log_mean2 = log(mean2**2 / sqrt(stddev2**2 + mean2**2))
    log_stddev2 = sqrt(log(1.0 + (stddev2**2 / mean2**2)))

    ! Initialize random number generator
    call random_seed()

    ! Generate correlated normal variables and transform to log-normal
    do i = 1, n
        ! Generate two independent uniform random numbers u1 and u2
        call random_number(u1)
        call random_number(u2)

        ! Convert to standard normal using Box-Muller transform
        z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793238 * u2)
        z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * 3.141592653589793238 * u2)

        ! Introduce the desired correlation
        z2 = corr * z1 + sqrt(1.0 - corr**2) * z2

        ! Scale and shift to create the two log-normal distributions
        x1(i) = exp(log_mean1 + log_stddev1 * z1)
        x2(i) = exp(log_mean2 + log_stddev2 * z2)
    end do


end subroutine generate_correlated_lognormal


subroutine correlated_angles(solza, satza, solaz,relaz)
    !use lapack95  ! Use LAPACK for Cholesky decomposition
    implicit none

    ! Define correlation matrix and standard deviation array
    real, dimension(4,4) :: R  ! Your correlation matrix
    real, dimension(4) :: mu  ! Mean for each variable (set to 0 for simplicity)
    real, dimension(4) :: sd  ! Standard deviations (adjust if needed)
    real,    intent(out)  :: solza, satza,solaz, relaz
    real                 :: sataz
    real, dimension(4) :: z        ! Array for uncorrelated random normal values
    real, dimension(4) :: x        ! Correlated values after applying Cholesky factor
    integer :: info                   ! Variable for LAPACK function status
    external :: dpotrf


    R = reshape([1.0d0, -0.47d0, -0.42d0, -0.62d0, &
                 -0.47d0, 1.0d0, 0.0d0, 0.0d0, &
                 -0.42d0, 0.0d0, 1.0d0, 0.32d0, &
                 -0.62d0, 0.0d0, 0.32d0, 1.0d0], [4,4])

    ! Set means and standard deviations
    mu = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    sd = [1.0d0, 1.0d0, 1.0d0, 1.0d0]


    ! Generate Cholesky decomposition of correlation matrix
    call dpotrf('U', 4, R, 4, info)
    !if (info /= 0) stop 'Cholesky decomposition failed'

    ! Generate uncorrelated standard normal random numbers
    call random_number(z)  ! z is a 4-element array of random values
    z = sqrt(-2.0d0 * log(z)) * cos(2.0d0 * 3.14159d0 * z)  ! Standard normal transformation

    ! Apply Cholesky factor to introduce correlation
    x = matmul(R, z)

    ! Map x values to angles
    solza = map_to_range(x(1), 0.01, 70.)
    satza = map_to_range(x(2), 0., 70.)
    solaz = map_to_range(x(3), 0., 360.)
    sataz = map_to_range(x(4), 0., 360.)
    relaz = abs(solaz - sataz)
    if (relaz > 180) relaz = 360 - relaz

contains

    ! Internal function to map a value from standard normal to a specified range
    real function map_to_range(value, min_val, max_val)
        real, intent(in) :: value, min_val, max_val
        map_to_range = min_val + (max_val - min_val) * (0.5d0 * (1.0d0 + erf(value / sqrt(2.0d0))))
    end function map_to_range

end subroutine correlated_angles
