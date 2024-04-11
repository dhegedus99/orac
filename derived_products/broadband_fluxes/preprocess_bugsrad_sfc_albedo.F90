!-------------------------------------------------------------------------------
! Name: preprocess_bugsrad_sfc_albedo
!
! Purpose:
! Interpolate spectral surface black and white sky albedos from MODIS to BUGSrad
! bands. For more accurate results pass 8 spectral channels into the program.
! For less accurate results use heritage channels (4). It will run using both.
!
! Inputs:
! nc_alb: number of input spectral bands
! rho_0d: black sky albedo for each spectral band
! rho_dd: white sky albedo for each spectral band
!
! Outputs:
! rho0d_bugsrad: interpolated albedo to bugsrad band
! rhodd_bugsrad: interpolated albedo to bugsrad band
!
! History:
! xxxx/xx/xx, MC: Initial implementation
!
! Bugs:
! None known.
!-------------------------------------------------------------------------------

subroutine preprocess_bugsrad_sfc_albedo(nc_alb, nc_conf, rho_0d,rho_dd,&
      msi_ch_swflag, msi_ch_modisref,rho0d_bugsrad, rhodd_bugsrad)
   
   implicit none

   ! Input arguments
   integer, intent(in) :: nc_alb, nc_conf
   real, intent(in) :: rho_0d(nc_alb)
   real, intent(in) :: rho_dd(nc_alb)
   real, intent(in) :: msi_ch_swflag(nc_conf) ! flags for sw channels used
   real, intent(in) :: msi_ch_modisref(nc_conf) ! modis ref. bands for albedo

   ! Output arguments
   real, intent(out) :: rho0d_bugsrad(6) ! rho_0d for each BUGSrad band
   real, intent(out) :: rhodd_bugsrad(6) ! rho_dd for each BUGSrad band
   
   
   ! Local variables    
   real :: modBand(nc_alb) ! center of modis bands used
   real :: modBandall(8) ! center of each possible modis band (um)
   real :: bugsBand(6) ! center of each BUGSrad band (um)
   real(kind=8) :: m,b
   integer :: i, tmploc(1)
   integer :: l, j, ii
   integer :: abs_ch_size ! nb. of modis channels to ignore in interpolation
   real :: bb
   real, allocatable :: masked_rho_0d(:), masked_rho_dd(:)
   integer :: modID1(6),modID2(6), firstmixedch
   real, allocatable :: abs_ch(:) ! modis channels to ignore in interpolation
   real, allocatable :: modBandcorr(:), rho_0dcorr(:), rho_ddcorr(:) ! corrected 
   ! modis bands, albedo values based on channels to ignore

   ! Center location of each BUGsrad band (um)
   data bugsBand/0.4445,0.994,1.602,1.8995,3.0045,3.7545/
   
   ! Center location of each possible modis band (um)
   data ModBandall/0.67,0.87,0.47,0.55,1.24,1.6,2.13,3.7/
   
   ! We expect the nb. of albedo channels to be equal to the nb. of sw channels
   ! see preprocessing/get_surface_reflectance.F90
   if (count(msi_ch_swflag/=0) .ne. nc_alb) then
   write(*,*) 'Error: Albedo file length and number of SW channels do not match'
   end if
   
   ! Check the nb. of mixed channels used for albedo 
   write(*,*) 'nb. of mixed channels', &
              (count(msi_ch_swflag/=0)) - (count(msi_ch_modisref>0))
   
   ! Get the center location of each MODIS reference band for albedo
   ! if no mixed channels are used
   if (count(msi_ch_modisref>0) .eq. nc_alb) then
      write(*,*) 'No mixed channels used for albedo'
      do l=1,nc_alb
         !if (msi_ch_modisref(l)>0) then
            modBand(l) = modBandall(int(msi_ch_modisref(l)))
         !end if
      end do
   ! If the length of the albedo file has one extra channel (other than the 
   ! sw channels with MODIS reference bands), we assume it comes from the 
   ! 3.7 mixed channel (see in preprocessing/get_surface_reflectance)
   else if (count(msi_ch_modisref>0) .eq. (nc_alb-1))  then
      write(*,*) 'One additional mixed channel is used for albedo'
      write(*,*) msi_ch_modisref, minloc(abs(msi_ch_modisref-0), dim=1)
      firstmixedch = minloc(abs(msi_ch_modisref-0), dim=1)
      do l=1,firstmixedch-1!nc_alb-1
         modBand(l) = modBandall(int(msi_ch_modisref(l)))
      end do
      modBand(nc_alb) = 3.7
   ! If more than 1 mixed channels are used, we don't know
   ! which channels those are.
   else
      write(*,*) 'Error: Unexpected number of mixed channels used for albedo'
   end if
   
   ! Select modis absorption channels to ignore 
   ! (uncomment and edit if needed)
   !abs_ch_size = 1
   !allocate(abs_ch(abs_ch_size))
   !abs_ch(1) = 0.87
   !abs_ch(2) = 3.7
   
   ! Ignore any absorption channels specified above
   if (allocated(abs_ch)) then
       allocate(modBandcorr(nc_alb-abs_ch_size))
       allocate(rho_0dcorr(nc_alb-abs_ch_size))
       allocate(rho_ddcorr(nc_alb-abs_ch_size))
       ii=1
       do i=1,size(modBand)
          if (any(modBand(i).eq.abs_ch)) then
              cycle
          else
              modBandcorr(ii) = modBand(i)
              rho_0dcorr(ii) = rho_0d(i)
              rho_ddcorr(ii) = rho_dd(i)
              ii = ii + 1
          end if
       end do
   ! If no channels are to be ignored, use all available channels
   else
       allocate(modBandcorr(nc_alb))
       modBandcorr = modBand
       allocate(rho_0dcorr(nc_alb))
       rho_0dcorr = rho_0d
       allocate(rho_ddcorr(nc_alb))
       rho_ddcorr = rho_dd
   end if
   
   ! Neighboring MODIS bands to interpolate to BUGSrad
   modID1 = 0
   modID2 = 0
   
   rho0d_bugsrad(:) = 0.
   rhodd_bugsrad(:) = 0.
   
   do j=1,6
      bb = float(int(bugsBand(j)*100))/100
      ! Nearest modis band < bugsrad band
      tmploc = minloc(abs(modBandcorr-bb), mask=(modBandcorr-bb).lt.0)
      if (tmploc(1).eq.0) then
          tmploc(1) = 1 ! no modis band less than bugsrad band
      elseif (tmploc(1).eq.size(modBandcorr)) then
          tmploc(1) = size(modBandcorr)-1 ! no modis band larger than bugsrad band
      endif
      modID1(j) = tmploc(1)
      ! nearest modis band => bugsrad band
      tmploc = minloc(abs(modBandcorr-bb), mask=(modBandcorr-bb).ge.0)
      if (tmploc(1).eq.0) then
          tmploc(1) = size(modBandcorr) ! no modis band larger than bugsrad band
      elseif (tmploc(1).eq.1) then
          tmploc(1) = 2 ! no modis band less than bugsrad band
      endif
      modID2(j) = tmploc(1)
      ! Interpolate to bugsrad bands
      if (j.eq.1) then
          !find channels =<0.67 and average them
          masked_rho_0d = merge(rho_0dcorr, 0., modBandcorr.le.0.67)
          masked_rho_dd = merge(rho_ddcorr, 0., modBandcorr.le.0.67)
          
          rho0d_bugsrad(1) = sum(masked_rho_0d, mask=(masked_rho_0d.ne.0)) / (max(1,count(masked_rho_0d/=0)))
          rhodd_bugsrad(1) = sum(masked_rho_dd, mask=(masked_rho_dd.ne.0)) / (max(1,count(masked_rho_dd/=0)))
          
      else
          ! rho_0d
          m = (rho_0dcorr(modID2(j)) - rho_0dcorr(modID1(j))) / &
             (modBandcorr(modID2(j)) - modBandcorr(modID1(j)))
          b = rho_0dcorr(modID1(j)) - m*modBandcorr(modID1(j))
          rho0d_bugsrad(j) = m*bugsBand(j) + b
          if (rho0d_bugsrad(j) .ge. 1.) rho0d_bugsrad(j) = 1.0

          ! rho_dd
          m = (rho_ddcorr(modID2(j)) - rho_ddcorr(modID1(j))) / &
             (modBandcorr(modID2(j)) - modBandcorr(modID1(j)))
          b = rho_ddcorr(modID1(j)) - m*modBandcorr(modID1(j))
          rhodd_bugsrad(j) = m*bugsBand(j) + b
          if (rhodd_bugsrad(j) .ge. 1.) rhodd_bugsrad(j) = 1.0
      end if

      
      
   end do

end subroutine preprocess_bugsrad_sfc_albedo
