$(OBJS)/aatsr_corrections.o: $(OBJS)/calender.o $(OBJS)/date_type_structure.o \
        $(OBJS)/preproc_constants.o
$(OBJS)/bright_m.o: $(OBJS)/preproc_constants.o fundamental_constants.inc
$(OBJS)/build_preproc_fields.o: $(OBJS)/imager_structures.o \
        $(OBJS)/preproc_constants.o $(OBJS)/preproc_structures.o
$(OBJS)/calender.o: $(OBJS)/date_type_structure.o $(OBJS)/preproc_constants.o
$(OBJS)/channel_structures.o: $(OBJS)/preproc_constants.o \
        allocate_channel_info.F90 deallocate_channel_info.F90
$(OBJS)/cloud_typing_pavolonis.o: $(OBJS)/constants_cloud_typing_pavolonis.o \
        $(OBJS)/ecmwf.o $(OBJS)/imager_structures.o $(OBJS)/neural_net_preproc.o \
        $(OBJS)/surface_structures.o pavolonis_fill_coefficients.inc
$(OBJS)/correct_for_ice_snow.o: $(OBJS)/channel_structures.o \
        $(OBJS)/constants_cloud_typing_pavolonis.o $(OBJS)/imager_structures.o \
        $(OBJS)/nise_def.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o $(OBJS)/surface_structures.o
$(OBJS)/cox_munk.o: $(OBJS)/nr.o $(OBJS)/preproc_constants.o
$(OBJS)/date_type_structure.o: $(OBJS)/preproc_constants.o
$(OBJS)/define_preprop_grid.o: $(OBJS)/imager_structures.o \
        $(OBJS)/preproc_constants.o $(OBJS)/preproc_structures.o
$(OBJS)/ecmwf.o: $(OBJS)/preproc_constants.o $(OBJS)/preproc_structures.o \
        compute_geopot_coordinate.F90 deallocate_ecmwf_structures.F90 \
        read_ecmwf_grib.F90 read_ecmwf_nc.F90 read_ecmwf_wind_badc.F90 \
        read_ecmwf_wind_grib.F90 read_ecmwf_wind_nc.F90 rearrange_ecmwf.F90
$(OBJS)/fill_grid.o: $(OBJS)/preproc_constants.o
$(OBJS)/get_USGS_data.o: $(OBJS)/imager_structures.o $(OBJS)/read_USGS_file.o
$(OBJS)/imager_structures.o: $(OBJS)/channel_structures.o \
        $(OBJS)/preproc_constants.o allocate_imager_structures.F90 \
        deallocate_imager_structures.F90
$(OBJS)/map_time_to_pixel.o: $(OBJS)/imager_structures.o \
        $(OBJS)/preproc_constants.o
$(OBJS)/mcd43c.o: $(OBJS)/preproc_constants.o read_mcd43c1.F90 read_mcd43c3.F90
$(OBJS)/netcdf_output.o: $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o $(OBJS)/surface_structures.o \
        netcdf_output_close.F90 netcdf_output_create.F90 \
        netcdf_output_create_file.F90 netcdf_output_write_swath.F90
$(OBJS)/neural_net_preproc.o: $(OBJS)/constants_cloud_typing_pavolonis.o \
        $(OBJS)/neural_net_constants.o
$(OBJS)/nise_def.o: $(OBJS)/preproc_constants.o read_nsidc_nise.F90
$(OBJS)/preparation.o: $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o
$(OBJS)/preproc_structures.o: $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        allocate_preproc_structures.F90 deallocate_preproc_structures.F90
$(OBJS)/read_aatsr.o: $(OBJS)/aatsr_corrections.o $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        read_aatsr_dimensions.F90 read_aatsr_l1b.F90
$(OBJS)/read_avhrr.o: $(OBJS)/calender.o $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        create_time_for_pixel.F90 read_avhrr_angles.F90 \
        read_avhrr_dimensions.F90 read_avhrr_l1b_radiances.F90 \
        read_avhrr_l1b_radiances_2.F90 read_avhrr_land_sea_mask.F90 \
        read_avhrr_land_sea_mask_2.F90 read_avhrr_lat_lon.F90 \
        read_avhrr_time.F90 read_avhrr_time_lat_lon_angles.F90
$(OBJS)/read_cimss_emissivity.o: $(OBJS)/preproc_constants.o
$(OBJS)/read_imager.o: $(OBJS)/channel_structures.o $(OBJS)/imager_structures.o \
        $(OBJS)/preproc_constants.o $(OBJS)/read_aatsr.o $(OBJS)/read_avhrr.o \
        $(OBJS)/read_modis.o
$(OBJS)/read_modis.o: $(OBJS)/bright_m.o $(OBJS)/calender.o \
        $(OBJS)/channel_structures.o $(OBJS)/imager_structures.o \
        $(OBJS)/preproc_constants.o get_modis_time.F90 map_time_to_pixel.F90 \
        modis_bright.F90 read_modis_angles.F90 read_modis_dimensions.F90 \
        read_modis_l1b_radiances.F90 read_modis_l1b_radiances_2.F90 \
        read_modis_land_sea_mask.F90 read_modis_lat_lon.F90 read_modis_time.F90 \
        read_modis_time_lat_lon_angles.F90
$(OBJS)/read_USGS_file.o: $(OBJS)/preproc_constants.o
$(OBJS)/ross_thick_li_sparse_r.o: $(OBJS)/nr.o $(OBJS)/preproc_constants.o
$(OBJS)/rttov_driver.o: $(OBJS)/channel_structures.o $(OBJS)/netcdf_output.o \
        $(OBJS)/preproc_constants.o $(OBJS)/preproc_structures.o \
        write_ir_rttov.F90 write_solar_rttov.F90
$(OBJS)/set_ecmwf.o: $(OBJS)/preproc_constants.o
$(OBJS)/setup.o: $(OBJS)/calender.o $(OBJS)/channel_structures.o \
        $(OBJS)/date_type_structure.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o
$(OBJS)/surface_emissivity.o: $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o $(OBJS)/read_cimss_emissivity.o \
        $(OBJS)/surface_structures.o get_surface_emissivity.F90 \
        select_modis_emiss_file.F90
$(OBJS)/surface_reflectance.o: $(OBJS)/channel_structures.o $(OBJS)/cox_munk.o \
        $(OBJS)/ecmwf.o $(OBJS)/fill_grid.o $(OBJS)/imager_structures.o \
        $(OBJS)/mcd43c.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o $(OBJS)/ross_thick_li_sparse_r.o \
        $(OBJS)/surface_structures.o get_surface_reflectance.F90 \
        select_modis_albedo_file.F90
$(OBJS)/surface_structures.o: $(OBJS)/channel_structures.o \
        $(OBJS)/imager_structures.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o allocate_surface_structures.F90 \
        deallocate_surface_structures.F90
$(OBJS)/test_rttov.o: $(OBJS)/preproc_constants.o $(OBJS)/preproc_structures.o
$(OBJS)/preprocessing_for_orac.o: $(OBJS)/channel_structures.o \
        $(OBJS)/cloud_typing_pavolonis.o $(OBJS)/correct_for_ice_snow.o \
        $(OBJS)/ecmwf.o $(OBJS)/imager_structures.o $(OBJS)/netcdf_output.o \
        $(OBJS)/preparation.o $(OBJS)/preproc_constants.o \
        $(OBJS)/preproc_structures.o $(OBJS)/read_USGS_file.o \
        $(OBJS)/read_aatsr.o $(OBJS)/read_avhrr.o $(OBJS)/read_imager.o \
        $(OBJS)/read_modis.o $(OBJS)/rttov_driver.o $(OBJS)/setup.o \
        $(OBJS)/surface_emissivity.o $(OBJS)/surface_reflectance.o \
        $(OBJS)/surface_structures.o