FUNCTION global_sfld_minmax_over_slices, i_sfld, dsc_lab, first_layer, last_layer, step

  EULER                         = 2.7182818284590452D

  i_k                           =  0 ; k-value
  i_pe                          =  1 ; kinetic        energy spectrum in layer ? of time step ?
  i_ae                          =  2 ; magnetic       energy spectrum in layer ? of time step ?
  i_ts                          =  3 ; total          energy spectrum in layer ? of time step ?
  i_zp                          =  4 ; elsasser+      energy spectrum in layer ? of time step ?
  i_zm                          =  5 ; elsasser-      energy spectrum in layer ? of time step ?
  i_tz                          =  6 ; total elsasser energy spectrum in layer ? of time step ?
  i_pef                         =  7 ;
  i_aef                         =  8 ;
  i_zpf                         =  9 ;
  i_zmf                         = 10 ;

  total_layers = last_layer - first_layer + 1

  sfld_minmax_over_slices = FLTARR(2, total_layers)

  FOR I = first_layer, last_layer DO BEGIN

    E             = open_spec_data_file( dsc_lab, I, step)

    IF (i_sfld LE i_tz) THEN BEGIN

      min_E_by_slice = MIN(E[*,i_sfld])
      max_E_by_slice = MAX(E[*,i_sfld])

    ENDIF ELSE BEGIN

      min_E_by_slice = MIN(EULER^E[*,i_sfld])
      max_E_by_slice = MAX(EULER^E[*,i_sfld])

    ENDELSE

    sfld_minmax_over_slices[0,I - first_layer] = min_E_by_slice 
    sfld_minmax_over_slices[1,I - first_layer] = max_E_by_slice 

  ENDFOR

  nz_mins      = WHERE(sfld_minmax_over_slices[0,*] NE 0)
  glb_sfld_min = MIN(sfld_minmax_over_slices[0,nz_mins])
  glb_sfld_max = MAX(sfld_minmax_over_slices[1,*])

  glb_minmax   = [glb_sfld_min, glb_sfld_max]

  RETURN, glb_minmax

END
