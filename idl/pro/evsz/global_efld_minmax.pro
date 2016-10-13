FUNCTION global_efld_minmax, i_efld, dsc_lab, first_step, last_step

  total_steps = last_step - first_step + 1

  efld_minmax_by_step = FLTARR(2, total_steps)

  FOR I = first_step, last_step DO BEGIN

    E             = open_evsz_data_file( dsc_lab, I)
    extE          = calc_extE( E, dsc_lab , I, i_efld)

    min_E_by_step = MIN(extE[*,i_efld])
    max_E_by_step = MAX(extE[*,i_efld])

    efld_minmax_by_step[0,I - first_step] = min_E_by_step 
    efld_minmax_by_step[1,I - first_step] = max_E_by_step 

  ENDFOR

  glb_efld_min = MIN(efld_minmax_by_step[0,*])
  glb_efld_max = MAX(efld_minmax_by_step[1,*])

  glb_minmax = [glb_efld_min, glb_efld_max]

  RETURN, glb_minmax

END
