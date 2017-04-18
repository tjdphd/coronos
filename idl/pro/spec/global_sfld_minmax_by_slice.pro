FUNCTION global_sfld_minmax_by_slice, i_sfld, dsc_lab, layer, first_step, last_step

  total_steps = last_step - first_step + 1

  sfld_minmax_by_step = FLTARR(2, total_steps)

  FOR I = first_step, last_step DO BEGIN

    E             = open_spec_data_file( dsc_lab, layer, I)

    min_E_by_step = MIN(E[*,i_sfld])
    max_E_by_step = MAX(E[*,i_sfld])

    sfld_minmax_by_step[0,I - first_step] = min_E_by_step 
    sfld_minmax_by_step[1,I - first_step] = max_E_by_step 

  ENDFOR

  nz_mins      = WHERE(sfld_minmax_by_step[0,*] NE 0)
  glb_sfld_min = MIN(sfld_minmax_by_step[0,nz_mins])
  glb_sfld_max = MAX(sfld_minmax_by_step[1,*])

  glb_minmax = [glb_sfld_min, glb_sfld_max]

  RETURN, glb_minmax

END
