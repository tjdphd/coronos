FUNCTION getCaseFileString, desc_label, sfld, layer, first_step, last_step


  res_str          = getResString(desc_label)
  str_layer        = STRTRIM(layer,2)

  str_first_step   = STRTRIM(first_step,2)
  str_last_step    = STRTRIM(last_step, 2)

  WHILE (STRLEN(str_layer) LT 3) DO BEGIN
    str_layer      = '0' + str_layer
  ENDWHILE
  WHILE (STRLEN(str_first_step) LT 3) DO BEGIN
    str_first_step = '0' + str_first_step
  ENDWHILE
  WHILE (STRLEN(str_last_step) LT 3) DO BEGIN
    str_last_step = '0' + str_last_step
  ENDWHILE

  n3               = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  zl               = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume

  str_zl           = STRTRIM(UINT(zl),2)

  IF (layer GT n3) THEN BEGIN
    is_layer_even  = (layer MOD 2) EQ 0
    IF (is_layer_even) THEN BEGIN
      n_proc       = CEIL(UINT(layer/n3)-1)
    ENDIF ELSE BEGIN
      n_proc       = FLOOR(UINT(layer/n3))
   ENDELSE
  ENDIF ELSE BEGIN
     n_proc        = 0
  ENDELSE

  str_proc         = STRTRIM(n_proc,   2)
  len_str_proc     = STRLEN(str_proc)

  WHILE (STRLEN(str_proc) LT 3) DO BEGIN
    str_proc       = '0' + str_proc
  ENDWHILE

  IF (first_step EQ last_step) THEN BEGIN
    str_srs_ave    = 'srun_' + str_first_step
  ENDIF ELSE BEGIN
    str_srs_ave    = 'srs_'  + str_first_step + '-' + str_last_step
  ENDELSE
 
  str_sfld         = getSfldString(sfld)
  case_file_str    =  str_sfld   + "ra_spec_" + res_str   + '_zl_'     + str_zl                       $
                   +  '_proc-'   + str_proc   + '_layer-' + str_layer                                 $
                   +  '_'        + str_srs_ave
  Return, case_file_str

END
