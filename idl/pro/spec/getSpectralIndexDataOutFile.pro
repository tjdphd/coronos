FUNCTION getSpectralIndexDataOutFile, plot_mode, first_step, last_step

  str_first_step                = STRTRIM(first_step,2)
  str_last_step                 = STRTRIM(last_step, 2)

  zl                            = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
  str_zl                        = STRTRIM(UINT(zl),2)

  WHILE (STRLEN(str_first_step) LT 3) DO BEGIN
    str_first_step              = '0' + str_first_step
  ENDWHILE

  WHILE (STRLEN(str_last_step) LT 3) DO BEGIN
    str_last_step               = '0' + str_last_step
  ENDWHILE

  IF (first_step EQ last_step) THEN BEGIN
    str_srs_ave                 = 'srun_' + str_first_step
  ENDIF ELSE BEGIN
    str_srs_ave                 = 'srs_'  + str_first_step + '-' + str_last_step
  ENDELSE

  data_out_dir                  = getDataOutDir(plot_mode)
  spec_idx_out_file             = data_out_dir + 'sp_idx_' + str_srs_ave + '_L-' + str_zl + '.dat'

  RETURN, spec_idx_out_file

END
