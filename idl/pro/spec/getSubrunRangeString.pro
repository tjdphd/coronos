FUNCTION getSubrunRangeString, first_step, last_step

  str_first_step   = STRTRIM(first_step,2)
  WHILE (STRLEN(str_first_step) LT 3) DO BEGIN
    str_first_step = '0' + str_first_step
  ENDWHILE
  str_last_step    = STRTRIM(last_step, 2)
  WHILE (STRLEN(str_last_step) LT 3) DO BEGIN
    str_last_step = '0' + str_last_step
  ENDWHILE

  IF (first_step EQ last_step) THEN BEGIN
    str_srs_ave    = 'srun_' + str_first_step
  ENDIF ELSE BEGIN
    str_srs_ave    = 'srs_'  + str_first_step + '-' + str_last_step
  ENDELSE

  RETURN, str_srs_ave

END
