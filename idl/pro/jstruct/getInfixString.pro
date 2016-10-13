FUNCTION getInfixString, desc_label, action, last_step, threshold

 str_last_stp             = STRTRIM(last_step,2)               ; prepare an infix string for file-naming
 len_steps                = STRLEN(str_last_stp)

 IF (len_steps LT 3) THEN len_steps = 3                        ; want uniform number fields on file names

  zl                      = scan_parameters('zl',0, desc_label)

  str_zl                  = STRTRIM(UINT(zl),   2)

  str_thr                 = STRTRIM(threshold,  2)
  str_thr                 = STRMID(str_thr, 0,  6)

  dec_loc                 = STRPOS(str_thr, ".")
  str_thr_p               = str_thr
  STRPUT, str_thr_p, 'p', dec_loc

  res_str                 = getResString(desc_label)
  IF (action EQ "vol" ) THEN BEGIN
  str_ifx                 = action + '_' + res_str + '_zl_' + str_zl +  '_thr_' + str_thr_p + '-'
  ENDIF ELSE BEGIN
    IF (action EQ "isf" ) THEN BEGIN
      str_ifx             = action + '_' + res_str + '_zl_' + str_zl +  '_width_' + str_thr_p + '-'
    ENDIF ELSE BEGIN
      str_ifx             = '_unknown_'
    ENDELSE
  ENDELSE


  RETURN, str_ifx

END
