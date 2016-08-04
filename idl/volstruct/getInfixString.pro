FUNCTION getInfixString, desc_label,            $
                         type,                  $
                         last_step,             $ 
                         THRESHOLD = threshold, $
                         WIDTH     = width

 str_last_stp             = STRTRIM(last_step,2)               ; prepare an infix string for file-naming
 len_steps                = STRLEN(str_last_stp)

 IF (len_steps LT 3) THEN len_steps = 3                        ; want uniform number fields on file names

  zl                      = scan_parameters('zl',0, desc_label)

  str_zl                  = STRTRIM(UINT(zl),   2)

  IF (KEYWORD_SET(THRESHOLD)) THEN BEGIN
    str_thr               = STRTRIM(threshold,  2)
    str_thr               = STRMID(str_thr, 0,  6)
  ENDIF ELSE BEGIN
    IF (KEYWORD_SET(WIDTH)) THEN BEGIN
      str_thr             = STRTRIM(width,      2)
      str_thr             = STRMID( str_thr, 0, 6)
    ENDIF ELSE BEGIN
      PRINT, 'Why the hell not?' & STOP
    ENDELSE
  ENDELSE

  dec_loc                 = STRPOS(str_thr, ".")

  str_thr_p               = str_thr
  STRPUT, str_thr_p, 'p', dec_loc

  len_thr                 = STRLEN(str_thr_p)
  count                   = 0
  WHILE (len_thr LT 6) DO BEGIN
    count                 = count + 1
    str_thr_p             = str_thr_p + '0'
    len_thr               = len_thr + 1
  ENDWHILE

  res_str                 = getResString(desc_label)
  IF (type EQ "vol" ) THEN BEGIN
  str_ifx                 = type + '_' + res_str + '_zl_' + str_zl +  '_thr_' + str_thr_p + '-'
  ENDIF ELSE BEGIN
    IF (type EQ "isf" ) THEN BEGIN
      str_ifx             = type + '_' + res_str + '_zl_' + str_zl +  '_width_' + str_thr_p + '-'
    ENDIF ELSE BEGIN
      str_ifx             = '_unknown_'
    ENDELSE
  ENDELSE

  RETURN, str_ifx

END
