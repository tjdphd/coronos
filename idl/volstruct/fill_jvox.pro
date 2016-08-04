FUNCTION fill_jvox, desc_label, n_step, threshold, first_slice, last_slice, global_j_minmax, pls_min, JVOX_P, JVOX_M, KK

  n3                                                 = scan_parameters('n3',  0, desc_label )

  global_j_min                                       = global_j_minmax[0]
  global_j_max                                       = global_j_minmax[1]

  PRINT, 'fill_jvox: global_j_max = ', global_j_max
  PRINT, 'fill_jvox: global_j_min = ', global_j_min

  IF (pls_min EQ 'p' OR pls_min EQ 'pm' ) THEN size_jvox = SIZE(JVOX_P, /DIMENSIONS)
  IF (pls_min EQ 'm' OR pls_min EQ 'pm' ) THEN size_jvox = SIZE(JVOX_M, /DIMENSIONS)

  sz1                                                = size_jvox[0] ; Z
  sz2                                                = size_jvox[1] ; X
  sz3                                                = size_jvox[2] ; Y

  JV_SQR                                             = FLTARR(     sz2, sz3)

  IF (pls_min EQ 'p' OR pls_min EQ 'pm' ) THEN JV_FLT_P = FLTARR(sz1, sz2, sz3)
  IF (pls_min EQ 'm' OR pls_min EQ 'pm' ) THEN JV_FLT_M = FLTARR(sz1, sz2, sz3)

; multi-thread from here can I make this a function?

  slice_range                                        = getSliceRange(first_slice, last_slice  )
  size_sr                                            = SIZE(slice_range, /DIMENSIONS)
  n_threads                                          = size_sr[0]

  IF (n_threads GT 1) THEN oBridge                   = OBJARR(n_threads - 1)

  FOR J = 0, n_threads - 1 DO BEGIN

    thread_sz1                                       = slice_range[J, 1] - slice_range[J, 0] + 1
    start_slice                                      = slice_range[J, 0]
    stop_slice                                       = slice_range[J, 1]

    thread_size_jvox                                 = [thread_sz1, sz2, sz3]

    IF (J EQ n_threads - 1 ) THEN BEGIN

      thread_JV_FLT_PM                               = $
      calc_flt_Jsqr(desc_label, n_step, threshold, start_slice, stop_slice, global_j_minmax, pls_min, thread_size_jvox, KK)

      IF (pls_min EQ 'p' OR pls_min EQ 'pm') THEN local_JV_FLT_P = FLTARR(thread_size_jvox)
      IF (pls_min EQ 'm' OR pls_min EQ 'pm') THEN local_JV_FLT_M = FLTARR(thread_size_jvox)

      IF (pls_min EQ 'pm' ) THEN BEGIN
        local_JV_FLT_P[0:(stop_slice-start_slice),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
        local_JV_FLT_M[0:(stop_slice-start_slice),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,1]
      ENDIF ELSE BEGIN
        IF (pls_min EQ 'p') THEN $
        local_JV_FLT_P[0:(stop_slice-start_slice),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
        IF (pls_min EQ 'm') THEN $
        local_JV_FLT_M[0:(stop_slice-start_slice),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
      ENDELSE

    ENDIF ELSE BEGIN

      oBridge[J]                                     = OBJ_NEW('IDL_IDLBRIDGE')
      oBridge[J]                                     -> SetVar, 'desc_label',         desc_label
      oBridge[J]                                     -> SetVar, 'n_step',             n_step
      oBridge[J]                                     -> SetVar, 'threshold',          threshold
      oBridge[J]                                     -> SetVar, 'start_slice',        slice_range[J,0]
      oBridge[J]                                     -> SetVar, 'stop_slice',         slice_range[J,1]
      oBridge[J]                                     -> SetVar, 'global_j_minmax',    global_j_minmax
      oBridge[J]                                     -> SetVar, 'pls_min',            pls_min
      oBridge[J]                                     -> SetVar, 'thread_size_jvox',   thread_size_jvox
      oBridge[J]                                     -> Execute, ".run calc_flt_Jsqr"
      oBridge[J]                                     -> Execute,                                                                $
        "thread_JV_FLT_PM = "                                                                                                   $
      + "calc_flt_Jsqr(desc_label, n_step, threshold,start_slice, stop_slice, global_j_minmax, pls_min, thread_size_jvox, KK)", $
      /nowait
    ENDELSE
  ENDFOR

  IF (n_threads GT 1) THEN BEGIN
    notdone                                          = 1
    WHILE notdone DO BEGIN
      done                                           = 0
      FOR J = 0, n_threads - 2 DO done               = done + oBridge[J] -> Status()
      IF (done EQ 0) THEN notdone                    = done
    ENDWHILE
  ENDIF

  FOR J = 0, n_threads - 1 DO BEGIN

    thread_sz1                                       = slice_range[J, 1] - slice_range[J, 0] + 1
    start_slice                                      = slice_range[J, 0]
    stop_slice                                       = slice_range[J, 1]

    IF (J EQ n_threads - 1) THEN BEGIN

      IF (pls_min EQ 'p' OR pls_min EQ 'pm') THEN $
      JV_FLT_P[(start_slice-1):(stop_slice-1),*,*]   = local_JV_FLT_P[0:(stop_slice-start_slice),*,*]
      IF (pls_min EQ 'm' OR pls_min EQ 'pm') THEN $
      JV_FLT_M[(start_slice-1):(stop_slice-1),*,*]   = local_JV_FLT_M[0:(stop_slice-start_slice),*,*]

    ENDIF ELSE BEGIN

        thread_JV_FLT_PM = oBridge[J] -> GetVar('thread_JV_FLT_PM')

        IF (pls_min EQ 'pm') THEN BEGIN
          JV_FLT_P[(start_slice-1):(stop_slice-1),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
          JV_FLT_M[(start_slice-1):(stop_slice-1),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,1]
        ENDIF ELSE BEGIN
          IF (pls_min EQ 'p') THEN $
          JV_FLT_P[(start_slice-1):(stop_slice-1),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
          IF (pls_min EQ 'm') THEN $
          JV_FLT_M[(start_slice-1):(stop_slice-1),*,*] = thread_JV_FLT_PM[0:(stop_slice-start_slice),*,*,0]
        ENDELSE

      obj_destroy, oBridge[J]

    ENDELSE
  ENDFOR

; multi-thread ends here

  IF (pls_min EQ 'p' OR pls_min EQ 'pm') THEN BEGIN
    max_jp_flt = MAX(JV_FLT_P)
    scale_p    = max_jp_flt / global_j_max
    JVOX_P     = scale_p * BYTSCL(JV_FLT_P[*,*,*])
    PRINT, 'fill_jvox: max_jp_flt = ', max_jp_flt
    PRINT, 'fill_jvox: scale_p    = ', scale_p
  ENDIF

  IF (pls_min EQ 'm' OR pls_min EQ 'pm') THEN BEGIN
    max_jm_flt = MAX(JV_FLT_M)
    scale_m    = max_jm_flt / ABS(global_j_min)
    JVOX_M     = scale_m * BYTSCL(JV_FLT_M[*,*,*])
    PRINT, 'fill_jvox: max_jm_flt = ', max_jm_flt
    PRINT, 'fill_jvox: scale_m    = ', scale_m
  ENDIF

RETURN, 0

END
