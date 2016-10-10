PRO cts,  desc_label,       $
          qty,              $
          n_slice,          $
          first_step,       $ 
          last_step,        $
          global_in_minmax, $
          PREFIX  = prefix, $
          N_CNTRS = n_cntrs

  total_steps        = last_step - first_step + 1
  cur_dir            = GETENV('PWD')
  out_dir            = cur_dir + '/cts/eps' 
      
  IF (NOT KEYWORD_SET(PREFIX) ) THEN prefix   = 'rmct2'
  IF (NOT KEYWORD_SET(N_CNTRS)) THEN n_cntrs  = 31
 
  ip1                = scan_parameters('ip1', 0, desc_label)
  ip2                = scan_parameters('ip2', 0, desc_label)
  n1                 = 2^ip1
  n2                 = 2^ip2
  KK                 = calcKK(n1,n2)

  IF (STRLEN(global_in_minmax) > 0 ) THEN BEGIN

    in_file          = cur_dir + '/glb_ext.out'

    OPENR, in_unit, in_file, /GET_LUN, ERROR = op_err

    IF ( op_err EQ 0 ) THEN BEGIN
    
      PRINT, 'glb_ext: reading previous output from ', in_file

      global_qty_minmax        = FLTARR(4,4)
      line                     = FLTARR(4)

      READF, in_unit,  FORMAT = '(4(e24.16,1x),:)',   line
      global_qty_minmax[*, 0]                    =    line
      READF, in_unit,  FORMAT = '(4(e24.16,1x),:)',   line
      global_qty_minmax[*, 1]                    =    line
      READF, in_unit,  FORMAT = '(4(e24.16,1x),:)',   line
      global_qty_minmax[*, 2]                     =   line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',    line
      global_qty_minmax[*, 3]                    =    line
    
      CLOSE, in_unit

    ENDIF ELSE BEGIN

;     PRINT, 'glb_ext: WARNING - could not open file ', in_file, '. Assuming no global minmax given...' & STOP
      PRINT, 'glb_ext: WARNING - could not open file ', in_file, '. Assuming no global minmax given...'

      ip1               = scan_parameters('ip1', 0, desc_label)
      ip2               = scan_parameters('ip2', 0, desc_label)

      n1                = 2^ip1
      n2                = 2^ip2

      global_qty_minmax = global_q_minmax('x',                                $
                                          first_step,                         $
                                          last_step,                          $
                                          n_slice,                            $
                                          n_slice,                            $
                                          desc_label,                         $
                                          KK,                                 $
                                          GLOBAL_IN_MINMAX = global_in_minmax $
                                         )

    ENDELSE

  ENDIF ELSE BEGIN

    PRINT, 'no global minmax given...'

    global_qty_minmax   = global_q_minmax('x',                                $
                                          first_step,                         $
                                          last_step,                          $
                                          n_slice,                            $
                                          n_slice,                            $
                                          desc_label,                         $
                                          KK,                                 $
                                          GLOBAL_IN_MINMAX = global_in_minmax $
                                         )

  ENDELSE

  FOR I = first_step, last_step DO BEGIN

    IF (I EQ first_step) THEN BEGIN
      IF (qty NE 'j') THEN BEGIN
        Q_CNTRS         = set_contour_levels(desc_label,qty,n_slice,n_cntrs,first_step,last_step,global_qty_minmax)
      ENDIF ELSE BEGIN
        Q_CNTRS         = set_contour_levels(desc_label,qty,n_slice,n_cntrs,first_step,last_step,global_qty_minmax)
        A_CNTRS         = set_contour_levels(desc_label,'a',n_slice,n_cntrs,first_step,last_step,global_qty_minmax)
      ENDELSE
    ENDIF

    dummy = plot_surface(n_slice, qty, I, 999,desc_label, global_qty_minmax, KK)
;   Pic                 = construct_contour( n_slice,           $
;                                            qty,               $
;                                            I,                 $
;                                            999,               $
;                                            desc_label,        $
;                                            global_qty_minmax, $
;                                            KK,                $
;                                            PREFIX=prefix,     $
;                                            QCNTRS=Q_CNTRS,    $
;                                            ACNTRS=A_CNTRS     $
;                                          )

;  IF (Pic EQ -1) THEN PRINT, FORMAT = '(a48,1x,i4,1x,a12)',                                 $
;                                        'single_slice_contour_time_series: WARNING - step', $
;                                         I, 'was skipped.'

  ENDFOR

END
