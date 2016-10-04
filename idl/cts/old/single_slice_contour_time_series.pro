PRO single_slice_contour_time_series , first_step, last_step,                  $
                                       preval   , INC_PREFIX    = inc_prefix,  $
                                       str_res  , INC_RES_STR   = inc_str_res, $
                                       lab_run  , INC_RUN_LAB   = inc_lab_run, $
                                       out_dir  , INC_OUT_DIR   = inc_out_dir, $ 
                                       n_slc    , INC_N_SLC     = inc_n_slc,   $
                                       qty      , INC_QTY       = inc_qty,     $
                                       n_cntrs  , INC_N_CNTRS   = inc_n_cntrs, $
                                       tot_steps, INC_TOT_STEPS = inc_tot_steps

COMMON save_r, Pre
COMMON First, first_call
COMMON fnc,   prefix, inc_res, nfld
COMMON loc,   eps_out_dir
COMMON CC,    Cntrs
COMMON JC,    J_Cntrs

IF (KEYWORD_SET(inc_qty)) THEN BEGIN

ENDIF ELSE BEGIN
          qty         = 'p'
      ENDELSE
      
IF (KEYWORD_SET(inc_out_dir)) THEN BEGIN
          eps_out_dir = out_dir
ENDIF ELSE BEGIN
          eps_out_dir = GETENV('PWD')
      ENDELSE
      
PRINT, 'postscript output will be directed to: ', eps_out_dir

IF (KEYWORD_SET(inc_prefix)) THEN BEGIN
          prefix      = preval
ENDIF ELSE BEGIN
          prefix      = 'rmct2'
      ENDELSE
 IF (KEYWORD_SET(inc_str_res)) THEN BEGIN
            inc_res   = str_res
 ENDIF ELSE BEGIN
            inc_res   = 'y'
       ENDELSE
 IF (KEYWORD_SET(inc_lab_run)) THEN BEGIN
           desc_label = lab_run
 ENDIF ELSE BEGIN
           desc_label = 'wn'
       ENDELSE
 IF (KEYWORD_SET(inc_n_slc)) THEN BEGIN
           n_slice    = n_slc
 ENDIF ELSE BEGIN
           n_slice    = 16
       ENDELSE
 IF (KEYWORD_SET(inc_n_cntrs)) THEN BEGIN
           n_contours = n_cntrs
 ENDIF ELSE BEGIN
           n_contours = 31
       ENDELSE
IF (KEYWORD_SET(inc_tot_steps)) THEN BEGIN

ENDIF ELSE BEGIN
          tot_steps   = 200
      ENDELSE

Pre                   = FLTARR(1,7)
Pre                   = TRANSPOSE(Pre)
first_call            = 0

out_dev               = 'X'

SET_PLOT, out_dev

FOR I = first_step, last_step DO BEGIN

   IF (I EQ first_step) THEN BEGIN
      IF (qty NE 'j') THEN BEGIN
         Cntrs        = set_jcontour_levels( n_slice, qty, 0, tot_steps, n_contours, desc_label )
      ENDIF ELSE BEGIN
              Cntrs   = set_jcontour_levels( n_slice, 'a', 0, tot_steps, n_contours, desc_label )
              J_Cntrs = set_jcontour_levels( n_slice, qty, 0, tot_steps, n_contours, desc_label )
            ENDELSE
   ENDIF

   Pic                = construct_jcontour(n_slice, qty, I, 999, desc_label)

;  IF (Pic EQ -1) THEN PRINT, FORMAT = '(a48,1x,i4,1x,a12)',                               $
;                                        'single_slice_contour_time_series: WARNING - step', $
;                                         I, 'was skipped.'

ENDFOR

END
