PRO single_time_contour_tomography_series , first_slice, last_slice,                $
                                            preval   , INC_PREFIX    = inc_prefix,  $
                                            str_res  , INC_RES_STR   = inc_str_res, $
                                            lab_run  , INC_RUN_LAB   = inc_lab_run, $
                                            out_dir  , INC_OUT_DIR   = inc_out_dir, $ 
                                            n_step   , INC_N_SLC     = inc_n_step,  $
                                            qty      , INC_QTY       = inc_qty,     $
                                            n_cntrs  , INC_N_CNTRS   = inc_n_cntrs, $
                                            tot_slc  , INC_TOT_SLCS  = inc_tot_slcs

COMMON save_r, Pre
COMMON First, first_call
COMMON fnc,   prefix, inc_res, nfld
COMMON loc,   eps_out_dir
COMMON CC,    Cntrs
COMMON JC,    J_Cntrs

IF (KEYWORD_SET(inc_qty)) THEN BEGIN
dummy = 1.0
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
 IF (KEYWORD_SET(inc_tot_slcs)) THEN BEGIN
           tot_slcs   = 200
 ENDIF ELSE BEGIN
           
       ENDELSE

Pre                   = FLTARR(1,7)
Pre                   = TRANSPOSE(Pre)
first_call            = 0

out_dev               = 'X'

tot_steps             = 200

SET_PLOT, out_dev

print, 'tomo: first_slice = ', first_slice
print, 'tomo: last_slice  = ', last_slice
print, 'tomo: preval      = ', preval
print, 'tomo: str_res     = ', str_res
print, 'tomo: lab_run     = ', lab_run
print, 'tomo: out_dir     = ', out_dir
print, 'tomo: n_step      = ', n_step
print, 'tomo: n_cntrs     = ', n_cntrs
print, 'tomo: tot_slc     = ', tot_slc

FOR I = first_slice, last_slice DO BEGIN

  print, 'single_time_contour_tomography_series: I = ', I

    IF (I EQ first_slice) THEN BEGIN
       IF (qty NE 'j') THEN BEGIN
          Cntrs   = set_levels_by_slice( n_step, qty, 1, tot_slc, n_contours, desc_label )
       ENDIF ELSE BEGIN
          Cntrs   = set_levels_by_slice( n_step, 'a', 1, tot_slc, n_contours, desc_label )
          J_Cntrs = set_levels_by_slice( n_step, qty, 1, tot_slc, n_contours, desc_label )
             ENDELSE
    ENDIF
 
    Pic                = construct_contour(I, qty, n_step, 999, desc_label)

;  IF (Pic EQ -1) THEN PRINT, FORMAT = '(a48,1x,i4,1x,a12)',                               $
;                                        'single_slice_contour_time_series: WARNING - step', $
;                                         I, 'was skipped.'

ENDFOR

END
