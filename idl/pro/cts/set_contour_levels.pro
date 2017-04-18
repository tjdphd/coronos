FUNCTION set_contour_levels,  $
         desc_label,          $
         qty,                 $
         n_slice,             $
         n_cntrs,             $
         first_step,          $
         last_step,           $
         global_in_minmax

  IF n_cntrs MOD 2 EQ 0 THEN n_cntrs = n_cntrs + 1

  IF (qty EQ 'j' ) THEN global_qty_minmax = global_in_minmax[*,3]
  IF (qty EQ 'o' ) THEN global_qty_minmax = global_in_minmax[*,2]
  IF (qty EQ 'a' ) THEN global_qty_minmax = global_in_minmax[*,1]
  IF (qty EQ 'p' ) THEN global_qty_minmax = global_in_minmax[*,0]

  n_time_steps = last_step - first_step + 1

  n_cntrs      = n_cntrs - 2
  CNTRS        = FLTARR(n_cntrs + 2)

  q_max        = global_qty_minmax[0] 
  q_min        = global_qty_minmax[1]

  PRINT, 'set_contour_levels: q_max = ', q_max
  PRINT, 'set_contour_levels: q_min = ', q_min

  IF  ABS(q_min) GT q_max THEN q_max =  ABS(q_min)
  IF -ABS(q_max) LT q_min THEN q_min = -ABS(q_max)
 
  inc_q                 = q_max / ((n_cntrs - 1)/2)
  inc_q                 = (q_max - q_min) / (n_cntrs - 1)

  CNTRS[1]              = q_min 
  CNTRS[n_cntrs]        = q_max

  q_max                 = q_max + (inc_q/3.0)
  q_min                 = q_min - (inc_q/3.0)

  CNTRS[0]              = q_min 
  CNTRS[n_cntrs+1]      = q_max
 
; PRINT, "CNTRS[", 0, "]  = ", CNTRS[0]

  FOR I = 2, n_cntrs - 1 DO BEGIN

    CNTRS[I]              = CNTRS[I-1] + inc_q
;   PRINT, "CNTRS[", I, "] = ", CNTRS[I]
  ENDFOR
; PRINT, "CNTRS[", n_cntrs+1, "] = ", CNTRS[n_cntrs+1]

  n_cntrs = n_cntrs + 2
  I_mid = (n_cntrs - 1)/2

  q_max        = global_qty_minmax[0] 
  IF ABS(CNTRS[I_mid]/q_max) LT 1.0e-06 THEN CNTRS[I_mid] = 0.0
; PRINT, "CNTRS[", I_mid, "] = ", CNTRS[I_mid]

  RETURN, CNTRS

END
