FUNCTION set_jcontour_levels, slice_n, qty, first_step, last_step, n_contours, desc_label

COMMON CC, Cntrs
COMMON fnc, prefix, inc_res, nfld
COMMON JC, J_Cntrs

IF n_contours MOD 2 EQ 0 THEN n_contours = n_contours + 1

n_time_steps = last_step - first_step + 1


Res_Cntrs    = FLTARR(1, n_contours) ; Resulting contour array
Res_Cntrs    = TRANSPOSE(Res_Cntrs)

q_max        = max_q(slice_n, qty, first_step, last_step, desc_label)
q_min        = min_q(slice_n, qty, first_step, last_step, desc_label)

IF  ABS(q_min) GT q_max THEN q_max =  ABS(q_min)
IF -ABS(q_max) LT q_min THEN q_min = -ABS(q_max)

  Res_Cntrs[0]              = q_min
  Res_Cntrs[n_contours - 1] = q_max

  inc_q                     = q_max / ((n_contours - 1)/2)
 
FOR I = 1, (n_contours - 1) DO BEGIN
  Res_Cntrs[I]              = Res_Cntrs[I-1] + inc_q
ENDFOR

I_mid = (n_contours - 1)/2

IF ABS(Res_Cntrs[I_mid]/q_max) LT 1.0e-08 THEN Res_Cntrs[I_mid] = 0.0

IF (qty NE 'j') THEN BEGIN
  Cntrs                     = Res_Cntrs
ENDIF ELSE BEGIN
  J_Cntrs                   = Res_Cntrs
      ENDELSE

RETURN, Res_Cntrs
END
