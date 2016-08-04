FUNCTION calc_flt_Jsqr, desc_label, n_step, threshold, start_slice, stop_slice, global_j_minmax, pls_min, size_jvox, KK

  global_j_min            = global_j_minmax[0]
  global_j_max            = global_j_minmax[1]

  sz1                     = size_jvox[0] ; Z
  sz2                     = size_jvox[1] ; X
  sz3                     = size_jvox[2] ; Y

  JVP                     = FLTARR(     sz2, sz3)
  JVM                     = FLTARR(     sz2, sz3)
  JV_SQR                  = FLTARR(     sz2, sz3)

  JV_FLT_PM               = FLTARR(sz1,sz2,sz3,2)

  FOR I_Slc = start_slice, stop_slice DO BEGIN

    Slice                 = fetch_layer(I_Slc, n_step, desc_label, KK)
    J_z                   = Slice[   *, 5]

    JV_SQR                = TRANSPOSE(REFORM(J_z, sz2, sz3))
    JVP[*,*]              = 0.0
    JVM[*,*]              = 0.0
    idx_p                 = WHERE(JV_SQR GT 0 AND JV_SQR GE threshold*global_j_max)
    idx_m                 = WHERE(JV_SQR LT 0 AND JV_SQR LE threshold*global_j_min)
    ndx                   = I_Slc - start_slice
    JVP[idx_p]            = ABS(JV_SQR[idx_p])
    JVM[idx_m]            = ABS(JV_SQR[idx_m])

    JV_FLT_P              = FLTARR(sz1, sz2, sz3)
    JV_FLT_M              = FLTARR(sz1, sz2, sz3)

    JV_FLT_P[ndx, *, *]   = JVP[*,*]
    JV_FLT_M[ndx, *, *]   = JVM[*,*]

    JV_FLT_PM[ndx, *,*,0] = JV_FLT_P[ndx,*,*]
    JV_FLT_PM[ndx, *,*,1] = JV_FLT_M[ndx,*,*]

  ENDFOR

  RETURN, JV_FLT_PM

END
