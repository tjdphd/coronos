FUNCTION max_diff, qty, I_step, first_slab, last_slab, str_res, lab_one, lab_two, tolerance

 n3                          = scan_parameters('p3', 0, lab_one )
 mp                          = scan_parameters('np', 0, lab_one )
 n_layers                    = n3 * mp

 MXDF_V                      = DBLARR(  n_layers)
 MXDF_L                      = LON64ARR(n_layers)

 MNDF_V                      = DBLARR(  n_layers)
 MNDF_L                      = LON64ARR(n_layers)

 glb_count                   = 0

 FOR J = first_slab, last_slab DO BEGIN

   SLCONE                    = fetch_layer( J, I_step, lab_one )
   SLCTWO                    = fetch_layer( J, I_step, lab_two )
   
   size_one                  = SIZE(SLCONE,/DIMENSIONS)
   size_two                  = SIZE(SLCTWO,/DIMENSIONS)

   IF (size_one[1] EQ size_two[1]) THEN BEGIN

     IF (STRCMP(qty, 'p')) THEN idx = 2
     IF (STRCMP(qty, 'a')) THEN idx = 3
     IF (STRCMP(qty, 'o')) THEN idx = 4
     IF (STRCMP(qty, 'j')) THEN idx = 5

     slc_one_mean            = MEAN(ABS(SLCONE[WHERE(ABS(SLCONE[*,idx]) GE 1.0E-14), idx]))
     slc_one_max             = MAX( ABS(SLCONE[WHERE(ABS(SLCONE[*,idx]) GE 1.0E-14), idx]))
     slc_one_min             = MIN( ABS(SLCONE[WHERE(ABS(SLCONE[*,idx]) GE 1.0E-14), idx]))

     slc_two_mean            = MEAN(ABS(SLCTWO[WHERE(ABS(SLCTWO[*,idx]) GE 1.0E-14), idx]))
     slc_two_max             = MAX( ABS(SLCTWO[WHERE(ABS(SLCTWO[*,idx]) GE 1.0E-14), idx]))
     slc_two_min             = MIN( ABS(SLCTWO[WHERE(ABS(SLCTWO[*,idx]) GE 1.0E-14), idx]))

     slc_one_mean_max_ratio  = slc_one_mean / slc_one_max
     slc_one_min_max_ratio   = slc_one_min  / slc_one_max
     slc_one_ratio_ratio     = slc_one_min_max_ratio / slc_one_mean_max_ratio

     slc_two_mean_max_ratio  = slc_two_mean / slc_two_max
     slc_two_min_max_ratio   = slc_two_min  / slc_two_max
     slc_two_ratio_ratio     = slc_two_min_max_ratio / slc_two_mean_max_ratio

     ave_mean                = 0.5*(slc_one_mean + slc_two_mean)
     ave_min                 = 0.5*(slc_one_min  + slc_two_min )

     SLCDIFF                 = ABS(SLCONE[*,*] - SLCTWO[*,*])

     nz_one_n_two            = WHERE(ABS(SLCONE[*, idx])  GE 1.0E-14 AND $
                                     ABS(SLCTWO[*, idx])  GE 1.0E-14 AND $
                                         SLCDIFF[*,idx]   GE ave_min     $
                                    )

     sz_nz_one_n_two         = SIZE(nz_one_n_two, /DIMENSIONS)

     IF (sz_nz_one_n_two[0] NE 0) THEN BEGIN

       SLCDIFF[nz_one_n_two[0],idx] = SLCDIFF[nz_one_n_two[0],idx] / ABS(SLCONE[nz_one_n_two[0],idx])
       diff_scale            = MAX(SLCDIFF[nz_one_n_two[0],idx])
       min_diff_scale        = MIN(SLCDIFF[nz_one_n_two[0],idx], SUBSCRIPT_MAX = loc_diff_scale )

     ENDIF ELSE BEGIN

       diff_scale            = MAX(SLCDIFF[*,idx])
       min_diff_scale        = MIN(SLCDIFF[*,idx], SUBSCRIPT_MAX = loc_diff_scale )

     ENDELSE

     PRINT, FORMAT = '(A20,E24.16)', "max_diff: diff_scale      = ", diff_scale

     tol_idx                 = WHERE(SLCDIFF[*,idx] GE tolerance, count_tol)

     size_tol_idx            = SIZE(tol_idx,/DIMENSIONS)

     IF (size_tol_idx[0] GT 0) THEN BEGIN
       FOR K = 0, size_tol_idx[0] - 1  DO BEGIN
         PRINT, 'tol_idx[', K, '] = ', tol_idx[K], ' SLCDIFF[',tol_idx[K],'] = ', SLCDIFF[tol_idx[K],idx]
       ENDFOR
     ENDIF

     glb_count               = glb_count + count_tol

     mxdf                    = MAX(SLCDIFF[*,idx], SUBSCRIPT_MIN = loc_mndf )
     mndf                    = MIN(SLCDIFF[*,idx], SUBSCRIPT_MAX = loc_mxdf )

   ENDIF ELSE BEGIN

     mxdf                    =  1.0e+10
     mndf                    = -1.0e+10

   ENDELSE

   MXDF_V[J-1]               = mxdf
   MXDF_L[J-1]               = loc_mxdf

   MNDF_V[J-1]               = mndf
   MNDF_L[J-1]               = loc_mndf

 ENDFOR

 DF_OUT                      = {df_out,        $
                                mxdf_v:MXDF_V, $
                                mxdf_l:MXDF_L, $
                                mndf_v:MNDF_V, $
                                mndf_l:MNDF_L, $
                                df_count:glb_count}
 RETURN, DF_OUT

END
