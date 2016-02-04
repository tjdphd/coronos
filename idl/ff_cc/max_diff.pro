FUNCTION max_diff, qty, I_step, first_slab, last_slab, str_res, lab_one, lab_two, tolerance

 n3              = scan_parameters('n3', 0, lab_one )
 mp              = scan_parameters('mp', 0, lab_one )
 n_layers        = n3 * mp

 MXDF_V          = DBLARR(  n_layers)
 MXDF_L          = LON64ARR(n_layers)

 MNDF_V          = DBLARR(  n_layers)
 MNDF_L          = LON64ARR(n_layers)

;FOR J = 1, n_layers DO BEGIN

glb_count        = 0

FOR J = first_slab, last_slab DO BEGIN

   SLCONE        = fetch_layer( J, I_step, lab_one )
   SLCTWO        = fetch_layer( J, I_step, lab_two )
   
   size_one      = SIZE(SLCONE)
   size_two      = SIZE(SLCTWO)


   IF (size_one[1] EQ size_two[1]) THEN BEGIN

   IF (STRCMP(qty, 'p')) THEN idx = 0
   IF (STRCMP(qty, 'a')) THEN idx = 1

   SLCDIFF     = ABS(SLCONE[*,idx] - SLCTWO[*,idx]) 

   too_small   = WHERE(SLCDIFF LE 1.0E-15, count_ts)

   IF (count_ts GT 0) THEN SLCDIFF[too_small] = 0.0
   
   FOR K = 0, size_one[1] - 1 DO BEGIN
     
     IF ( SLCDIFF[K] NE 0.0 ) THEN BEGIN

       IF ( SLCONE[K,idx] NE 0.0 ) THEN BEGIN

         SLCDIFF[K] = SLCDIFF[K] / ABS(SLCONE[K,idx])

       ENDIF ELSE BEGIN

         IF (SLCTWO[K, idx] NE 0.0) THEN BEGIN

           SLCDIFF[K] = SLCDIFF[K] / ABS(SLCTWO[K,idx])

         ENDIF ELSE BEGIN

           PRINT, "max_diff: WARNING - SLCDIFF is zero with both slices non-zero:"
           PRINT, "I_step = ", I_step
           PRINT, "J      = ", J
           PRINT, "K      = ", K

         ENDELSE

       ENDELSE

     ENDIF

   ENDFOR   

     tol_idx     = WHERE(SLCDIFF GE tolerance, count_tol)

;    IF ( count_tol GT 0 ) THEN PRINT, "tol_idx = ", tol_idx

     glb_count   = glb_count + count_tol

     mxdf        = MAX(SLCDIFF, SUBSCRIPT_MIN = loc_mndf )
     mndf        = MIN(SLCDIFF, SUBSCRIPT_MAX = loc_mxdf )


   ENDIF ELSE BEGIN

     mxdf        =  1.0e+10
     mndf        = -1.0e+10

   ENDELSE

   MXDF_V[J-1]   = mxdf
   MXDF_L[J-1]   = loc_mxdf

   MNDF_V[J-1]   = mndf
   MNDF_L[J-1]   = loc_mndf

 ENDFOR

 DF_OUT         = {df_out, mxdf_v:MXDF_V, mxdf_l:MXDF_L, mndf_v:MNDF_V, mndf_l:MNDF_L, df_count:glb_count}
 RETURN, DF_OUT

END
