PRO cross_check, qty       , $
                 first_step, $ 
                 last_step , $
                 first_slab, $ 
                 last_slab , $
                 str_res   , $
                 lab_one   , $
                 lab_two   , $
                 tolerance

FOR I = first_step, last_step DO BEGIN

  MDF_OUT    = max_diff(qty, I, first_slab, last_slab, str_res, lab_one, lab_two, tolerance)

  glb_max_df = MAX(MDF_OUT.mxdf_v)
  glb_min_df = MIN(MDF_OUT.mndf_v)

; IF (MDF_OUT.df_count NE 0) THEN BEGIN

  PRINT, " "
  PRINT, " ************************************************** "
  PRINT, " "
  PRINT, "glb_max_df                = ", glb_max_df
  PRINT, "glb_min_df                = ", glb_min_df
  PRINT, "count                     = ", MDF_OUT.df_count
  PRINT, " "
  PRINT, " ************************************************** "
  PRINT, " "

; ENDIF

ENDFOR

END
