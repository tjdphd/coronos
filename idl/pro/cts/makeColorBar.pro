FUNCTION  makeColorBar, qty, r_values, g_values, b_values, $
                        QCNTRS = Q_CNTRS

  bar_dims         = [0.05, 1.0 ]  

  CASE qty OF

  'p'  : qty_label = '!9F!X'
  'a'  : qty_label = 'A'
  'o'  : qty_label = '!9W!X'
  'j'  : qty_label = 'J'
  ELSE : qty_label = 'wtf'
  ENDCASE
  
  Bar_TTLfont      = OBJ_NEW('IDLgrFont', 'helvetica*bold', SIZE=16)

  BarTitle         = OBJ_NEW( 'IDLgrText',                       $
                              ENABLE_FORMATTING = 1,             $
                              FONT              = Bar_TTLfont,   $
                              qty_label                          $
                            )

  ContColorBar     = OBJ_NEW( 'IDLgrColorBar',                   $
                               r_values,                         $
                               g_values,                         $
                               b_values,                         $
                               MAJOR        = 5,                 $
                               TITLE        = BarTitle,          $
                               DIMENSIONS   = bar_dims,          $
                               SHOW_AXIS    = 1,                 $
                               SHOW_OUTLINE = 1                  $
                             )

  tick_vals        = FLTARR(5)
  tick_labs        = OBJARR(5)

  size_tv          = SIZE(tick_vals, /DIMENSIONS)

  qmax             = MAX(Q_CNTRS)
  qmin             = MIN(Q_CNTRS)

  q_inc            = (qmax - qmin) / (size_tv[0] - 1)

 tick_vals[0]      = qmin
 IF (tick_vals[0] LT 0.0) THEN BEGIN
   tick_labs[0]    = OBJ_NEW('IDLgrText',STRMID(STRTRIM(tick_vals[0],2),0, 6))
 ENDIF ELSE BEGIN
   tick_labs[0]    = OBJ_NEW('IDLgrText',STRMID(STRTRIM(tick_vals[0],2),0, 5))
 ENDELSE

 FOR I = 1, size_tv[0] - 1 DO BEGIN

   tick_vals[I]    = tick_vals[I - 1] + q_inc 

   IF (tick_vals[I] LT 0.0) THEN BEGIN
     tick_labs[I]  = OBJ_NEW('IDLgrText', STRMID(STRTRIM(tick_vals[I],2),0, 6))
   ENDIF ELSE BEGIN
     tick_labs[I]  = OBJ_NEW('IDLgrText', STRMID(STRTRIM(tick_vals[I],2),0, 5))
   ENDELSE


 ENDFOR

  ContColorBar -> SetProperty, TICKTEXT=tick_labs

  RETURN, ContColorBar

END
