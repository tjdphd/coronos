FUNCTION makeAxes, cc, size_jvox

  n_xtics         = 5
  n_ytics         = 5

  JAxisFont       = OBJ_NEW('IDLgrFont', 'helvetica', SIZE=12)
  JAxisLargeFont  = OBJ_NEW('IDLgrFont', 'helvetica', SIZE=24)

  JxTitle         = OBJ_NEW('IDLgrText', 'X',                                   FONT=JAxisFont     )
  JyTitle         = OBJ_NEW('IDLgrText', 'Y', UPDIR=[0,0,-1], BASELINE=[0,1,0], FONT=JAxisFont     )

  cx              = cc[*,0]
  cy              = cc[*,1]

  x_rng           = [0, size_jvox[0]]

  x_rng[0]        = x_rng[0] * cx[0]
  x_rng[1]        = x_rng[1] * cx[1]

  y_rng           = [0, size_jvox[1]]

  y_rng[0]        = y_rng[0] * cy[0]
  y_rng[1]        = y_rng[1] * cy[1]
 
  jx_tic_len      =cc[1,0]*cc[1,0] * size_jvox[0]
  jy_tic_len      =cc[1,1]*cc[1,1] * size_jvox[1]

  x_tics          = FLTARR(n_xtics)
  y_tics          = FLTARR(n_ytics)

  str_x_tics      = STRARR(n_xtics)
  str_y_tics      = STRARR(n_ytics)

  inc_xtic        = 1.0 * (x_rng[1] - x_rng[0])/(n_xtics - 1)
  inc_ytic        = 1.0 * (y_rng[1] - y_rng[0])/(n_ytics - 1)

  ITx             = INDGEN(n_xtics)
  ITy             = INDGEN(n_ytics)

  x_tics[ITx]     = (x_rng[0] + ITx * inc_xtic)/x_rng[1]
  y_tics[ITy]     = (y_rng[0] + ITy * inc_ytic)/y_rng[1]

  str_x_tics[ITx] = STRTRIM(x_tics[ITx], 2)
  str_y_tics[ITy] = STRTRIM(y_tics[ITy], 2)

  str_x_tics[ITx] = STRMID(str_x_tics[ITx], 0, 4)
  str_y_tics[ITy] = STRMID(str_y_tics[ITy], 0, 4)

  xtickLabels     = OBJ_NEW('IDLgrText', str_x_tics, FONT=JAxisFont)
  ytickLabels     = OBJ_NEW('IDLgrText', str_y_tics, FONT=JAxisFont)

  jVolAxes        = OBJARR(2,2)

;; X-Axes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  ;JxAxis

    jVolAxes[0,0]  = OBJ_NEW('IDLgrAxis',                             $
                              COLOR        = [0,0,0],                 $
                              DIRECTION    = 0,                       $
                              /EXACT,                                 $
                              GRIDSTYLE    = 0,                       $
                              LOCATION     = [0, 0],                  $
                              MAJOR        = 5,                       $
                              MINOR        = 4,                       $
                              RANGE        = x_rng,                   $
                              THICK        = 2,                       $
                              TICKLEN      = jx_tic_len,              $
                              TITLE        = JxTitle,                 $
                              DEPTH_TEST_FUNCTION = 8                 $
                             )
    ;JxAxisp1
     jVolAxes[0,1] = OBJ_NEW('IDLgrAxis',                             $
                              COLOR        = [0,0,0],                 $
                              DIRECTION    = 0,                       $
                              /EXACT,                                 $
                              GRIDSTYLE    = 0,                       $
                              LOCATION     = [0, y_rng[1]],           $
                              MAJOR        = 5,                       $
                              MINOR        = 4,                       $
                              /NOTEXT,                                $
                              TICKLEN      = jx_tic_len,              $
                              TICKDIR      = 1,                       $
                              RANGE        = x_rng,                   $
                              THICK        = 2,                       $
                              DEPTH_TEST_FUNCTION = 8                 $
                             )
;; Y-Axes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; ;JyAxis
   jVolAxes[1,0]   = OBJ_NEW('IDLgrAxis',                       $
                              COLOR        = [0,0,0],           $
                              DIRECTION   = 1,                  $
                              /EXACT,                           $
                              GRIDSTYLE   = 0,                  $
                              LOCATION    = [0,0],              $
                              MAJOR       = 5,                  $
                              MINOR       = 4,                  $
                              RANGE       = y_rng,              $
                              THICK       = 2,                  $
                              TICKLEN     = jy_tic_len,         $
                              TITLE       = JyTitle,            $
                              DEPTH_TEST_FUNCTION = 8           $
                             )
  ;JyAxisp1
   jVolAxes[1,1]   = OBJ_NEW('IDLgrAxis',                       $
                              COLOR        = [0,0,0],           $
                              DIRECTION   = 1,                  $
                              /EXACT,                           $
                              GRIDSTYLE   = 0,                  $
                              LOCATION=[x_rng[1],0],            $
                              MAJOR       = 5,                  $
                              MINOR       = 4,                  $
                              /NOTEXT,                          $
                              RANGE       = y_rng,              $
                              THICK       = 2,                  $
                              TICKDIR     = 1,                  $
                              TICKLEN     = jy_tic_len,         $
                              DEPTH_TEST_FUNCTION = 8           $
                             )

  RETURN, jVolAxes

END
