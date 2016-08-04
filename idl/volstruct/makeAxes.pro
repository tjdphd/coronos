FUNCTION makeAxes, cc, size_jvox

  n_ztics         = 5
  n_xtics         = 5
  n_ytics         = 5

  JAxisFont       = OBJ_NEW('IDLgrFont', 'helvetica', SIZE=12)
  JAxisLargeFont  = OBJ_NEW('IDLgrFont', 'helvetica', SIZE=24)

  JzTitle         = OBJ_NEW('IDLgrText', 'Z',                                   FONT=JAxisLargeFont)
  JxTitle         = OBJ_NEW('IDLgrText', 'X',                                   FONT=JAxisFont     )
  JyTitle         = OBJ_NEW('IDLgrText', 'Y', UPDIR=[0,0,-1], BASELINE=[0,1,0], FONT=JAxisFont     )

  cx              = cc[*,0]
  cy              = cc[*,1]
  cz              = cc[*,2]

  z_rng           = [0, size_jvox[0]]
  x_rng           = [0, size_jvox[1]]
  y_rng           = [0, size_jvox[2]]

  jx_tic_len      = 0.16*size_jvox[1]
  jy_tic_len      = 0.16*size_jvox[2]
  jz_tic_len      = 0.02*size_jvox[0]

  zL              = z_rng[1] / 32

  z_tics          = FLTARR(n_ztics)
  x_tics          = FLTARR(n_xtics)
  y_tics          = FLTARR(n_ytics)

  str_z_tics      = STRARR(n_ztics)
  str_x_tics      = STRARR(n_xtics)
  str_y_tics      = STRARR(n_ytics)

  inc_ztic        = zL*(z_rng[1] - z_rng[0])/(n_ztics - 1)
  inc_xtic        = 1.0 * (x_rng[1] - x_rng[0])/(n_xtics - 1)
  inc_ytic        = 1.0 * (y_rng[1] - y_rng[0])/(n_ytics - 1)

  ITz             = INDGEN(n_ztics)
  ITx             = INDGEN(n_xtics)
  ITy             = INDGEN(n_ytics)

  z_tics[ITz]     = (z_rng[0] + ITz * inc_ztic)/z_rng[1]
  x_tics[ITx]     = (x_rng[0] + ITx * inc_xtic)/x_rng[1]
  y_tics[ITy]     = (y_rng[0] + ITy * inc_ytic)/y_rng[1]

  str_z_tics[ITz] = STRTRIM(z_tics[ITz], 2)
  str_x_tics[ITx] = STRTRIM(x_tics[ITx], 2)
  str_y_tics[ITy] = STRTRIM(y_tics[ITy], 2)

  str_z_tics[ITz] = STRMID(str_z_tics[ITz], 0, 4)
  str_x_tics[ITx] = STRMID(str_x_tics[ITx], 0, 4)
  str_y_tics[ITy] = STRMID(str_y_tics[ITy], 0, 4)

  zrs             = WHERE(y_tics EQ 0.0)

  str_y_tics[zrs] = ''


  ztickLabels     = OBJ_NEW('IDLgrText', str_z_tics, FONT=JAxisLargeFont)
  xtickLabels     = OBJ_NEW('IDLgrText', str_x_tics, FONT=JAxisFont)
  ytickLabels     = OBJ_NEW('IDLgrText', str_y_tics, FONT=JAxisFont)

  jVolAxes        = OBJARR(4,3)

; Z-Axes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; JzAxis
  jVolAxes[0,0]   = OBJ_NEW('IDLgrAxis',                           $
                             COLOR       = [0,0,0],                $
                             DIRECTION   = 0,                      $
                             /EXACT,                               $
                             GRIDSTYLE   = 0,                      $
                             LOCATION    = [0, y_rng[0], 0],       $
                             MAJOR       = n_ztics,                $
                             MINOR       = 4,                      $
                             /NOTEXT,                              $
                             RANGE       = z_rng,                  $
                             THICK       = 1,                      $
                             TICKDIR     = 0,                      $
                             TICKLEN     = 0,                      $
                             XCOORD_CONV = cz,                     $
                             YCOORD_CONV = cx,                     $
                             ZCOORD_CONV = cy                      $
                           )
  ; JzAxisp1
  jVolAxes[1,0]   = OBJ_NEW('IDLgrAxis',                           $ 
                             COLOR       = [0,0,0],                $
                             DIRECTION   = 0,                      $
                             /EXACT,                               $
                             GRIDSTYLE   = 0,                      $
                             LOCATION    = [0, y_rng[1], 0],       $
                             MAJOR       = n_ztics,                $
                             MINOR       = 4,                      $
                             /NOTEXT,                              $
                             RANGE       = z_rng,                  $
                             THICK       = 2,                      $
                             TICKDIR     = 0,                      $
                             TICKLEN     = 0,                      $
                             XCOORD_CONV = cz,                     $
                             YCOORD_CONV = cx,                     $
                             ZCOORD_CONV = cy                      $
                           )

  ; JzAxisp2
  jVolAxes[2,0]   = OBJ_NEW('IDLgrAxis',                           $
                             COLOR       = [0,0,0],                $
                             DIRECTION   = 0,                      $
                             /EXACT,                               $
                             GRIDSTYLE   = 0,                      $
                             LOCATION    = [0, y_rng[1], x_rng[1]],$
                             MAJOR       = n_ztics,                $
                             MINOR       = 4,                      $
                             /NOTEXT,                              $
                             RANGE       = z_rng,                  $
                             THICK       = 2,                      $
                             TICKDIR     = 1,                      $
                             TICKLEN     = 0,                      $
                             XCOORD_CONV = cz,                     $
                             YCOORD_CONV = cx,                     $
                             ZCOORD_CONV = cy                      $
                           )

; JzAxisp3
  jVolAxes[3,0]   = OBJ_NEW('IDLgrAxis',                           $
                             COLOR        = [0,0,0],               $
                             DIRECTION    = 0,                     $
                             /EXACT,                               $
                             GRIDSTYLE    = 0,                     $
                             LOCATION     = [0, 0, x_rng[1]],      $
                             MAJOR        = n_ztics,               $
                             MINOR        = 4,                     $
                             RANGE        = z_rng,                 $
                             TEXTBASELINE =  [1, 0,  0],           $
                             TEXTPOS      = 0,                     $
                             TEXTUPDIR    =  [0, 1,  0],           $
                             THICK        = 2,                     $
                             TICKLEN      = jz_tic_len,            $
                             TITLE        = JzTitle,               $
                             TICKTEXT     = ztickLabels,           $
                             XCOORD_CONV  = cz,                    $
                             YCOORD_CONV  = cx,                    $
                             ZCOORD_CONV  = cy                     $
                           )

; X-Axes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;JxAxis
   jVolAxes[0,1]  = OBJ_NEW('IDLgrAxis',                             $
                             COLOR        = [0,0,0],                 $
                             DIRECTION    = 2,                       $
                             /EXACT,                                 $
                             GRIDSTYLE    = 0,                       $
                             LOCATION     = [0, 0,999],              $
                             MAJOR        = 5,                       $
                             MINOR        = 4,                       $
                             RANGE        = x_rng,                   $
                             TEXTBASELINE = [1, 0,  0],              $   
                             TEXTUPDIR    = [0, 0, -1],              $
                             THICK        = 2,                       $
                             TICKLEN      = jx_tic_len,              $
                             TITLE        = JxTitle,                 $
                             TICKTEXT     = xtickLabels,             $
                             XCOORD_CONV  = cz,                      $
                             YCOORD_CONV  = cx,                      $
                             ZCOORD_CONV  = cy                       $
                            )
   ;JxAxisp1
    jVolAxes[1,1] = OBJ_NEW('IDLgrAxis',                             $
                             COLOR        = [0,0,0],                 $
                             DIRECTION    = 2,                       $
                             /EXACT,                                 $
                             GRIDSTYLE    = 0,                       $
                             LOCATION     = [0, y_rng[1], 0],        $
                             MAJOR        = 5,                       $
                             MINOR        = 4,                       $
                             /NOTEXT,                                $
                             RANGE        = x_rng,                   $
                             THICK        = 2,                       $
                             TICKLEN      = 0.0,                     $
                             XCOORD_CONV  = cz,                      $
                             YCOORD_CONV  = cx,                      $
                             ZCOORD_CONV  = cy                       $
                            )

   ;JxAxisp2
    jVolAxes[2,1] = OBJ_NEW('IDLgrAxis',                             $
                             COLOR        = [0,0,0],                 $
                             DIRECTION    = 2,                       $
                             /EXACT,                                 $
                             GRIDSTYLE    = 0,                       $
                             LOCATION     = [z_rng[1], y_rng[1], 0], $
                             MAJOR        = 5,                       $
                             MINOR        = 4,                       $
                             /NOTEXT,                                $
                             RANGE        = x_rng,                   $
                             THICK        = 2,                       $
                             TICKLEN      = 0.0,                     $
                             XCOORD_CONV  = cz,                      $
                             YCOORD_CONV  = cx,                      $
                             ZCOORD_CONV  = cy                       $
                            )

   ;JxAxisp3
    jVolAxes[3,1] = OBJ_NEW('IDLgrAxis',                             $
                             COLOR        = [0,0,0],                 $
                             DIRECTION    = 2,                       $
                             /EXACT,                                 $
                             GRIDSTYLE    = 0,                       $
                             LOCATION     = [z_rng[1], 0, 0],        $
                             MAJOR        = 5,                       $
                             MINOR        = 4,                       $
                             /NOTEXT,                                $
                             RANGE        = x_rng,                   $
                             THICK        = 2,                       $
                             TICKLEN      = 0.0,                     $
                             XCOORD_CONV  = cz,                      $
                             YCOORD_CONV  = cx,                      $
                             ZCOORD_CONV  = cy                       $
                            )

; Y-Axes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;JyAxis
  jVolAxes[0,2]   = OBJ_NEW('IDLgrAxis',                       $
                             COLOR        = [0,0,0],           $
                             DEPTH_TEST_FUNCTION=8,            $
                             DIRECTION   = 1,                  $
                             /EXACT,                           $
                             GRIDSTYLE   = 0,                  $
                             LOCATION    = [0,0,0],            $
                             MAJOR       = 5,                  $
                             MINOR       = 4,                  $
                             RANGE       = y_rng,              $
                             THICK       = 2,                  $
                             TICKLEN     = jy_tic_len,         $
                             TICKTEXT    = ytickLabels,        $
                             TITLE       = JyTitle,            $
                             XCOORD_CONV = cz,                 $
                             YCOORD_CONV = cx,                 $
                             ZCOORD_CONV = cy                  $
                            )
 ;JyAxisp1
  jVolAxes[1,2]   = OBJ_NEW('IDLgrAxis',                       $
                             COLOR        = [0,0,0],           $
                             DIRECTION   = 1,                  $
                             /EXACT,                           $
                             GRIDSTYLE   = 0,                  $
                             LOCATION=[0, 0, x_rng[1]],        $
                             MAJOR       = 5,                  $
                             MINOR       = 4,                  $
                             /NOTEXT,                          $
                             RANGE       = y_rng,              $
                             THICK       = 2,                  $
                             TICKLEN     = 0.0,                $
                             XCOORD_CONV = cz,                 $
                             YCOORD_CONV = cx,                 $
                             ZCOORD_CONV = cy                  $
                            )
 ;JyAxisp2
  jVolAxes[2,2]   = OBJ_NEW('IDLgrAxis',                       $
                             COLOR        = [0,0,0],           $
                             DIRECTION   = 1,                  $
                             /EXACT,                           $
                             GRIDSTYLE   = 0,                  $
                             LOCATION=[z_rng[1], 0, x_rng[1]], $
                             MAJOR       = 5,                  $
                             MINOR       = 4,                  $
                             /NOTEXT,                          $
                             RANGE       = y_rng,              $
                             THICK       = 2,                  $
                             TICKLEN     = 0.0,                $
                             XCOORD_CONV = cz,                 $
                             YCOORD_CONV = cx,                 $
                             ZCOORD_CONV = cy                  $
                            )
 ;JyAxisp3
  jVolAxes[3,2]   = OBJ_NEW('IDLgrAxis',                       $
                             COLOR        = [0,0,0],           $
                             DEPTH_TEST_FUNCTION=8,            $
                             DIRECTION   = 1,                  $
                             /EXACT,                           $
                             GRIDSTYLE   = 0,                  $
                             LOCATION=[z_rng[1], 0, 0],        $
                             MAJOR       = 5,                  $
                             MINOR       = 4,                  $
                             /NOTEXT,                          $
                             RANGE       = y_rng,              $
                             THICK       = 2,                  $
                             TICKLEN     = 0.0,                $
                             XCOORD_CONV = cz,                 $
                             YCOORD_CONV = cx,                 $
                             ZCOORD_CONV = cy                  $
                            )

  RETURN, jVolAxes

END
