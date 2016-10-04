FUNCTION  plot_surface, n_slice,           $
                        qty,               $
                        n_step,            $
                        tot_steps,         $
                        desc_label,        $
                        global_qty_minmax, $
                        KK,                $
                        QCNTRS=Q_CNTRS

  IF (NOT KEYWORD_SET(Q_CNTRS) ) THEN BEGIN
    n_cntrs       = 15
  ENDIF
  !EXCEPT         =  2

  i_x             =  6
  i_j             =  5
  i_o             =  4
  i_a             =  3
  i_p             =  2

  i_qty           = -1

  CASE qty OF
  'j' : i_qty     =  i_j
  'o' : i_qty     =  i_o
  'a' : i_qty     =  i_a
  'p' : i_qty     =  i_p
  ELSE: i_qty     = -1
  ENDCASE

  ip1             = scan_parameters('ip1', 0, desc_label)
  ip2             = scan_parameters('ip2', 0, desc_label)
  
  n1              = 2^ip1
  n2              = 2^ip2   
  
  i_x_res         = LONG64(n1)
  i_y_res         = LONG64(n2)
  
  dat_file        = fetch_datafile(n_slice, n_step, desc_label     )
  SLB             = fetch_layer(   n_slice, n_step, desc_label, KK )
  
  Q               = REFORM(SLB[*,i_qty], n1, n2)
  XX              = REFORM(SLB[*,0    ], n1, n2); <- works
  YY              = REFORM(SLB[*,1    ], n1, n2); <- works

  X               = REFORM(XX[0,*])
  Y               = REFORM(YY[*,0])

  xmin            = MIN(X)
  xmax            = MAX(X)

  x_range         = [xmin, xmax]

  ymin            = MIN(Y)
  ymax            = MAX(Y)

  y_range         = [ymin, ymax]

  PRINT, "size Q:", SIZE(Q,/DIMENSIONS)
  PRINT, "size X:", SIZE(X,/DIMENSIONS)
  PRINT, "size Y:", SIZE(Y,/DIMENSIONS)
  PRINT, ""

  qmax             = MAX(Q)
  qmin             = MIN(Q)

  z_range         = [qmin, qmax]

  PRINT, 'qmax = ', qmax
  PRINT, 'qmin = ', qmin

  sfc_title = "THIS IS THE TITLE"

  surf = SURFACE(Q, X, Y,                                                          $
                 STYLE        = 2,                                                 $
                 PERSPECTIVE  = 1,                                                 $
                 TITLE        = sfc_title,                                         $
                 DEPTH_CUE    = [0,2],                                             $
                 AXIS_STYLE   = 1,                                                 $
                 FONT_SIZE    = 16,                                                $
                 COLOR        = 'blue',                                            $
                 XRANGE       = x_range,                                           $
                 YRANGE       = y_range,                                           $
                 ZRANGE       = z_range,                                           $
                 TRANSPARENCY = 00,                                                $
                 CLIP         = 0                                                  $
               )
;  QBuffer          = OBJ_NEW('IDLgrBuffer')
;  QScene           = OBJ_NEW('IDLgrScene')
;  SurfView         = OBJ_NEW('IDLgrView', LOCATION=[0.125, 0.0], VIEWPLANE_RECT=[-0.25, -0.25, 2.0, 1.5], $
;                                              DIMENSIONS=[0.9, 0.9], UNITS=3, COLOR=[255,255,255])
;                                              
;  QScene->Add, SurfView
;           
;  QPalette         = OBJ_NEW('IDLgrPalette')
;  QPalette         -> LoadCT, 33
;           
;  QPalette         -> GetProperty, N_COLORS     = n_colors
;  QPalette         -> GetProperty, BLUE_VALUES  = b_values
;  QPalette         -> GetProperty, RED_VALUES   = r_values
;  QPalette         -> GetProperty, GREEN_VALUES = g_values
;           
;  color_skip       = n_colors / n_cntrs
;  rgb_colors       = INTARR(3, n_cntrs)
;           
;  col_idx          = 0
;           
;  FOR I=0,n_cntrs-1 DO BEGIN
;
;    rgb_colors[0,I]   = r_values[col_idx]
;     rgb_colors[1,I]  = g_values[col_idx]
;     rgb_colors[2,I]  = b_values[col_idx]
;     col_idx          = col_idx + color_skip
;
;  ENDFOR
;
;; QSurf                = OBJ_NEW('IDLgrSurface', Q, X, Y, $
;  QSurf                = OBJ_NEW('IDLgrSurface', Q, $
;                                 STYLE = 0          $
;                                )
;
;  SurfModel         = OBJ_NEW('IDLgrModel')
;  SurfModel->Add, Qsurf
;  SurfView->Add,  SurfModel
;
;  QBuffer->Draw, QScene
;  QBuffer->SetProperty, GRAPHICS_TREE = QScene
;  QBuffer->GetProperty, DIMENSIONS    = windowSize
;  QBuffer->GetProperty, RESOLUTION    = screenResolution
;  QClipboard=OBJ_NEW('IDLgrClipboard',  QUALITY       = 2,                 $
;                                        DIMENSIONS    = windowSize,        $
;                                        RESOLUTION    = screenResolution,  $
;                                        GRAPHICS_TREE = QScene)
;   
;  eps_prefix                             = qty + '_contour'
;  eps_infix                              = '-ff-spec'
;  postfile             = construct_eps_outfile_name(qty, eps_prefix, eps_infix, n_slice, n_step, tot_steps) 
;
;  QClipboard->Draw,     FILENAME         = postfile, /POSTSCRIPT, /VECTOR 
;  QBuffer->SetProperty, QUALITY          = 2
;  QBuffer->SetProperty, PALETTE          = QPalette
;  QBuffer->SetProperty, GRAPHICS_TREE    = QScene
;  QBuffer->GetProperty, ZBUFFER_DATA     = Pic 
;  
;  QClipboard->SetProperty, GRAPHICS_TREE = OBJ_NEW()
;  OBJ_DESTROY, QClipboard
; RETURN, Pic
  RETURN, 0

END
