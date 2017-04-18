FUNCTION  makeContourPlot, n_slice,           $
                           qty,               $
                           n_step,            $
                           tot_steps,         $
                           desc_label,        $
                           global_qty_minmax, $
                           KK,                $
                           QCNTRS = Q_CNTRS,  $
                           PREFIX = prefix

  IF (NOT KEYWORD_SET(Q_CNTRS) ) THEN BEGIN
    n_cntrs            = 31
  ENDIF ELSE BEGIN
    size_q_cntrs       = SIZE(Q_CNTRS,/DIMENSIONS)
    n_cntrs            = size_q_cntrs[0]
  ENDELSE
  IF (NOT KEYWORD_SET(PREFIX) ) THEN BEGIN
    prefix             = 'rmct2'
  ENDIF

  !EXCEPT              =  2

  i_x                  =  6
  i_j                  =  5
  i_o                  =  4
  i_a                  =  3
  i_p                  =  2

  i_qty                = -1

  CASE qty OF
  'j' : i_qty          =  i_j
  'o' : i_qty          =  i_o
  'a' : i_qty          =  i_a
  'p' : i_qty          =  i_p
  ELSE: i_qty          = -1
  ENDCASE

  ip1                  = scan_parameters('p1', 0, desc_label)
  ip2                  = scan_parameters('p2', 0, desc_label)
  
  n1                   = 2^ip1
  n2                   = 2^ip2   
  
  i_x_res              = LONG64(n1)
  i_y_res              = LONG64(n2)
  
  dat_file             = fetch_datafile(n_slice, n_step, desc_label     )
  SLB                  = fetch_layer(   n_slice, n_step, desc_label, KK )
  
  Q                    = REFORM(SLB[*,i_qty], n1, n2)
  XX                   = REFORM(SLB[*,0    ], n1, n2)
  YY                   = REFORM(SLB[*,1    ], n1, n2)

  xmin                 = MIN(XX)
  xmax                 = MAX(XX)

  X                    = REFORM(XX[0,*])
  Y                    = REFORM(YY[*,0])

  XX                   = FLTARR(n1+1)
  XX[0:n1-1]           = X[*]
  XX[n1]               = X[n1 - 1] + (X[n1-1] - X[n1-2])

  YY                   = FLTARR(n2+1)
  YY[0:n2-1]           = Y[*]
  YY[n2]               = Y[n2 - 1] + (Y[n2-1] - Y[n2-2])

  QQ                   = FLTARR(n1+1,n2+1)
  QQ[0:n1-1,0:n2-1]    = Q[ *,*]
  QQ[n1,* ]            = QQ[0,*]
  QQ[* ,n2]            = QQ[*,0]

  size_qq              = SIZE(QQ, /DIMENSIONS)
  size_xx              = SIZE(XX, /DIMENSIONS)
  size_yy              = SIZE(YY, /DIMENSIONS)

  X                    = XX
  Y                    = YY
  Q                    = QQ

  ymin                 = MIN(Y)
  ymax                 = MAX(Y)

  size_q               = SIZE(Q,/DIMENSIONS)

  QBuffer              = OBJ_NEW('IDLgrBuffer')
  QScene               = OBJ_NEW('IDLgrScene')

  ContView             = OBJ_NEW('IDLgrView',                                $
                                  LOCATION       = [-0.5, -0.625],           $
                                  VIEWPLANE_RECT = [-1.0 , -1.0 , 2.0, 2.0], $
                                  DIMENSIONS     = [1.5, 1.5],               $
                                  UNITS          = 3,                        $
                                  COLOR          = [180,180,180]             $
                                )
                                               
  QScene               ->Add, ContView
            
  QPalette             = OBJ_NEW('IDLgrPalette')
  QPalette             -> LoadCT, 33
            
  QPalette             -> GetProperty, N_COLORS     = n_colors
  QPalette             -> GetProperty, BLUE_VALUES  = b_values
  QPalette             -> GetProperty, RED_VALUES   = r_values
  QPalette             -> GetProperty, GREEN_VALUES = g_values
            
  color_skip           = n_colors / n_cntrs
  rgb_colors           = INTARR(3,  n_cntrs)
            
  col_idx              = 0
            
  FOR I = 0, n_cntrs - 1 DO BEGIN
 
    rgb_colors[0,I]    = r_values[col_idx]
    rgb_colors[1,I]    = g_values[col_idx]
    rgb_colors[2,I]    = b_values[col_idx]

    col_idx            = col_idx + color_skip

  ENDFOR

  ctr_thick            = INTARR(size_q_cntrs)
  ctr_thick[*]         =  1.5

  qmax                 = MAX(Q)
  qmin                 = MIN(Q)
  xmax                 = MAX(X)
  xmin                 = MIN(X)
  ymax                 = MAX(Y)
  ymin                 = MIN(Y)
  cmax                 = MAX(Q_CNTRS)
  cmin                 = MIN(Q_CNTRS)
  
  QCont                = OBJ_NEW( 'IDLgrContour',            $
                                   DATA_VALUES = Q,          $
                                   GEOMX       = X,          $
                                   GEOMY       = Y,          $
                                   GEOMZ       = Q,          $
                                   C_VALUE     = Q_CNTRS,    $
                                   C_COLOR     = [0,0,0],    $
                                   PLANAR      = 1,          $
                                   C_THICK     = ctr_thick   $
                                )

  QcMap                = OBJ_NEW( 'IDLgrContour',            $
                                   DATA_VALUES = Q,          $
                                   GEOMX       = X,          $
                                   GEOMY       = Y,          $
                                   C_VALUE     = Q_CNTRS,    $
                                   C_COLOR     = rgb_colors, $
                                   FILL        = 1,          $
                                   PLANAR      = 1           $
                                ) 

   ContModel           =  OBJ_NEW('IDLgrModel')

   cc                  = getCoordConvs(desc_label, size_q)
   ContAxes            = makeAxes(cc, size_q)

   str_title           = getTitle(qty, desc_label, n_step, n_slice, prefix, global_qty_minmax)

   CTR_TTLfont         = OBJ_NEW('IDLgrFont', 'helvetica*bold', SIZE=16)

   ContTitle           = OBJ_NEW( 'IDLgrText',                          $
                                   str_title,                           $
                                   LOCATIONS         = [-0.15 , 1.05],  $
                                   DEPTH_TEST_FUNCTION = 8,             $
                                   ALIGNMENT         = 0.0,             $
                                   FONT              = CTR_TTLfont,     $
                                   ENABLE_FORMATTING = 1                $
                                )

   BarModel            = OBJ_NEW('IDLgrModel')
   QColorBar           = makeColorBar(qty, r_values, g_values, b_values, QCNTRS=Q_CNTRS)

   ContModel           -> Add, Qcont
   ContModel           -> Add, QcMap

   ContModel           -> Add, ContTitle

   ContModel           -> Add, ContAxes[0,0] ; bottom x axis
   ContModel           -> Add, ContAxes[0,1] ; top    x axis
   ContModel           -> Add, ContAxes[1,0] ; left   y axis
   ContModel           -> Add, ContAxes[1,1] ; right  y axis


   BarModel            -> Add, QColorBar
   BarModel            -> Translate, -0.20,0.0, 0.0 

   ContView            -> Add, ContModel
   ContView            -> Add, BarModel

   QBuffer             -> Draw, QScene

   str_ifx             = getInfixString(desc_label, "ctr",  tot_steps, 0.0)

   sx = 0.87 & sy = 0.87 & sz = 0.87

   ContModel           -> Scale, sx, sy, sz
   BarModel            -> Scale, sx, sy, sz

   IF (0 EQ 1) THEN BEGIN

     ; this doesn't work yet
     output_path       = getOutputPath(n_step, qty, str_ifx, tot_steps, 'tif')

     QImage            =  QBuffer->Read()
     QImage            -> GetProperty, DATA = image
     image             =  REVERSE(image,3)

     WRITE_TIFF, output_path, image

   ENDIF ELSE BEGIN

     output_path       = getOutputPath(n_step, qty, str_ifx, tot_steps, 'eps')



     QBuffer           -> GetProperty, Resolution = screenResolution
     QBuffer           -> GetProperty, Dimensions = windowSize

     QClipboard        = OBJ_NEW('IDLgrClipboard',                         $
                                  QUALITY         =  2,                    $
                                  DIMENSIONS      =  2.0*windowSize,       $
                                  RESOLUTION      =  0.5*screenResolution, $
                                  GRAPHICS_TREE   =      QScene)

     QClipboard        ->         Draw,                                    $
                                  VECTOR          = 0,                     $
                                  POSTSCRIPT      = 1,                     $
                                  FILENAME        = output_path

     QClipboard        ->         SetProperty,                             $
                                  GRAPHICS_TREE   = OBJ_NEW()

     OBJ_DESTROY, QClipboard

   ENDELSE

    OBJ_DESTROY, QBuffer


  RETURN, 0

END
