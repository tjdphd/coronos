FUNCTION construct_jcontour, n_slice, qty, n_step, tot_steps, desc_label

COMMON save_r, Pre
COMMON First, first_call
COMMON fnc, prefix, inc_res, nfld
COMMON CC, Cntrs
COMMON JC, J_Cntrs
COMMON loc, eps_out_dir

ip1      = scan_parameters('ip1', 0, desc_label)
ip2      = scan_parameters('ip2', 0, desc_label)

x_res    = 2^ip1
y_res    = 2^ip2   

i_x_res  = LONG64(x_res)
i_y_res  = LONG64(y_res)

dat_file = fetch_datafile(n_slice, n_step, desc_label)            ; Get name of data file
Slc      = fetch_slices(n_slice, n_slice, n_step, desc_label)     ; prepare data array for slice

size_slc = SIZE(Slc)
rank     = size_slc[0]
nlines   = LONG64(0)

IF (rank EQ 2) THEN BEGIN

  nlines = size_slc[1]
  cols   = size_slc[2]

ENDIF ELSE BEGIN
          IF (rank EQ 3) THEN BEGIN
             nlines = size_slc[2]
             cols   = size_slc[3]
          ENDIF
      ENDELSE

IF (nlines NE (i_x_res * i_y_res)) THEN BEGIN

  PRINT, FORMAT = '(a38,1x,i4,1x,a11,1x,i4,1x,a18)', 'construct_jcontour: WARNING - skipping slice', $
                                                      n_slice, 'for sub-run', n_step,          $
                                                     'data not available'
  RETURN, -1

ENDIF

n_cntrs  = SIZE(Cntrs)                                            ; get size/rank info for Cntrs

CLABS    = INTARR(1,n_cntrs[1])                                   ; create contour label array
CLABS    = TRANSPOSE(CLABS)

   FOR I = 1,n_cntrs[1] DO BEGIN                                  ; initialize CLABS

     IF ( I MOD 6 EQ 0 ) THEN BEGIN
               CLABS[I-1] = 1
     ENDIF ELSE BEGIN
               CLABS[I-1] = 0
           ENDELSE
              CLABS[I-1]  = 0
   ENDFOR

A             = FLTARR(cols - 2, nlines)                           ; Allocate space for plotting
X             = FLTARR(1,        (i_x_res + 1)*(i_y_res + 1) )                           ; arrays
Y             = FLTARR(1,        (i_x_res + 1)*(i_y_res + 1) )

A             = TRANSPOSE(A)  
X             = TRANSPOSE(X)
Y             = TRANSPOSE(Y)  

dx = 1.0 / x_res
dy = 1.0 / y_res

PRINT, 'construct_jcontour: dx = ', dx
PRINT, 'construct_jcontour: dy = ', dy

IF (rank EQ 2) THEN BEGIN

  FOR idx = 0, i_x_res DO BEGIN
    FOR jdx = 0, i_y_res DO BEGIN
  
      cur_ndx    = LONG64((idx * i_x_res) + jdx)

      X[cur_ndx] = LONG64(idx)*dx
      Y[cur_ndx] = LONG64(jdx)*dy

    ENDFOR
  ENDFOR
  
  A[*,0:cols-3]  = Slc[*,2:cols-1]

ENDIF ELSE BEGIN

          FOR idx = 0, i_x_res DO BEGIN
            FOR jdx = 0, i_y_res DO BEGIN
          
              cur_ndx    = LONG64((idx * (i_x_res+1)) + jdx)
          
              X[cur_ndx] = LONG64(idx)*dx
              Y[cur_ndx] = LONG64(jdx)*dy

            ENDFOR
          ENDFOR
  
          A[*,0:cols-3] = Slc[0, *, 2:cols-1]

      ENDELSE
      
time     = scan_parameters('tstart', n_step, desc_label)
str_time = STRTRIM(time,2)
str_time = STRMID(str_time,0,5)

IF (qty NE 'j' ) THEN BEGIN

          CASE qty OF
          'vz': i_qty     =  3
          'bz': i_qty     =  2
          'a' : i_qty     =  1
          'p' : i_qty     =  0
          ELSE: i_qty     = -1
          ENDCASE

          n3              = scan_parameters('n3' , 0, desc_label)
          mp              = scan_parameters('mp' , 0, desc_label)

          tot_slc         = mp * n3

          str_tot_slc     = STRTRIM(tot_slc,2)
          str_slc         = STRTRIM(n_slice,2)

          str_layer       = '(layer ' + str_slc + ' of ' + str_tot_slc + ')'
          
          CASE qty OF
          'vz': str_title = 'V_z ' + str_layer + ': t = ' + str_time
          'bz': str_title = 'Z   ' + str_layer + ': t = ' + str_time
          'a' : str_title = 'A   ' + str_layer + ': t = ' + str_time
          'p' : str_title = 'P   ' + str_layer + ': t = ' + str_time
          ELSE: str_title = 'default title string'
          ENDCASE

          Q_sqr             = FLTARR(i_x_res + 1, i_y_res + 1)

          FOR I=0, i_x_res - 1 DO BEGIN
             FOR J=0, i_y_res - 1 DO BEGIN
          
                idx         = (I * i_x_res) + J
                Q_sqr[I, J] = A[idx, i_qty]
          
             ENDFOR
          ENDFOR

          Q_sqr[*,       i_y_res ] = Q_sqr[*, 0 ]
          Q_sqr[i_x_res, *       ] = Q_sqr[0, * ]

          QBuffer     = OBJ_NEW('IDLgrBuffer')

          QScene      = OBJ_NEW('IDLgrScene')
          CntrView    = OBJ_NEW('IDLgrView', LOCATION=[0.125, 0.0], VIEWPLANE_RECT=[-0.25, -0.25, 2.0, 1.5], $
                                             DIMENSIONS=[0.9, 0.9], UNITS=3, COLOR=[255,255,255])
                                             
          BarView     = OBJ_NEW('IDLgrView', LOCATION=[0.75,  0.0], VIEWPLANE_RECT=[-0.25, -0.25, 2.0, 1.5], $
                                             DIMENSIONS=[0.75,0.9], UNITS=3, COLOR=[255,255,255])
          QScene->Add, CntrView
          QScene->Add, barView
          
          QPalette   = OBJ_NEW('IDLgrPalette')
          QPalette->LoadCT, 33
          
          QPalette->GetProperty, N_COLORS     = n_colors
          QPalette->GetProperty, BLUE_VALUES  = b_values
          QPalette->GetProperty, RED_VALUES   = r_values
          QPalette->GetProperty, GREEN_VALUES = g_values
          
          b_values    = TRANSPOSE(b_values)
          r_values    = TRANSPOSE(r_values)
          g_values    = TRANSPOSE(g_values)
          
          color_skip  = n_colors / n_cntrs[1]
          rgb_colors  = INTARR(3,n_cntrs[1])
          
          col_idx     = 0
          
          FOR I=0,n_cntrs[1]-1 DO BEGIN
             rgb_colors[0,I] = r_values[col_idx]
             rgb_colors[1,I] = g_values[col_idx]
             rgb_colors[2,I] = b_values[col_idx]
             col_idx         = col_idx + color_skip
          ENDFOR
          
          c_label_orient     = INTARR(1,n_cntrs[1])
          c_label_orient     = 1
          
          c_colors           = IndGen(n_cntrs[1]) + 1
          QContourTitle      = OBJ_NEW('IDLgrText',str_title, LOCATIONS=[0.5, 1.15], ALIGNMENT=0.5)
          QContourXlabel     = OBJ_NEW('IDLgrText', 'X'     , LOCATIONS=[0.5, -0.2], ALIGNMENT=0.5)
          QContourYlabel     = OBJ_NEW('IDLgrText', 'Y'     , LOCATIONS=[-0.2, 0.5], ALIGNMENT=0.5)

          QFContour          = OBJ_NEW('IDLgrContour', Q_sqr,                                   $
                                                      /FILL,                                    $
                                                       C_VALUE                 = Cntrs,         $
                                                      /PLANAR,                                  $
                                                       GEOMZ                   = 0,             $
                                                       GEOMX                   = X,             $
                                                       GEOMY                   = Y,             $
                                                       PALETTE                 = QPalette,      $
                                                       C_COLOR                 = rgb_colors       )
          
          QLContour          = OBJ_NEW('IDLgrContour', Q_sqr,                                   $
                                                      /PLANAR,                                  $
                                                       C_VALUE                 = Cntrs,         $
                                                       GEOMZ                   = 0,             $
                                                       GEOMX                   = X,             $
                                                       GEOMY                   = Y,             $
                                                       PALETTE                 = QPalette,      $
                                                       COLOR                   = 0,             $
                                                       C_LABEL_SHOW            = CLABS,         $
                                                       C_USE_LABEL_ORIENTATION = c_label_orient   )
          
          QFContour->GetProperty, XRANGE=x_rng
          QFContour->GetProperty, YRANGE=y_rng
          
          x_min = x_rng[0]
          x_max = x_rng[1]

          y_min = y_rng[0]
          y_max = y_rng[1]

          QLeftYaxis        = OBJ_NEW('IDLgrAxis', DIRECTION=1, THICK=2, TICKLEN=0.02, TICKDIR=1, TEXTPOS=0, RANGE=y_rng)
          QBotXaxis         = OBJ_NEW('IDLgrAxis', DIRECTION=0, THICK=2, TICKLEN=0.02, TICKDIR=1, TEXTPOS=0, RANGE=x_rng)
          
          QRightYaxis       = OBJ_NEW('IDLgrAxis',  DIRECTION=1, THICK=2, TICKLEN=0.02, TICKDIR=0, TEXTPOS=1, $
                                       RANGE=y_rng, LOCATION=[x_max,y_min] )
          QTopXaxis         = OBJ_NEW('IDLgrAxis',  DIRECTION=0, THICK=2, TICKLEN=0.02, TICKDIR=0, TEXTPOS=1, $
                                       RANGE=x_rng, LOCATION=[x_min,y_max] )
          
          CntrModel         = OBJ_NEW('IDLgrModel')
          CntrModel->Add, QLContour
          CntrModel->Add, QFContour
          CntrModel->Add, QContourTitle
          CntrModel->Add, QContourXlabel
          CntrModel->Add, QContourYlabel
          CntrModel->Add, QBotXaxis
          CntrModel->Add, QLeftYaxis
          CntrModel->Add, QTopXaxis
          CntrModel->Add, QRightYaxis
          CntrModel->Scale, 1.0, 1.0, 1.0
          
          CntrView->Add,  CntrModel

          c_max              = MAX(Cntrs)
          c_min              = MIN(Cntrs)
          
          n_tics             = 5
          inc_c              = (c_max - c_min) / (n_tics-1)
          
          c_tics             = FLTARR(1,n_tics)
          c_tics             = TRANSPOSE(c_tics)
          str_c_tics         = STRARR(1, n_tics)
          str_c_tics         = TRANSPOSE(str_c_tics)
          
          FOR I = 0, n_tics - 1 DO BEGIN
          
             c_tics[I]       = c_min + (I * inc_c)
             
             IF (I EQ (n_tics - 1)/2) THEN BEGIN
               IF (ABS(c_tics[I]/c_max) LT 1.0e-06) THEN c_tics[I] = 0.0
             ENDIF

             str_c_tics[I]   = STRTRIM(c_tics[I],2)                             ; Convert this tic label to a string
             dot_loc         = STRPOS(str_c_tics[I],'.')                        ; find the location of its decimal point
             tic_len         = STRLEN(str_c_tics[I])                            ; get the length of the string
             nz_found        = 0                                                ; set to 1 if non-zero digit is found
             nz_loc          = dot_loc                                          ; begin search for non-zero digit at location 
                                                                                ; of decimal point

             IF (ABS(c_tics[I]) LT 1.0) THEN BEGIN                              ; If the tic label is smaller than 1.0...
               FOR J=0, tic_len - 1 DO BEGIN                                    ; ...then scan to the right for non-zero digit
                 nz_loc      = nz_loc + 1                                       ;  -> start just to the right of the decimal point
                 next_char   = STRMID(str_c_tics[I],nz_loc,1)                   ;  -> get the next character in the string
                 IF (next_char NE '0') THEN BEGIN                               ;  -> if it's not a zero find out what it is
                   digits    = ['1','2','3','4','5','6','7','8','9']            ;  -> we're looking for one of these
                   char_scan = STRMATCH(digits,next_char)                       ;  -> if next_char is a match for any of the above 
                                                                                ;     this array will have a non-zero component
                   nz_iprod  = TRANSPOSE(char_scan) # char_scan                 ;  -> and this dot product will be non-zero
                   IF (nz_iprod GE 1) THEN nz_found = 1                         ;  -> and therefore a non-zero digit has been found
                   BREAK                                                        ;  -> so we've done what we've needed to do
                ENDIF
               ENDFOR
               IF (nz_found EQ 1) THEN BEGIN                                    ; If we found a non-zero digit
                   dot_shift = nz_loc - dot_loc                                 ;  -> figure out how many negative powers 
                                                                                ;     of ten we need
                   ten_pow   = 'e' + STRTRIM(-dot_shift,2)                      ;  -> and construct the appropriate 
                                                                                ;     base/exponent factor 
                   mantbeg   = STRMID(str_c_tics[I],nz_loc,1)                   ;  -> First digit of the mantissa 
                   IF (c_tics[I] LT 0.0) THEN mantbeg = '-' + mantbeg           ;  -> tack on either a '+' sign...
                   IF (c_tics[I] GT 0.0) THEN mantbeg = '+' + mantbeg           ;  -> .. or a minus sign
                   mantend   = STRMID(str_c_tics[I],nz_loc+1, tic_len - nz_loc) ;  -> the part of the mantissa to the right 
                                                                                ;     of the decimal
                   mantissa  = mantbeg + '.' + mantend                          ;  -> construct the mantissa
                   mantissa  = STRMID(mantissa,0,6)                             ;  -> truncate to desired length - will want 
                                                                                ;     to make this more flexible
                   str_c_tics[I] = mantissa + ten_pow                           ;  -> now combine mantissa and base/exponent 
                                                                                ;     factor for complete tic string
               ENDIF ELSE BEGIN
                   str_c_tics[I]   = STRMID(STRTRIM(c_tics[I],2),0,5)           ; If nz_found is 0 then the value of c_tic[I] 
                                                                                ; must be zero. So...
                     ENDELSE                                                    ; truncate to length one less than non-zero 
                                                                                ; values (no need for '+' or '-')
             ENDIF ELSE BEGIN

                   print, 'construct_jcontour: dot_loc = ', dot_loc
                   IF (c_tics[I] GE 0) THEN BEGIN                               ; use is_neg to keep track of any minus signs
                     is_neg = 0                                                 ; (see usage below)
                   ENDIF ELSE BEGIN
                     is_neg = 1
                         ENDELSE
                   dot_shift = dot_loc - (1 + is_neg)                           ; powers of ten

                   mantbeg  = STRMID(str_c_tics[I],0, 1 + is_neg)               ;
                   mantmid  = STRMID(str_c_tics[I],is_neg + 1, dot_shift)       ;
                   mantend  = STRMID(str_c_tics[I],dot_loc+1,tic_len-(dot_loc+1))
                   mantissa = mantbeg + '.' + mantmid + mantend

                   IF (is_neg EQ 0) THEN BEGIN
                      IF (c_tics[I] NE 0) THEN mantissa = '+' + mantissa        ; -> For symmetric appearance put '+' on 
                                                                                ;    positive mantissas
                   ENDIF
                   mantissa = STRMID(mantissa,0,6)                              ; -> Trim the result

                   print, 'mantbeg  = ', mantbeg
                   print, 'mantmid  = ', mantmid
                   print, 'mantend  = ', mantend
                   print, 'mantissa = ', mantissa

                   IF (dot_shift GT 0) THEN BEGIN
                     ten_pow = 'e+' + STRTRIM(dot_shift , 2)
                     str_c_tics[I] = mantissa + ten_pow
                   ENDIF ELSE BEGIN
                     str_c_tics[I] = mantissa
                         ENDELSE
                   ENDELSE
          ENDFOR
          
          QtickLabels          = OBJ_NEW('IDLgrText', str_c_tics)
          barDims              = [0.05, 1.0 ]

          CASE qty OF
          'vz': Q_cb_title_str = 'V_z'
          'bz': Q_cb_title_str = 'Z'
          'a' : Q_cb_title_str = 'A'
          'p' : Q_cb_title_str = 'Stream Function Contours '
          ELSE: str_title = 'default title string'
          ENDCASE

          QColorBarTitle       = OBJ_NEW('IDLgrText', Q_cb_title_str)
          QColorBar            = OBJ_NEW('IDLgrColorbar', r_values,                         $
                                                          g_values,                         $
                                                          b_values,                         $
                                                           TITLE         = QColorBarTitle,  $
                                                           DIMENSIONS    = barDims,         $
                                                           PALETTE       = QPalette,        $
                                                           SHOW_AXIS     = 2,               $
                                                           MAJOR         = 5,               $
                                                           MINOR         = 5,               $
                                                           TICKTEXT      = QtickLabels,     $
                                                          /SHOW_OUTLINE                       )
          
          barModel             = OBJ_NEW('IDLgrModel')
          barModel->Scale, 1.0, 1.0, 1.0
          barView->Add, barModel
          barModel->Add, QColorBar
          barPlusTextDims=QColorBar->ComputeDimensions(QBuffer)
          QBuffer->Draw, QScene
          QBuffer->SetProperty, GRAPHICS_TREE = QScene
          QBuffer->GetProperty, DIMENSIONS    = windowSize
          QBuffer->GetProperty, RESOLUTION    = screenResolution
          QClipboard=OBJ_NEW('IDLgrClipboard',  QUALITY       = 2,                 $
                                                DIMENSIONS    = windowSize,        $
                                                RESOLUTION    = screenResolution,  $
                                                GRAPHICS_TREE = QScene)
           
          eps_prefix                             = qty + '_contour'
          eps_infix                              = '-ff-spec'
          postfile                               = construct_eps_outfile_name(eps_prefix, eps_infix, n_slice, n_step, tot_steps) 

          QClipboard->Draw,     FILENAME         = postfile, /POSTSCRIPT, /VECTOR 
          QBuffer->SetProperty, QUALITY          = 2
          QBuffer->SetProperty, PALETTE          = QPalette
          QBuffer->SetProperty, GRAPHICS_TREE    = QScene
          QBuffer->GetProperty, ZBUFFER_DATA     = Pic 
          
          QClipboard->SetProperty, GRAPHICS_TREE = OBJ_NEW()
          OBJ_DESTROY, QClipboard

ENDIF ELSE BEGIN

          J_z             = calc_jz(Slc, n_slice, n_step, desc_label )
          i_qty           = '1'
          str_title       = 'A (contours) and J (color map): t = ' + str_time

          J_Sqr           = FLTARR(i_x_res + 1, i_y_res + 1)
          A_Sqr           = FLTARR(i_x_res + 1, i_y_res + 1)

          FOR I=0, i_x_res - 1 DO BEGIN
             FOR J=0, i_y_res - 1 DO BEGIN
          
                idx         = (I * i_x_res) + J
                A_Sqr[I, J] = A[idx, i_qty]
                J_Sqr[I, J] = J_z[idx]
          
             ENDFOR
          ENDFOR

          A_sqr[*,       i_y_res ] = A_sqr[*, 0 ]
          A_sqr[i_x_res, *       ] = A_sqr[0, * ]

          J_sqr[*,       i_y_res ] = J_sqr[*, 0 ]
          J_sqr[i_x_res, *       ] = J_sqr[0, * ]

          myBuffer    = OBJ_NEW('IDLgrBuffer')
          
          myScene     = OBJ_NEW('IDLgrScene')
          CntrView    = OBJ_NEW('IDLgrView', LOCATION=[0.125, 0.0], VIEWPLANE_RECT=[-0.25, -0.25, 2.0, 1.5], $
                                             DIMENSIONS=[0.9, 0.9], UNITS=3, COLOR=[255,255,255])
                                             
          BarView     = OBJ_NEW('IDLgrView', LOCATION=[0.75, 0.0 ], VIEWPLANE_RECT=[-0.25, -0.25, 2.0, 1.5], $
                                             DIMENSIONS=[0.75,0.9], UNITS=3, COLOR=[255,255,255])
          myScene->Add, CntrView
          myScene->Add, barView
          
          myPalette   = OBJ_NEW('IDLgrPalette')
          myPalette->LoadCT, 33
          
          myPalette->GetProperty, N_COLORS     = n_colors
          myPalette->GetProperty, BLUE_VALUES  = b_values
          myPalette->GetProperty, RED_VALUES   = r_values
          myPalette->GetProperty, GREEN_VALUES = g_values
          
          b_values    = TRANSPOSE(b_values)
          r_values    = TRANSPOSE(r_values)
          g_values    = TRANSPOSE(g_values)
          
          color_skip  = n_colors / n_cntrs[1]
          rgb_colors  = INTARR(3,n_cntrs[1])
          
          col_idx     = 0
          
          print, "construct_jcontours: n_colors   = ", n_colors
          print, "construct_jcontours: color_skip = ", color_skip

          FOR I=0,n_cntrs[1]-1 DO BEGIN
             rgb_colors[0,I] = r_values[col_idx]
             rgb_colors[1,I] = g_values[col_idx]
             rgb_colors[2,I] = b_values[col_idx]
             col_idx         = col_idx + color_skip
          ENDFOR
          
          c_label_orient     = INTARR(1,n_cntrs[1])
          c_label_orient     = 1
          
          c_colors           = IndGen(n_cntrs[1]) + 1
          myContourTitle     = OBJ_NEW('IDLgrText',str_title, LOCATIONS=[0.5, 1.15], ALIGNMENT=0.5)
          myContourXlabel    = OBJ_NEW('IDLgrText', 'X'     , LOCATIONS=[0.5, -0.2], ALIGNMENT=0.5)
          myContourYlabel    = OBJ_NEW('IDLgrText', 'Y'     , LOCATIONS=[-0.2, 0.5], ALIGNMENT=0.5)
          
          myContour          = OBJ_NEW('IDLgrContour', J_sqr,                                   $
                                                      /FILL,                                    $
                                                       C_VALUE                 = J_Cntrs,       $
                                                      /PLANAR,                                  $
                                                       GEOMZ                   = 0,             $
                                                       GEOMX                   = X,             $
                                                       GEOMY                   = Y,             $
                                                       PALETTE                 = myPalette,     $
                                                       C_COLOR                 = rgb_colors       )
          
          AContour           = OBJ_NEW('IDLgrContour', A_sqr,                                   $
                                                      /PLANAR,                                  $
                                                       C_VALUE                 = Cntrs,         $
                                                       GEOMZ                   = 0,             $
                                                       GEOMX                   = X,             $
                                                       GEOMY                   = Y,             $
                                                       PALETTE                 = myPalette,     $
                                                       COLOR                   = 0,             $
                                                       C_LABEL_SHOW            = CLABS,         $
                                                       C_USE_LABEL_ORIENTATION = c_label_orient   )

          myContour->GetProperty, XRANGE=x_rng
          myContour->GetProperty, YRANGE=y_rng

          x_min = x_rng[0]
          x_max = x_rng[1]

          y_min = y_rng[0]
          y_max = y_rng[1]
          
          myLeftYaxis        = OBJ_NEW('IDLgrAxis', DIRECTION=1, THICK=2, TICKLEN=0.02, TICKDIR=1, TEXTPOS=0, RANGE=y_rng)
          myBotXaxis         = OBJ_NEW('IDLgrAxis', DIRECTION=0, THICK=2, TICKLEN=0.02, TICKDIR=1, TEXTPOS=0, RANGE=y_rng)
          
          myRightYaxis       = OBJ_NEW('IDLgrAxis', DIRECTION=1, THICK=2, TICKLEN=0.02, TICKDIR=0, TEXTPOS=1, $
                                        RANGE=y_rng, LOCATION=[x_max,y_min])
          myTopXaxis         = OBJ_NEW('IDLgrAxis', DIRECTION=0, THICK=2, TICKLEN=0.02, TICKDIR=0, TEXTPOS=1, $
                                        RANGE=x_rng, LOCATION=[x_min,y_max])
          
          CntrModel          = OBJ_NEW('IDLgrModel')
          CntrModel->Add, AContour
          CntrModel->Add, myContour
          CntrModel->Add, myContourTitle
          CntrModel->Add, myContourXlabel
          CntrModel->Add, myContourYlabel
          CntrModel->Add, myBotXaxis
          CntrModel->Add, myLeftYaxis
          CntrModel->Add, myTopXaxis
          CntrModel->Add, myRightYaxis
          CntrModel->Scale, 1.0, 1.0, 1.0
          
          CntrView->Add,  CntrModel
          
          j_max              = MAX(J_Cntrs)
          j_min              = MIN(J_Cntrs)
          n_tics             = 5
          inc_j              = (j_max - j_min) / (n_tics-1)
          
          j_tics             = FLTARR(1,n_tics)
          j_tics             = TRANSPOSE(j_tics)
          str_j_tics         = STRARR(1, n_tics)
          str_j_tics         = TRANSPOSE(str_j_tics)
          
          FOR I = 0, n_tics - 1 DO BEGIN
          
             j_tics[I]       = j_min + (I * inc_j)
             
             IF (I EQ (n_tics - 1)/2) THEN BEGIN
               IF (ABS(j_tics[I]/j_max) LT 1.0e-06) THEN j_tics[I] = 0.0
             ENDIF

             str_j_tics[I]   = STRTRIM(j_tics[I],2)                             ; Convert this tic label to a string
             dot_loc         = STRPOS(str_j_tics[I],'.')                        ; find the location of its decimal point
             tic_len         = STRLEN(str_j_tics[I])                            ; get the length of the string
             nz_found        = 0                                                ; set to 1 if non-zero digit is found
             nz_loc          = dot_loc                                          ; begin search for non-zero digit at location 
                                                                                ; of decimal point

             IF (ABS(j_tics[I]) LT 1.0) THEN BEGIN                              ; If the tic label is smaller than 1.0...
               FOR K=0, tic_len - 1 DO BEGIN                                    ; ...then scan to the right for non-zero digit
                 nz_loc      = nz_loc + 1                                       ;  -> start just to the right of the decimal point
                 next_char   = STRMID(str_j_tics[I],nz_loc,1)                   ;  -> get the next character in the string
                 IF (next_char NE '0') THEN BEGIN                               ;  -> if it's not a zero find out what it is
                   digits    = ['1','2','3','4','5','6','7','8','9']            ;  -> we're looking for one of these
                   char_scan = STRMATCH(digits,next_char)                       ;  -> if next_char is a match for any of the above 
                                                                                ;     this array will have a non-zero component
                   nz_iprod  = TRANSPOSE(char_scan) # char_scan                 ;  -> and this dot product will be non-zero
                   IF (nz_iprod GE 1) THEN nz_found = 1                         ;  -> and therefore a non-zero digit has been found
                   BREAK                                                        ;  -> so we've done what we've needed to do
                ENDIF
               ENDFOR
               IF (nz_found EQ 1) THEN BEGIN                                    ; If we found a non-zero digit
                   dot_shift = nz_loc - dot_loc                                 ;  -> figure out how many negative powers 
                                                                                ;     of ten we need
                   ten_pow   = 'e' + STRTRIM(-dot_shift,2)                      ;  -> and construct the appropriate 
                                                                                ;     base/exponent factor 
                   mantbeg   = STRMID(str_j_tics[I],nz_loc,1)                   ;  -> First digit of the mantissa 
                   IF (j_tics[I] LT 0.0) THEN mantbeg = '-' + mantbeg           ;  -> tack on either a '+' sign...
                   IF (j_tics[I] GT 0.0) THEN mantbeg = '+' + mantbeg           ;  -> .. or a minus sign
                   mantend   = STRMID(str_j_tics[I],nz_loc+1, tic_len - nz_loc) ;  -> the part of the mantissa to the right 
                                                                                ;     of the decimal
                   mantissa  = mantbeg + '.' + mantend                          ;  -> construct the mantissa
                   mantissa  = STRMID(mantissa,0,6)                             ;  -> truncate to desired length - will want 
                                                                                ;     to make this more flexible
                   str_j_tics[I] = mantissa + ten_pow                           ;  -> now combine mantissa and base/exponent 
                                                                                ;     factor for complete tic string
               ENDIF ELSE BEGIN
                   str_j_tics[I]   = STRMID(STRTRIM(j_tics[I],2),0,5)           ; If nz_found is 0 then the value of c_tic[I] 
                                                                                ; must be zero. So...
                     ENDELSE                                                    ; truncate to length one less than non-zero 
                                                                                ; values (no need for '+' or '-')
             ENDIF ELSE BEGIN

                   print, 'construct_jcontour: dot_loc = ', dot_loc
                   IF (j_tics[I] GE 0) THEN BEGIN                               ; use is_neg to keep track of any minus signs
                     is_neg = 0                                                 ; (see usage below)
                   ENDIF ELSE BEGIN
                     is_neg = 1
                         ENDELSE
                   dot_shift = dot_loc - (1 + is_neg)                           ; powers of ten

                   mantbeg  = STRMID(str_j_tics[I],0, 1 + is_neg)               ;
                   mantmid  = STRMID(str_j_tics[I],is_neg + 1, dot_shift)       ;
                   mantend  = STRMID(str_j_tics[I],dot_loc+1,tic_len-(dot_loc+1))
                   mantissa = mantbeg + '.' + mantmid + mantend

                   IF (is_neg EQ 0) THEN BEGIN
                      IF (j_tics[I] NE 0) THEN mantissa = '+' + mantissa        ; -> For symmetric appearance put '+' on 
                                                                                ;    positive mantissas
                   ENDIF
                   mantissa = STRMID(mantissa,0,6)                              ; -> Trim the result

                   print, 'mantbeg  = ', mantbeg
                   print, 'mantmid  = ', mantmid
                   print, 'mantend  = ', mantend
                   print, 'mantissa = ', mantissa

                   IF (dot_shift GT 0) THEN BEGIN
                     ten_pow = 'e+' + STRTRIM(dot_shift , 2)
                     str_j_tics[I] = mantissa + ten_pow
                   ENDIF ELSE BEGIN
                     str_j_tics[I] = mantissa
                         ENDELSE
                   ENDELSE
          ENDFOR
          
          mytickLabels       = OBJ_NEW('IDLgrText', str_j_tics)
          barDims            = [0.05, 1.0 ]
          
          myColorBarTitle    = OBJ_NEW('IDLgrText', 'Z-Component of Current Density (J) ')
          myColorBar         = OBJ_NEW('IDLgrColorbar', r_values,                        $
                                                        g_values,                        $
                                                        b_values,                        $
                                                        TITLE         = myColorBarTitle, $
                                                        DIMENSIONS    = barDims,         $
                                                        PALETTE       = myPalette,       $
                                                        SHOW_AXIS     = 2,               $
                                                        MAJOR         = 5,               $
                                                        MINOR         = 5,               $
                                                        TICKTEXT      = mytickLabels,    $
                                                       /SHOW_OUTLINE                       )
          
          barModel           = OBJ_NEW('IDLgrModel')
          barModel->Scale, 1.0, 1.0, 1.0
          barView->Add, barModel
          barModel->Add, myColorBar
          barPlusTextDims=myColorBar->ComputeDimensions(myBuffer)
          myBuffer->Draw, myScene
          myBuffer->SetProperty, GRAPHICS_TREE=myScene
          myBuffer->GetProperty, DIMENSIONS=windowSize
          myBuffer->GetProperty, RESOLUTION=screenResolution
          myClipboard=OBJ_NEW('IDLgrClipboard', QUALITY       = 2,                $
                                                DIMENSIONS    = windowSize,       $
                                                RESOLUTION    = screenResolution, $
                                                GRAPHICS_TREE = myScene)
           
          
          eps_prefix       = qty + '_contour'
          eps_infix        = '-ff-spec'
          postfile         = construct_eps_outfile_name(eps_prefix, eps_infix, n_slice, n_step, tot_steps) 
          myClipboard->Draw, FILENAME = postfile, /POSTSCRIPT, /VECTOR 
          
          myBuffer->SetProperty, QUALITY        = 2
          myBuffer->SetProperty, PALETTE        = myPalette
          myBuffer->SetProperty, GRAPHICS_TREE  = myScene
          myBuffer->GetProperty, ZBUFFER_DATA   = Pic 
          
          myClipboard->SetProperty, GRAPHICS_TREE = OBJ_NEW()
          OBJ_DESTROY, myClipboard

ENDELSE

RETURN, Pic
END
