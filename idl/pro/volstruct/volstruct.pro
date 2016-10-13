PRO volstruct, qty,                        $
               desc_label,                 $
               type,                       $
               action,                     $
               global_in_minmax,           $
               PREFIX      = prefix,       $
               PLS_MIN     = pls_min,      $
               SLICE_RANGE = slice_range,  $
               THRESHOLD   = threshold,    $
               WIDTH       = width,        $
               N_CNTRS     = n_cntrs,      $
               INC_CONT    = inc_cont,     $
               FIRST_STEP  = first_step,   $
               LAST_STEP   = last_step,    $
               N_STEP      = n_step,       $
               YROT        = yrot,         $
               START_ANGLE = start_angle,  $
               STOP_ANGLE  = stop_angle,   $
               N_FRAMES    = n_frames,     $
               GAP         = gap,          $
               WINDOW_SIZE = window_size

ZERO                           = 0.0D

PRINT, 'volstruct: first_step = ', first_step
PRINT, 'volstruct: last_step  = ', last_step
PRINT, 'volstruct: width      = ', width

; note need a default for yrot and slice_range

; type        - volume with threshold
;             - n_cntrs isosurfaces with width
;action       ~ time-series 'tseries' (fixed yrot, fixed threshold or width)
;             ~ yrotate     'rotate'  (fixed step, fixed threshold or width)
;             ~ scan        'scan'    (fixed step, fixed yrot, fixed width, n_cntrs = 1, varying contour value) 
;features     ~ gap         (specified by z values)
;             ~ contour     (location at far end of gap, or specified)

;req'd by all ~ qty
;             ~ pm (plus, minus, or plus and minus)
;             ~ desc_label
;             ~ type
;             ~ action
;             ~ first and last slice default value or as specified by slice_range
;             ~ prefix
;             ~ global_in_minmax

;req'd by type:
;             ~ vol -> threshold
;             ~ isf -> width
;                   -> n_cntrs
;req'd by action:
;             ~ time series -> first_step, last_step
;             ~ yrotate     -> step
;                           -> start_angle, stop_angle
;                           -> number of plots (n_frames)
;             ~ scan        -> step
;                           -> number of plots (n_frames)
;req'd by feature:
;             ~ gap         -> beginning and ending z-coordinates for gap
;             ~ contour     -> layer to be contoured
;                           -> quantity to be contoured
;                           -> number of contours

  !EXCEPT                      = 2

  ZERO                    = 0.0D

  IF (NOT KEYWORD_SET(N_STEP) ) THEN BEGIN
    IF (NOT (KEYWORD_SET(FIRST_STEP) AND KEYWORD_SET(LAST_STEP)) ) THEN BEGIN
      PRINT, 'volstruct: ERROR - no step specified, must set either n_step or first_step and last_step'
    ENDIF
  ENDIF ELSE BEGIN
    IF ((KEYWORD_SET(FIRST_STEP) OR KEYWORD_SET(LAST_STEP)) ) THEN BEGIN
      PRINT, 'volstruct: WARNING - first_step and/or last_step specified with n_step also specified, !C defaulting to single time step.'
    ENDIF
    first_step                 = n_step
    last_step                  = n_step
  ENDELSE

  IF (type EQ 'vol') THEN BEGIN
    IF (NOT KEYWORD_SET(THRESHOLD) ) THEN BEGIN
      PRINT, 'volstruct: WARNING - threshold for volume plot not specified, using default value of 0.0'
      threshold                = ZERO
    ENDIF
    IF (KEYWORD_SET(WIDTH)) THEN BEGIN
      PRINT, 'volstruct: WARNING - ignoring specification of width for volume plot'
      width                    = 0.1D
    ENDIF
  ENDIF

  IF (type EQ 'isf') THEN BEGIN
    IF (KEYWORD_SET(THRESHOLD) AND threshold GT ZERO ) THEN BEGIN
      PRINT, 'volstruct: WARNING - ignoring non-zero value of threshold in isosurface plot'
      threshold                = ZERO
    ENDIF
    IF (NOT KEYWORD_SET(WIDTH)) THEN BEGIN
      PRINT, 'volstruct: WARNING - isosurface width not specified, using default value of 0.1'
      width                    = 0.1D
    ENDIF
  ENDIF

  IF (NOT KEYWORD_SET(PLS_MIN)) THEN pls_min = 'pm'
  IF (NOT KEYWORD_SET(PREFIX) ) THEN prefix  = 'rmct2'

  IF (KEYWORD_SET(SLICE_RANGE)) THEN BEGIN

    first_slice                = slice_range[0]
    last_slice                 = slice_range[1]

  ENDIF ELSE BEGIN

    n3                         = scan_parameters('n3',  0, desc_label )
    mp                         = scan_parameters('mp',  0, desc_label )

    first_slice                = 1
    last_slice                 = n3*mp

  ENDELSE

  ip1                          = scan_parameters('ip1', 0, desc_label )
  ip2                          = scan_parameters('ip2', 0, desc_label )
  n1                           = 2^ip1
  n2                           = 2^ip2
  KK                           = calcKK(n1,n2)   
  
  IF (STRLEN(global_in_minmax) > 0 ) THEN BEGIN

    cur_dir                    = GETENV('PWD')
    in_file                    = cur_dir + '/' + 'glb_ext.out'

    OPENR, in_unit, in_file, /GET_LUN, ERROR = op_err

    IF ( op_err EQ 0 ) THEN BEGIN
    
      PRINT, 'glb_ext: reading previous output from ', in_file

      global_in_minmax         = FLTARR(4,4)
      line                     = FLTARR(4)

      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,0]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,1]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,2]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,3]                     =    line
    
      CLOSE, in_unit

      IF (qty EQ 'j' ) THEN global_j_minmax = global_in_minmax[*,3]
      IF (qty EQ 'o' ) THEN global_j_minmax = global_in_minmax[*,2]
      IF (qty EQ 'a' ) THEN global_j_minmax = global_in_minmax[*,1]
      IF (qty EQ 'p' ) THEN global_j_minmax = global_in_minmax[*,0]

    ENDIF ELSE BEGIN

      PRINT, 'glb_ext: WARNING - could not open file ', in_file, '. Assuming no global minmax given...'

      global_j_minmax          = global_q_minmax(qty, first_step, last_step, first_slice, last_slice, desc_label, KK)

    ENDELSE

  ENDIF ELSE BEGIN

    PRINT, 'no global minmax given...'
    global_j_minmax            = global_q_minmax(qty, first_step, last_step, first_slice, last_slice, desc_label, KK)

  ENDELSE

  IF (action EQ 'tseries') THEN BEGIN

    IF (pls_min EQ 'p' || pls_min EQ 'pm') THEN JVOX_P = init_jvox(first_slice, last_slice) ;  create BYTE-array space for positive J
    IF (pls_min EQ 'm' || pls_min EQ 'pm') THEN JVOX_M = init_jvox(first_slice, last_slice) ;  create BYTE-array space for negative J

    IF (pls_min EQ 'p' || pls_min EQ 'pm') THEN BEGIN
      size_jvox                  = SIZE(JVOX_P, /DIMENSIONS) ; upper bnds provides xyz-range
    ENDIF ELSE BEGIN
      size_jvox                  = SIZE(JVOX_M, /DIMENSIONS) ; upper bnds provides xyz-range
    ENDELSE

    cc                           = getCoordConvs(desc_label, size_jvox) ; need before axes

    jVolAxes                     = makeAxes(cc, size_jvox)

    coltabs                      = getColorTables(qty)
    coltabP                      = coltabs[*,*,0]
    coltabM                      = coltabs[*,*,1]

    IF (qty EQ 'j') THEN str_pfx = 'jstruct'
    IF (qty EQ 'o') THEN str_pfx = 'ostruct'
    IF (qty EQ 'a') THEN str_pfx = 'astruct'
    IF (qty EQ 'p') THEN str_pfx = 'pstruct'


    IF (type EQ 'vol' ) THEN str_ifx = getInfixString(desc_label, type, last_step, THRESHOLD = threshold )
    IF (type EQ 'isf' ) THEN str_ifx = getInfixString(desc_label, type, last_step, WIDTH     = width     )

    FOR I = first_step, last_step DO BEGIN

      output_path              = getOutputPath(I, str_pfx, str_ifx, last_step)

      IF (type EQ 'vol' ) THEN BEGIN

        PRINT, 'jstruct: making volume plot...'

        dummy                  = fill_jvox(desc_label, I, threshold, first_slice, last_slice, global_j_minmax, pls_min, JVOX_P, JVOX_M, KK)

      ENDIF ELSE BEGIN

        IF (type EQ 'isf' ) THEN  BEGIN

          PRINT, 'jstruct: making isosurfaces...'

          dummy                = fill_jvox(desc_label, I, threshold, first_slice, last_slice, global_j_minmax, pls_min, JVOX_P, JVOX_M, KK)

          IF (pls_min EQ 'p' || pls_min EQ 'pm') THEN BEGIN

            APContour          = makeContourPlot(n_cntrs, 384, I, desc_label, width, global_j_minmax, JVOX_P, cc, 'p')
            APContour          -> GetProperty, C_VALUE = CNTRS_P

            MASK_P             = jMask(CNTRS_P, width, JVOX_P)
            JVOX_P             = MASK_P
          
          ENDIF

          IF (pls_min EQ 'm' || pls_min EQ 'pm') THEN BEGIN

            AMContour          = makeContourPlot(n_cntrs, 384, I, desc_label, width, global_j_minmax, JVOX_M, cc, 'm')
            AMContour          -> GetProperty, C_VALUE = CNTRS_M

            MASK_M             = jMask(CNTRS_M, width, JVOX_M)
            JVOX_M             = MASK_M

          ENDIF

        ENDIF
      ENDELSE

      IF (KEYWORD_SET(WINDOW_SIZE)) THEN BEGIN
        windowSize          = window_size
      ENDIF ELSE BEGIN
        windowSize          = [800, 600]
      ENDELSE

      JBuffer               =  OBJ_NEW('IDLgrBuffer', DIMENSIONS = windowSize)
      JViewgroup            =  OBJ_NEW('IDLgrViewgroup'                      )

      JBuffer               -> SetProperty, GRAPHICS_TREE = JViewgroup

      JView                 =  OBJ_NEW('IDLgrView',                              $
                                        PROJECTION     =   2,                    $
                                        VIEWPLANE_RECT = [-2.0, -2.0, 7.5, 7.5], $
                                        ZCLIP          = [4.9, -4.9],            $
                                        EYE            =  5.0                    $
                                      )

      JViewgroup            -> Add, JView
      JModel                =  OBJ_NEW('IDLgrModel', /LIGHTING )
      JView                 -> Add, JModel

      IF (pls_min EQ 'p' || pls_min EQ 'pm' ) THEN JVolumeP = makeVolume(JVOX_P, coltabP, cc)
      IF (pls_min EQ 'm' || pls_min EQ 'pm' ) THEN JVolumeM = makeVolume(JVOX_M, coltabM, cc)

      str_title             = getTitle(qty, type, action, desc_label, I, threshold, width, prefix, global_j_minmax)

      JTTLfont              = OBJ_NEW('IDLgrFont', 'helvetica*bold', SIZE=16)
      JTitle                = OBJ_NEW('IDLgrText',              $
                                       str_title,               $
                                       LOCATIONS = [0.0, 1.15], $
                                       ALIGNMENT = 0.0,         $
                                       FONT      = JTTLfont,    $
                                      /ENABLE_FORMATTING)

      JModel                -> Add, JTitle

      JzTitle               = OBJ_NEW('IDLgrText', 'Z')
      JxTitle               = OBJ_NEW('IDLgrText', 'X')
      JyTitle               = OBJ_NEW('IDLgrText', 'Y', UPDIR=[0,0,-1], BASELINE=[0,1,0])

      IF (pls_min EQ 'p' || pls_min EQ 'pm') THEN JModel -> Add, JVolumeP
      IF (pls_min EQ 'm' || pls_min EQ 'pm') THEN JModel -> Add, JVolumeM

      JModel                -> Add, jVolAxes[0,1]
      JModel                -> Add, jVolAxes[1,1]
      JModel                -> Add, jVolAxes[2,1]
      JModel                -> Add, jVolAxes[3,1]

      JModel                -> Add, jVolAxes[0,2]
      JModel                -> Add, jVolAxes[1,2]
      JModel                -> Add, jVolAxes[2,2]
      JModel                -> Add, jVolAxes[3,2]

      JModel                -> Add, jVolAxes[0,0]
      JModel                -> Add, jVolAxes[1,0]
      JModel                -> Add, jVolAxes[2,0]
      JModel                -> Add, jVolAxes[3,0]

      IF (type EQ 'isf' AND KEYWORD_SET(inc_cont) ) THEN BEGIN

        IF (pls_min EQ 'p' || pls_min EQ 'pm' ) THEN JModel -> Add, APContour
        IF (pls_min EQ 'm' || pls_min EQ 'pm' ) THEN JModel -> Add, AMContour

      ENDIF

      sx = 1.8 & sy = 1.8 & sz = 1.8

      JModel                -> Scale, sx, sy, sz
      JModel                -> rotate, [0,0,1],  00
      JModel                -> rotate, [0,1,0],  yrot
      JModel                -> rotate, [1,0,0],  30

      JBuffer               -> GetProperty, Resolution = screenResolution
      JBuffer               -> GetProperty, Dimensions = windowSize
      JClipboard            = OBJ_NEW('IDLgrClipboard',                         $
                                       QUALITY         =  2,                    $
                                       DIMENSIONS      =  2.0*windowSize,       $
                                       RESOLUTION      =  0.5*screenResolution, $
                                       GRAPHICS_TREE   =      JViewgroup)

      JClipboard            -> Draw, VECTOR=0, POSTSCRIPT=1, FILENAME = output_path
      JClipboard            -> SetProperty, GRAPHICS_TREE = OBJ_NEW()

      OBJ_DESTROY, JClipboard
      OBJ_DESTROY, JBuffer

    ENDFOR
  ENDIF
  
END
