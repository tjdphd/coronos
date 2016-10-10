PRO jstruct, qty,                 $
             desc_label,          $
             action,              $
             first_step,          $
             last_step,           $
             first_slice,         $
             last_slice,          $
             threshold,           $
             prefix,              $
             yrot,                $
             n_cntrs,             $
             global_in_minmax

  !EXCEPT                 = 2

  ZERO                    = 0.0D

  JVOX_P                  = init_jvox(first_slice, last_slice) ;  create BYTE-array space for positive J
  JVOX_M                  = init_jvox(first_slice, last_slice) ;  create BYTE-array space for negative J

  str_ifx                 = getInfixString(desc_label, action, last_step, threshold)
  IF (qty EQ 'j') THEN BEGIN
    str_pfx               = 'jstruct'
  ENDIF ELSE BEGIN
    IF (qty EQ 'o') THEN BEGIN
    str_pfx               = 'ostruct'
    ENDIF
  ENDELSE
  out_dir                 = GETENV('PWD') + '/' + str_pfx + '/eps'     ;  prepare an output directory

  ip1                     = scan_parameters('ip1', 0, desc_label )
  ip2                     = scan_parameters('ip2', 0, desc_label )
  n1                      = 2^ip1
  n2                      = 2^ip2
  KK                      = calcKK(n1,n2)   

  IF (STRLEN(global_in_minmax) > 0 ) THEN BEGIN

    cur_dir               = GETENV('PWD')
    in_file               = cur_dir + '/' + 'glb_ext.out'

    OPENR, in_unit, in_file, /GET_LUN, ERROR = op_err


    IF ( op_err EQ 0 ) THEN BEGIN
    
      PRINT, 'jstruct: reading previous output from ', in_file

      global_in_minmax    = FLTARR(4,4)
      line                = FLTARR(4)
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,0]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,1]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,2]                     =    line
      READF, in_unit, FORMAT = '(4(e24.16,1x),:)',   line
      global_in_minmax[*,3]                     =    line
    
      CLOSE, in_unit
      IF (qty EQ 'j' ) THEN BEGIN
        global_j_minmax     = global_in_minmax[*,3]
      ENDIF ELSE BEGIN
        IF (qty EQ 'o') THEN BEGIN
          global_j_minmax     = global_in_minmax[*,2]
        ENDIF
      ENDELSE

    ENDIF ELSE BEGIN

        PRINT, 'glb_ext: WARNING - could not open file ', in_file, '. Assuming no global minmax given...'
        global_j_minmax   = global_q_minmax(qty, first_step, last_step, first_slice, last_slice, desc_label)

    ENDELSE

  ENDIF ELSE BEGIN

    PRINT, 'no global minmax given...'
    global_j_minmax       = global_q_minmax(qty, first_step, last_step, first_slice, last_slice, desc_label, KK)

  ENDELSE

  size_jvox               = SIZE(JVOX_P, /DIMENSIONS)            ; upper bnds provides xyz-range

  cc                      = getCoordConvs(desc_label, size_jvox) ; need before axes

; jVolAxes                = makeAxes(cc, size_jvox)

  coltabs                 = getColorTables(qty)
  coltabP                 = coltabs[*,*,0]
  coltabM                 = coltabs[*,*,1]

  FOR I = first_step, last_step DO BEGIN

    output_path           = getOutputPath(I, str_pfx, str_ifx, last_step)

    IF (action EQ 'vol' ) THEN  BEGIN

      PRINT, 'jstruct: making volume plot...'

      dummy               = fill_jvox(desc_label, I, threshold, first_slice, last_slice, global_j_minmax, JVOX_P, JVOX_M)

    ENDIF ELSE BEGIN

      IF (action EQ 'isf' ) THEN  BEGIN

      PRINT, 'jstruct: making isosurfaces...'

      width               = threshold

      dummy               = fill_jvox(desc_label, I, ZERO, first_slice, last_slice, global_j_minmax, JVOX_P, JVOX_M)

      APContour           = makeContourPlot(n_cntrs, 384, I, desc_label, width, global_j_minmax, JVOX_P, cc, 'p')
      APContour           -> GetProperty, C_VALUE = CNTRS_P

      AMContour           = makeContourPlot(n_cntrs, 384, I, desc_label, width, global_j_minmax, JVOX_M, cc, 'm')
      AMContour           -> GetProperty, C_VALUE = CNTRS_M

      MASK_P              = jMask(CNTRS_P, width, JVOX_P)
      JVOX_P              = MASK_P

      MASK_M              = jMask(CNTRS_M, width, JVOX_M)
      JVOX_M              = MASK_M

      size_cntrs_p        = SIZE(CNTRS_P, /DIMENSIONS)

      ENDIF
    ENDELSE

;   windowSize            = [640, 480]
    windowSize            = [800, 600]

    JBuffer               =  OBJ_NEW('IDLgrBuffer', DIMENSIONS = windowSize)
    JViewgroup            =  OBJ_NEW('IDLgrViewgroup'                      )

    JBuffer               -> SetProperty, GRAPHICS_TREE = JViewgroup

    JView                 =  OBJ_NEW('IDLgrView',                              $
                                      PROJECTION     =   2,                    $
                                      VIEWPLANE_RECT = [-2.0, -2.0, 7.5, 7.5], $
                                      ZCLIP          = [4.5, -4.5],            $
                                      EYE            =  4.6                    $
                                    )

    JViewgroup            -> Add, JView
    JModel                =  OBJ_NEW('IDLgrModel', /LIGHTING )
    JView                 -> Add, JModel

    JVolumeP              = makeVolume(JVOX_P, coltabP, cc)
    JVolumeM              = makeVolume(JVOX_M, coltabM, cc)

    str_title             = getTitle(qty, action, desc_label, I, threshold, width, prefix, global_j_minmax)

    JTTLfont              = OBJ_NEW('IDLgrFont', 'helvetica*bold', SIZE=16)
    JTitle                = OBJ_NEW('IDLgrText', str_title, LOCATIONS=[0.0, 1.15], ALIGNMENT = 0.0, FONT=JTTLfont, /ENABLE_FORMATTING)

    JModel                -> Add, JTitle

    JzTitle               = OBJ_NEW('IDLgrText', 'Z')
    JxTitle               = OBJ_NEW('IDLgrText', 'X')
    JyTitle               = OBJ_NEW('IDLgrText', 'Y', UPDIR=[0,0,-1], BASELINE=[0,1,0])

    JModel                -> Add, JVolumeP
    JModel                -> Add, JVolumeM

    jVolAxes              = makeAxes(cc, size_jvox)

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

    IF (action EQ 'isf' AND KEYWORD_SET(inc_cont) ) THEN BEGIN

;     JModel              -> Add, APContour
;     JModel              -> Add, AMContour

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
  
END
