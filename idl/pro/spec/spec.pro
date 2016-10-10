PRO spec, dsc_lab, sfld, layer, first_step, last_step, minmax_mode

COMMON step, time

; i_t   = 0 ; time
; i_z   = 1 ; z - coordinate

  i_k    = 0 ; k-value
  i_pe   = 1 ; kinetic        energy spectrum in layer ? of time step ?
  i_ae   = 2 ; magnetic       energy spectrum in layer ? of time step ?
  i_ts   = 3 ; total          energy spectrum in layer ? of time step ?
  i_szp  = 4 ; elsasser+      energy spectrum in layer ? of time step ?
  i_szm  = 5 ; elsasser-      energy spectrum in layer ? of time step ?
  i_tz   = 6 ; total elsasser energy spectrum in layer ? of time step ?

  CASE sfld OF
  'pe' : i_sfld   = i_pe
  'ae' : i_sfld   = i_ae
  'ts' : i_sfld   = i_ts
  'zp' : i_sfld   = i_zp
  'zm' : i_sfld   = i_zm
  'tz' : i_sfld   = i_tz
  ELSE:  i_sfld   = -1

  ENDCASE
  CASE i_sfld OF
  '0' :  str_sfld = 'spec_kk_'
  '1' :  str_sfld = 'spec_pe_'
  '2' :  str_sfld = 'spec_ae_'
  '3' :  str_sfld = 'spec_ts_'
  '4' :  str_sfld = 'spec_zp_'
  '5' :  str_sfld = 'spec_zm_'
  '6' :  str_sfld = 'spec_tz_'
  ELSE:  str_sfld = 'spec_xx_'
  ENDCASE

  ip1                = scan_parameters('ip1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2                = scan_parameters('ip2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n3                 = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp                 = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl                 = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res              = 2^ip1                                                   ; resolution in x
  y_res              = 2^ip2                                                   ; resolution in y
  z_res              = n3 * mp                                                 ; resolution in z
  
  i_x_res            = UINT(x_res)
  i_y_res            = UINT(y_res)
  i_z_res            = UINT(z_res)
   
  str_x_res          = STRTRIM(i_x_res, 2)
  str_y_res          = STRTRIM(i_y_res, 2)
  str_z_res          = STRTRIM(i_z_res, 2)
  
  str_zl             = STRTRIM(UINT(zl),2)

  cur_dir            = GETENV('PWD')
  eps_out_dir        = cur_dir + '/spec/' + sfld + '/eps'

  res_str            = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res
; case_res_str       = res_str   + ' (zl = ' + str_zl    + ')'
  case_res_str       = + ' (L = ' + str_zl    + ')'

  IF (minmax_mode EQ  'global') THEN BEGIN
    glb_minmax_y     = FLTARR(2)
    glb_minmax_y     = global_sfld_minmax(i_sfld, dsc_lab, layer, first_step, last_step)
  ENDIF

; FOR I = first_step, last_step DO BEGIN

;      E             = open_spec_data_file( dsc_lab, I)

;      str_time      = STRTRIM(time,2)

;      str_step      = STRTRIM(I,2)
;      str_last_step = STRTRIM(last_step,2)
;      max_z_digits  = STRLEN(str_last_step)
;      IF (max_z_digits LT 3) THEN max_z_digits = 3
;      z_digits      = STRLEN(str_step)
;      zero_pad      = max_z_digits - z_digits
;      zero_str      = ''
;      IF(zero_pad NE 0) THEN BEGIN
;        FOR J = 1, zero_pad DO BEGIN
;          zero_str  = zero_str + '0'
;        ENDFOR
;      ENDIF
;      str_step      = zero_str + str_step

;      case_file_str = str_sfld + str_x_res + '_' + str_z_res + '_zl_' + str_zl + '-' + str_step
;      eps_out       = eps_out_dir + '/' + case_file_str + '.eps'

;      aumlaut = STRING(228B)
;      CASE i_sfld OF
;      '1' :  BEGIN 
;             str_title   = 'Magnetic Energy vs z at t = '
;             str_y_title = 'E!DM,!9x!X!N(z)'
;             END
;      '2' :  BEGIN 
;             str_title   = 'Kinetic Energy vs z at t = '
;             str_y_title = 'E!DK,!9x!X!N(z)'
;             END
;      '3' :  BEGIN 
;             str_title   = 'Cross Helicity vs z at t = '
;             str_y_title = 'H!IC!N(z)'
;             END
;      '4' :  BEGIN 
;             str_title   = 'Z!E+!N Els' + aumlaut + 'sser Energy vs z at t = '
;             str_y_title = 'E!DZ!E+!N(z)'
;             END
;      '5' :  BEGIN 
;             str_title   = 'Z!E-!N Els' + aumlaut + 'sser Energy vs z at t = '
;             str_y_title = 'E!DZ!E-!N(z)'
;             END
;      '6' :  BEGIN 
;             str_title   = 'Square-Integrated Current Density vs z at t = '
;             str_y_title = '<!9!!G!Dx!U!X2!NA!9!!!X!U2!N>(z)'
;             END
;      '7' :  BEGIN 
;             str_title   = 'Square-Integrated Vorticity vs z at t = '
;             str_y_title = '<!9!!G!Dx!U!X2!N!7u!9!!!X!U2!N>(z)'
;             END
;      '8' :  BEGIN 
;             str_title   = 'Square-Integrated Current/Vorticity Sum vs z at t = '
;             str_y_title = '<!9!!!XJ+!7X!9!!!X>(z)'
;             END
;      '9' :  BEGIN 
;             str_title   = 'Square-Integrated Current/Vorticity Difference vs z at t = '
;             str_y_title = '<!9!!!XJ-!7X!9!!!X!U2!N>(z)'
;             END
;      '10':  BEGIN 
;             str_title   = 'Normalized Cross-Helicity vs z at t = '
;             str_y_title = '!7r!X!DC!N(z)'
;             END
;      '11':  BEGIN 
;             str_title   = 'Inverse Square Magnetic Energy Length Scale vs z at t = '
;             str_y_title = 'L!DE!IM!U-2!N'
;             END
;      '12':  BEGIN 
;             str_title   = 'Inverse Square Kinetic Energy Length Scale vs z at t = '
;             str_y_title = 'L!DE!IK!U-2!N'
;             END
;      '13': BEGIN 
;              str_title   = 'Inverse Square Z!E+!N Els' + aumlaut + 'sser Length Scale vs z at t = '
;              str_y_title = 'L!DE!IZ!E+!U-2!N'
;             END
;      '14': BEGIN 
;              str_title   = 'Inverse Square Z!E-!N Els' + aumlaut + 'sser Length Scale vs z at t = '
;              str_y_title = 'L!DE!IZ!E-!U-2!N'
;             END
;      ELSE: BEGIN 
;             str_title   = 'Something is terribly wrong at t =  '
;             str_y_title = ''
;             END
;      ENDCASE

;      str_title          = str_title + str_time + case_res_str

;      size_E = SIZE(E)

;      min_x  = MIN(E[*, i_z   ])
;      max_x  = MAX(E[*, i_z   ])

;      IF (min_x GT 0.0) THEN min_x = 0.0

;      x_rng  = [min_x, max_x]

;      IF (minmax_mode EQ 'global') THEN BEGIN
;        y_rng  = glb_minmax_y
;      ENDIF ELSE BEGIN
;        min_y  = MIN(E[*, i_sfld])
;        max_y  = MAX(E[*, i_sfld])
;        y_rng  = [min_y, max_y]
;      ENDELSE

;        SET_PLOT, 'PS' 
;        DEVICE, /ENCAPSULATED
;        DEVICE, FILENAME = eps_out

;        PLOT, E[*,i_z], E[*, i_sfld],    $
;              CHARSIZE    = 1.0,         $
;              LINESTYLE   = 0,           $
;              YTICKFORMAT = '(E8.1)',    $
;              XMARGIN     =[12,4],       $
;              YMARGIN     =[4,4],        $
;              XRANGE      = x_rng,       $
;              YRANGE      = y_rng,       $
;              XTITLE      = 'z',         $
;              YTITLE      = str_y_title, $
;              THICK       = 2,           $
;              TITLE       = str_title

;       
;        DEVICE, /CLOSE
;       
; ENDFOR

; SET_PLOT, 'X'

END
