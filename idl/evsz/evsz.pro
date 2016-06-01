PRO evsz, dsc_lab, efld, first_step, last_step, minmax_mode
  COMPILE_OPT IDL2

  COMMON step, time

  i_z    = 0 ; z - coordinate
  i_me   = 1 ; megnetic energy in slab at z
  i_pe   = 2 ; kinetic energy in slab at z
  i_ch   = 3 ; cross helicity in slab at z
  i_ep   = 4 ; elsasser energy
  i_em   = 5 ; elsasser energy
  i_ce   = 6 ; int (J^2 dz)
  i_oe   = 7 ; int (Omega^2 dz)
  i_zp   = 8 ; int( ( J + Omega)^2 dz)
  i_zm   = 9 ; int( ( J - Omega)^2 dz)
  i_nch  = 10 ; normalized cross helicity
  i_km   = 11 ; square - "magnetic" wave number 
  i_kp   = 12 ; square - "kinetic"  wave number 
  i_kzp  = 13 ; square - "Z+" wave number
  i_kzm  = 14 ; square - "Z+" wave number
  i_te   = 15 ; total energy at z
  i_re   = 16 ; residual energy at z
  i_ne   = 17 ; normalized residual energy
  i_zpm  = 18 ; Z+Z-  where  Z+ = |(Z+)^2|^(1/2) = |ep^2|^(1/2) and similarly for Z- (em).
  i_lkop = 19 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom = 20 ; Kolmogorov length scale lambda_Kol,- (see notes)
  i_likp = 21 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm = 22 ; IK length length scale  lambda_IK,-  (see notes)
  
  CASE efld OF
  'zz'  : i_efld   = i_z
  'me'  : i_efld   = i_me
  'pe'  : i_efld   = i_pe
  'ch'  : i_efld   = i_ch
  'ep'  : i_efld   = i_ep
  'em'  : i_efld   = i_em
  'ce'  : i_efld   = i_ce
  'oe'  : i_efld   = i_oe
  'zp'  : i_efld   = i_zp
  'zm'  : i_efld   = i_zm
  'nch' : i_efld   = i_nch
  'km'  : i_efld   = i_km
  'kp'  : i_efld   = i_kp
  'kzp' : i_efld   = i_kzp
  'kzm' : i_efld   = i_kzm
  'te'  : i_efld   = i_te
  're'  : i_efld   = i_re
  'ne'  : i_efld   = i_ne
  'zpm' : i_efld   = i_zpm
  'lkop': i_efld   = i_lkop
  'lkom': i_efld   = i_lkom
  'likp': i_efld   = i_likp
  'likm': i_efld   = i_likm
  ELSE: i_efld   = -1

  ENDCASE
  CASE i_efld OF
  '0' :  str_efld = 'evsz_zz_'
  '1' :  str_efld = 'evsz_me_'
  '2' :  str_efld = 'evsz_pe_'
  '3' :  str_efld = 'evsz_ch_'
  '4' :  str_efld = 'evsz_ep_'
  '5' :  str_efld = 'evsz_em_'
  '6' :  str_efld = 'evsz_ce_'
  '7' :  str_efld = 'evsz_oe_'
  '8' :  str_efld = 'evsz_zp_'
  '9' :  str_efld = 'evsz_zm_'
  '10':  str_efld = 'evsz_nch_'
  '11':  str_efld = 'evsz_km_'
  '12':  str_efld = 'evsz_kp_'
  '13':  str_efld = 'evsz_kzp_'
  '14':  str_efld = 'evsz_kzm_'
  '15':  str_efld = 'evsz_te_'
  '16':  str_efld = 'evsz_re_'
  '17':  str_efld = 'evsz_ne_'
  '18':  str_efld = 'evsz_zpm_'
  '19':  str_efld = 'evsz_lkop_'
  '20':  str_efld = 'evsz_lkom_'
  '21':  str_efld = 'evsz_likp_'
  '22':  str_efld = 'evsz_likm_'
  ELSE: str_efld = 'evsz_xx_'
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
  eps_out_dir        = cur_dir + '/evsz/' + efld + '/eps'

  res_str            = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res

  case_res_str       = + ' (L = ' + str_zl    + ')'

  IF (minmax_mode EQ  'global') THEN BEGIN
    glb_minmax_y     = FLTARR(2)
    glb_minmax_y     = global_efld_minmax(i_efld, dsc_lab, first_step, last_step)
  ENDIF

  FOR I = first_step, last_step DO BEGIN

       E             = open_evsz_data_file( dsc_lab, I)
       extE          = calc_extE(E, desc_lab)

       size_extE     = SIZE(extE)

       ineq          = WHERE(E[*,0:14] - extE[*,0:14], count)

       IF (count NE 0) THEN PRINT, ' evsz: WARNING - E and extE inconsistent. count = ', count

       E             = extE

       str_time      = STRTRIM(time,2)

       str_step      = STRTRIM(I,2)
       str_last_step = STRTRIM(last_step,2)
       max_z_digits  = STRLEN(str_last_step)
       IF (max_z_digits LT 3) THEN max_z_digits = 3
       z_digits      = STRLEN(str_step)
       zero_pad      = max_z_digits - z_digits
       zero_str      = ''
       IF(zero_pad NE 0) THEN BEGIN
         FOR J = 1, zero_pad DO BEGIN
           zero_str  = zero_str + '0'
         ENDFOR
       ENDIF
       str_step      = zero_str + str_step

       case_file_str = str_efld + str_x_res + '_' + str_z_res + '_zl_' + str_zl + '-' + str_step
       eps_out       = eps_out_dir + '/' + case_file_str + '.eps'
       tif_out       = eps_out_dir + '/' + case_file_str + '.tif'

       aumlaut = STRING(228B)
       CASE i_efld OF
       '1' :  BEGIN 
              str_title   = 'Magnetic Energy vs z at t = '
              str_y_title = 'E!DM,!9x!X!N(z)'
              END
       '2' :  BEGIN 
              str_title   = 'Kinetic Energy vs z at t = '
              str_y_title = 'E!DK,!9x!X!N(z)'
              END
       '3' :  BEGIN 
              str_title   = 'Cross Helicity vs z at t = '
              str_y_title = 'H!IC!N(z)'
              END
       '4' :  BEGIN 
              str_title   = 'Z!E+!N Els' + aumlaut + 'sser Energy vs z at t = '
              str_y_title = 'E!DZ!E+!N(z)'
              END
       '5' :  BEGIN 
              str_title   = 'Z!E-!N Els' + aumlaut + 'sser Energy vs z at t = '
              str_y_title = 'E!DZ!E-!N(z)'
              END
       '6' :  BEGIN 
              str_title   = 'Square-Integrated Current Density vs z at t = '
              str_y_title = '<!9!!G!Dx!U!X2!NA!9!!!X!U2!N>(z)'
              END
       '7' :  BEGIN 
              str_title   = 'Square-Integrated Vorticity vs z at t = '
              str_y_title = '<!9!!G!Dx!U!X2!N!7u!9!!!X!U2!N>(z)'
              END
       '8' :  BEGIN 
              str_title   = 'Square-Integrated Current/Vorticity Sum vs z at t = '
              str_y_title = '<!9!!!XJ+!7X!9!!!X>(z)'
              END
       '9' :  BEGIN 
              str_title   = 'Square-Integrated Current/Vorticity Difference vs z at t = '
              str_y_title = '<!9!!!XJ-!7X!9!!!X!U2!N>(z)'
              END
       '10':  BEGIN 
              str_title   = 'Normalized Cross-Helicity vs z at t = '
              str_y_title = '!7r!X!DC!N(z)'
              END
       '11':  BEGIN 
              str_title   = 'Inverse Square Magnetic Energy Length Scale vs z at t = '
              str_y_title = 'L!DE!IM!U-2!N'
              END
       '12':  BEGIN 
              str_title   = 'Inverse Square Kinetic Energy Length Scale vs z at t = '
              str_y_title = 'L!DE!IK!U-2!N'
              END
       '13': BEGIN 
               str_title   = 'Inverse Square Z!E+!N Els' + aumlaut + 'sser Length Scale vs z at t = '
               str_y_title = 'L!DE!IZ!E+!U-2!N'
              END
       '14': BEGIN
               str_title   = 'Inverse Square Z!E-!N Els' + aumlaut + 'sser Length Scale vs z at t = '
               str_y_title = 'L!DE!IZ!E-!U-2!N'
              END
       '15':  BEGIN
               str_title   = 'Total Energy vs z at t = '
               str_y_title = 'E!DT!N(z)'
              END
       '16':  BEGIN
               str_title   = 'Residual Energy vs z at t = '
               str_y_title = 'E!Dr!N(z)'
              END
       '17':  BEGIN
               str_title   = 'Normalized Residual Energy vs z at t = '
               str_y_title = '!7r!X!DD!N(z)'
              END
       '18':  BEGIN
               str_title   =  'Z!D+!NZ!D-!N vs z at t = '
               str_y_title =  'Z!D+!NZ!D-!N(z)'
              END
       '19':  BEGIN
               str_title   = 'Kolmogorov Length Scale !7k!X!DKol,+!N vs z at t = '
               str_y_title = '!7k!X!DKol,+!N(z)'
              END
       '20':  BEGIN
               str_title   = 'Kolmogorov Length Scale !7k!X!DKol,-!N vs z at t = '
               str_y_title = '!7k!X!DKol,-!N(z)'
              END
       '21':  BEGIN
               str_title   = 'IK Length Scale !7k!X!DIK,+!N vs z at t = '
               str_y_title = '!7k!X!DIK,+!N(z)'
              END
       '22':  BEGIN
               str_title   = 'IK Length Scale !7k!X!DIK,-!N vs z at t = '
               str_y_title = '!7k!X!DIK,-!N(z)'
              END
       ELSE: BEGIN 
              str_title    = 'Something is terribly wrong at t =  '
              str_y_title  = ''
              END
       ENDCASE

       str_title          = str_title + str_time + case_res_str

;      size_E = SIZE(E)

       min_x  = MIN(E[*, i_z   ])
       max_x  = MAX(E[*, i_z   ])

       IF (min_x GT 0.0) THEN min_x = 0.0

       x_rng  = [min_x, max_x]

       IF (minmax_mode EQ 'global') THEN BEGIN
         y_rng  = glb_minmax_y
       ENDIF ELSE BEGIN
         min_y  = MIN(E[*, i_efld])
         max_y  = MAX(E[*, i_efld])
         y_rng  = [min_y, max_y]
       ENDELSE

;      X[*] = E[*, i_z   ]
;      Y[*] = E[*, i_efld]

;      the_plot = PLOT( E[*,i_z], E[*,i_efld]             $
;                      AXIS_STYLE = 2,                    $
;                      MARGIN     = [0.0, 0.0, 0.0, 0.0], $
;                      XRANGE     = x_rng,                $
;                      YRANGE     = y_rng,                $
;                      TITLE      = str_title,            $
;                      XTITLE     = 'z',                  $
;                      YTITLE     = str_y_title,          $
;                      FONT_SIZE  = 12,                   $
;                      COLOR      = 'red',                $
;                      HIDE       = 1                     $
;                     )


;      image    = the_plot.copyWindow()
;      image    = REVERSE(image,3)

;      WRITE_TIFF, tif_out, image

        SET_PLOT, 'PS' 
        DEVICE, /ENCAPSULATED
        DEVICE, FILENAME = eps_out

        PLOT, E[*,i_z], E[*, i_efld],    $
              CHARSIZE    = 1.2,         $
              LINESTYLE   = 0,           $
              YTICKFORMAT = '(E8.1)',    $
              XMARGIN     =[12,4],       $
              YMARGIN     =[4,4],        $
              XRANGE      = x_rng,       $
              YRANGE      = y_rng,       $
              XTITLE      = 'z',         $
              YTITLE      = str_y_title, $
              THICK       = 2,           $
              TITLE       = str_title

        DEVICE, /CLOSE
        
  ENDFOR

; SET_PLOT, 'X'

END
