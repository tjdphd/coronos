PRO run_ave_evsz, dsc_lab, efld, first_step, last_step

  ZERO                  = 0.0D                                      ; Define some constants for convenience
  ONE                   = 1.0D                                      ; everywhere
  TWO                   = 2.0D
  THREE                 = 3.0D
  PI                    = 3.141592653589793D
  TWO_THIRDS            = TWO / THREE
  TWO_PI                = TWO*PI

  CZERO                 = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  CI                    = DCOMPLEX(ZERO,ONE)                       ; and to make sure that DOUBLE's are used

  i_z                   = 0 ; z - coordinate
  i_me                  = 1 ; megnetic energy in slab at z
  i_pe                  = 2 ; kinetic energy in slab at z
  i_ch                  = 3 ; cross helicity in slab at z
  i_ep                  = 4 ; elsasser energy
  i_em                  = 5 ; elsasser energy
  i_ce                  = 6 ; int (J^2 dz)
  i_oe                  = 7 ; int (Omega^2 dz)
  i_zp                  = 8 ; int( ( J + Omega)^2 dz)
  i_zm                  = 9 ; int( ( J - Omega)^2 dz)
  i_nch                 = 10 ; normalized cross helicity
  i_km                  = 11 ; square - "magnetic" wave number 
  i_kp                  = 12 ; square - "kinetic"  wave number 
  i_kzp                 = 13 ; square - "Z+" wave number
  i_kzm                 = 14 ; square - "Z+" wave number

  i_te                  = 15 ; total energy at z
  i_re                  = 16 ; residual energy at z
  i_ne                  = 17 ; normalized residual energy
  i_rzp                 = 18 ; (Before Kolmogorov length scale calculation)
  i_rzm                 = 19 ; (Before Kolmogorov length scale calculation)

  i_zpm                 = 20 ; Z+Z-  where  Z+ = |(Z+)^2|^(1/2) = |ep^2|^(1/2) and similarly for Z- (em).

  i_likp                = 21 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm                = 22 ; IK length length scale  lambda_IK,-  (see notes)

  i_lkop                = 23 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom                = 24 ; Kolmogorov length scale lambda_Kol,- (see notes)
  i_lpl                 = 25 ; length scale from F.T. of z^+
  i_lmi                 = 26 ; length scale from F.T. of z^+
  i_lpla                = 27 ; alternative definition for lpl
  i_lmia                = 28 ; alternative definition for lmi
  i_lsfp                = 29 ; length scale from structure function for lambda +
  i_lsfm                = 30 ; length scale from structure function for lambda -
  
  aumlaut               = STRING(228B)
  CASE efld OF
  'zz'  : BEGIN 
          i_efld        = i_z
          str_efld      = 'evsz_zz_'
          END
  'me'  : BEGIN 
          i_efld        = i_me
          str_efld      = 'evsz_me_'
          str_title     = 'Magnetic Energy vs z'
          str_y_title   = 'E!DM,!9x!X!N(z)'
          END
  'pe'  : BEGIN 
          i_efld        = i_pe
          str_efld      = 'evsz_pe_'
          str_title     = 'Kinetic Energy vs z'
          str_y_title   = 'E!DK,!9x!X!N(z)'
          END
  'ch'  : BEGIN 
          i_efld        = i_ch
          str_efld      = 'evsz_ch_'
          str_title     = 'Cross Helicity vs z'
          str_y_title   = 'H!IC!N(z)'
          END
  'ep'  : BEGIN 
          i_efld        = i_ep
          str_efld      = 'evsz_ep_'
          str_title     = 'Z!E+!N Els' + aumlaut + 'sser Energy vs z'
          str_y_title   = 'E!DZ!E+!N(z)'
          END
  'em'  : BEGIN 
          i_efld        = i_em
          str_efld      = 'evsz_em_'
          str_title     = 'Z!E-!N Els' + aumlaut + 'sser Energy vs z'
          str_y_title   = 'E!DZ!E-!N(z)'
          END
  'ce'  : BEGIN 
          i_efld        = i_ce
          str_efld      = 'evsz_ce_'
          str_title     = 'Square-Integrated Current vs z'
          str_y_title   = '<!9!!G!Dx!U!X2!NA!9!!!X!U2!N>(z)'
          END
  'oe'  : BEGIN 
          i_efld        = i_oe
          str_efld      = 'evsz_oe_'
          str_title     = 'Square-Integrated Vorticity vs z'
          str_y_title   = '<!9!!G!Dx!U!X2!N!7u!9!!!X!U2!N>(z)'
          END
  'zp'  : BEGIN 
          i_efld        = i_zp
          str_efld      = 'evsz_zp_'
          str_title     = 'Square-Integrated Sum of Current and Vorticity vs z'
          str_y_title   = '<!9!!!XJ+!7X!9!!!X>(z)'
          END
  'zm'  : BEGIN 
          i_efld        = i_zm
          str_efld      = 'evsz_zm_'
          str_title     = 'Square-Integrated Difference of Current and Vorticity vs z'
          str_y_title   = '<!9!!!XJ-!7X!9!!!X!U2!N>(z)'
          END
  'nch' : BEGIN 
          i_efld        = i_nch
          str_efld      = 'evsz_nch_'
          str_title     = 'Normalized Cross-Helicity vs z'
          str_y_title   = '!7r!X!DC!N(z)'
          END
  'km'  : BEGIN 
          i_efld        = i_km
          str_efld      = 'evsz_km_'
          str_title     = 'Inverse Square Magnetic Energy Length Scale vs z'
          str_y_title   = 'L!DE!IM!U-2!N'
          END
  'kp'  : BEGIN 
          i_efld        = i_kp
          str_efld      = 'evsz_kp_'
          str_title     = 'Inverse Square Kinetic Energy Length Scale vs z'
          str_y_title   = 'L!DE!IK!U-2!N'
          END
  'kzp' : BEGIN 
          i_efld        = i_kzp
          str_efld      = 'evsz_kzp_'
          str_title     = 'Inverse Square Z!E+!N Els' + aumlaut + 'sser Length Scale vs z'
          str_y_title   = 'L!DE!IZ!E+!U-2!N'
          END
  'kzm' : BEGIN 
          i_efld        = i_kzm
          str_efld      = 'evsz_kzm_'
          str_title     = 'Inverse Square Z!E-!N Els' + aumlaut + 'sser Length Scale vs z'
          str_y_title   = 'L!DE!IZ!E-!U-2!N'
          END
  'te'  : BEGIN 
          i_efld        = i_te
          str_efld      = 'evsz_te_'
          str_title     = 'Total Energy vs z'
          str_y_title   = 'E!DT!N(z)'
          END
  're'  : BEGIN 
          i_efld        = i_re
          str_efld      = 'evsz_re_'
          str_title     = 'Residual Energy vs z'
          str_y_title   = 'E!Dr!N(z)'
          END
  'ne'  : BEGIN 
          i_efld        = i_ne
          str_efld      = 'evsz_ne_'
          str_title     = 'Normalized Residual Energy vs z'
          str_y_title   = '!7r!X!DD!N(z)'
          END
  'rzp' : BEGIN 
          i_efld        = i_rzp
          str_efld      = 'evsz_rzp_'
          str_title     =  '(Z!D+!N!U2!N)!U1/2!N vs z'
          str_y_title   =  'Z!D+!NZ!D-!N(z)'
          END
  'rzm' : BEGIN 
          i_efld        = i_rzm
          str_efld      = 'evsz_rzm_'
          str_title     =  '(Z!D-!N!U2!N)!U1/2!N vs z'
          END
  'zpm' : BEGIN 
          i_efld        = i_zpm
          str_efld      = 'evsz_zpm_'
          str_title     =  'Z!D+!NZ!D-!N vs z'
          str_y_title   =  'Z!D+!NZ!D-!N(z)'
          END
  'likp': BEGIN 
          i_efld        = i_likp
          str_efld      = 'evsz_likp_'
          str_title     = 'IK Length Scale !7k!X!DIK,+!N vs z'
          str_y_title   = '!7k!X!DIK,+!N(z)'
          END
  'likm': BEGIN 
          i_efld        = i_likm
          str_efld      = 'evsz_likm_'
          str_title     = 'IK Length Scale !7k!X!DIK,-!N vs z'
          str_y_title   = '!7k!X!DIK,-!N(z)'
          END
  'lkop': BEGIN 
          i_efld        = i_lkop
          str_efld      = 'evsz_lkop_'
          str_title     = 'Kolmogorov Length Scale !7k!X!DKol,+!N vs z'
          str_y_title   = '!7k!X!DKol,+!N(z)'
          END
  'lkom': BEGIN 
          i_efld        = i_lkom
          str_efld      = 'evsz_lkom_'
          str_title     = 'Kolmogorov Length Scale !7k!X!DKol,-!N vs z'
          str_y_title   = '!7k!X!DKol,-!N(z)'
          END
  'lkom': BEGIN 
          i_efld        = i_lkom
          str_efld      = 'evsz_lkom_'
          str_title     = 'Kolmogorov Length Scale !7k!X!DKol,-!N vs z'
          str_y_title   = '!7k!X!DKol,-!N(z)'
          END
  'lpl' : BEGIN 
          i_efld        = i_lpl
          str_efld      = 'evsz_lpl_'
          str_title     = 'Length Scale From Z!U+!N(k); !7k!X!D+!N vs z'
          str_y_title   = '!7k!X!D+!N(z)'
          END
  'lmi' : BEGIN 
          i_efld        = i_lmi
          str_efld      = 'evsz_lmi_'
          str_title     = 'Length Scale From Z!U-!N(k); !7k!X!D-!N vs z'
          str_y_title   = '!7k!X!D-!N(z)'
          END
  'lpla': BEGIN 
          i_efld        = i_lpla
          str_efld      = 'evsz_lpla_'
          str_title     = 'Length Scale From Z!U+!N(k); !7k!X!D+,alt!N vs z '
          str_y_title   = '!7k!X!D+,alt!N(z)'
          END
  'lmia': BEGIN 
          i_efld        = i_lmia
          str_efld      = 'evsz_lmia_'
          str_title     = 'Length Scale From Z!U-!N(k); !7k!X!D-,alt!N vs z '
          str_y_title   = '!7k!X!D-,alt!N(z)'
          END
  'lsfp': BEGIN 
          i_efld        = i_lsfp
          str_efld      = 'evsz_lsfp_'
          str_title     = 'Length Scale From Z!U+!N(z,x!D!9^!X!N);  vs z'
          str_y_title   = '!7k!X!DC,+!N(z)'
          END
  'lsfm': BEGIN 
          i_efld        = i_lsfm
          str_efld      = 'evsz_lsfm_'
          str_title     = 'Length Scale From Z!U-!N(z,x!D!9^!X!N);  vs z'
          str_y_title   = '!7k!X!DC,-!N(z)'
          END
  ELSE: BEGIN 
        i_efld          = -1
        str_efld        = 'evsz_xx_' 
        str_title       = 'Something is terribly wrong'
        str_y_title     = 'This is just not right'
        END
  ENDCASE

  ip1                   = scan_parameters('ip1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2                   = scan_parameters('ip2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n1                    = 2^ip1
  n2                    = 2^ip2
  n3                    = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp                    = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl                    = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res                 = n1                                                      ; resolution in x
  y_res                 = n2                                                      ; resolution in y
  z_res                 = n3 * mp                                                 ; resolution in z

  dz                    = zl / z_res
  dz_m1                 = z_res / zl                                              ; dz^(-1)

; i_x_res               = UINT(x_res)
; i_y_res               = UINT(y_res)
; i_z_res               = UINT(z_res)
   
; str_x_res             = STRTRIM(i_x_res, 2)
; str_y_res             = STRTRIM(i_y_res, 2)
; str_z_res             = STRTRIM(i_z_res, 2)
  

; cur_dir               = GETENV('PWD')
; eps_out_dir           = cur_dir + '/ra_evsz/' + efld + '/eps'

; res_str               = str_x_res + '_' + str_z_res
; str_run_span          = 'srs_' + STRTRIM(first_step,2) + '-' + STRTRIM(last_step,2) + '_'
; data_out_file         = cur_dir + '/ra_evsz/' + 'ra_' + efld + '_' + str_run_span + res_str + '_L-' + str_zl + '.dat'


; tfirst                = scan_parameters('tstart', first_step, dsc_lab)
; tlast                 = scan_parameters('tstart', last_step,  dsc_lab)
; str_dt                = STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)
; str_zl                = STRTRIM(UINT(zl),2)
; case_res_str          = + ' (L = ' + str_zl    + ')'
; str_title             = str_title + '!C Averaged Over ' + str_dt + case_res_str

  IF (i_efld GE i_lpl && i_efld LE i_lsfm) THEN BEGIN

      IF (i_efld EQ i_lsfm) THEN PRINT, 'run_ae_esz: calculating KK'
      KK                = calcKK(n1,n2)

      IF (i_efld EQ i_lsfp || i_efld EQ i_lsfm) THEN BEGIN

        IF (i_efld EQ i_lsfm) THEN PRINT, 'run_ae_esz: calculating GOFX'
          GOFX          = calcGofx(KK, dsc_lab) & STOP

      ENDIF ELSE BEGIN

        GOFX            = FLTARR(n1,n2)
        GOFX[*,*]       = ONE

      ENDELSE

  ENDIF ELSE BEGIN

      KK                = FLTARR(n1,n2)
      GOFX              = FLTARR(n1,n2)
      KK[*,*]           = ONE
      GOF[*,*]          = ONE

  ENDELSE

  E                     = open_evsz_data_file( dsc_lab, first_step)
  extE                  = calc_extE( E, dsc_lab, first_step, i_efld, KK, GOFX )

  ineq                  = WHERE(E[*,0:14] - extE[*,0:14], count)

  IF (count NE 0) THEN PRINT, ' ra_evsz: WARNING - E and extE inconsistent. count = ', count

  E                     = extE

  size_E                = SIZE(extE)

  n_lines               = size_E[1]
  n_cols                = size_E[2]


  FOR I = first_step + 1, last_step DO BEGIN

       E_next           = open_evsz_data_file( dsc_lab, I)
       extE_next        = calc_extE(E_next, dsc_lab, I, i_efld, KK, GOFX)
       E_next           = extE_next

       FOR J = 0, i_lsfm DO E[*,J] = E[*, J] + E_next[*, J] 

  ENDFOR

  tot_steps             = FLOAT(last_step - first_step + 1)

  E[*,*]                = E[*,*] / tot_steps

; ~ extended portion ~;

  dzpdz                 = FLTARR(n_lines)
  dzmdz                 = FLTARR(n_lines)

 FOR I = 0, n_lines - 2 DO BEGIN

   delta_z             = E[I+1,i_z] - E[I, i_z]

   IF (delta_z NE dz) THEN BEGIN
     PRINT, 'ra_evsz: WARNING - value of dz inconsistent with deltaE[*,0] = ',  delta_z
   ENDIF

   dzpdz[I]            = (E[I+1,i_ep] - E[I, i_ep]) * dz_m1
   dzmdz[I]            = (E[I+1,i_em] - E[I, i_em]) * dz_m1

 ENDFOR

 dzpdz[n_lines - 1]    = dzpdz[n_lines - 2]
 dzmdz[n_lines - 1]    = dzmdz[n_lines - 2]

 ; Kolmogorov: 

 E[*,i_lkop]           = ABS((E[*, i_ep] * E[*,i_rzm]) / dzpdz[*])
 E[*,i_lkom]           = ABS((E[*, i_em] * E[*,i_rzp]) / dzmdz[*])


 ; (ep^1/2 * em^1/2)

 E[*, i_zpm]           = E[*, i_rzp] * E[*, i_rzm]

 ; Iroshnikov-Kraichnan

 E[*, i_likp]          = ABS(E[*,i_zpm]^2 / dzpdz[*])
 E[*, i_likm]          = ABS(E[*,i_zpm]^2 / dzmdz[*])

 ;~ extended portion ~;

  i_x_res               = UINT(x_res)
  i_y_res               = UINT(y_res)
  i_z_res               = UINT(z_res)
   
  str_x_res             = STRTRIM(i_x_res, 2)
  str_y_res             = STRTRIM(i_y_res, 2)
  str_z_res             = STRTRIM(i_z_res, 2)

  tfirst                = scan_parameters('tstart', first_step, dsc_lab)
  tlast                 = scan_parameters('tstart', last_step,  dsc_lab)

  str_dt                = STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)

  str_zl                = STRTRIM(UINT(zl),2)
  case_res_str          = ' (L = ' + str_zl    + ')'
  str_title             = str_title + '!C Averaged Over ' + str_dt + case_res_str
  cur_dir               = GETENV('PWD')

  eps_out_dir           = cur_dir + '/ra_evsz/' + efld + '/eps'
  res_str               = str_x_res + '_' + str_z_res
  str_run_span          = 'srs_' + STRTRIM(first_step,2) + '-' + STRTRIM(last_step,2) + '_'
  data_out_file         = cur_dir + '/ra_evsz/' + 'ra_' + efld + '_' + str_run_span + res_str + '_L-' + str_zl + '.dat'

  OpenW, data_unit, data_out_file, /GET_LUN

  FOR K = 0, n_lines - 1 DO BEGIN

     PRINTF, data_unit,  E[K, i_z:i_lsfm], FORMAT='(31(E16.8,1x),:/)'

  ENDFOR

  FREE_LUN, data_unit

  case_file_str         = "ra_" + str_efld + str_run_span + res_str + '_zl_' + str_zl
  eps_out               = eps_out_dir + '/' + case_file_str + '.eps'

  min_x                 = MIN(E[*, i_z] )
  max_x                 = MAX(E[*, i_z] )

  IF (min_x GT 0.0) THEN min_x = 0.0

  min_y                 = MIN(E[*, i_efld])
  max_y                 = MAX(E[*, i_efld])

  x_rng                 = [min_x, max_x]
  y_rng                 = [min_y, max_y]

  SET_PLOT, 'PS'
  DEVICE, /ENCAPSULATED
  DEVICE, FILENAME      = eps_out

  PLOT, E[*,i_z], E[*, i_efld],        $
        CHARSIZE        = 1.2,         $
        LINESTYLE       = 0,           $
        YTICKFORMAT     = '(E8.1)',    $
        XMARGIN         = [12,4],      $
        YMARGIN         = [4,4],       $
        XRANGE          = x_rng,       $
        YRANGE          = y_rng,       $
        XTITLE          = 'z',         $
        YTITLE          = str_y_title, $
        THICK           = 2,           $
        TITLE           = str_title

  DEVICE, /CLOSE
        
  SET_PLOT, 'X'

END
