PRO run_ave_evsz, dsc_lab, efld, first_step, last_step, minmax_mode

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

  i_rzp  = 18 ; (Before Kolmogorov length scale calculation)
  i_rzm  = 19 ; (Before Kolmogorov length scale calculation)

  i_zpm  = 18 ; Z+Z-  where  Z+ = |(Z+)^2|^(1/2) = |ep^2|^(1/2) and similarly for Z- (em).

  i_likp = 19 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm = 20 ; IK length length scale  lambda_IK,-  (see notes)

  i_lkop = 21 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom = 22 ; Kolmogorov length scale lambda_Kol,- (see notes)
  i_lpl  = 23 ; length scale from F.T. of z^+
  i_lmi  = 24 ; length scale from F.T. of z^+
  

  aumlaut             = STRING(228B)
  CASE efld OF
  'zz'  : BEGIN 
          i_efld      = i_z
          str_efld    = 'evsz_zz_'
          END
  'me'  : BEGIN 
          i_efld      = i_me
          str_efld    = 'evsz_me_'
          str_title   = 'Magnetic Energy vs z'
          str_y_title = 'E!DM,!9x!X!N(z)'
          END
  'pe'  : BEGIN 
          i_efld      = i_pe
          str_efld    = 'evsz_pe_'
          str_title   = 'Kinetic Energy vs z'
          str_y_title = 'E!DK,!9x!X!N(z)'
          END
  'ch'  : BEGIN 
          i_efld      = i_ch
          str_efld    = 'evsz_ch_'
          str_title   = 'Cross Helicity vs z'
          str_y_title = 'H!IC!N(z)'
          END
  'ep'  : BEGIN 
          i_efld      = i_ep
          str_efld    = 'evsz_ep_'
          str_title   = 'Z!E+!N Els' + aumlaut + 'sser Energy vs z'
          str_y_title = 'E!DZ!E+!N(z)'
          END
  'em'  : BEGIN 
          i_efld      = i_em
          str_efld    = 'evsz_em_'
          str_title   = 'Z!E-!N Els' + aumlaut + 'sser Energy vs z'
          str_y_title = 'E!DZ!E-!N(z)'
          END
  'ce'  : BEGIN 
          i_efld      = i_ce
          str_efld    = 'evsz_ce_'
          str_title   = 'Square-Integrated Current vs z'
          str_y_title = '<!9!!G!Dx!U!X2!NA!9!!!X!U2!N>(z)'
          END
  'oe'  : BEGIN 
          i_efld      = i_oe
          str_efld    = 'evsz_oe_'
          str_title   = 'Square-Integrated Vorticity vs z'
          str_y_title = '<!9!!G!Dx!U!X2!N!7u!9!!!X!U2!N>(z)'
          END
  'zp'  : BEGIN 
          i_efld      = i_zp
          str_efld    = 'evsz_zp_'
          str_title   = 'Square-Integrated Sum of Current and Vorticity vs z'
          str_y_title = '<!9!!!XJ+!7X!9!!!X>(z)'
          END
  'zm'  : BEGIN 
          i_efld      = i_zm
          str_efld    = 'evsz_zm_'
          str_title   = 'Square-Integrated Difference of Current and Vorticity vs z'
          str_y_title = '<!9!!!XJ-!7X!9!!!X!U2!N>(z)'
          END
  'nch' : BEGIN 
          i_efld      = i_nch
          str_efld    = 'evsz_nch_'
          str_title   = 'Normalized Cross-Helicity vs z'
          str_y_title = '!7r!X!DC!N(z)'
          END
  'km'  : BEGIN 
          i_efld      = i_km
          str_efld    = 'evsz_km_'
          str_title   = 'Inverse Square Magnetic Energy Length Scale vs z'
          str_y_title = 'L!DE!IM!U-2!N'
          END
  'kp'  : BEGIN 
          i_efld      = i_kp
          str_efld    = 'evsz_kp_'
          str_title   = 'Inverse Square Kinetic Energy Length Scale vs z'
          str_y_title = 'L!DE!IK!U-2!N'
          END
  'kzp' : BEGIN 
          i_efld      = i_kzp
          str_efld    = 'evsz_kzp_'
          str_title   = 'Inverse Square Z!E+!N Els' + aumlaut + 'sser Length Scale vs z'
          str_y_title = 'L!DE!IZ!E+!U-2!N'
          END
  'kzm' : BEGIN 
          i_efld      = i_kzm
          str_efld    = 'evsz_kzm_'
          str_title   = 'Inverse Square Z!E-!N Els' + aumlaut + 'sser Length Scale vs z'
          str_y_title = 'L!DE!IZ!E-!U-2!N'
          END
  'te'  : BEGIN 
          i_efld      = i_te
          str_efld    = 'evsz_te_'
          str_title   = 'Total Energy vs z'
          str_y_title = 'E!DT!N(z)'
          END
  're'  : BEGIN 
          i_efld      = i_re
          str_efld    = 'evsz_re_'
          str_title   = 'Residual Energy vs z'
          str_y_title = 'E!Dr!N(z)'
          END
  'ne'  : BEGIN 
          i_efld      = i_ne
          str_efld    = 'evsz_ne_'
          str_title   = 'Normalized Residual Energy vs z'
          str_y_title = '!7r!X!DD!N(z)'
          END
  'zpm' : BEGIN 
          i_efld      = i_zpm
          str_efld    = 'evsz_zpm_'
          str_title   =  'Z!D+!NZ!D-!N vs z'
          str_y_title =  'Z!D+!NZ!D-!N(z)'
          END
  'likp': BEGIN 
          i_efld      = i_likp
          str_efld    = 'evsz_likp_'
          str_title   = 'IK Length Scale !7k!X!DIK,+!N vs z'
          str_y_title = '!7k!X!DIK,+!N(z)'
          END
  'likm': BEGIN 
          i_efld      = i_likm
          str_efld    = 'evsz_likm_'
          str_title   = 'IK Length Scale !7k!X!DIK,-!N vs z'
          str_y_title = '!7k!X!DIK,-!N(z)'
          END
  'lkop': BEGIN 
          i_efld      = i_lkop
          str_efld    = 'evsz_lkop_'
          str_title   = 'Kolmogorov Length Scale !7k!X!DKol,+!N vs z'
          str_y_title = '!7k!X!DKol,+!N(z)'
          END
  'lkom': BEGIN 
          i_efld      = i_lkom
          str_efld    = 'evsz_lkom_'
          str_title   = 'Kolmogorov Length Scale !7k!X!DKol,-!N vs z'
          str_y_title = '!7k!X!DKol,-!N(z)'
          END
  'lkom': BEGIN 
          i_efld      = i_lkom
          str_efld    = 'evsz_lkom_'
          str_title   = 'Kolmogorov Length Scale !7k!X!DKol,-!N vs z'
          str_y_title = '!7k!X!DKol,-!N(z)'
          END
  'lpl' : BEGIN 
          i_efld      = i_lpl
          str_efld    = 'evsz_lpl_'
          str_title   = 'Length Scale From Z!U+!N(k); !7k!X!D+!N vs z'
          str_y_title = '!7k!X!D+!N(z)'
          END
  'lmi' : BEGIN 
          i_efld      = i_lmi
          str_efld    = 'evsz_lpl_'
          str_title   = 'Length Scale From Z!U-!N(k); !7k!X!D+!N vs z'
          str_y_title = '!7k!X!D+!N(z)'
          END
  ELSE: BEGIN 
        i_efld        = -1
        str_efld      = 'evsz_xx_' 
        str_title     = 'Something is terribly wrong'
        str_y_title   = 'This is just not right'
        END
  ENDCASE

  ip1                = scan_parameters('ip1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2                = scan_parameters('ip2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n3                 = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp                 = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl                 = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res              = 2^ip1                                                   ; resolution in x
  y_res              = 2^ip2                                                   ; resolution in y
  z_res              = n3 * mp                                                 ; resolution in z

  dz                 = zl / z_res
  dz_m1              = z_res / zl                                              ; dz^(-1)

  i_x_res            = UINT(x_res)
  i_y_res            = UINT(y_res)
  i_z_res            = UINT(z_res)
   
  str_x_res          = STRTRIM(i_x_res, 2)
  str_y_res          = STRTRIM(i_y_res, 2)
  str_z_res          = STRTRIM(i_z_res, 2)
  
  str_zl             = STRTRIM(UINT(zl),2)

  cur_dir            = GETENV('PWD')
  eps_out_dir        = cur_dir + '/ra_evsz/' + efld + '/eps'

  res_str            = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res
  case_res_str       = + ' (L = ' + str_zl    + ')'

  tfirst             = scan_parameters('tstart',first_step, dsc_lab)
  tlast              = scan_parameters('tstart',last_step,  dsc_lab)

  str_dt             = STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)

  str_title          = str_title + '!C Averaged Over ' + str_dt + case_res_str

  IF (minmax_mode EQ  'global') THEN BEGIN
    glb_minmax_y     = FLTARR(2)
    glb_minmax_y     = global_efld_minmax(i_efld, dsc_lab, first_step, last_step)
  ENDIF


  E                  = open_evsz_data_file( dsc_lab, first_step)
  extE               = calc_extE( E, dsc_lab, first_step, i_efld )

  ineq               = WHERE(E[*,0:14] - extE[*,0:14], count)
  IF (count NE 0) THEN PRINT, ' ra_evsz: WARNING - E and extE inconsistent. count = ', count

  E                  = extE

  size_E             = SIZE(extE)

  n_lines            = size_E[1]
  n_cols             = size_E[2]

  FOR I = first_step + 1, last_step DO BEGIN

       E_next        = open_evsz_data_file( dsc_lab, I)
       extE_next     = calc_extE(E_next, desc_lab, I, i_efld)
       E_next        = extE_next

       FOR J = 0, 24 DO BEGIN
         E[*,J]      = E[*, J] + E_next[*, J] 
       ENDFOR

  ENDFOR

  tot_steps          = FLOAT(last_step - first_step + 1)

  E[*,*]             = E[*,*] / tot_steps

; ~ extended portion ~;

  dzpdz              = FLTARR(n_lines)
  dzmdz              = FLTARR(n_lines)

  FOR I = 0, n_lines - 2 DO BEGIN

    delta_z          = E[I+1,i_z] - E[I, i_z]

    IF (delta_z NE dz) THEN BEGIN
      PRINT, 'ra_evsz: WARNING - value of dz inconsistent with deltaE[*,0] = ',  delta_z
    ENDIF

    dzpdz[I]         = (E[I+1,i_ep] - E[I, i_ep]) * dz_m1
    dzmdz[I]         = (E[I+1,i_em] - E[I, i_em]) * dz_m1

  ENDFOR

  dzpdz[n_lines - 1] = dzpdz[n_lines - 2]
  dzmdz[n_lines - 1] = dzmdz[n_lines - 2]

  ; Kolmogorov: 

  E[*,i_lkop]        = ABS((E[*, i_ep] * E[*,i_rzm]) / dzpdz[*])
  E[*,i_lkom]        = ABS((E[*, i_em] * E[*,i_rzp]) / dzmdz[*])


  ; lambda_(+/-) still need to be normalized by Z+^2 and Z-^2 time averages respectively

  E[*, i_lpl]        = E[*,i_lpl] / (E[*, i_rzp]^2)
  E[*, i_lmi]        = E[*,i_lmi] / (E[*, i_rzm]^2)

  ; (ep^1/2 * em^1/2)

  E[*, i_rzp]        = E[*, i_rzp] * E[*, i_rzm]

  ; now use i_zpm = 18 and i_likp =19 

  ; Iroshnikov-Kraichnan

  E[*, i_likp]       = ABS(E[*,i_zpm]^2 / dzpdz[*])
  E[*, i_likm]       = ABS(E[*,i_zpm]^2 / dzmdz[*])

; ~ extended portion ~;

  res_str            = str_x_res + '_' + str_z_res
  str_run_span       = 'srs_' + STRTRIM(first_step,2) + '-' + STRTRIM(last_step,2) + '_'

  data_out_file      = cur_dir + '/ra_evsz/' + 'ra_' + efld + '_' + str_run_span + res_str + '_L-' + str_zl + '.dat'

  OpenW, data_unit, data_out_file, /GET_LUN

  FOR K = 0, n_lines - 1 DO BEGIN

;    PRINTF, data_unit,  E[K, 0:24], FORMAT='(23(E24.16,1x),:/)'
     PRINTF, data_unit,  E[K, 0:24], FORMAT='(23(E16.8,1x),:/)'

  ENDFOR

  FREE_LUN, data_unit

  case_file_str      = "ra_" + str_efld + str_run_span + res_str + '_zl_' + str_zl
  eps_out            = eps_out_dir + '/' + case_file_str + '.eps'

  min_x              = MIN(E[*, i_z] )
  max_x              = MAX(E[*, i_z] )

  IF (min_x GT 0.0) THEN min_x = 0.0

  min_y              = MIN(E[*, i_efld])
  max_y              = MAX(E[*, i_efld])

  x_rng              = [min_x, max_x]
  y_rng              = [min_y, max_y]

  SET_PLOT, 'PS' 
  DEVICE, /ENCAPSULATED
  DEVICE, FILENAME   = eps_out

  PLOT, E[*,i_z], E[*, i_efld],     $
        CHARSIZE     = 1.2,         $
        LINESTYLE    = 0,           $
        YTICKFORMAT  = '(E8.1)',    $
        XMARGIN      =[12,4],       $
        YMARGIN      =[4,4],        $
        XRANGE       = x_rng,       $
        YRANGE       = y_rng,       $
        XTITLE       = 'z',         $
        YTITLE       = str_y_title, $
        THICK        = 2,           $
        TITLE        = str_title

  
  DEVICE, /CLOSE
        
  SET_PLOT, 'X'

END
