PRO spec, dsc_lab, sfld, first_layer, last_layer, first_step, last_step, minmax_mode, plot_mode

COMMON step, time

LOADCT, 13

ZERO                   = 0.0D
HALF                   = 0.5D
five_thirds            = 5.0D/3.0D
three_halves           = 3.0D/2.0D

  i_k                  = 0 ; k-value
  i_pe                 = 1 ; kinetic        energy spectrum in layer ? of time step ?
  i_ae                 = 2 ; magnetic       energy spectrum in layer ? of time step ?
  i_ts                 = 3 ; total          energy spectrum in layer ? of time step ?
  i_zp                 = 4 ; elsasser+      energy spectrum in layer ? of time step ?
  i_zm                 = 5 ; elsasser-      energy spectrum in layer ? of time step ?
  i_tz                 = 6 ; total elsasser energy spectrum in layer ? of time step ?

  CASE sfld OF
  'pe' : i_sfld        = i_pe
  'ae' : i_sfld        = i_ae
  'ts' : i_sfld        = i_ts
  'zp' : i_sfld        = i_zp
  'zm' : i_sfld        = i_zm
  'tz' : i_sfld        = i_tz
  ELSE:  i_sfld        = -1

  ENDCASE
  CASE i_sfld OF
  '0' :  str_sfld      = 'kk_'
  '1' :  str_sfld      = 'pe_'
  '2' :  str_sfld      = 'ae_'
  '3' :  str_sfld      = 'ts_'
  '4' :  str_sfld      = 'zp_'
  '5' :  str_sfld      = 'zm_'
  '6' :  str_sfld      = 'tz_'
  ELSE:  str_sfld      = 'xx_'
  ENDCASE

  tot_steps            = FLOAT(last_step  - first_step  + 1)
  tot_layers           = FLOAT(last_layer - first_layer + 1)

  slc_minmax_y         = DBLARR(tot_layers, 2) 

  FOR I_L = first_layer, last_layer DO BEGIN
    IF (minmax_mode EQ  'global') THEN BEGIN
      slc_minmax_y[I_L - first_layer, *] = global_sfld_minmax_by_slice(i_sfld, dsc_lab, I_L, first_step, last_step)
    ENDIF
  ENDFOR

  y_global_max         = MAX(slc_minmax_y[*,1])
  y_global_min         = MIN(slc_minmax_y[*,0])

  PRINT, "y_global_max = ", y_global_max
  PRINT, "y_global_min = ", y_global_min

  ip1                  = scan_parameters('ip1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2                  = scan_parameters('ip2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n3                   = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp                   = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl                   = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res                = 2^ip1                                                   ; resolution in x
  y_res                = 2^ip2                                                   ; resolution in y
  z_res                = n3 * mp                                                 ; resolution in z

  dz                   = zl / z_res

  i_x_res              = UINT(x_res)
  i_y_res              = UINT(y_res)
  i_z_res              = UINT(z_res)
   
  str_x_res            = STRTRIM(i_x_res, 2)
  str_y_res            = STRTRIM(i_y_res, 2)
  str_z_res            = STRTRIM(i_z_res, 2)
  
  str_zl               = STRTRIM(UINT(zl),2)
; str_zl               = '_L-' + str_zl

  cur_dir              = GETENV('PWD')

  data_out_dir         = cur_dir      + '/spec/'
  eps_out_dir          = data_out_dir + sfld + '/eps'
  res_str              = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res
;
;  EVSZ                 = open_evsz_data_file(dsc_lab, first_step)
;
;  size_evsz            = SIZE(EVSZ, /DIMENSIONS)
;
;  NCHVSZ               = DBLARR(size_evsz[0],2)
;  
;  i_z                  = 0
;  i_nch                = 10
;
;  NCHVSZ[*,0]          = EVSZ[*,i_z]
;  NCHVSZ[*,1]          = EVSZ[*,i_nch]
;
;  FOR S = first_step+1, last_step DO BEGIN
;    EVSZ               = open_evsz_data_file(dsc_lab, S)
;    NCHVSZ[*,1]        = NCHVSZ[*,1] + EVSZ[*,i_nch]
;  ENDFOR
;
;  NCHVSZ[*,1]          = NCHVSZ[*,1] / tot_steps
;
;; layer dependencies
;
;  tot_layers           = FLOAT(last_layer - first_layer + 1)
;
  str_tot_layers       = str_z_res
  str_first_step       = STRTRIM(first_step,2)
  str_last_step        = STRTRIM(last_step,2)


  WHILE (STRLEN(str_first_step) LT 3) DO BEGIN
    str_first_step     = '0' + str_first_step
  ENDWHILE

  WHILE (STRLEN(str_last_step) LT 3) DO BEGIN
    str_last_step      = '0' + str_last_step
  ENDWHILE

  WHILE (STRLEN(str_tot_layers) LT 3) DO BEGIN
    str_tot_layers     = '0' + str_tot_layers

  ENDWHILE

  str_srs_ave          = 'srs_' + str_first_step + '-' + str_last_step

; spec_idx_out_file    = data_out_dir + 'sp_idx_' + str_srs_ave + '_L-' + str_zl + '.dat'

;; OpenW, spidx_unit, spec_idx_out_file, /GET_LUN
;
;  SPIDX                = DBLARR(tot_layers, 5)
;
   FOR I_T = first_step, last_step DO BEGIN
     str_step             = STRTRIM(I_T,2)
     WHILE (STRLEN(str_step) LT 3) DO BEGIN
       str_step     = '0' + str_step
     ENDWHILE
     FOR L = first_layer, last_layer DO BEGIN
       
       PRINT, "spec: I_T = ", I_T, " L = ", L
   
       IF (L GT n3) THEN BEGIN
         is_layer_even    = (L MOD 2) EQ 0
         IF (is_layer_even) THEN BEGIN
           n_proc         = CEIL(UINT(L/n3)-1)
         ENDIF ELSE BEGIN
           n_proc         = FLOOR(UINT(L/n3))
         ENDELSE
       ENDIF ELSE BEGIN
         n_proc           = 0
       ENDELSE
   
       loc_layer          = UINT(L / ( n_proc+1 ))
       str_proc           = STRTRIM(n_proc,   2)
       len_str_proc       = STRLEN(str_proc)
   
       WHILE (STRLEN(str_proc) LT 3) DO BEGIN
         str_proc         = '0' + str_proc
       ENDWHILE
   
       str_layer          = STRTRIM(L,2)
       WHILE (STRLEN(str_layer) LT 3) DO BEGIN
         str_layer        = '0' + str_layer
       ENDWHILE
       
       str_loc_layer      = STRTRIM(loc_layer,2)
       WHILE (STRLEN(str_loc_layer) LT 3) DO BEGIN
         str_loc_layer    = '0' + str_loc_layer
       ENDWHILE
       
       ELOG               = open_spec_data_file( dsc_lab, L, I_T)
       size_ELOG          = SIZE(ELOG)
   
       n_lines            = size_ELOG[1]
       n_cols             = size_ELOG[2]
   
       ENLG               = DBLARR(n_lines,n_cols)
       ENLG_next          = DBLARR(n_lines,n_cols)
   
       ENLG[*,*]          = ELOG[*,*]
   
       FOR J = 1, 6 DO BEGIN
   
          e_zz            = WHERE(ELOG[*,J] EQ ZERO, zz_cnt)
          e_nz            = WHERE(ELOG[*,J] NE ZERO, nz_cnt)
   
   
          IF (zz_cnt NE 0) THEN BEGIN
             ELOG[e_zz,J]    = -30.0
          ENDIF
          IF (nz_cnt NE 0) THEN BEGIN
             ELOG[e_nz,J]    = ALOG10(ELOG[e_nz,J])
          ENDIF
          
          e_nz_max        = MAX(ELOG[*,J])
          e_nz_min        = MIN(ELOG[*,J])
   
       ENDFOR
 
     this_time            = time
     str_time             = STRTRIM(this_time,2)
;
;    FOR I = first_step + 1, last_step DO BEGIN
;
;      ELOG_next        = open_spec_data_file( dsc_lab, L, I)
;      ENLG_next[*,*]   = ELOG_next[*,*]
;
;         FOR J = 1, 6 DO BEGIN
;
;            e_zz       = WHERE(ELOG_next[*,J] EQ ZERO, zz_cnt)
;            e_nz       = WHERE(ELOG_next[*,J] NE ZERO, nz_cnt)
;
;            IF (zz_cnt NE 0) THEN BEGIN
;               ELOG_next[e_zz,J] = -30.0
;            ENDIF
;            IF (nz_cnt NE 0) THEN BEGIN
;               ELOG_next[e_nz,J] = ALOG10(ELOG_next[e_nz,J])
;            ENDIF
;
;         ENDFOR
;
;         FOR J = 0, 6 DO BEGIN
;           ELOG[*,J]    = ELOG[*, J] + ELOG_next[*, J]
;           ENLG[*,J]    = ENLG[*, J] + ENLG_next[*, J]
;         ENDFOR
;
;    ENDFOR
;
;    tlast              = time
;
;    str_dt             = STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)
;
;    ELOG[*,*]             = ELOG[*,*] / tot_steps
;    ENLG[*,*]             = ENLG[*,*] / tot_steps
;
      FOR J = 1, 6 DO BEGIN
        ELOG[*,J]          = 10^ELOG[*,J]
      ENDFOR
;
;;   data_out_file      = cur_dir + '/ra_spec/' + 'ra_' + str_srs_ave + '_' + sfld + '_L-' + str_zl + '_lyr-' + str_layer + '.dat'
;;   OpenW, data_unit, data_out_file, /GET_LUN
;;   FOR K = 0, n_lines - 1 DO BEGIN
;;      PRINTF, data_unit,  E[K, 0:6], FORMAT='(7(E16.8,1x),:/)'
;;   ENDFOR
;;   FREE_LUN, data_unit
;
;;;; merged from sp_vs_z ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
     DK                 = DBLARR(n_lines-1)
;
     ELOG               = TRANSPOSE(ELOG)
;    ENLG               = TRANSPOSE(ENLG)
;
     FOR I = 0, n_lines - 2 DO BEGIN
 
       DK[I]            = ELOG[i_k,I+1] - ELOG[i_k, I]
 
     ENDFOR
 
     dk                 = MEAN(DK)
 
     min_x              = MIN(ELOG[i_k,   *])
     max_x              = MAX(ELOG[i_k,   *])
;
     min_y              = MIN(ELOG[i_sfld,*])
     max_y              = MAX(ELOG[i_sfld,*])
;
     IF (L GT n3) THEN BEGIN
       is_layer_even    = (L MOD 2) EQ 0
       IF (is_layer_even) THEN BEGIN
         n_proc         = CEIL(UINT(L/n3)-1)
       ENDIF ELSE BEGIN
         n_proc         = FLOOR(UINT(L/n3))
      ENDELSE
     ENDIF ELSE BEGIN
      n_proc            = 0
     ENDELSE
;
     zcoord             =    dz*L
;
     str_zpos           = STRTRIM(zcoord,2)
     str_zpos           = STRMID(str_zpos,0,5)
;    loc_layer          = UINT(L/(n_proc+1))
;    str_proc           = STRTRIM(n_proc,2)
;    len_str_proc       = STRLEN(str_proc)
;
     WHILE (STRLEN(str_proc) LT 3) DO BEGIN
       str_proc         = '0' + str_proc
     ENDWHILE
;
;    str_loc_layer      = STRTRIM(loc_layer,2)
;    WHILE (STRLEN(str_loc_layer) LT 3) DO BEGIN
;      str_loc_layer    = '0' + str_loc_layer
;    ENDWHILE
;
     case_file_str      =  str_sfld   + "spec_" + str_x_res + '_' + str_z_res + '_zl_' + str_zl  $
                         + '_proc-'   + str_proc   + '_layer-' + str_layer + '_step-' + str_step
;
     eps_out            = eps_out_dir + '/'        + case_file_str   + '.eps'
     PRINT, 'eps_out = ', eps_out
;    str_ave            = ' Averaged over '        + str_dt
     str_ps             = ' Power Spectrum for layer ' + str_layer + ' of ' + str_tot_layers
     str_line_two       = 'at z = ' + str_zpos + ' (L = ' + str_zl + '),' + ' and t = ' + str_time
     aumlaut            = STRING(228B)
     str_els            = ' Els' + aumlaut + 'sser '
     str_elsp           = 'Z!E+!N Els' + aumlaut + 'sser' 
     str_elsm           = 'Z!E-!N Els' + aumlaut + 'sser'
  
     CASE i_sfld OF
     '1' : str_title    = 'Kinetic Energy'             + str_ps + '!C' + str_line_two
     '2' : str_title    = 'Magnetic Energy'            + str_ps + '!C' + str_line_two
     '3' : str_title    = 'Total Energy'               + str_ps + '!C' + str_line_two
     '4' : str_title    = str_elsp                     + str_ps + '!C' + str_line_two
     '5' : str_title    = str_elsm                     + str_ps + '!C' + str_line_two
     '6' : str_title    = 'Total' + str_els + 'Energy' + str_ps + '!C' + str_line_two
     ELSE: str_title    = 'Something is terribly wrong'
     ENDCASE
  
     k_drive            = 12.0
     x_rng              = [k_drive, max_x]
;    y_rng              = [min_y, max_y]
     y_rng              = [y_global_min, y_global_max]
;    y_rng              = [min_y, y_global_max]
 
     sy_offset_a        = 000000.0D
     sy_offset_b        = 000000.0D
 
     PLAWA              = DBLARR(n_lines)
;    PLAWB              = DBLARR(n_lines)
  
     PLAWA[*]           = max_y * ( ELOG[i_k,*]^(-five_thirds) )
  
     !P.FONT = -1
 
     k_low              = 20.0
     k_high             = 2.0 * !pi * 128.0/12.0
     
     k_drive_idx        = WHERE((ELOG[i_k, *] GE k_drive - half*dk) AND (ELOG[i_k, *] LE k_drive + half*dk),  low_count )
     k_drive            = ELOG[i_k, k_drive_idx]
 
     k_low_idx          = WHERE((ELOG[i_k, *] GE k_low   - half*dk) AND (ELOG[i_k, *] LE k_low   + half*dk),  low_count )
     k_high_idx         = WHERE((ELOG[i_k, *] GE k_high  - half*dk) AND (ELOG[i_k, *] LE k_high  + half*dk), high_count )
     
;;   PRINT, "k_drive_idx          = ", k_drive_idx
;;   PRINT, "k_low_idx            = ", k_low_idx
;;   PRINT, "k_high_idx           = ", k_high_idx
;
;;   PRINT, "K[", k_drive_idx, "] = ", ELOG[i_k, k_drive_idx]
;;   PRINT, "K[", k_drive_idx, "] = ", ENLG[i_k, k_drive_idx]
;;   PRINT, "K[", k_low_idx,   "] = ", ELOG[i_k, k_low_idx  ]
;;   PRINT, "K[", k_low_idx,   "] = ", ENLG[i_k, k_low_idx  ]
;
     fit_range          = k_high_idx[0] - k_low_idx[0] + 1
 
;    CHKEDIF            = DBLARR(fit_range)
;    CHKEDIF            = ((ENLG[i_sfld, k_low_idx[0]:k_high_idx[0]] - ELOG[i_sfld, k_low_idx[0]:k_high_idx[0]])/ELOG[i_sfld,k_low_idx[0]:k_high_idx[0]])
;
;    FOR I_I = k_low_idx[0], k_high_idx[0] DO BEGIN
;      PRINT, "ENLG[", i_sfld, ",", I_I,"] = "  , ENLG[i_sfld, I_I], " ELOG[", i_sfld, ",", I_I,"] = ", ELOG[i_sfld, I_I] , " CHKEDIF[",I_I,"] = ", CHKEDIF[I_I - k_low_idx[0]]
;    ENDFOR
;
;    max_edif           = MAX(CHKEDIF, SUBSCRIPT_MIN=loc_max_edif )
;    min_edif           = MIN(CHKEDIF, SUBSCRIPT_MAX=loc_min_edif )
;
;    PRINT, 'max_edif     = ', max_edif
;    PRINT, 'loc_max_edif = ', loc_max_edif
;    PRINT, 'min_edif     = ', min_edif
;    PRINT, 'loc_min_edif = ', loc_min_edif
;
;    KXNLG              = DBLARR(fit_range)
;    KYNLG              = DBLARR(fit_range)
     KXLOG              = DBLARR(fit_range)
     KYLOG              = DBLARR(fit_range)
;    KYERR              = DBLARR(fit_range)
;    
;    KXNLG[*]              = ALOG10(ENLG[i_k,   k_low_idx:k_high_idx])
;    KYNLG[*]              = ALOG10(ENLG[i_sfld,k_low_idx:k_high_idx])
;
;    FOR I_I = 0, fit_range - 1 DO BEGIN
;      PRINT, "KYNLG[", I_I, "] = " , 10^KYNLG[I_I]
;    ENDFOR
;
;    lfit_nlg           = LINFIT(KXNLG, KYNLG, CHISQR = chi_sqr)
;
     KXLOG[*]           = ALOG10(ELOG[i_k,   k_low_idx:k_high_idx])
     KYLOG[*]           = ALOG10(ELOG[i_sfld,k_low_idx:k_high_idx])
;
;    FOR I_I = 0, fit_range - 1 DO BEGIN
;      PRINT, "KYLOG[", I_I,"] = ", 10^KYLOG[I_I]
;    ENDFOR
;
     lfit               = LINFIT(KXLOG, KYLOG, CHISQR = chi_sqr)
 
     PRINT, "using elog: index = ", lfit[1]
 ;   PRINT, "using enlg: index = ", lfit_nlg[1]
;
;    z_coord_idx        = WHERE((NCHVSZ[*,0] GE zcoord  - half*dz) AND (NCHVSZ[*,0] LE zcoord  + half*dz),  z_count )
;    size_zci           = SIZE(z_coord_idx,/DIMENSIONS)
;    
;    SPIDX[L-first_layer,0] = zcoord    ; z coordinate for layer L
;    SPIDX[L-first_layer,1] = lfit[1]   ; spectral index
;    SPIDX[L-first_layer,2] = lfit[0]   ; intercept
;    SPIDX[L-first_layer,3] = chi_sqr   ; chi square for fit
;
;    IF (size_zci[0] EQ 1) THEN BEGIN
;      SPIDX[L-first_layer,4] = NCHVSZ[z_coord_idx[0],1]      ; will contain normalized cross-helicity.
;    ENDIF ELSE BEGIN
;      SPIDX[L-first_layer,4] = ZERO
;      PRINT, 'ra_spec: too many z coordinates for layer ', L
;    ENDELSE
;
;;   PRINTF,  spidx_unit, SPIDX[L-first_layer,0:4], FORMAT='(5(E16.8,1x),:/)'
;    
;; in-loop plotting section begins here ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
     IF (STRCMP(plot_mode, "spectra", 6, /FOLD_CASE) EQ 1) THEN BEGIN
 
        PRINT, "plotting spectrum to file ", eps_out
 
        SET_PLOT, 'PS' 
        DEVICE, /COLOR, /ENCAPSULATED
        DEVICE, FILENAME   = eps_out
 
        PLOT, ELOG[i_k,*], ELOG[i_sfld,*],                          $
              CHARSIZE     = 1.0,                                   $
              LINESTYLE    = 0,                                     $
              YTICKFORMAT  = '(E8.1)',                              $
              XMARGIN      = [12,4],                                $
              YMARGIN      = [4,4],                                 $
              XRANGE       = x_rng,                                 $
              YRANGE       = y_rng,                                 $
              XTITLE       = 'k',                                   $
              YTITLE       = str_y_title,                           $
              THICK        = 2,                                     $
              CLIP         = [k_drive, y_global_min, max_x, y_global_max],        $
              /XLOG,                                                $
              /YLOG,                                                $
              TITLE        = str_title
  
        OPLOT, ELOG[i_k,*], PLAWA[*],                               $
              LINESTYLE    = 2,                                     $
              THICK        = 3
  
        k_high             = ELOG[i_k, k_high_idx]
        k_low              = ELOG[i_k, k_low_idx ]
 
        OPLOT, [k_low,  k_low],  [y_global_min, y_global_max],                    $
                LINESTYLE  = 2,                                     $
                THICK      = 4
        
        OPLOT, [k_high, k_high], [y_global_min, y_global_max],                    $
                LINESTYLE  = 2,                                     $
                THICK      = 4
 
        LOGKXX             = DBLARR(n_lines)
        KXX                = DBLARR(n_lines)
        KYY                = DBLARR(n_lines)
 
        LOGKXX             = ALOG10(ELOG[i_k,*])
        KYY[*]             = lfit[0] + lfit[1]*LOGKXX[*]
        
        KXX[*]             = 10^(LOGKXX[*])
        KYY[*]             = 10^(KYY[*])
        
        KXLOG[*]           = 10^KXLOG[*]
        KYLOG[*]           = 10^KYLOG[*]
 
        OPLOT, KXLOG, KYLOG,                      $
               COLOR       = 254      ,           $
               SYMSIZE     = 1.5,                 $
;              SYMSIZE     = 3.0,                 $
               PSYM        = 4
;              PSYM        = 6
 
        OPLOT, KXX, KYY,                          $
               COLOR       = 254      ,           $
               LINESTYLE   = 1,                   $
               THICK       = 3
;
;       LOGKXX             = ALOG10(ENLG[i_k,*])
;       KYY[*]             = lfit_nlg[0] + lfit_nlg[1]*LOGKXX[*]
;       KYY[*]             = 10^(KYY[*])
;
;       KXNLG[*]           = 10^KXNLG[*]
;       KYNLG[*]           = 10^KYNLG[*]
;
;       OPLOT, KXNLG, KYNLG,                      $
;              COLOR       = 80       ,           $
;              SYMSIZE     = 1.5,                 $
;              PSYM        = 4
;
;       OPLOT, KXX, KYY,                          $
;              COLOR       = 80       ,           $
;              LINESTYLE   = 1,                   $
;              THICK       = 3
;
        str_adx            = STRTRIM(lfit[1], 2)
        str_adx            = STRMID(str_adx,0,5)
        str_mn_five_thirds = '-5/3'
;
        XS = DOUBLE(2)
        YS = DOUBLE(2)
        XS = [ ELOG[i_k, 80], ELOG[i_k,135] ]
        YS = [ KYY[0],  KYY[0] ]

        OPLOT, XS, YS,       $
        COLOR      = 254,    $
        THICK      = 3,      $
        LINESTYLE  = 1

        XYOUTS, ELOG[i_k,140], KYY[0],   'E(k)!9?!Xk!U' + str_adx + '!N', CHARSIZE=1.2

        XK = DOUBLE(2)
        YK = DOUBLE(2)
        XK = [ ELOG[i_k, 80], ELOG[i_k,135] ]
        YK = [ PLAWA[3],  PLAWA[3] ]

        OPLOT, XK, YK,       $
        COLOR      = 0,      $
        THICK      = 3,      $
        LINESTYLE  = 2

        XYOUTS, ELOG[i_k,140], PLAWA[3], 'E(k)!9?!Xk!U' + str_mn_five_thirds + '!N', CHARSIZE=1.2
;       
        DEVICE, /CLOSE
             
        SET_PLOT, 'X'
   
     ENDIF 
 
    ENDFOR  ; layer
  ENDFOR    ; step
;; in-loop plotting section ends here ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;  IF ((STRCMP(plot_mode, "idxvsz", 6, /FOLD_CASE) EQ 1) OR (STRCMP(plot_mode, "idxvsnch", 6, /FOLD_CASE) EQ 1)) THEN BEGIN
;
;    i_y   = 1
;    IF (STRCMP(plot_mode, "idxvsz",     6, /FOLD_CASE) EQ 1) THEN BEGIN
;      sp_vs_z_eps_out = data_out_dir +  sfld + '/' + 'spidx_vs_z_' +  sfld + '_' + str_srs_ave + '.eps'
;      PRINT, "plotting spectral index vs z to file ", sp_vs_z_eps_out
;      i_x = 0
;      str_x_title = 'z'
;      str_title   = 'Spectral Index (a) vs. z'
;      
;    ENDIF ELSE BEGIN
;      IF (STRCMP(plot_mode, "idxvsnch", 6, /FOLD_CASE) EQ 1) THEN BEGIN
;      sp_vs_z_eps_out = data_out_dir +  sfld + '/' + 'spidx_vs_nch_' +  sfld + '_' + str_srs_ave + '.eps'
;      PRINT, "plotting spectral index vs normalized cross helicity to file ", sp_vs_z_eps_out
;      i_x = 4
;      str_x_title = '!7r!X'
;      str_title   = 'Spectral Index (a) vs. Normalized Cross-Helicity (!7r!X)'
;      ENDIF
;    ENDELSE
;
;    str_y_title  = 'a'
;
;
;    x_max = MAX(SPIDX[*,i_x])
;    x_min = MIN(SPIDX[*,i_x])
;
;    y_max = MAX(SPIDX[*,i_y])
;    y_min = MIN(SPIDX[*,i_y])
;
;    x_rng = [x_min, x_max]
;    y_rng = [y_min, y_max]
;
;    SET_PLOT, 'PS'
;    DEVICE, /COLOR, /ENCAPSULATED
;    DEVICE, FILENAME   = sp_vs_z_eps_out
;
;    PLOT, SPIDX[*,i_x], SPIDX[*,i_y],         $
;          CHARSIZE     = 1.2,                 $
;          LINESTYLE    = 0,                   $
;          YTICKFORMAT  = '(E8.1)',            $
;          XMARGIN      = [12,4],              $
;          YMARGIN      = [4,4],               $
;          XRANGE       = x_rng,               $
;          YRANGE       = y_rng,               $
;          XTITLE       = str_x_title,         $
;          YTITLE       = str_y_title,         $
;          THICK        = 2,                   $
;          TITLE        = str_title
;
;    DEVICE, /CLOSE
;            
;    SET_PLOT, 'X'
;
;  ENDIF
;
;; FREE_LUN, spidx_unit
;
  PRINT, "y_global_max = ", y_global_max
  PRINT, "y_global_min = ", y_global_min
 END
