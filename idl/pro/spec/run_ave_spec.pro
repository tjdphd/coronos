PRO run_ave_spec, dsc_lab, sfld, first_layer, last_layer, first_step, last_step, minmax_mode, plot_mode

  COMMON step, time

  ZERO                          = 0.0D
  HALF                          = 0.5D
  five_thirds                   = 5.0D/3.0D
  three_halves                  = 3.0D/2.0D
  EULER                         = 2.7182818284590452D

  i_k                           =  0 ; k-value
  i_pe                          =  1 ; kinetic        energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_ae                          =  2 ; magnetic       energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_ts                          =  3 ; total          energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_zp                          =  4 ; elsasser+      energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_zm                          =  5 ; elsasser-      energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_tz                          =  6 ; total elsasser energy spectrum in layer ? of time step ? averaged in code before logarithm
  i_pef                         =  7 ; kinetic        energy spectrum in layer ? of time step ? averaged in code after  logarithm
  i_aef                         =  8 ; magnetic       energy spectrum in layer ? of time step ? averaged in code after  logarithm
  i_zpf                         =  9 ; elsasser+      energy spectrum in layer ? of time step ? averaged in code after  logarithm
  i_zmf                         = 10 ; elsasser-      energy spectrum in layer ? of time step ? averaged in code after  logarithm

  i_sfld                        = getSfldIndex(sfld)
  str_sfld                      = getSfldString(sfld)

  n3                            = scan_parameters('p3' , 0, dsc_lab )                     ; number of slices per data file
  mp                            = scan_parameters('np' , 0, dsc_lab )                     ; number of processors used in run
  zl                            = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  dz                            = zl / (n3 * mp)
  cur_dir                       = GETENV('PWD')

  tot_steps                     = FLOAT(last_step - first_step + 1)

  EVSZ                          = open_evsz_data_file(dsc_lab, first_step)

  size_evsz                     = SIZE(EVSZ, /DIMENSIONS)

  NCHVSZ                        = DBLARR(size_evsz[0],2)
  
  i_z                           = 0
  i_nch                         = 10

  NCHVSZ[*,0]                   = EVSZ[*,i_z]
  NCHVSZ[*,1]                   = EVSZ[*,i_nch]

  FOR S = first_step + 1, last_step DO BEGIN
    EVSZ                        = open_evsz_data_file(dsc_lab, S)
    NCHVSZ[*,1]                 = NCHVSZ[*,1] + EVSZ[*,i_nch]
  ENDFOR

  NCHVSZ[*,1]                   = NCHVSZ[*,1] / tot_steps

; layer dependencies

  SPIDX                         = DBLARR(last_layer - first_layer + 1, 5)
  data_out_dir                  = getDataOutDir(plot_mode)
  spec_idx_out_file             = getSpectralIndexDataOutFile(plot_mode, first_step, last_step)

; OpenW, spidx_unit, spec_idx_out_file, /GET_LUN

  IF (minmax_mode EQ  'global') THEN BEGIN
    slc_minmax_y                = global_sfld_minmax_over_slices(i_sfld, dsc_lab, first_layer, last_layer, last_step)
  ENDIF

  LOADCT, 13

  FOR L = first_layer, last_layer DO BEGIN
    
    ELOG                        = open_spec_data_file( dsc_lab, L, first_step)
    size_ELOG                   = SIZE(ELOG,/DIMENSIONS)

    n_lines                     = size_ELOG[0]
    n_cols                      = size_ELOG[1]

    ENLG                        = DBLARR(n_lines,n_cols)
    ENLG_next                   = DBLARR(n_lines,n_cols)

    ENLG[*,*]                   = ELOG[*,*]

    FOR J = 1, i_tz DO BEGIN

       e_zz                     = WHERE(ELOG[*,J] EQ ZERO, zz_cnt)
       e_nz                     = WHERE(ELOG[*,J] NE ZERO, nz_cnt)

       IF (zz_cnt NE 0) THEN BEGIN
          ELOG[e_zz,J]          = -30.0
       ENDIF
       IF (nz_cnt NE 0) THEN BEGIN
          ELOG[e_nz,J]          = ALOG10(ELOG[e_nz,J])
       ENDIF
       
       e_nz_max                 = MAX(ELOG[*,J])
       e_nz_min                 = MIN(ELOG[*,J])

    ENDFOR

    FOR J = i_tz + 1, i_zmf DO BEGIN

       ELOG[*,J]                = ALOG10(EULER) * ELOG[*,J]      ; adjust elogs in ELOG to base-ten logs
       ENLG[*,J]                = EULER^ENLG[*,J]                ; adjust elogs to "un" logged values in ENLG

    ENDFOR



    tfirst                      = time

    FOR I = first_step + 1, last_step DO BEGIN                   ; getting ELOG and ENLG starts here

      ELOG_next                 = open_spec_data_file( dsc_lab, L, I)
      ENLG_next[*,*]            = ELOG_next[*,*]

         FOR J = 1, i_tz DO BEGIN

            e_zz                = WHERE(ELOG_next[*,J] EQ ZERO, zz_cnt)
            e_nz                = WHERE(ELOG_next[*,J] NE ZERO, nz_cnt)

            IF (zz_cnt NE 0) THEN BEGIN
              ELOG_next[e_zz,J] = -30.0
            ENDIF
            IF (nz_cnt NE 0) THEN BEGIN
              ELOG_next[e_nz,J] = ALOG10(ELOG_next[e_nz,J])
            ENDIF

         ENDFOR

         FOR J = i_tz + 1, i_zmf DO BEGIN

           ELOG_next[*,J]       = ALOG10(EULER) * ELOG_next[*,J] ; adjust elogs in ELOG to base-ten logs
           ENLG_next[*,J]       = EULER^ENLG_next[*,J]           ; adjust elogs to "un" logged values in ENLG

         ENDFOR

         FOR J = 0, n_cols - 1 DO BEGIN

           ELOG[*,J]            = ELOG[*, J] + ELOG_next[*, J]
           ENLG[*,J]            = ENLG[*, J] + ENLG_next[*, J]

         ENDFOR


    ENDFOR

    tlast                       = time

    IF (tfirst EQ tlast) THEN BEGIN
      str_dt                    = ' for t = '      + STRTRIM(tfirst,2)
    ENDIF ELSE BEGIN
      str_dt                    = ' Average over ' + STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)
    ENDELSE

    ELOG[*,*]                   = ELOG[*,*] / tot_steps
    ENLG[*,*]                   = ENLG[*,*] / tot_steps

    FOR J = 1, n_cols - 1 DO BEGIN
      ELOG[*,J]                 = 10^ELOG[*,J]
    ENDFOR

    ELOG                        = TRANSPOSE(ELOG)
    ENLG                        = TRANSPOSE(ENLG)

; getting ELOG and ENLG ends here

    !P.FONT = -1

    k_data                      = setKData( 12.0, ELOG)

    k_drive                     = k_data.k_drive
    k_high_idx                  = k_data.k_high_idx
    k_low_idx                   = k_data.k_low_idx
    fit_range                   = k_high_idx[0] - k_low_idx[0] + 1

;   CHKEDIF                     = (                                                                                         $
;                                   ( ENLG[i_sfld, k_low_idx[0]:k_high_idx[0]] - ELOG[i_sfld, k_low_idx[0]:k_high_idx[0]] ) $
;                                   /                                                                                       $
;                                     ELOG[i_sfld,k_low_idx[0]:k_high_idx[0]]                                               $
;                                 )

;   max_edif                    = MAX(CHKEDIF, SUBSCRIPT_MIN=loc_max_edif )
;   min_edif                    = MIN(CHKEDIF, SUBSCRIPT_MAX=loc_min_edif )

    KXNLG                       = DBLARR(fit_range)
    KYNLG                       = DBLARR(fit_range)

    KXNLG[*]                    = ALOG10(ENLG[i_k,   k_low_idx:k_high_idx])
    KYNLG[*]                    = ALOG10(ENLG[i_sfld,k_low_idx:k_high_idx])

    lfit_nlg                    = LINFIT(KXNLG, KYNLG, CHISQR = chi_sqr)

    KXLOG                       = DBLARR(fit_range)
    KYLOG                       = DBLARR(fit_range)

    KXLOG[*]                    = ALOG10(ELOG[i_k,   k_low_idx:k_high_idx])
    KYLOG[*]                    = ALOG10(ELOG[i_sfld,k_low_idx:k_high_idx])

    lfit                        = LINFIT(KXLOG, KYLOG, CHISQR = chi_sqr)

    zcoord                      = dz*L
    z_coord_idx                 = WHERE((NCHVSZ[*,0] GE zcoord  - half*dz) AND (NCHVSZ[*,0] LE zcoord  + half*dz),  z_count )
    size_zci                    = SIZE(z_coord_idx,/DIMENSIONS)
    
    SPIDX[L-first_layer,0]      = zcoord                        ; z coordinate for layer L
    SPIDX[L-first_layer,1]      = lfit[1]                       ; spectral index
    SPIDX[L-first_layer,2]      = lfit[0]                       ; intercept
    SPIDX[L-first_layer,3]      = chi_sqr                       ; chi square for fit

    IF (size_zci[0] EQ 1) THEN BEGIN
      SPIDX[L-first_layer,4]    = NCHVSZ[z_coord_idx[0],1]      ; will contain normalized cross-helicity.
    ENDIF ELSE BEGIN
      SPIDX[L-first_layer,4]    = ZERO
      PRINT, 'ra_spec: too many z coordinates for layer ', L
    ENDELSE

    IF (STRCMP(plot_mode, "spectra", 6, /FOLD_CASE) EQ 1) THEN BEGIN

      dummy                     = plotSpectrum( desc_label,                                                                 $
                                                sfld,                                                                       $
                                                str_dt,                                                                     $
                                                k_data,                                                                     $
                                                ENLG, KXNLG, KYNLG, lfit_nlg,                                               $
                                                ELOG, KXLOG, KYLOG, lfit,                                                   $
                                                L,                                                                          $
                                                slc_minmax_y,                                                               $
                                                first_step, last_step,                                                      $
                                                plot_mode                                                                   $
                                              )

    ENDIF 

  ENDFOR

  IF ( (STRCMP(plot_mode, "idxvsz",   6, /FOLD_CASE) EQ 1)                                                                  $
  OR   (STRCMP(plot_mode, "idxvsnch", 6, /FOLD_CASE) EQ 1)                                                                  $
     )                                                                                                                      $
  THEN BEGIN

    dummy                         = plotSpVsZ( desc_label,                                                                  $
                                               sfld,                                                                        $
                                               str_dt,                                                                      $
                                               SPIDX,                                                                       $
                                               first_step, last_step,                                                       $
                                               plot_mode                                                                    $
                                             )

  ENDIF

; FREE_LUN, spidx_unit

END
