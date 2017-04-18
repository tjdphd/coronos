FUNCTION plotSpectrum, desc_label,                   $
                       sfld,                         $
                       str_dt,                       $
                       k_data,                       $
                       ENLG, KXNLG, KYNLG, lfit_nlg, $
                       ELOG, KXLOG, KYLOG, lfit,     $
                       L,                            $
                       slc_minmax_y,                 $
                       first_step, last_step,        $
                       plot_mode

  ZERO               = 0.0D
  HALF               = 0.5D
  five_thirds        = 5.0D/3.0D
  three_halves       = 3.0D/2.0D

  case_file_str      = getCaseFileString(desc_label, sfld, L, first_step, last_step)
  data_out_dir       = getDataOutDir(plot_mode)
  eps_out            = data_out_dir + sfld + '/eps/' + case_file_str + '.eps'
  str_title          = getSpectrumTitle(sfld, L, str_dt)

  PRINT, "creating ", eps_out

  SET_PLOT, 'PS'
  DEVICE, /COLOR, /ENCAPSULATED
  DEVICE, FILENAME   = eps_out

  i_k                = 0
  i_sfld             = getSfldIndex(sfld)

  k_drive            = k_data.k_drive

  min_x              = k_drive
  min_y              = 1.0e-15
  max_x              = MAX(ELOG[i_k,    *])
  max_y              = slc_minmax_y[1]

  x_rng              = [min_x, max_x]
  y_rng              = [min_y, max_y]


  PLOT, ENLG[i_k,*], ENLG[i_sfld,*],                                                                              $
        CHARSIZE     = 0.8,                                                                                       $
        LINESTYLE    = 0,                                                                                         $
        YTICKFORMAT  = '(E8.1)',                                                                                  $
        XMARGIN      = [12,4],                                                                                    $
        YMARGIN      = [4, 6],                                                                                    $
        XRANGE       = x_rng,                                                                                     $
        YRANGE       = y_rng,                                                                                     $
        XTITLE       = 'k',                                                                                       $
        YTITLE       = str_y_title,                                                                               $
        THICK        = 2,                                                                                         $
        CLIP         = [k_drive, min_y, max_x, max_y],                                                            $
        /XLOG,                                                                                                    $
        /YLOG,                                                                                                    $
        TITLE        = str_title

  size_ELOG          = SIZE(ELOG,/DIMENSIONS)
  n_lines            = size_ELOG[1]
  PLAWA              = DBLARR(n_lines)
  PLAWA[*]           = max_y * ( ELOG[i_k,*]^(-five_thirds) )

  OPLOT, ELOG[i_k,*], PLAWA[*],                                                                                   $
          COLOR      = 254,                                                                                       $
         LINESTYLE   = 2,                                                                                         $
         THICK       = 3

  k_high             = k_data.k_high
  k_low              = k_data.k_low 

  OPLOT, [k_low,  k_low],  [min_y, max_y],                                                                        $
          LINESTYLE  = 2,                                                                                         $
          THICK      = 4
      
  OPLOT, [k_high, k_high], [min_y, max_y],                                                                        $
          LINESTYLE  = 2,                                                                                         $
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

  OPLOT, KXLOG, KYLOG,                                                                                            $
         COLOR       = 80       ,                                                                                 $
         SYMSIZE     = 3.0,                                                                                       $
         PSYM        = 6

  OPLOT, KXX, KYY,                                                                                                $
         COLOR       = 254      ,                                                                                 $
         LINESTYLE   = 1,                                                                                         $
         THICK       = 3

  LOGKXX             = ALOG10(ENLG[i_k,*])
  KYY[*]             = lfit_nlg[0] + lfit_nlg[1]*LOGKXX[*]
  KYY[*]             = 10^(KYY[*])

  KXNLG[*]           = 10^KXNLG[*]
  KYNLG[*]           = 10^KYNLG[*]

  OPLOT, KXNLG, KYNLG,                                                                                            $
         COLOR       = 80       ,                                                                                 $
         SYMSIZE     = 1.5,                                                                                       $
         PSYM        = 4

  OPLOT, KXX, KYY,                                                                                                $
         COLOR       = 80       ,                                                                                 $
         LINESTYLE   = 1,                                                                                         $
         THICK       = 3

  str_adx            = STRTRIM(lfit[1], 2)
  str_adx            = STRMID(str_adx,0,5)
  str_mn_five_thirds = '-5/3'

  LGND_empX          = [600.0,   1000.0 ]
  LGND_empY          = [1.0e-04, 1.0e-04]

  OPLOT, LGND_empX, LGND_empY,                                                                                    $
         COLOR       = 80,                                                                                        $
         LINESTYLE   = 1,                                                                                         $
         THICK       = 3

  XYOUTS, 1100.0, 1.0e-04, 'E(k)!9?!Xk!U' + str_adx + '!N'

  LGND_kolX          = [600.0,   1000.0 ]
  LGND_kolY          = [1.0e-05, 1.0e-05]

  OPLOT, LGND_kolX, LGND_kolY,                                                                                    $
         COLOR       = 254,                                                                                       $
         LINESTYLE   = 2,                                                                                         $
         THICK       = 4

  XYOUTS, 1100.0, 1.0e-05, 'E(k)!9?!Xk!U' + str_mn_five_thirds + '!N'
  
  DEVICE, /CLOSE
        
  SET_PLOT, 'X'

  RETURN, 0
END
