FUNCTION calc_extE, E, desc_lab, step, i_efld, KK, GOFX

  ZERO                         = 0.0D                                      ; Define some constants for convenience
  ONE                          = 1.0D                                      ; everywhere
  TWO                          = 2.0D
  THREE                        = 3.0D
  PI                           = 3.141592653589793D
  TWO_THIRDS                   = TWO / THREE
  TWO_PI                       = TWO*PI

  CZERO                        = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  CI                           = DCOMPLEX(ZERO,ONE)                        ; and to make sure that DOUBLE's are used

  i_z                          = 0 ; z - coordinate
  i_me                         = 1 ; megnetic energy in slab at z
  i_pe                         = 2 ; kinetic energy in slab at z
  i_ch                         = 3 ; cross helicity in slab at z
  i_ep                         = 4 ; elsasser energy
  i_em                         = 5 ; elsasser energy
  i_ce                         = 6 ; int (J^2 dz)
  i_oe                         = 7 ; int (Omega^2 dz)
  i_zp                         = 8 ; int( ( J + Omega)^2 dz)
  i_zm                         = 9 ; int( ( J - Omega)^2 dz)
  i_nch                        = 10 ; normalized cross helicity
  i_km                         = 11 ; square - "magnetic" wave number 
  i_kp                         = 12 ; square - "kinetic"  wave number 
  i_kzp                        = 13 ; square - "Z+" wave number
  i_kzm                        = 14 ; square - "Z+" wave number

  i_te                         = 15 ; total energy at z
  i_re                         = 16 ; residual energy at z
  i_ne                         = 17 ; normalized residual energy

  i_rzp                        = 18 
  i_rzm                        = 19 

  i_zpm                        = 20 ; Z+Z-  where  Z+ = |(Z+)^2|^(1/2) = |ep^2|^(1/2) and similarly for Z- (em).

  i_likp                       = 21 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm                       = 22 ; IK length length scale  lambda_IK,-  (see notes)

  i_lkop                       = 23 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom                       = 24 ; Kolmogorov length scale lambda_Kol,- (see notes)
  i_lpl                        = 25 ; length scale from F.T. of z^+
  i_lmi                        = 26 ; length scale from F.T. of z^+
  i_lpla                       = 27 ; alternative definition for lpl
  i_lmia                       = 28 ; alternative definition for lmi
  i_lsfp                       = 29 ; length scale from structure function for lambda +
  i_lsfm                       = 30 ; length scale from structure function for lambda -
  
  ip1                          = scan_parameters('ip1', 0, desc_lab )      ; power of 2 giving x-resolution
  ip2                          = scan_parameters('ip2', 0, desc_lab )      ; power of 2 giving y-resolution
  n3                           = scan_parameters('n3' , 0, desc_lab ) 
  mp                           = scan_parameters('mp' , 0, desc_lab ) 
  zl                           = scan_parameters('zl' , 0, desc_lab ) 
  
  n1                           = 2^ip1 
  n2                           = 2^ip2
  
  dx                           = ONE / FLOAT(n1)
  dy                           = ONE / FLOAT(n2)
  dz                           = zl / (n3 * mp)

  size_e                       = SIZE(E)

  extE                         = FLTARR(size_e[1], size_e[2] + (i_lsfm - i_kzm) )
  extE[*, i_z:i_kzm]           =  E[*,i_z:i_kzm]                           ; copy E to 0:14 of Eext
  
  extE[*,   i_te ]             =  E[*, 1] + E[*, 2]                        ; total energy               i_te  = 15
  extE[*,   i_re ]             =  E[*, 2] - E[*, 1]                        ; residual energy            i_re  = 16
  extE[*,   i_ne ]             = (E[*, 2] - E[*, 1]) / (E[*, 1] + E[*, 2]) ; normalize residual energy  i_ne  = 17

  extE[*,   i_rzp]             = SQRT(E[*,4])                              ; sqrt(Z+^2)                 i_rzp = 18
  extE[*,   i_rzm]             = SQRT(E[*,5])                              ; sqrt(Z-^2)                 i_rzm = 19

  extE[*,   i_zpm]             = extE[*,i_rzp] * extE[*,i_rzm]             ; sqrt(Z+^2Z-^2)             i_zpm = 20
  
  extE[*, i_likp:i_lsfm]       = ZERO                                     ; set the rest to zero for now


  IF ( (i_efld GE i_lpl) && (i_efld LE i_lsfm) ) THEN BEGIN

    n_layers                   = size_e[1]
    multithread                = 1

    IF (multithread) THEN BEGIN

      slice_range              = getSliceRange(1,   n_layers  )
      size_sr                  = SIZE(slice_range, /DIMENSIONS)

      n_threads                = size_sr[0]
      
      PRINT, 'calc_extE: n_threads = ', n_threads
      oBridge                  = OBJARR(n_threads - 1)

      FOR J = 0, n_threads - 1 DO BEGIN

        first_layer            = slice_range[J,0]
        last_layer             = slice_range[J,1]
        n_layers_this_thread   = last_layer - first_layer + 1

        IF (J EQ n_threads - 1 ) THEN BEGIN

          thread_extE_for_lpm  = calc_lpm(desc_lab, step, n1, n2, first_layer, last_layer, KK, GOFX) 
;         local_extE_for_lpm   = FLTARR(n_layers_this_thread, 2) ; should not be necessary
          local_extE_for_lpm   = thread_extE_for_lpm

        ENDIF ELSE BEGIN

          oBridge[J]           = OBJ_NEW('IDL_IDLBRIDGE')
          oBridge[J]           -> setVar, 'desc_lab',    desc_lab
          oBridge[J]           -> setVar, 'step',        step
          oBridge[J]           -> setVar, 'n1',          n1
          oBridge[J]           -> setVar, 'n2',          n2
          oBridge[J]           -> setVar, 'first_layer', slice_range[J,0]
          oBridge[J]           -> setVar, 'last_layer',  slice_range[J,1]
          oBridge[J]           -> setVar, 'KK',           KK
          oBridge[J]           -> setVar, 'GOFX',         GOFX
          oBridge[J]           -> Execute, ".run calc_lpm"
          oBridge[J]           -> Execute, $
   "thread_extE_for_lpm  = calc_lpm(desc_lab, step, n1, n2, first_layer, last_layer, KK,GOFX)",/nowait

        ENDELSE
      ENDFOR

      IF (n_threads GT 1) THEN BEGIN
        notdone                            = 1
        WHILE notdone DO BEGIN
          done                             = 0
          FOR J = 0, n_threads - 2 DO done = done + oBridge[J] -> Status()
          IF (done EQ 0 ) THEN notdone     = done
        ENDWHILE
      ENDIF

      FOR J = 0, n_threads - 1 DO BEGIN

        first_layer            = slice_range[J,0]
        last_layer             = slice_range[J,1]
        n_layers_this_thread   = last_layer - first_layer + 1

        IF (J EQ n_threads - 1) THEN BEGIN

          extE[(first_layer-1):(last_layer-1), i_rzp  : i_rzm]  = local_extE_for_lpm[*,0:1]
          extE[(first_layer-1):(last_layer-1), i_lpl  : i_lmi]  = local_extE_for_lpm[*,2:3]
          extE[(first_layer-1):(last_layer-1), i_lpla : i_lmia] = local_extE_for_lpm[*,4:5]
          extE[(first_layer-1):(last_layer-1), i_lsfp : i_lsfm] = local_extE_for_lpm[*,6:7]

        ENDIF ELSE BEGIN
       
          thread_extE_for_lpm                                   = oBridge[J] -> GetVar('thread_extE_for_lpm')
          extE[(first_layer-1):(last_layer-1), i_rzp  : i_rzm]  = thread_extE_for_lpm[*,0:1]
          extE[(first_layer-1):(last_layer-1), i_lpl  : i_lmi]  = thread_extE_for_lpm[*,2:3]
          extE[(first_layer-1):(last_layer-1), i_lpla : i_lmia] = thread_extE_for_lpm[*,4:5]
          extE[(first_layer-1):(last_layer-1), i_lsfp : i_lsfm] = thread_extE_for_lpm[*,6:7]
    
          obj_destroy, oBridge[J]

        ENDELSE

      ENDFOR 

    ENDIF ELSE BEGIN

      FOR  I = 0, n_layers - 1 DO BEGIN
      
        SLB                    = fetch_layer( I+1, step, desc_lab)
      
        P                      = REFORM(SLB[*,0], n1, n2)
        A                      = REFORM(SLB[*,1], n1, n2)
      
        FTP                    = FFT(P)
        FTA                    = FFT(A)

        FT_P_PLS_A             = FTP  + FTA
        FT_P_MIN_A             = FTP  - FTA

        FTZ_PLS                = FT_P_PLS_A * KK
        FTZ_MNS                = FT_P_MIN_A * KK

        SQR_MOD_PPA            = (ABS(FT_P_PLS_A))^2
        SQR_MOD_PMA            = (ABS(FT_P_MIN_A))^2

        SQR_MOD_FTZ_PLS        = (ABS(FTZ_PLS))^2
        SQR_MOD_FTZ_MNS        = (ABS(FTZ_MNS))^2

        KK[0,0]                = ONE

        extE[I, i_rzp]         = SQRT(TOTAL( SQR_MOD_FTZ_PLS[*,*]) )
        extE[I, i_rzm]         = SQRT(TOTAL( SQR_MOD_FTZ_MNS[*,*]) )

        extE[I, i_lpl]         = TOTAL( SQR_MOD_FTZ_PLS[*,*] / KK[*,*])
        extE[I, i_lmi]         = TOTAL( SQR_MOD_FTZ_MNS[*,*] / KK[*,*])

        extE[I, i_lpla]        = SQRT(TOTAL(SQR_MOD_PPA[*,*]) / (extE[I,i_rzp]^2))
        extE[I, i_lmia]        = SQRT(TOTAL(SQR_MOD_PMA[*,*]) / (extE[I,i_rzm]^2))
        
        extE[I, i_lsfp]        = TOTAL( (SQR_MOD_FTZ_PLS[*,*] * GOFX[*,*] ) / KK[*,*])
        extE[I, i_lsfm]        = TOTAL( (SQR_MOD_FTZ_MNS[*,*] * GOFX[*,*] ) / KK[*,*])

      ENDFOR

    ENDELSE

    extE[*,i_lpl]              =  PI * extE[*,i_lpl ] / (extE[*,i_rzp]^2)
    extE[*,i_lmi]              =  PI * extE[*,i_lmi ] / (extE[*,i_rzm]^2)

    extE[*,i_lsfp]             =  PI * extE[*,i_lsfp] / (extE[*,i_rzp]^2)
    extE[*,i_lsfm]             =  PI * extE[*,i_lsfm] / (extE[*,i_rzm]^2)

  ENDIF
  
  RETURN, extE

END
