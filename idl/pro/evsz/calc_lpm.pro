FUNCTION calc_lpm, desc_lab, step, n1, n2, first_layer, last_layer, KK, GOFX
  
  ZERO                       = 0.0D                                      ; Define some constants for convenience
  CZERO                      = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  ONE                        = 1.0D                                      ; everywhere
  TWO                        = 2.0D
  THREE                      = 3.0D
  PI                         = 3.141592653589793D
  TWO_THIRDS                 = TWO / THREE
  
  TWO_PI                     = TWO*PI

  n_layers_this_thread       = last_layer - first_layer + 1

  thread_extE_for_lpm        = FLTARR(n_layers_this_thread, 8)

  FOR  I = 0, n_layers_this_thread - 1 DO BEGIN
   
    SLB                      = fetch_layer( I + first_layer, step, desc_lab)
   
    P                        = REFORM(SLB[*,0], n1, n2)
    A                        = REFORM(SLB[*,1], n1, n2)
   
    FTP                      = FFT(P)
    FTA                      = FFT(A)

    FT_P_PLS_A               = FTP  + FTA
    FT_P_MIN_A               = FTP  - FTA

    FTZ_PLS                  = FT_P_PLS_A * KK
    FTZ_MNS                  = FT_P_MIN_A * KK

    SQR_MOD_PPA              = (ABS(FT_P_PLS_A))^2
    SQR_MOD_PMA              = (ABS(FT_P_MIN_A))^2

    SQR_MOD_FTZ_PLS          = (ABS(FTZ_PLS))^2
    SQR_MOD_FTZ_MNS          = (ABS(FTZ_MNS))^2

    KK[0,0]                  = ONE

    thread_extE_for_lpm[I,0] = SQRT(TOTAL( SQR_MOD_FTZ_PLS[*,*]) )                            ; placed in extE[I,i_rzp]
    thread_extE_for_lpm[I,1] = SQRT(TOTAL( SQR_MOD_FTZ_MNS[*,*]) )                            ; placed in extE[I,i_rzm]

    thread_extE_for_lpm[I,2] = TOTAL( SQR_MOD_FTZ_PLS[*,*] / KK[*,*])                         ; placed in extE[I,i_lpl]
    thread_extE_for_lpm[I,3] = TOTAL( SQR_MOD_FTZ_MNS[*,*] / KK[*,*])                         ; placed in extE[I,i_lmi]

    thread_extE_for_lpm[I,4] = SQRT(TOTAL(SQR_MOD_PPA[*,*]) / (thread_extE_for_lpm[I,0]^2))   ; placed in extE[I,i_lpla]
    thread_extE_for_lpm[I,5] = SQRT(TOTAL(SQR_MOD_PMA[*,*]) / (thread_extE_for_lpm[I,1]^2))   ; placed in extE[I,i_lmia]

    thread_extE_for_lpm[I,6] = TOTAL( (SQR_MOD_FTZ_PLS[*,*] * GOFX[*,*]) / KK[*,*])           ; placed in extE[I,i_lsfp]
    thread_extE_for_lpm[I,7] = TOTAL( (SQR_MOD_FTZ_MNS[*,*] * GOFX[*,*]) / KK[*,*])           ; placed in extE[I,i_lsfm]

  ENDFOR

    RETURN, thread_extE_for_lpm

END
