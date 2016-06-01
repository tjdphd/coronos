FUNCTION calc_extE, E, desc_lab, step, i_efld

  ZERO                     = 0.0D                                      ; Define some constants for convenience
  CZERO                    = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  ONE                      = 1.0D                                      ; everywhere
  TWO                      = 2.0D
  THREE                    = 3.0D
  PI                       = 3.141592653589793D
  TWO_THIRDS               = TWO / THREE
  
  TWO_PI                   = TWO*PI
  
  ip1                      = scan_parameters('ip1', 0, desc_lab )      ; power of 2 giving x-resolution
  ip2                      = scan_parameters('ip2', 0, desc_lab )      ; power of 2 giving y-resolution
  
  n1                       = 2^ip1 
  n2                       = 2^ip2
  
  
  dx                       = ONE / FLOAT(n1)
  dy                       = ONE / FLOAT(n2)
  
  size_e                   = SIZE(E)

  extE                     = FLTARR(size_e[1], size_e[2] + 10)
  
  extE[*, 0:size_e[2] - 1] =  E[*,0:size_e[2]-1]
  
  extE[*,   size_e[2] + 0] =  E[*, 1] + E[*, 2]                        ; total energy
  extE[*,   size_e[2] + 1] =  E[*, 2] - E[*, 1]                        ; residual energy
  extE[*,   size_e[2] + 2] = (E[*, 2] - E[*, 1]) / (E[*, 1] + E[*, 2]) ; normalize residual energy
  
  extE[*,   size_e[2] + 3] = SQRT(E[*,4])                              ; sqrt(Z+^2)
  extE[*,   size_e[2] + 4] = SQRT(E[*,5])                              ; sqrt(Z-^2)
  
  i_ext_i                  = size_e[2] + 5
  i_ext_f                  = size_e[2] + 9

  extE[*, i_ext_i:i_ext_f] = ZERO                                     ;

; Loop over layers (z's)
; fetch_layer
; calculate zp and zm
; calculate FFT's of zp and zm
; calculate square modulus sums over k-vectors
; store in extE at i_lpl and i_lmi
; average over time along with other qty's
; multiply by pi / Z_{+/-} appropriately for final result

  IF (i_efld EQ 23 || i_efld EQ 24) THEN BEGIN

    n_layers                 = size_e[1]
    
    Xf                       = FINDGEN((n1-1)/2) + ONE
    is_n1_even               = (n1 MOD 2) EQ 0

    IF (is_n1_even) THEN BEGIN
      KX                     = TWO_PI * [ ZERO, Xf,   FLOAT(n1)/TWO, FLOAT(-n1)/TWO  + Xf ] / ( FLOAT(n1)*dx)
    ENDIF ELSE BEGIN 
      KX                     = TWO_PI * [ ZERO, Xf, -(FLOAT(n1)/TWO + ONE) + Xf ] / (FLOAT(n1)*dx)
    ENDELSE

    Yf                       = FINDGEN((n2-1)/2) + ONE
    is_n2_even               = (n2 MOD 2) EQ 0
    
    IF (is_n2_even) THEN BEGIN 
      KY                     = TWO_PI * [ ZERO, Yf,   FLOAT(n2)/TWO, FLOAT(-n2)/TWO  + Yf ] / ( FLOAT(n2)*dy)
    ENDIF ELSE BEGIN
      KY                     = TWO_PI * [ ZERO, Yf, -(FLOAT(n2)/TWO + ONE) + Yf ] / (FLOAT(n2)*dy)
    ENDELSE

    KK                       = SQRT((KX^2)#(KY^2))

    size_kk                  = SIZE(KK, /DIMENSIONS)    ; dims
    nz_kk                    = WHERE(KK NE ZERO, count) ;location
    ind_nz_kk                = ARRAY_INDICES(size_kk, nz_kk, /DIMENSIONS)


    FOR  I = 0, n_layers - 1 DO BEGIN
;   FOR  I = 0, 0 DO BEGIN
    
      SLB                    = fetch_layer( I+1, step, desc_lab)
    
      P                      = REFORM(SLB[*,0], n1, n2)
      A                      = REFORM(SLB[*,1], n1, n2)
    
      FTP                    = FFT(P)
      FTA                    = FFT(A)

      FTZ_PLS                = FTP + FTA
      FTZ_MNS                = FTP - FTA

      SQR_MOD_FTZ_PLS        = (ABS(FTZ_PLS))^2
      SQR_MOD_FTZ_MNS        = (ABS(FTZ_MNS))^2

      extE[I, 23]            = TOTAL( SQR_MOD_FTZ_PLS[ind_nz_kk[0,*], ind_nz_kk[1,*]] / KK[ind_nz_kk[0,*], ind_nz_kk[1,*]] )
      extE[I, 24]            = TOTAL( SQR_MOD_FTZ_MNS[ind_nz_kk[0,*], ind_nz_kk[1,*]] / KK[ind_nz_kk[0,*], ind_nz_kk[1,*]] )
      
    ENDFOR

    extE[*,23:24]            =  PI * extE[*,23:24]

  ENDIF
  
  RETURN, extE

END
