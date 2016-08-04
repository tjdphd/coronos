FUNCTION calcKK, n1, n2

  ZERO                  = 0.0D                                      ; Define some constants for convenience
  ONE                   = 1.0D                                      ; everywhere
  TWO                   = 2.0D
  THREE                 = 3.0D
  PI                    = 3.141592653589793D
  TWO_THIRDS            = TWO / THREE
  TWO_PI                = TWO*PI

  CZERO                 = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  CI                    = DCOMPLEX(ZERO,ONE)                       ; and to make sure that DOUBLE's are used

    Xf                       = FINDGEN((n1-1)/2) + ONE
    is_n1_even               = (n1 MOD 2) EQ 0

    IF (is_n1_even) THEN BEGIN
      KX                     = TWO_PI * [ ZERO, Xf,   FLOAT(n1)/TWO, FLOAT(-n1)/TWO  + Xf ] ; / ( FLOAT(n1)*dx)
    ENDIF ELSE BEGIN 
      KX                     = TWO_PI * [ ZERO, Xf, -(FLOAT(n1)/TWO + ONE) + Xf ] ; / (FLOAT(n1)*dx)
    ENDELSE

    Yf                       = FINDGEN((n2-1)/2) + ONE
    is_n2_even               = (n2 MOD 2) EQ 0
    
    IF (is_n2_even) THEN BEGIN 
      KY                     = TWO_PI * [ ZERO, Yf,   FLOAT(n2)/TWO, FLOAT(-n2)/TWO  + Yf ] ; / ( FLOAT(n2)*dy)
    ENDIF ELSE BEGIN
      KY                     = TWO_PI * [ ZERO, Yf, -(FLOAT(n2)/TWO + ONE) + Yf ] ; / (FLOAT(n2)*dy)
    ENDELSE

    size_kx                  = SIZE(KX, /DIMENSIONS)
    size_ky                  = SIZE(KY, /DIMENSIONS)

    KK                       = FLTARR(size_kx[0],size_ky[0])

    FOR I = 0, n1-1 DO BEGIN
      FOR J = 0, n2-1 DO BEGIN
        KK[I,J]              = SQRT(KX[I]^2 + KY[J]^2)
      ENDFOR
    ENDFOR

    PRINT, 'calcKK: completed...'
    RETURN, KK

END
