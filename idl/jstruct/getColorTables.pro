FUNCTION getColorTables, qty

  inttab                = indgen(256)
  fintab                = FLTARR(256)

; base                  = SQRT(SQRT(2.0))
  base                  = SQRT(2.0)
; base                  = SQRT(3.0)
; base                  = SQRT(4.0)
; base                  = SQRT(9.0)

  FOR I = 0, 255 DO BEGIN

    fintab[I]           = (256.0/base) * (base^(FLOAT(I+1) / 256.0))
;   PRINT, "getColorTables: inttab = ", fintab[I]

  ENDFOR
   
  coltab0               = BYTARR(256,3)
  coltab1               = BYTARR(256,3)

  IF (qty EQ 'j') THEN BEGIN
    coltab0[*,0]        = BYTE(fintab[*])                 ; red
    coltab1[*,2]        = BYTE(fintab[*])                 ; blue
  ENDIF ELSE BEGIN
    IF (qty EQ 'o') THEN BEGIN
      coltab0[*,0]      = BYTE(fintab[*])                 ; orange
      coltab0[*,1]      = BYTE(fintab[*]/2.0)             ; 

      coltab1[*,1]      = BYTE(fintab[*])                 ; green
    ENDIF
  ENDELSE

  coltabs               = BYTARR(256,3,2)

  coltabs[*,*,0]        = coltab0[*,*]
  coltabs[*,*,1]        = coltab1[*,*]

  RETURN, coltabs

END
