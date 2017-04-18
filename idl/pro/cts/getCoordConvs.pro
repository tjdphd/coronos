FUNCTION getCoordConvs, desc_label, size_q

  size_size_q               = SIZE(size_q, /DIMENSIONS)

  zl                        = scan_parameters('zl',0, desc_label)

  IF (size_size_q[0] EQ 3) THEN BEGIN

    cx1                     =  0.0
    cx2                     =  1.0 / size_q[1]
                           
    cy1                     =  0.0
    cy2                     =  1.0 / size_q[2]
                           
    cz1                     =  0.0
    cz2                     = (0.125 * zl) / size_q[0]
                           
    cx                      = [cx1, cx2]
    cy                      = [cy1, cy2]
    cz                      = [cz1, cz2]

    cc                      = FLTARR(2,3)

    cc[*,0]                 = cx[*]
    cc[*,1]                 = cy[*]
    cc[*,2]                 = cz[*]

  ENDIF ELSE BEGIN

    IF (size_size_q[0] EQ 2) THEN BEGIN

      cx1                   =  0.0
      cx2                   =  1.0 / size_q[0]
                            
      cy1                   =  0.0
      cy2                   =  1.0 / size_q[1]
                            
      cx                    = [cx1, cx2]
      cy                    = [cy1, cy2]

      cc                    = FLTARR(2,2)

      cc[*,0]               = cx[*]
      cc[*,1]               = cy[*]

    ENDIF
  ENDELSE

  RETURN, cc

END
