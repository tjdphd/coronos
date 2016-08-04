FUNCTION getCoordConvs, desc_label, size_jvox

  zl                      = scan_parameters('zl',0, desc_label)

  cx1                     =  0.0
  cx2                     =  1.0 / size_jvox[1]
                         
  cy1                     =  0.0
  cy2                     =  1.0 / size_jvox[2]
                         
  cz1                     =  0.0
  cz2                     = (0.125 * zl) / size_jvox[0]
                         
  cx                      = [cx1, cx2]
  cy                      = [cy1, cy2]
  cz                      = [cz1, cz2]

  cc                      = FLTARR(2,3)

  cc[*,0]                 = cx[*]
  cc[*,1]                 = cy[*]
  cc[*,2]                 = cz[*]

  RETURN, cc

END
