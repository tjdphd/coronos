FUNCTION makeContourPlot, n_cntrs, n_slice, n_step, desc_label, width, global_j_minmax, JVOX, cc, pom

  ip1                   = scan_parameters('ip1', 0, desc_label )
  ip2                   = scan_parameters('ip2', 0, desc_label )
  n1                    = 2^ip1
  n2                    = 2^ip2
  n3                    = scan_parameters('n3', 0, desc_label)
  mp                    = scan_parameters('mp', 0, desc_label)

  SLB                   = fetch_layer(n_slice, n_step, desc_label )

  X                     = TRANSPOSE(REFORM(SLB[*,0], n1, n2))
  Y                     = TRANSPOSE(REFORM(SLB[*,1], n1, n2))

  x_max                 = MAX(X)
  y_max                 = MAX(Y)

  X                     = X * (n1 / x_max)
  Y                     = Y * (n2 / y_max)

  A                     = BYTARR(n1,n2)
  A[*,*]                = JVOX[n3*mp/2, *, *]

  Z                     = FLTARR(n1,n2)
  Z[*,*]                = 0.5*n3*mp

  CNTRS                 = FLTARR(n_cntrs)

  max_jvox              = MAX(JVOX) ; should be 255 or 256

  PRINT, 'makeContourPlot: max_jvox = ', max_jvox

  a_max                 = (255.0) / (1.0 + width)
  a_min                 = 0.01 * a_max

  inc_a                 = (a_max - a_min) / ((n_cntrs - 1))

  CNTRS[0]              = a_min
  CNTRS[n_cntrs -1]     = a_max

  FOR cdx = 1, n_cntrs - 1 DO BEGIN
    IF (cdx LT n_cntrs - 1) THEN $
    CNTRS[cdx]          = CNTRS[cdx-1] + inc_a
  ENDFOR

  AContour              = OBJ_NEW('IDLgrContour', A,   $
                          COLOR       = [0, 0, 255],   $
                          C_VALUE     = CNTRS,         $
                          PLANAR      = 0,             $
                          XCOORD_CONV = cc[*,2],       $
                          YCOORD_CONV = cc[*,0],       $
                          ZCOORD_CONV = cc[*,1],       $
                          GEOMX       = Z,             $
                          GEOMY       = X,             $
                          GEOMZ       = Y              $
                          ) 

RETURN, AContour

END
