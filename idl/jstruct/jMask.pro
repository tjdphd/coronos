FUNCTION jMask, CNTRS, width, JVOX

  size_jvox          = SIZE(JVOX, /DIMENSIONS)
  MASK               = BYTARR(size_jvox)

  size_cntrs         = SIZE(CNTRS, /DIMENSIONS)

  n_cntrs            = size_cntrs[0]

  FOR K = 0, n_cntrs - 1 DO BEGIN

    this_contour     = CNTRS[K]

    cntr_lb          = (1.0 - width) * this_contour
    cntr_ub          = (1.0 + width) * this_contour

    IF  ((this_contour - cntr_lb) LE 0.75) THEN  cntr_lb = this_contour - 0.75
    IF  ((cntr_ub - this_contour) LE 0.75) THEN  cntr_ub = this_contour + 0.75

    PRINT, 'jMask: K            = ', K
    PRINT, 'jMask: this_contour = ', this_contour
    PRINT, 'jMask: cntr_lb      = ', cntr_lb
    PRINT, 'jMask: cntr_ub      = ', cntr_ub

    IF (K LE n_cntrs - 1) THEN BEGIN
      iso_surf       = WHERE( (JVOX GE cntr_lb) AND (JVOX LE cntr_ub), count)
      IF (count EQ 0) THEN PRINT, 'jMask: count is zero for contour K = ', K
      MASK[iso_surf] = JVOX[iso_surf]
    ENDIF

  ENDFOR

RETURN, MASK

END
