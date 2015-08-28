FUNCTION cfdp, X, Y, M

COMMON loc, eps_out_dir

  ip1         = scan_parameters('ip1', 0, desc_label )
  ip2         = scan_parameters('ip2', 0, desc_label )
  x_res       = 2^ip1
  y_res       = 2^ip2
  dx          = 1.0 / (x_res)
  dy          = 1.0 / (y_res)
  i_x_res     = LONG64(x_res)
  i_y_res     = LONG64(y_res)
  slc_lines   = LONG64(i_x_res   * i_y_res)
  dxyM        = FLTARR(2, slc_lines)
  dxyM        = TRANSPOSE(dxyM)

      FOR I   = 0, slc_lines - 1 DO BEGIN
         IF ( ( I  MOD i_y_res) EQ 0 ) THEN BEGIN
            IF ( I EQ 0 ) THEN BEGIN 
              im1_idx = I + (i_y_res * (i_x_res - 1))
              i_idx   = I                                                                  ; This is the bottom boundary x
              ip1_idx = I + i_y_res
            ENDIF ELSE BEGIN
                      IF (I GE (i_x_res * (i_y_res - 1))) THEN BEGIN
                        im1_idx = I - i_y_res 
                        i_idx   = I                                                        ; This is the top boundary in x
                        ip1_idx = I - ( i_y_res * (i_x_res - 1) )
                      ENDIF ELSE BEGIN
                                 im1_idx = I - i_y_res 
                                 i_idx   = I                                               ; Here is where we do the interior 
                                 ip1_idx = I + i_y_res                                     ; points in x
                            ENDELSE
                  ENDELSE
            FOR J = 0, i_y_res - 1 DO BEGIN
               IF ( J EQ 0 ) THEN BEGIN
                                    ip1_idy = i_idx + J + 1                                ; y = 0 boundary case
                                    im1_idy = i_idx + i_y_res - 1
               ENDIF ELSE BEGIN
                          IF ( J EQ i_y_res - 1 ) THEN BEGIN
                             ip1_idy = i_idx                                               ; y = 1 boundary case
                             im1_idy = i_idx + J - 1
                          ENDIF ELSE BEGIN
                                    ip1_idy = i_idx + J + 1                                ; interior y case
                                    im1_idy = i_idx + J - 1
                                ENDELSE
                     ENDELSE
                dxyM[i_idx + J, 0] = (M[ip1_idx + J,  1] - M[im1_idx + J, 1]) / ( 2 * dx ) ; partial with respect to x
                dxyM[i_idx + J, 1] = (M[ip1_idy,      1] - M[im1_idy,     1]) / ( 2 * dy ) ; partial with respect to y
            ENDFOR
         ENDIF
      ENDFOR
  RETURN, dxyM
END
