FUNCTION setKData, k_drive_nom, ELOG

  HALF                        = 0.5D
  i_k                         = 0
  k_drive                     = k_drive_nom

  size_ELOG                   = SIZE(ELOG,/DIMENSIONS)
  n_lines                     = size_ELOG[0]

  DK                          = DBLARR(n_lines-1)
  FOR I = 0, n_lines - 2 DO BEGIN
    DK[I]                     = ELOG[i_k,I+1] - ELOG[i_k, I]
  ENDFOR
  dk                          = MEAN(DK)
  k_drive_idx                 = WHERE((ELOG[i_k, *] GE k_drive - half*dk) AND (ELOG[i_k, *] LE k_drive + half*dk),  low_count )
  k_drive                     = ELOG[i_k, k_drive_idx]

  k_low                       = 20.0
  k_high                      = 4.0 * !pi * 128.0/12.0

  k_low_idx                   = WHERE((ELOG[i_k, *] GE k_low   - half*dk) AND (ELOG[i_k, *] LE k_low   + half*dk),  low_count )
  k_high_idx                  = WHERE((ELOG[i_k, *] GE k_high  - half*dk) AND (ELOG[i_k, *] LE k_high  + half*dk), high_count )

  k_high                      = ELOG[i_k, k_high_idx]
  k_low                       = ELOG[i_k, k_low_idx ]

  k_data                      = {k_drive:k_drive, k_low:k_low, k_low_idx:k_low_idx, k_high:k_high, k_high_idx:k_high_idx}

  RETURN, k_data

END
