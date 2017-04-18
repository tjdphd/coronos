FUNCTION getSliceRange, first_slice, last_slice

  n_layers                                 = last_slice - first_slice + 1

  nCPUs                                    = !CPU.HW_NCPU
  slices_per_cpu                           = CEIL(n_layers / nCPUs)

  IF (nCPUs GT 1) THEN BEGIN
    IF (slices_per_cpu GE 1) THEN BEGIN
      n_threads                            = nCPUs
    ENDIF ELSE BEGIN
      n_threads                            = n_layers
    ENDELSE
  ENDIF ELSE BEGIN
    n_threads                              = 1
  ENDELSE
  IF (n_threads GT 12) THEN n_threads = 12
  stagger                                  = 0
  IF ((slices_per_cpu GT 1) AND (last_slice  - (n_threads * slices_per_cpu)) GE n_threads / 2 ) THEN stagger = 1
  IF (slices_per_cpu EQ 0 ) THEN slices_per_cpu = 1

  slice_range                              = INTARR(n_threads,  2)

  bump                                     = 0

  FOR I = 0, n_threads - 1 DO BEGIN

    IF (I MOD 2 EQ 0) THEN  BEGIN

      slice_range[I, 0]                    = first_slice  + (I       * slices_per_cpu)      +   bump * stagger
      slice_range[I, 1]                    = first_slice  + ((I + 1) * slices_per_cpu) - 1  +   bump * stagger

    ENDIF ELSE BEGIN

      slice_range[I, 0]                    = first_slice  + (I       * slices_per_cpu)      +   bump       * stagger
      slice_range[I, 1]                    = first_slice  + ((I + 1) * slices_per_cpu) - 1  + ( bump + 1 ) * stagger

      bump = bump + 1

    ENDELSE

    IF ( (I EQ n_threads - 1) AND ( slice_range[I,1] LT last_slice)) THEN slice_range[I,1] = last_slice

  ENDFOR

  RETURN, slice_range

END
