PRO glb_ext, desc_label,         $
             qty,                $
             first_step,         $
             last_step,          $
             first_slice,        $
             last_slice


  !EXCEPT =2

  i_x                            =  6
  i_j                            =  5
  i_o                            =  4
  i_a                            =  3
  i_p                            =  2

  i_qty                          = -1

  CASE qty OF
  'x' : i_qty                    =  i_x
  'j' : i_qty                    =  i_j
  'o' : i_qty                    =  i_o
  'a' : i_qty                    =  i_a
  'p' : i_qty                    =  i_p
  ELSE: i_qty                    = -1
  ENDCASE

  cur_dir                        = GETENV('PWD')
  out_file                       = cur_dir + '/' + 'glb_ext.out'
  OPENR, out_unit, out_file, /GET_LUN, ERROR = op_err

  IF ( op_err EQ 0) THEN BEGIN

    PRINT, 'glb_ext: reading previous output from ', out_file
    global_out_minmax = FLTARR(4,4)
    line                         = FLTARR(4)
    FOR i_q = 0, 3 DO BEGIN
      READF, out_unit, FORMAT = '(4(e24.16,1x),:)', line
      global_out_minmax[*,i_q]   = line
    ENDFOR
    CLOSE, out_unit
    global_out_minmax            = global_q_minmax(i_qty, first_step, last_step, first_slice, last_slice, desc_label, $
                                                   GLOBAL_IN_MINMAX = global_out_minmax                               $
                                                  )
  ENDIF ELSE BEGIN

    PRINT, 'glb_ext: no previous output'
    global_out_minmax            = global_q_minmax(i_qty, first_step, last_step, first_slice, last_slice, desc_label)

    PRINT, 'glb_ext: found minmax'
  ENDELSE

  IF (i_qty EQ i_x) THEN BEGIN

    PRINT, FORMAT='(A19,I1,A4,4(e12.4,1x),:)', 'global_out_minmax[*',0,'] = ', global_out_minmax[*,0]
    PRINT, FORMAT='(A19,I1,A4,4(e12.4,1x),:)', 'global_out_minmax[*',1,'] = ', global_out_minmax[*,1]
    PRINT, FORMAT='(A19,I1,A4,4(e12.4,1x),:)', 'global_out_minmax[*',2,'] = ', global_out_minmax[*,2]
    PRINT, FORMAT='(A19,I1,A4,4(e12.4,1x),:)', 'global_out_minmax[*',3,'] = ', global_out_minmax[*,3]

    cur_dir                      = GETENV('PWD')
    out_file                     = cur_dir + '/' + 'glb_ext.out'

    OPENW, out_unit, out_file, /GET_LUN

    FOR qdx = 0, 3 DO BEGIN
      PRINTF, out_unit, FORMAT = '(4(e24.16,1x),:)', global_out_minmax[*, qdx]
    ENDFOR

    CLOSE, out_unit

  ENDIF ELSE BEGIN

    PRINT, FORMAT='(A20,4(e12.4,1x),:)', 'global_out_minmax = ', global_out_minmax[*]

  ENDELSE

END
