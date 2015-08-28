FUNCTION calc_jz, Slc, n_slice, n_step, desc_label

COMMON save_r, Pre
COMMON First, first_call
COMMON fnc, prefix, inc_res, nfld
COMMON CC, Cntrs
COMMON loc, eps_out_dir

size_slc     = SIZE(Slc)
rank         = size_slc[0]

nlines       = LONG64(0)

IF (rank EQ 2) THEN nlines = size_slc[1]
IF (rank EQ 3) THEN nlines = size_slc[2]

A            = FLTARR(2, nlines)                                       ; Allocate space for plotting
X            = FLTARR(1, nlines)                                       ; arrays
Y            = FLTARR(1, nlines)  

A            = TRANSPOSE(A)  
X            = TRANSPOSE(X)
Y            = TRANSPOSE(Y)  

IF (rank EQ 2) THEN BEGIN

  X          = Slc[*,0]                                                ; initialize plotting arrays
  Y          = Slc[*,1]                                                ; with data from Slc
  A[*,0]     = Slc[*,2]
  A[*,1]     = Slc[*,3]

ENDIF ELSE BEGIN

      X[*]   = Slc[0, *, 0]
      Y[*]   = Slc[0, *, 1]
      A[*,0] = Slc[0, *, 2]
      A[*,1] = Slc[0, *, 3]
         
      ENDELSE

dxyM         = cfdp(X, Y, A)

d2xyM        = cfdp(X,Y,dxyM)
d2AdY2       = d2xyM[*, 1]

size_d2xyM   = SIZE(d2xyM)
rows         = size_d2xyM[1]
cols         = size_d2xyM[2]

dyxM         = FLTARR(cols,rows)
dyxM         = TRANSPOSE(dyxM)

dyxM[*, 1]   = dxyM[*,  0]
dyxM[*, 0]   = dxyM[*,  1]

d2yxM        = cfdp(X,Y,dyxM)

d2AdX2       = d2yxM[*, 0]

J_z          = -(d2AdX2 + d2AdY2)

RETURN, J_z

END
