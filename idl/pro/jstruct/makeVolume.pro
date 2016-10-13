FUNCTION makeVolume, JVOX, coltab, cc

  JVolume =  OBJ_NEW('IDLgrVolume', JVOX, /ZERO_OPACITY_SKIP)
  JVolume -> SetProperty, RGB_TABLE0  = coltab
  JVolume -> SetProperty, XCOORD_CONV = cc[*,2], YCOORD_CONV = cc[*,0], ZCOORD_CONV = cc[*,1]

  RETURN, JVolume

END
