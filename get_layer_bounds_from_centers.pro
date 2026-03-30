FUNCTION GET_LAYER_BOUNDS_FROM_CENTERS, z
  nlayer = N_ELEMENTS(z)
  IF nlayer EQ 0 THEN RETURN, DBLARR(0)

  zloc = DOUBLE(z)
  bounds = DBLARR(nlayer + 1L)

  IF nlayer EQ 1 THEN BEGIN
    bounds(0) = zloc(0) > 0D
    bounds(1) = zloc(0) + 1D
    RETURN, bounds
  ENDIF

  bounds(0) = zloc(0) - 0.5D * (zloc(1) - zloc(0))
  IF bounds(0) LT 0D THEN bounds(0) = 0D

  FOR i = 1L, nlayer - 1L DO bounds(i) = 0.5D * (zloc(i - 1L) + zloc(i))

  bounds(nlayer) = zloc(nlayer - 1L) + 0.5D * (zloc(nlayer - 1L) - zloc(nlayer - 2L))

  RETURN, bounds
END
