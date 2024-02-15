FUNCTION GAUSSIANFILTER, x, y, s
  x = x(*)
  y = y(*)
  nx = N_ELEMENTS(x)
  dx = x(1) - x(0)
  dy = TRANSPOSE(y(1 : nx - 1) + y(0 : nx - 2))
  nG = ROUND(6 * s / dx)
  xG = DINDGEN(nG, INCREMENT = dx, START = dx)
  xG = [REVERSE(-xG), 0, xG]
  yG = 1 / (s * SQRT(2 * !dpi)) * EXP(-.5 * (xG / s)^2)
  ynew = DBLARR(nx)
  FOR i = 0, nx - 1 DO ynew(i) = dy() # yG(MAX([0, nG + 1 - i]) : MIN([2 * nG, nx - 1 - i + nG]))
  ynew *= .25 * dx
  RETURN, ynew
end