FUNCTION CONVOL_PRECISE, x, y, ils
  RETURN, CONVOL(y, ILS, /CENTER, /EDGE_WRAP) * (x(1) - x(0))
END
