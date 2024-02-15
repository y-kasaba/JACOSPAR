FUNCTION GAUSSIAN, x, A, fwmh, x0

  sigma = fwmh * 0.424660900144009534340483469350D ; %1/(2*sqrt(2*log(2)))
  sigmaP = 0.199471140200716351431609041356D / sigma ; %1/sqrt(2*pi)
  sigmaPP = 0.5D / sigma^2;
  
  RETURN, A * sigmaP * EXP(-(x - x0)^2 * sigmaPP)
end