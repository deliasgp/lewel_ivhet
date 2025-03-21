sim_mtriangular_ivlewel <- function(n){ #Modelo triangular
  x1  <- rnorm(n); U  <- rnorm(n)
  s1 <- rnorm(n); s2 <- rnorm(n)
  eps1 <- U + exp(x1)*s1 #Residuo con heterocedasticidad
  eps2 <- U + exp(-x1)*s2 #Residuo con heterocedasticidad
  b10 <- b11 <- b20 <- b21 <- gamma1 <- 1
  y2 <- b10 + b11*x1 + eps2 #Primera etapa
  y1 <- b20 + b21*x1 + gamma1*y2 + eps1 #Segunda etapa
  z <- x1
  #----------------------------------------------------------------------------*
  m1 <-  lm(y2~x1); m1 <- summary(m1)  #Modelo de etapa 1
  eps2_hat <- m1$residuals #Residuos (heterocedasticos del model) del modelo
  IV_HET <- (z - mean(z))*eps2_hat
  df_temp <- data.frame(y1=y1, y2=y2, x1=x1, z=z, inst_extra=IV_HET)
  ivreg(y1 ~ x1 + y2 | z + IV_HET, data=df_temp)
  m2 <- ivreg(y1 ~ x1 + y2 | z + IV_HET, data=df_temp)
  theta <- c(m1$coefficients[,1],m2$coefficients)
  return(theta)
}