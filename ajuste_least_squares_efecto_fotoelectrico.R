library(readxl)
library(ggplot2)

data_fe <- read_excel("Escritorio/f3/tp_efecto_fotoelectrico/EF 1C2021 datos exprimentales - Grupo 2.xlsx")
data_fe
c <- 2.99792458e8 # m/s 
wavelength <- data_fe[[1]]
voltage <- data_fe[[2]]

frequency = c / wavelength # v = c/lambda -> no hay error proveniente de los filtros

frequency
df <- data.frame(frequency, voltage)
df
plt.1 <- ggplot(df, aes(x=frequency, y=voltage)) + 
  geom_point(colour="darkblue") + 
  xlab("Frecuencia (Hz)") + ylab("Tension estimada (V)")

ggsave("recta_ajuste.png", plot = last_plot(), device = "png", scale = 1, width = NA, height = NA, dpi = 300)


## CALCULO DE ESTIMADORES DE LA RECTA (weighted least squares method)
# valor absoluto de la incertidumbre de medicion en la tension
calc_err <- function(voltage){
  return(voltage * (0.8 / 100) + 0.001)
}

data.frame(voltage)

voltage_errors <- function(voltages){
  return(sapply(data.frame(voltages), calc_err))
}

sigma <- voltage_errors(voltage)
sigma
# incertidumbre en las tensiones medidas
sigma <- t(sigma)

ggplot(df, aes(x=frequency, y=voltage)) + 
  geom_errorbar(aes(ymin=voltage - sigma, ymax=voltage + sigma), width=1) +
  geom_point() + 
  geom_point(colour = "darkblue") + 
  xlab("Frecuencia (Hz)") + ylab("Tension estimada (V)")
ggsave("recta_ajuste.png", plot = last_plot(), device = "png", scale = 1, width = NA, height = NA, dpi = 300)
S <- sum(1/sigma^2) # (1/V^2)

sigmaSqr <- sigma^2

sumX <- sum(frequency / sigmaSqr) #sigmaX(Hz / V^2)
sumY <- sum(voltage / sigmaSqr) #sigmaY (1 / V)
sumXsqr <- sum((frequency^2 / sigmaSqr)) #sigma(X^2) (Hz^2 / V^2)
sumXY <- sum((voltage * frequency) / sigmaSqr) #sigmaXY (Hz / V)

a1 <- ((sumY * sumXsqr) - (sumX * sumXY)) / (S * sumXsqr - (sumX)^2) # ordenada al origen (volts)
a2 <- (S * sumXY - (sumX * sumY)) / (S * sumXsqr - (sumX)^2) # pendiente de la recta (V/Hz)

a2
a1
plot(frequency, frequency * a2 + a1)

## INCERTIDUMBRE PARA LOS ESTIMADORES DE LA PENDIENTE Y LA ORDENADA AL ORIGEN

sigma.a1 <- sqrt(sumXsqr / (S*sumXsqr - (sumX)^2))
sigma.a2 <- sqrt(S / (S*sumXsqr - (sumX)^2))

sigma.a1 ## (Volts)
sigma.a2 ## (V/Hz)

## ESTIMADOR INSESGADO DE LA VARIANZA DE LA MUESTRA (para las barras de incertidumbre)

mse <- sum((voltage - (a2 * frequency + a1))^2)
mse
s.sqrt <- (1 / (length(voltage) - 2)) * mse
s.sqrt
std.dev <- sqrt(s.sqrt)
std.dev # (Volts)

data <- data.frame("Frecuencia_(Hz)" = frequency, "Tension_estimada_(V)" = a2 * frequency + a1)
data

## Plot con las barras de +- 1 desvio estandar
ggplot(data, aes(x=Frecuencia_.Hz., y=Tension_estimada_.V.)) + 
  geom_errorbar(aes(ymin=Tension_estimada_.V. - std.dev, ymax=Tension_estimada_.V. + std.dev), width=.5) +
  geom_line(colour = "royalblue") +
  geom_point() + 
  geom_point(data = data.frame(voltage), 
             mapping = aes(x = frequency, y = voltage), colour = "seagreen", size=2) + 
  xlab("Frecuencia (Hz)") + ylab("Tension estimada (V)")

## Calculo de la cte. de Planck
q <- 1.602176487e-19 # Coulomb
h <- q * a2 # Joule * sec
h
sigma.h <- h * (sigma.a2 / a2)
sigma.h
sigma.h / h # error relativo 0.9%

h.real <- 6.62606896e-34
(h - h.real) / h.real # 3% de error aprox considerando el valor de h.real como exacto

## Calculo de la funcion trabajo del material
w <- - a1 * q # volt * coulomb = joule
w
sigma.w <- w * (sigma.a1 / a1)
sigma.w
abs(sigma.w / w) #error relativo 1.17% 
w
## pasamos de joules a electronVolt
w.ev <- w/q
w.ev
sigma.w.ev <- abs(sigma.w / q)
sigma.w.ev
w.ev - sigma.w.ev
