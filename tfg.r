# Para borrar la consola y sea más fácil la lectura
cat("\014")

#install.packages("ggplot2")
#install.packages("pillar")
#install.packages("pscl")   
#install.packages("performance")

# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(MASS) 
library(pscl)
library(performance)
library(lmtest)

##############################################
#
# Carga de datos
#
##############################################
datos <- read.csv("G:/Mi unidad/UNED/5. TFG/arsenic.csv", header=T)
# Transformar a factor algunas variables
datos$gender <- factor(datos$gender, labels = c("Mujer", "Hombre"))
datos$type   <- factor(datos$type, labels = c("Vejiga", "Pulmón"))


##############################################
#
# Exploracion de datos conc
#
##############################################

summary_stats_conc <- c(
  Media   = mean(datos$conc),
  Min     = min(datos$conc),
  Q1      = quantile(datos$conc, 0.25),
  Mediana = median(datos$conc),
  Q3      = quantile(datos$conc, 0.75),
  Max     = max(datos$conc)
)

print(summary_stats_conc)

hist(datos$conc,
     breaks = 40,
     col = "steelblue",
     border = "black",
     main = "Distribución de la concentración de arsénico (conc)",
     xlab = "Concentración de arsénico",
     ylab = "Frecuencia")

tabla_events <- c(
  events_0      = count(filter(datos, events == 0)),
  events_1      = count(filter(datos, events == 1)),
  events_2      = count(filter(datos, events == 2)),
  events_3      = count(filter(datos, events == 3)),
  events_4      = count(filter(datos, events == 4)),
  events_5o_mas = count(filter(datos, events >= 5))
)

print(tabla_events)

tabla_casos <-c(
  poblacion_total = sum(datos$at.risk),
  casos_total = sum(datos$events),
  poblacion_no_expuesta = sum(filter(datos, conc == 0)$at.risk),
  casos_no_expuesta = sum(filter(datos, conc == 0)$events),
  poblacion_expuesta = sum(filter(datos, conc > 0)$at.risk),
  casos_expuesta = sum(filter(datos, conc > 0)$events)
)

print(tabla_casos)
print("Porcentaje casos total")
print(tabla_casos["casos_total"]/tabla_casos["poblacion_total"]*100)
print("Porcentaje casos no expuesta")
print(tabla_casos["casos_no_expuesta"]/tabla_casos["poblacion_no_expuesta"]*100)
print("Porcentaje casos expuesta")
print(tabla_casos["casos_expuesta"]/tabla_casos["poblacion_expuesta"]*100)

##############################################
#
# Analisis de datos
#
##############################################

# --------------------------------------------
# MODELO ALL DATA: TODOS LOS DATOS
# --------------------------------------------

cat("\n\nModelo all data \n=================================\n")

modelo_all_data <- glm(events ~ conc + gender + type, 
                       offset = log(at.risk), 
                       family = poisson(link = "log"), 
                       data = datos)

print(summary(modelo_all_data))


# Diagnóstico de sobredispersión
print(check_overdispersion(modelo_all_data))

# Calculo de los ceros observados y esperados
ceros_observados <- sum(datos$events == 0)
mu_hat <- predict(modelo_all_data, type = "response") 
prob_cero_poisson <- exp(-mu_hat)
ceros_esperados <- sum(prob_cero_poisson)
cat("Ceros observados:", ceros_observados, "\n")
cat("Ceros esperados por modelo Poisson:", round(ceros_esperados, 2), "\n")

# --------------------------------------------
# Datos filtrados
# --------------------------------------------
cat("\n\nFiltrado conc > 0\n=================================\n")
datos_filtrados <- filter(datos, conc > 0)


# --------------------------------------------
# MODELO: Poisson events ~ conc + gender + type
# --------------------------------------------
cat("\n\nModelo events ~ conc + gender + type\n=================================\n")

modelo_Poisson_CGT <- glm(events ~ conc + gender + type, 
                          offset = log(at.risk), 
                          family = poisson(link = "log"), 
                          data = datos_filtrados)

print(summary(modelo_Poisson_CGT))

# Diagnóstico de sobredispersión
print(check_overdispersion(modelo_Poisson_CGT))

# --------------------------------------------
# MODELO: Poisson events ~ conc + type
# --------------------------------------------
cat("\n\nModelo Poisson events ~ conc + type\n=================================\n")

modelo_Poisson_CT <- glm(events ~ conc + type, 
                         offset = log(at.risk), 
                         family = poisson(link = "log"), 
                         data = datos_filtrados)

print(summary(modelo_Poisson_CT))

# Diagnóstico de sobredispersión
print(check_overdispersion(modelo_Poisson_CT))

# Calculo de los ceros observados y esperados
ceros_observados <- sum(datos_filtrados$events == 0)
mu_hat <- predict(modelo_Poisson_CT, type = "response") 
prob_cero_poisson <- exp(-mu_hat)
ceros_esperados <- sum(prob_cero_poisson)
cat("Ceros observados:", ceros_observados, "\n")
cat("Ceros esperados por modelo Poisson:", round(ceros_esperados, 2), "\n")

# --------------------------------------------
# MODELO: QuasiPoisson
# --------------------------------------------
cat("\n\nModelo quasipoisson \n=================================\n")

modelo_qp <- glm(events ~ conc + gender + type, 
                 offset = log(at.risk), 
                 family = quasipoisson(link = "log"), 
                 data = datos_filtrados)


print(summary(modelo_qp))

# Diagnóstico de sobredispersión
print(check_overdispersion(modelo_qp))

# Diagnóstico de sobredispersión
m = modelo_qp
sum_res = sum(residuals(m, type = "deviance")^2)
df_res = m$df.residual  
cat("Factor de sobredispersión:", sum_res / df_res, "\n")
cat("Pearson's Chi-squared:", sum_res, df_res, "\n")
cat("Contraste:", 1-pchisq(sum_res, df_res), "\n")


# --------------------------------------------
# MODELO: Binomial Negativa events ~ conc + gender + type + offset(log(at.risk))
# --------------------------------------------
cat("\n\nModelo NB_CGT \n=================================\n")
modelo_NB_CGT <- glm.nb(events ~ conc + gender + type + offset(log(at.risk)), 
                        data = datos_filtrados)


print(summary(modelo_NB_CGT))

print(check_overdispersion(modelo_NB_CGT))
# Diagnóstico de sobredispersión
m = modelo_NB_CGT
sum_res = sum(residuals(m, type = "deviance")^2)
df_res = m$df.residual  
cat("Factor de sobredispersión:", sum_res / df_res, "\n")
cat("Pearson's Chi-squared:", sum_res, df_res, "\n")
cat("Contraste:", 1-pchisq(sum_res, df_res), "\n")

# --------------------------------------------
# MODELO: Binomial Negativa
# --------------------------------------------
cat("\n\nModelo NB_CT events ~ conc + type + offset(log(at.risk)) \n=================================\n")
modelo_NB_CT <- glm.nb(events ~ conc + type + offset(log(at.risk)), 
                       data = datos_filtrados)

print(summary(modelo_NB_CT))

print(check_overdispersion(modelo_NB_CT))

# --------------------------------------------
# MODELO: Poisson Cero Inflado conc + gender + type | 1
# --------------------------------------------
cat("\n\nModelo ZIP conc + gender + type | 1 \n=================================\n")

modelo_ZIP_CGT_1 <- zeroinfl(events ~ conc + gender + type | 1,
                             offset = log(at.risk),
                             data = datos_filtrados,
                             dist = "poisson", link="logit")

print(summary(modelo_ZIP_CGT_1))

dispersion_modelo_ZIP_CGT_1 <- sum(residuals(modelo_ZIP_CGT_1, type = "pearson")^2) / modelo_ZIP_CGT_1$df.residual
cat("Factor de sobredispersión:", round(dispersion_modelo_ZIP_CGT_1, 2), "\n")

mean(datos_filtrados$events == 0)  # proporción de ceros observados


pred_poisson <- predict(modelo_ZIP_CGT_1, type = "response")
esperados_ceros <- mean(dpois(0, lambda = pred_poisson))
cat("Proporción esperada de ceros bajo Poisson:", esperados_ceros, "\n")


# --------------------------------------------
# MODELO: Poisson Cero Inflado conc + gender + type | gender + type
# --------------------------------------------
cat("\n\nModelo ZIP conc + gender + type | gender + type \n=================================\n")
modelo_ZIP_CGT_GT <- zeroinfl(events ~ conc + gender + type | gender + type,
                              offset = log(at.risk),
                              data = datos_filtrados,
                              dist = "poisson")

print(summary(modelo_ZIP_CGT_GT))


# --------------------------------------------
# MODELO: Poisson Cero Inflado conc | 1
# --------------------------------------------
cat("\n\nModelo ZIP conc | 1 \n=================================\n")
modelo_ZIP_C_1 <- zeroinfl(events ~ conc | 1,
                           offset = log(at.risk),
                           data = datos_filtrados,
                           dist = "poisson")

print(summary(modelo_ZIP_C_1))


# --------------------------------------------
# MODELO: Poisson Cero Inflado conc | type
# --------------------------------------------
cat("\n\nModelo ZIP conc | type\n=================================\n")
modelo_ZIP_C_T <- zeroinfl(events ~ conc | type,
                           offset = log(at.risk),
                           data = datos_filtrados,
                           dist = "poisson")

print(summary(modelo_ZIP_C_T))

# Calculo de los ceros observados y esperados
ceros_observados <- sum(datos_filtrados$events == 0)
prob_cero_estimado <- predict(modelo_ZIP_C_T, type = "prob")[, 1]  # columna 1 = P(Y = 0)
ceros_esperados <- sum(prob_cero_estimado)
cat("Ceros observados:", ceros_observados, "\n")
cat("Ceros esperados por modelo ZIP_C_T:", round(ceros_esperados, 2), "\n")


# --------------------------------------------
# MODELO: Binomial Negativo Cero Inflado conc + gender + type | 1
# --------------------------------------------
cat("\n\nModelo ZINB conc + gender + type | 1\n=================================\n")
modelo_ZINB_CGT_1 <- zeroinfl(events ~ conc + gender + type | 1,
                              offset = log(at.risk),
                              data = datos_filtrados,
                              dist = "negbin")

print(summary(modelo_ZINB_CGT_1))


# --------------------------------------------
# MODELO: Binomial Negativo Cero Inflado conc + gender + type  | gender + type
# --------------------------------------------
cat("\n\nModelo ZINB conc + gender + type | gender + type \n=================================\n")
modelo_ZINB_CGT_GT <- zeroinfl(events ~ conc + gender + type | gender + type,
                               offset = log(at.risk),
                               data = datos_filtrados,
                               dist = "negbin")

print(summary(modelo_ZINB_CGT_GT))

# --------------------------------------------
# MODELO: Binomial Negativo Cero Inflado conc | 1
# --------------------------------------------
cat("\n\nModelo ZINB conc | 1\n=================================\n")
modelo_ZINB_C_1 <- zeroinfl(events ~ conc | 1,
                            offset = log(at.risk),
                            data = datos_filtrados,
                            dist = "negbin")

print(summary(modelo_ZINB_C_1))


# --------------------------------------------
# MODELO: Binomial Negativo Cero Inflado conc |  type
# --------------------------------------------
cat("\n\nModelo ZINB conc | type \n=================================\n")
modelo_ZINB_C_T <- zeroinfl(events ~ conc | type,
                            offset = log(at.risk),
                            data = datos_filtrados,
                            dist = "negbin")

print(summary(modelo_ZINB_C_T))

# Calculo de los ceros observados y esperados
ceros_observados <- sum(datos_filtrados$events == 0)
prob_cero_estimado <- predict(modelo_ZINB_C_T, type = "prob")[, 1]  # columna 1 = P(Y = 0)
ceros_esperados <- sum(prob_cero_estimado)
cat("Ceros observados:", ceros_observados, "\n")
cat("Ceros esperados por modelo ZINB_C_T:", round(ceros_esperados, 2), "\n")


# --------------------------------------------
# MODELO: Comparativa de modelos
# --------------------------------------------
cat("\n\nAIC\n=================================\n")

res = AIC(modelo_all_data, 
          modelo_Poisson_CGT, modelo_Poisson_CT, 
          modelo_ZIP_CGT_1, modelo_ZIP_CGT_GT, modelo_ZIP_C_1, modelo_ZIP_C_T,
          modelo_NB_CGT, modelo_NB_CT, 
          modelo_ZINB_CGT_1, modelo_ZINB_CGT_GT, modelo_ZINB_C_1, modelo_ZINB_C_T)
print(res)

cat("\n\nBIC\n=================================\n")

res = BIC(modelo_all_data, 
          modelo_Poisson_CGT, modelo_Poisson_CT, 
          modelo_ZIP_CGT_1, modelo_ZIP_CGT_GT, modelo_ZIP_C_1, modelo_ZIP_C_T,
          modelo_NB_CGT, modelo_NB_CT,
          modelo_ZINB_CGT_1, modelo_ZINB_CGT_GT, modelo_ZINB_C_1, modelo_ZINB_C_T)
print(res)

cat("\n\nMAE\n=================================\n")
y_real <- datos_filtrados$events
print(mean(abs(y_real - predict(modelo_Poisson_CGT, type = "response"))))    
print(mean(abs(y_real - predict(modelo_ZIP_CGT_1, type = "response")))) 
print(mean(abs(y_real - predict(modelo_ZIP_CGT_GT, type = "response"))))
print(mean(abs(y_real - predict(modelo_NB_CGT, type = "response"))))  
print(mean(abs(y_real - predict(modelo_ZINB_CGT_1, type = "response"))))
print(mean(abs(y_real - predict(modelo_ZINB_CGT_GT, type = "response")))) 
print(mean(abs(y_real - predict(modelo_rob, type = "response"))))   

cat("\n\nAJUSTE RESIDUOS PEARSON\n=================================\n")
print(mean(abs(residuals(modelo_Poisson_CGT, type = "pearson"))))
print(mean(abs(residuals(modelo_ZIP_CGT_1, type = "pearson"))))
print(mean(abs(residuals(modelo_ZIP_CGT_GT, type = "pearson"))))
print(mean(abs(residuals(modelo_NB_CGT, type = "pearson"))))
print(mean(abs(residuals(modelo_ZINB_CGT_1, type = "pearson"))))
print(mean(abs(residuals(modelo_ZINB_CGT_GT, type = "pearson"))))
print(mean(abs(residuals(modelo_rob, type = "pearson"))))


cat("\n\nAJUSTE RESIDUOS PREDICTIVOS\n=================================\n")
print(mean(abs(residuals(modelo_Poisson_CGT, type = "response"))))
print(mean(abs(residuals(modelo_ZIP_CGT_1, type = "response"))))
print(mean(abs(residuals(modelo_ZIP_CGT_GT, type = "response"))))
print(mean(abs(residuals(modelo_NB_CGT, type = "response"))))
print(mean(abs(residuals(modelo_NB_CT, type = "response"))))
print(mean(abs(residuals(modelo_ZINB_CGT_1, type = "response"))))
print(mean(abs(residuals(modelo_ZINB_CGT_GT, type = "response"))))
print(mean(abs(residuals(modelo_rob, type = "response"))))

##############################################
#
# Visualización de residuos
#
##############################################

lista_modelos_nzi <- list(
  AllData = modelo_all_data,
  Poisson_CGT = modelo_Poisson_CGT,
  NB_CGT = modelo_NB_CGT,
  NB_CT = modelo_NB_CT
)

for (nombre in names(lista_modelos_nzi)) {
  print(nombre)
  m  <- lista_modelos_nzi[[nombre]]
  
  disp <- sum(residuals(m, type = "pearson")^2) / m$df.residual
  cat("Factor de sobredispersión:", disp, "\n")
  cat("Pearson's Chi-squared", sum(residuals(m, type = "pearson")^2), m$df.residual, "\n")
}

lista_modelos <- list(
  AllData = modelo_all_data,
  Poisson_CGT = modelo_Poisson_CGT,
  Poisson_CT = modelo_Poisson_CT,
  ZIP_CGT_1 = modelo_ZIP_CGT_1,
  ZIP_CGT_GT = modelo_ZIP_CGT_GT,
  ZIP_C_1 = modelo_ZIP_C_1,
  ZIP_C_T = modelo_ZIP_C_T,
  NB_CGT = modelo_NB_CGT,
  NB_CT = modelo_NB_CT,
  ZINB_CGT_1 = modelo_ZINB_CGT_1,
  ZINB_CGT_GT = modelo_ZINB_CGT_GT,
  ZINB_C_1 = modelo_ZINB_C_1,
  ZINB_C_T = modelo_ZINB_C_T
)

for (nombre in names(lista_modelos)) {
  print(nombre)
  m  <- lista_modelos[[nombre]]
  
  
  # obtener residuos y valores ajustados
  res <- residuals(m, type = "response")
  
  y <- m$y
  mu <- predict(m, type = "response")
  
  res <- ifelse(y == 0,
                -sqrt(2 * mu),
                sign(y - mu) * sqrt(2 * (y * log(y / mu) - (y - mu))))
  
  valores_eje_X = fitted(m)
  #valores_eje_X = datos_filtrados$events
  #valores_eje_X = predict(m, type="response")
  #valores_eje_X = 1:length(res)
  
  plot(x = valores_eje_X, 
       y = res,
       main = paste("Residuos -", nombre),
       xlab = "Valores ajustados",
       ylab = "Residuos de Deviance",
       pch = 20, col = "blue")
  abline(h = 0, col = "red", lty = 2)
  # 
  # res_pearson <- residuals(m, type = "pearson")
  # res_response <- residuals(m, type = "response")
  # fitted_vals <- fitted(m)
  # 
  # # Gráfico 1: Residuos de Pearson vs valores ajustados
  # plot(fitted_vals, res_pearson,
  #      xlab = "Valores ajustados",
  #      ylab = "Residuos de Pearson",
  #      main = "Residuos de Pearson vs Valores ajustados",
  #      pch = 20, col = "blue")
  # abline(h = 0, col = "red", lty = 2)
  # 
  # # Gráfico 2: Residuos de desviancia vs valores ajustados
  # plot(fitted_vals, res_response,
  #      xlab = "Valores ajustados",
  #      ylab = "Residuos de Desviancia",
  #      main = "Residuos de Desviancia vs Valores ajustados",
  #      pch = 20, col = "darkgreen")
  # abline(h = 0, col = "red", lty = 2)
}


# Obtener coeficientes del componente de conteo del modelo ZIP
coefs <- modelo_ZIP_C_T$coefficients$count
intercepto <- coefs["(Intercept)"]
beta_conc <- coefs["conc"]

# Generar valores de concentración
conc <- seq(0, 1000, by = 10)

# Calcular valores de la función mu = exp(intercepto + beta * conc)
mu <- exp(intercepto + beta_conc * conc)*10000

# Graficar
plot(conc, mu,
     type = "l",
     col = "darkgreen",
     lwd = 2,
     xlab = "Concentración de arsénico (conc)",
     ylab = expression(mu), # == 10000 * e^(beta[0] + beta[1] * conc)),
     main = "Tasa esperada según concentración por 10000 habitantes\nZIP_C_T")
