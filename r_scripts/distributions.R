library(gamlss)
library(gamlss.dist)
library(gamlss.add)

distributions <- c("norm", "exp", "gamma", "logis", "beta", "cauchy", "t", "pois", "weibull", "chisq", "geom", "lnorm", "invgamma")
# Definir funciones para verificar requisitos y aplicar transformaciones
dlogistic <- function(x, location, scale) {
  1 / (scale * exp(-((x - location) / scale)) * (1 + exp(-((x - location) / scale)))^2)
}


# Definir la función de densidad para la distribución invgamma
dinvgamma <- function(x, shape, scale) {
  (scale^shape / gamma(shape)) * x^(-shape - 1) * exp(-scale / x)
}
pinvgamma <- function(q, shape, rate) {
  if (length(q) == 0) {
    numeric(0)  # Devolver un vector vacío si q tiene longitud cero
  } else {
    pgamma(1/q, shape, scale = 1/rate)
  }
}

dpoisson<-function(k,lambda){
    dpois(k, lambda)
}

cumple_requisitos <- function(datos, dist_name) {
  switch(
    dist_name,
    norm = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos))
    },
    exp = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos >= 0)
    },
    gamma = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos >= 0)
    },
    logis = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos))
    },
    beta = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos >= 0) && all(datos <= 1)
    },
    cauchy = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos))
    },
    t = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos))
    },
    poisson = {
      is.integer(datos) && !any(is.na(datos)) && all(datos >= 0)
    },
    pois = {
      is.integer(datos) && !any(is.na(datos)) && all(datos >= 0) && abs(mean(datos) - sd(datos)) < 0.1 * mean(datos)
    },
    weibull = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos > 0)
    },
    chisq = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos >= 0)
    },
    logistic = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos))
    },
    geom = {
      is.integer(datos) && !any(is.na(datos)) && all(datos > 0)
    },
    lnorm = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos >= 0)
    },
    weibull2 = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos > 0)
    },
    invgamma = {
      is.numeric(datos) && !any(is.na(datos)) && !any(is.infinite(datos)) && all(datos > 0)
    },
    stop("Distribución no reconocida")
  )
}
coth = function(x) cosh(x) / sinh(x)
sec = function(x) 1 / cos(x)
csc = function(x) 1 / sin(x)
log_squared = function(x) (log(x))^2  # Square of Logarithm
transform_functions <- list(
  original = function(x) x,  # Original
  scaled = function(x) x / sum(x),  # Scaled
  naturallog = function(x) log(x),  # Natural Logarithm
  raiz_cuadrada = function(x) sqrt(x),  # Square Root
  exponencial = function(x) exp(x),  # Exponential
  tanh = function(x) tanh(x),  # Hyperbolic Tangent
  cuadratica = function(x) x^2,  # Quadratic
  arcsin = function(x) asin(sqrt(x)),  # Arcsine
  exp_negativa = function(x) exp(-x),  # Negative Exponential
  arcotangente = function(x) atan(x),  # Arctangent
  logistica_inversa = function(x) log(1 / x - 1),  # Inverse Logistic
  asinh = function(x) asinh(x),  # Inverse Hyperbolic Sine
  cuadratica_inversa = function(x) sqrt(1 / x),  # Inverse Quadratic
  logistica_cuadrada = function(x) log(x^2 / (1 - x^2)),  # Square Logistic
  hip_sin = function(x) sinh(x),  # Hyperbolic Sine
  hip_tanh = function(x) tanh(x),  # Hyperbolic Tangent
  log_doble = function(x) log(log(x + 1) + 1),  # Double Logarithm
  reciproco = function(x) 1 / x,  # Reciprocal
  cubic_root = function(x) x^(1/3),  # Cubic Root
  sinh_inversa = function(x) asinh(1/x),  # Inverse Hyperbolic Sine
  logistica_doble = function(x) log(1 / (1 - x)),  # Double Logistic
  arcoseno_tangente = function(x) asin(tanh(x)),  # Arcsine of Hyperbolic Tangent
  seno_hiperbolico = function(x) sinh(asin(x)),  # Hyperbolic Sine of Arcsine
  raiz_cuarta = function(x) sqrt(sqrt(x)),  # Fourth Root
  hip_cuadrada = function(x) sinh(x)^2,  # Square of Hyperbolic Sine
  cosecante = function(x) 1/sin(x),  # Cosecant
  cuadrado_negativo = function(x) -x^2,  # Negative Square
  tanh_inversa = function(x) atanh(x),  # Inverse Hyperbolic Tangent
  raiz_quinta = function(x) x^(1/5),  # Fifth Root
  exp_doble = function(x) exp(exp(x)),  # Double Exponential
  hip_cubica = function(x) sinh(x)^3,  # Cubic Hyperbolic Sine
  cube = function(x) x^3,  # Cube
  sqrt_abs = function(x) sqrt(abs(x)),  # Square Root of Absolute Value
  log10 = function(x) log10(x),  # Logarithm Base 10
  inverse_hyperbolic_cosine = function(x) acosh(x),  # Inverse Hyperbolic Cosine
  reciprocal_sqrt = function(x) 1 / sqrt(x),  # Reciprocal Square Root
  hyperbolic_tangent_inverse = function(x) atanh(x),  # Inverse Hyperbolic Tangent
  inverse_sine_hyperbolic = function(x) 1 / sinh(x),  # Inverse Hyperbolic Sine
  exp_sqrt = function(x) exp(sqrt(x)),  # Exponential of Square Root
  cosecanthyperbolic = function(x) 1 / sinh(x),  # Cosecant Hyperbolic
  asinh_squared = function(x) asinh(x)^2,  # Square of Inverse Hyperbolic Sine
  exp_abs_squared = function(x) exp((abs(x))^2),  # Exponential of Absolute Value Squared
  cbrt_abs = function(x) sign(x) * abs(x)^(1/3),  # Cube Root of Absolute Value
  coth = function(x) cosh(x) / sinh(x),  # Hyperbolic Cotangent
  acoth = function(x) atanh(1/x),  # Inverse Hyperbolic Cotangent
  reciprocal_squared = function(x) (1 / x)^2,  # Squared Reciprocal
  tan_squared = function(x) tan(x)^2,  # Tangent Squared
  coth_squared = function(x) coth(x)^2,  # Square of Hyperbolic Cotangent
  atanh_inverse = function(x) atanh(x),  # Inverse Hyperbolic Tangent
  atan_abs = function(x) atan(abs(x)),  # Inverse Tangent of Absolute Value
  squared_cubert = function(x) (x^(1/3))^2,  # Squared Cube Root
  log_abs = function(x) log(abs(x) + 1),  # Logarithm of Absolute Value
  tan_log = function(x) tan(log(abs(x) + 1)),  # Tangent of Logarithm
  sqrt_exp = function(x) sqrt(exp(x)),  # Square Root of
  inv_log_squared = function(x) 1 / (log(abs(x) + 1))^2,  # Inverse Logarithm Squared
  log_acosh = function(x) log(acosh(x)),  # Logarithm of Inverse Hyperbolic Cosine
  sqrt_abs_quartic = function(x) sqrt(abs(x)^4),  # Square Root of Absolute Value to the Fourth Power
  tanh_inv_abs = function(x) tanh(1/abs(x)),  # Hyperbolic Tangent of the Inverse Absolute Value
  exp_inv_abs_quartic = function(x) exp(1/abs(x)^4),  # Exponential of the Inverse Absolute Value to the Fourth Power
  log_abs_cubed = function(x) log(abs(x)^3),  # Logarithm of the Absolute Value to the Third Power
  tan_sqrt_abs = function(x) tan(sqrt(abs(x))),  # Tangent of the Square Root of the Absolute Value
  csc_abs = function(x) 1/sin(abs(x)),  # Cosecant of the Absolute Value
  sec_squared = function(x) sec(x)^2,  # Square of Secant
  log_abs_fifth = function(x) log(abs(x)^5),  # Logarithm of the Absolute Value to the Fifth Power
  inv_csc = function(x) 1/csc(x),  # Inverse Cosecant
  sqrt_cube = function(x) sqrt(x^3),  # Square Root of Cube
  exp_div_square = function(x) exp(x / x^2),  # Exponential divided by Square
  atan_cube_root = function(x) atan(x^(1/3)),  # Arctangent of Cube Root
  cos_log = function(x) cos(log(abs(x) + 1)),  # Cosine of Logarithm
  exp_tanh = function(x) exp(tanh(x)),  # Exponential of Hyperbolic Tangent
  sinh_logistic = function(x) sinh(log(1 / x - 1)),  # Hyperbolic Sine of Inverse Logistic
  atan_cubic = function(x) atan(x^3),  # Arctangent of Cube
  tanh_log_inverse = function(x) tanh(1 / log(x + 1)),  # Hyperbolic Tangent of Inverse Logarithm
  sqrt_exp_squared = function(x) sqrt(exp(x^2)),  # Square Root of Exponential Squared

  exp_log_inv = function(x) exp(log(1/x)),  # Exponential of Inverse Logarithm
  tanh_exp = function(x) tanh(exp(x)),  # Hyperbolic Tangent of Exponential
  atan_log_exp = function(x) atan(log(exp(x))),  # Arctangent of Logarithm of Exponential
  log_exp_cube = function(x) log(exp(x)^3),  # Logarithm of Exponential Cubed
  sinh_inv_tanh = function(x) sinh(1 / tanh(x)),  # Hyperbolic Sine of Inverse Hyperbolic Tangent
  atan_cosh_squared = function(x) atan(cosh(x)^2),  # Arctangent of Square of Hyperbolic Cosine
  log_sqrt_exp = function(x) log(sqrt(exp(x))),  # Logarithm of Square Root of Exponential
  exp_tanh_squared = function(x) exp(tanh(x)^2),  # Exponential of Square of Hyperbolic Tangent
  asin_cubic_root = function(x) asin(x^(1/3)),  # Arcsine of Cube Root

  cos_inv_tan = function(x) cos(1 / tan(x)),  # Cosine of Inverse Tangent
  atan_tanh_inv_exp = function(x) atan(tanh(1 / exp(x))),  # Arctangent of Inverse Hyperbolic Tangent of Exponential
  exp_cos_inv_tanh = function(x) exp(cos(1 / tanh(x))),  # Exponential of Cosine of Inverse Hyperbolic Tangent
  atan_coth_exp = function(x) atan(coth(exp(x))),  # Arctangent of Hyperbolic Cotangent of Exponential
  sinh_inv_cosh_squared = function(x) sinh(1 / cosh(x)^2),  # Hyperbolic Sine of Inverse Square of Hyperbolic Cosine
  atan_log_inv_cosh = function(x) atan(log(1 / cosh(x))),  # Arctangent of Inverse Hyperbolic Cosine
  sinh_cubed_exp = function(x) sinh(x)^3 * exp(x),  # Cubic Hyperbolic Sine times Exponential
  atan_tanh_log_inv = function(x) atan(tanh(log(1/x))),  # Arctangent of Hyperbolic Tangent of Inverse Logarithm
  exp_sqrt_tan = function(x) exp(sqrt(tan(x))),  # Exponential of Square Root of Tangent
  asin_exp_cubert = function(x) asin(exp(x)^(1/3)),  # Arcsine of Exponential Cube Root

  atan_tanh_exp_sqrt = function(x) atan(tanh(exp(sqrt(x)))),  # Arctangent of Hyperbolic Tangent of Exponential Square Root
  sinh_inv_tan_cosh = function(x) sinh(1 / tan(cosh(x))),  # Hyperbolic Sine of Inverse Tangent of Hyperbolic Cosine
  atan_cosec = function(x) atan(1 / sin(x)),  # Arctangent of Cosecant
  sinh_cube_log = function(x) sinh(x)^3 * log(x + 1),  # Cubic Hyperbolic Sine times Logarithm
  exp_tan_squared_cubert = function(x) exp(tan(x)^2)^(1/3),  # Exponential of Square of Tangent Cube Root
  log_exp_cos = function(x) log(exp(cos(x))),  # Logarithm of Exponential of Cosine
  sinh_inv_tan_log = function(x) sinh(1 / tan(log(x))),  # Hyperbolic Sine of Inverse Tangent of Logarithm
  atan_cubert_exp = function(x) atan(x^(1/3) * exp(x)),  # Arctangent of Cube Root times Exponential
  sinh_inv_tan_exp = function(x) sinh(1 / tan(exp(x))),  # Hyperbolic Sine of Inverse Tangent of Exponential
  atan_tan_cube = function(x) atan(tan(x)^3),  # Arctangent of Cubic Tangent
  log_cosec_exp = function(x) log(1 / sin(x) * exp(x)),  # Logarithm of Cosecant times Exponential
  
  sqrt_cos_abs = function(x) sqrt(cos(abs(x))),  # Square Root of Cosine of Absolute Values
  exp_tanh_cubert = function(x) exp(tanh(x)^(1/3)),  # Exponential of Cube Root of Hyperbolic Tangent
  atan_log_inv_cube = function(x) atan(log(1 / x)^3),  # Arctangent of Inverse Logarithm Cubed
  coth_sqrt_abs = function(x) coth(sqrt(abs(x))),  # Hyperbolic Cotangent of Square Root of Absolute Values
  log_inv_sqrt_exp = function(x) log(1 / sqrt(exp(x))),  # Logarithm of Inverse Square Root of Exponential
  sinh_abs_log = function(x) sinh(abs(log(x + 1))),  # Hyperbolic Sine of Absolute Logarithm
  atan_cosh_inv = function(x) atan(cosh(1 / x)),  # Arctangent of Hyperbolic Cosine of Inverse
  tanh_cube_inv_log = function(x) tanh(x^3) / log(x + 1),  # Hyperbolic Tangent of Cube times Inverse Logarithm
  sqrt_cubert_exp = function(x) sqrt(x^(1/3) * exp(x)),  # Square Root of Cube Root times Exponential
  

  exp_inv_cosec = function(x) exp(1 / sin(x)),  # Exponential of Inverse Cosecant
  coth_inv_cubert = function(x) coth(1 / x^(1/3)),  # Hyperbolic Cotangent of Inverse Cube Root
  sinh_tan_exp = function(x) sinh(tan(exp(x))),  # Hyperbolic Sine of Tangent of Exponential
  exp_cos_cube = function(x) exp(cos(x)^3),  # Exponential of Cube of Cosine
  atan_log_sqrt_abs = function(x) atan(log(sqrt(abs(x)))),  # Arctangent of Logarithm of Square Root of Absolute Values
  sinh_tan_inv_cosec = function(x) sinh(tan(1 / sin(x))),  # Hyperbolic Sine of Inverse Tangent of Inverse Cosecant
  atan_exp_cubert = function(x) atan(exp(x)^(1/3)),  # Arctangent of Exponential Cube Root
  sinh_log_cube = function(x) sinh(log(x + 1)^3),  # Hyperbolic Sine of Logarithm Cubed
  atan_cosec_squared = function(x) atan(1 / sin(x)^2),  # Arctangent of Square of Cosecant
  coth_sqrt_inv_exp = function(x) coth(sqrt(1 / exp(x))),  # Hyperbolic Cotangent of Square Root of Inverse Exponential
  

  sinh_inv_cubert_exp = function(x) sinh(1 / x^(1/3) * exp(x)),  # Hyperbolic Sine of Inverse Cube Root times Exponential
  atan_exp_tan_inv = function(x) atan(exp(tan(1 / x))),  # Arctangent of Exponential times Inverse Tangent
  sinh_cube_tanh = function(x) sinh(x)^3 * tanh(x),  # Cubic Hyperbolic Sine times Hyperbolic Tangent
  atan_tan_cosec_inv = function(x) atan(tan(1 / sin(x))),  # Arctangent of Tangent times Inverse Cosecant
  coth_inv_cosec = function(x) coth(1 / sin(x)),  # Hyperbolic Cotangent of Inverse Cosecant
  exp_inv_cosec_squared = function(x) exp(1 / sin(x)^2),  # Exponential of Inverse Square of Cosecant
  sinh_tan_cosec_inv = function(x) sinh(tan(1 / sin(x))),  # Hyperbolic Sine of Tangent times Inverse Cosecant
  atan_cosec_cube = function(x) atan(1 / sin(x)^3),  # Arctangent of Cube of Cosecant
  exp_tan_cube = function(x) exp(tan(x)^3),  # Exponential of Cube of Tangent
  sinh_log_sqrt_inv_exp = function(x) sinh(log(sqrt(1 / exp(x)))),  # Hyperbolic Sine of Logarithm of Square Root of Inverse Exponential
  atan_tanh_cube_root = function(x) atan(tanh(x)^(1/3)),  # Arctangent of Cube Root of Hyperbolic Tangent
  

  sinh_inv_tanh_exp = function(x) sinh(1 / tanh(exp(x))),  # Hyperbolic Sine of Inverse Hyperbolic Tangent of Exponential
  exp_log_cubert = function(x) exp(log(x + 1)^(1/3)),  # Exponential of Logarithm Cubed
  tanh_cubert_inv_exp = function(x) tanh(x^(1/3)) / exp(x),  # Hyperbolic Tangent of Cube Root times Inverse Exponential
  atan_inv_tanh_exp = function(x) atan(1 / tanh(exp(x))),  # Arctangent of Inverse Hyperbolic Tangent of Exponential
  sinh_exp_cube_root = function(x) sinh(exp(x)^(1/3)),  # Hyperbolic Sine of Exponential Cube Root
  atan_tanh_inv_cosec = function(x) atan(tanh(1 / sin(x))),  # Arctangent of Hyperbolic Tangent of Inverse Cosecant
  sinh_inv_tan_cube = function(x) sinh(1 / tan(x)^3),  # Hyperbolic Sine of Inverse Cube of Tangent
  atan_log_inv_cosec = function(x) atan(log(1 / sin(x))),  # Arctangent of Inverse Logarithm of Cosecant
  exp_inv_cosec_cube = function(x) exp(1 / sin(x)^3),  # Exponential of Inverse Cube of Cosecant
  tanh_cosec_cube = function(x) tanh(1 / sin(x)^3),  # Hyperbolic Tangent of Inverse Cube of Cosecant
  sqrt_log_abs = function(x) sqrt(log(abs(x))),                  # Square Root of the Logarithm of the Absolute Value
  exp_tanh_inv_sqrt = function(x) exp(tanh(1 / sqrt(x))),        # Exponential of the Inverse Hyperbolic Tangent of the Square Root
  cos_log_inv_cube = function(x) cos(log(1 / x)^3),              # Cosine of the Logarithm of the Inverse Cube
  sqrt_exp_cube = function(x) sqrt(exp(x)^3),                    # Square Root of the Cube of the Exponential
  log_inv_tanh_cube = function(x) log(1 / tanh(x)^3),            # Logarithm of the Inverse Hyperbolic Tangent Cubed
  sqrt_tanh_inv_exp = function(x) sqrt(tanh(1 / exp(x))),        # Square Root of the Inverse Hyperbolic Tangent of the Exponential
  cos_cube_inv_sqrt = function(x) cos(x^3) / sqrt(x),            # Cosine Cubed divided by the Square Root
  exp_cube_sqrt = function(x) exp(x^3)^(1/2),                    # Exponential of the Cube Square Root
  sqrt_tanh_exp = function(x) sqrt(tanh(exp(x))),                # Square Root of the Hyperbolic Tangent of the Exponential
  log_inv_cosec_inv_sqrt = function(x) log(1 / (1 / sin(x))^(1/2)), # Logarithm of the Inverse Cosecant divided by the Square Root
  
  sinh_exp_sqrt = function(x) sinh(exp(sqrt(x))),               # Hyperbolic Sine of the Exponential Square Root
  tan_inv_cosec_exp = function(x) tan(1 / sin(exp(x))),          # Tangent of the Inverse Cosecant of the Exponential
  log_tanh_cube_inv = function(x) log(tanh(x)^3 / sqrt(1 / x)),  # Logarithm of the Hyperbolic Tangent Cubed divided by the Square Root
  sinh_inv_cubert_exp = function(x) sinh(1 / x^(1/3) * exp(x)),  # Hyperbolic Sine of the Inverse Cube Root times Exponential
  atan_tan_cosec = function(x) atan(tan(1 / sin(x))),            # Arctangent of the Tangent of the Inverse Cosecant
  log_inv_tanh_cube = function(x) log(1 / tanh(x)^3),            # Logarithm of the Inverse Hyperbolic Tangent Cubed
  sinh_cosec_sqrt = function(x) sinh(1 / sin(x))^(1/2),          # Hyperbolic Sine of the Inverse Cosecant divided by the Square Root
  atan_exp_cosec_inv = function(x) atan(exp(1 / sin(x))),        # Arctangent of the Exponential of the Inverse Cosecant
  sqrt_tan_cube = function(x) sqrt(tan(x)^3),                    # Square Root of the Tangent Cubed
  sinh_cubert_exp = function(x) sinh(x^(1/3) * exp(x)),          # Hyperbolic Sine of the Cube Root times Exponential
  
  exp_tan_inv_cosec = function(x) exp(tan(1 / sin(x))),         # Exponential of the Tangent of the Inverse Cosecant
  log_inv_tan_cube = function(x) log(1 / tan(x)^3),             # Logarithm of the Inverse Tangent Cubed
  tanh_exp_cube = function(x) tanh(exp(x)^3),                   # Hyperbolic Tangent of the Exponential Cubed
  sqrt_tan_inv_exp = function(x) sqrt(tan(1 / exp(x))),         # Square Root of the Inverse Tangent of the Exponential
  cos_exp_tanh = function(x) cos(exp(tanh(x))),                 # Cosine of the Exponential of the Hyperbolic Tangent
  exp_tanh_inv_cosec = function(x) exp(tanh(1 / sin(x))),       # Exponential of the Inverse Hyperbolic Tangent of the Cosecant
  log_inv_cosec_cube = function(x) log(1 / sin(x)^3),           # Logarithm of the Inverse Cosecant Cubed
  sinh_exp_cosec_inv = function(x) sinh(exp(1 / sin(x))),       # Hyperbolic Sine of the Exponential of the Inverse Cosecant
  exp_tanh_cube = function(x) exp(tanh(x)^3),                   # Exponential of the Hyperbolic Tangent Cubed
  cos_inv_tan_exp = function(x) cos(1 / tan(exp(x))),           # Cosine of the Inverse Tangent of the Exponential
  
  exp_cube_cosec = function(x) exp(x^3) / sin(x),               # Exponential of the Cube divided by the Cosecant
  sinh_cube_tanh = function(x) sinh(x)^3 * tanh(x),             # Hyperbolic Sine Cubed times the Hyperbolic Tangent
  tanh_inv_cube_sqrt = function(x) tanh(1 / x^3)^(1/2),         # Hyperbolic Tangent of the Inverse Cube divided by the Square Root
  sinh_exp_tanh_cube = function(x) sinh(exp(tanh(x)^3)),        # Hyperbolic Sine of the Exponential of the Hyperbolic Tangent Cubed
  exp_inv_tanh_cube = function(x) exp(1 / tanh(x)^3),           # Exponential of the Inverse Hyperbolic Tangent Cubed
  log_cosec_cube = function(x) log(1 / sin(x)^3),               # Logarithm of the Inverse Cosecant Cubed
  sinh_cube_tan_inv = function(x) sinh(x)^3 / tan(1 / x),       # Hyperbolic Sine Cubed divided by the Inverse Tangent
  tanh_inv_cube_exp = function(x) tanh(1 / x^3)^exp(x),         # Hyperbolic Tangent of the Inverse Cube raised to the Exponential
  exp_tan_cube_sqrt = function(x) exp(tan(x)^3)^(1/2),         # Exponential of the Tangent Cubed Square Root
  tanh_cubert_cosec = function(x) tanh(x^(1/3)) / sin(x),    
  
  sinh_exp_tanh = function(x) sinh(exp(tanh(x))),                     # Hyperbolic Sine of the Exponential of the Hyperbolic Tangent
  log_inv_cube_tanh = function(x) log(1 / x^3) / tanh(x),             # Logarithm of the Inverse Cube divided by the Hyperbolic Tangent
  exp_cube_inv_tanh = function(x) exp(x^3) / tanh(1 / x),             # Exponential of the Cube divided by the Inverse Hyperbolic Tangent
  tanh_cube_inv_sqrt = function(x) tanh(x^3)^(1/2),                   # Hyperbolic Tangent of the Cube divided by the Square Root
  sqrt_cubert_exp = function(x) sqrt(x^(1/3) * exp(x)),               # Square Root of the Cube Root times the Exponential
  exp_inv_cube_cosec = function(x) exp(1 / x^3) / sin(x),             # Exponential of the Inverse Cube divided by the Cosecant
  tanh_exp_cubert = function(x) tanh(exp(x)^(1/3)),                   # Hyperbolic Tangent of the Exponential Cube Root
  sinh_cubert_inv_exp = function(x) sinh(x^(1/3)) / exp(x),           # Hyperbolic Sine of the Cube Root divided by the Exponential
  sqrt_tan_inv_cosec = function(x) sqrt(tan(1 / sin(x))),             # Square Root of the Inverse Tangent of the Cosecant
  log_exp_inv_tanh = function(x) log(exp(1 / tanh(x))),               # Logarithm of the Exponential of the Inverse Hyperbolic Tangent
  
  exp_cube_inv_cubert = function(x) exp(x^3) / x^(1/3),               # Exponential of the Cube divided by the Cube Root
  sinh_cube_tan_inv_sqrt = function(x) sinh(x)^3 / sqrt(1 / tan(x)),  # Hyperbolic Sine Cubed divided by the Square Root of the Inverse Tangent
  exp_sqrt_cubert = function(x) exp(sqrt(x)^(1/3)),                   # Exponential of the Square Root of the Cube Root
  tanh_inv_cubert_exp = function(x) tanh(1 / x^(1/3))^exp(x),         # Hyperbolic Tangent of the Inverse Cube Root raised to the Exponential
  exp_inv_cosec_inv_sqrt = function(x) exp(1 / sin(x))^(1/2),         # Exponential of the Inverse Cosecant divided by the Square Root
  sinh_cosec_cube = function(x) sinh(1 / sin(x)^3),                   # Hyperbolic Sine of the Inverse Cosecant Cubed
  tan_inv_exp_cosec = function(x) tan(1 / exp(sin(x))),               # Tangent of the Inverse Exponential of the Cosecant
  exp_tan_cube_inv_sqrt = function(x) exp(tan(x)^3)^(1/2),           # Exponential of the Tangent Cubed divided by the Square Root
  tanh_cosec_cube_inv = function(x) tanh(1 / sin(x))^3 / sqrt(x),    # Hyperbolic Tangent of the Inverse Cosecant Cubed divided by the Square Root
  exp_cube_inv_cosec_inv_sqrt = function(x) exp(x^3) / (1 / sin(x))^(1/2), # Exponential of the Inverse Cube divided by the Inverse Cosecant divided by the Square Root
  
  sinh_cube_tan_cube_inv = function(x) sinh(x)^3 / tan(x)^3,         # Hyperbolic Sine Cubed divided by the Tangent Cubed
  exp_tanh_cosec_inv_sqrt = function(x) exp(tanh(1 / sin(x)))^(1/2),  # Exponential of the Hyperbolic Tangent of the Inverse Cosecant divided by the Square Root
  log_tan_cube_inv_cosec = function(x) log(tan(x)^3 / 1 / sin(x)),   # Logarithm of the Tangent Cubed divided by the Inverse Cosecant
  tanh_inv_exp_cosec_inv = function(x) tanh(1 / exp(1 / sin(x))),    # Hyperbolic Tangent of the Inverse Exponential of the Inverse Cosecant
  sqrt_cosec_exp_inv_tan = function(x) sqrt(1 / sin(x)) / exp(1 / tan(x)), # Square Root of the Inverse Cosecant divided by the Exponential of the Inverse Tangent
  tanh_cubert_inv_exp = function(x) tanh(x^(1/3)) / exp(x),           # Hyperbolic Tangent of the Cube Root divided by the Exponential
  exp_cube_sqrt_tan_inv = function(x) exp(x^3)^(1/2) / tan(1 / x),    # Exponential of the Cube Square Root divided by the Inverse Tangent
  sinh_cosec_cube_inv_sqrt = function(x) sinh(1 / sin(x))^3 / sqrt(x), # Hyperbolic Sine of the Inverse Cosecant Cubed divided by the Square Root
  log_inv_tan_exp_cube = function(x) log(1 / tan(exp(x))^3),         # Logarithm of the Inverse Tangent of the Exponential Cubed
  tan_cube_sqrt_cosec = function(x) tan(x)^3^(1/2) / sin(x),          # Tangent Cubed Square Root divided by the Cosecant
  
  log_squared_inv_tanh = function(x) log_squared(1 / tanh(x)),                # log_squared of the Inverse Hyperbolic Tangent
  exp_log_squared = function(x) exp(log_squared(x)),                          # Exponential of log_squared
  log_squared_cube_inv_sqrt = function(x) log_squared(x^3) / sqrt(x),         # log_squared of the Cube divided by the Square Root
  sqrt_log_squared_exp = function(x) sqrt(log_squared(exp(x))),              # Square Root of log_squared of the Exponential
  tanh_log_squared_inv_exp = function(x) tanh(log_squared(1 / exp(x))),       # Hyperbolic Tangent of log_squared of the Inverse Exponential
  log_squared_cubert_inv_tanh = function(x) log_squared(x^(1/3)) / tanh(x),   # log_squared of the Cube Root divided by the Hyperbolic Tangent
  exp_inv_sqrt_log_squared = function(x) exp(1 / sqrt(log_squared(x))),      # Exponential of the Inverse Square Root of log_squared
  log_squared_inv_cosec_sqrt = function(x) log_squared(1 / sin(x))^(1/2),     # log_squared of the Inverse Cosecant divided by the Square Root
  sinh_log_squared_cube = function(x) sinh(log_squared(x)^3),                # Hyperbolic Sine of log_squared cubed
  log_squared_exp_inv_cosec = function(x) log_squared(exp(1 / sin(x))),        # log_squared of the Exponential of the Inverse Cosecant
  
  log_squared_inv_exp_sqrt = function(x) log_squared(1 / exp(sqrt(x))),           # log_squared of the Inverse Exponential Square Root
  tanh_log_squared_cube = function(x) tanh(log_squared(x)^3),                   # Hyperbolic Tangent of log_squared cubed
  exp_sqrt_log_squared_inv_cosec = function(x) exp(sqrt(log_squared(1 / sin(x)))), # Exponential of the Square Root of log_squared of the Inverse Cosecant
  log_squared_cosec_exp_inv_tanh = function(x) log_squared(1 / sin(exp(1 / tanh(x)))), # log_squared of the Inverse Cosecant Exponential of the Inverse Hyperbolic Tangent
  log_squared_inv_tanh_cube_sqrt = function(x) log_squared(1 / tanh(x)^3)^(1/2),   # log_squared of the Inverse Hyperbolic Tangent Cubed divided by the Square Root
  exp_tanh_log_squared = function(x) exp(tanh(log_squared(x))),                 # Exponential of the Hyperbolic Tangent of log_squared
  sqrt_exp_log_squared_inv = function(x) sqrt(exp(log_squared(1 / x))),         # Square Root of the Exponential of log_squared of the Inverse
  log_squared_inv_cube_sqrt = function(x) log_squared(1 / x^3)^(1/2),           # log_squared of the Inverse Cube divided by the Square Root
  tanh_exp_log_squared_cube = function(x) tanh(exp(log_squared(x)^3)),          # Hyperbolic Tangent of the Exponential of log_squared cubed
  exp_tanh_inv_log_squared = function(x) exp(tanh(1 / log_squared(x))),         # Exponential of the Inverse Hyperbolic Tangent of log_squared
  log_squared_inv_tan_cube_sqrt = function(x) log_squared(1 / tan(x)^3)^(1/2),   # log_squared of the Inverse Tangent Cubed divided by the Square Root
  sqrt_log_squared_exp_inv_cosec = function(x) sqrt(log_squared(exp(1 / sin(x)))), # Square Root of log_squared of the Exponential of the Inverse Cosecant
  tanh_inv_exp_log_squared = function(x) tanh(1 / exp(log_squared(x))),         # Hyperbolic Tangent of the Inverse Exponential of log_squared
  log_squared_cosec_cube_inv_exp = function(x) log_squared(1 / sin(x)^3) / exp(x), # log_squared of the Inverse Cosecant Cubed divided by the Exponential
  exp_inv_tanh_log_squared_cube = function(x) exp(1 / tanh(log_squared(x)^3)),  # Exponential of the Inverse Hyperbolic Tangent of log_squared cubed
  log_squared_cube_inv_sqrt_exp = function(x) log_squared(x^3) / sqrt(exp(x)),   # log_squared of the Cube divided by the Square Root of the Exponential
  sqrt_log_squared_inv_exp_cosec = function(x) sqrt(log_squared(1 / exp(1 / sin(x)))), # Square Root of log_squared of the Inverse Exponential of the Inverse Cosecant
  tanh_inv_log_squared_cube_exp = function(x) tanh(1 / log_squared(x)^(1/3))^exp(x), # Hyperbolic Tangent of the Inverse log_squared Cube Root raised to the Exponential
  exp_log_squared_inv_cosec_cube = function(x) exp(log_squared(1 / sin(x)))^3      # Exponential of log_squared of the Inverse Cosecant cubed

)
#Kullback-Leibler (KL) Divergence
kl_divergence <- function(data) {
  p_data <- density(data)
  p_norm <- dnorm(p_data$x, mean = mean(data), sd = sd(data))
  kl_div <- sum(p_data$y * log(p_data$y / p_norm))
  return(kl_div)
}


# Ajustar distribución y aplicar transformación si es necesario
ajustar_distribucion <- function(datos, dist_name) {
    val<-tryCatch(
      {
        # Ajustar a la distribución especificada
        if(dist_name=="t" || dist_name=="chisq"){
          dist_fit <- fitdist(datos, dist_name, start=list(df=2), discrete = F)
        }else if(dist_name=="pois"){
          dist_fit <- fitdist(datos, "pois", start=list(lambda=mean(datos)), method = "mse", discrete = F)
        }else if(dist_name=="weibull"){
          start  <- list(shape=1, scale=1)
          dist_fit <- fitdist(datos, "weibull", start=start, method="mle", discrete = F)
        }
        else if(dist_name=="geom"){
          dist_fit <- fitdist(datos, "geom", discrete = F, method = "mme")
        }else{
          dist_fit <- fitdist(datos, dist_name, discrete = F)
        }
        return(dist_fit)
      },
      error = function(e) {
        print(e)
        return(100000)
      },
      finally = function(){
      return(100000)
      }
    )
    return(val)
}

# Crear una función para ajustar y obtener el AIC de cada distribución
obtener_fit <- function(dist_name, datos=NA) {
  
  print(dist_name)
  transformed_data<-datos
  dist_fit<-10000
  for(i in 1:length(transform_functions)){
    print(i)
    tfunction<-transform_functions[[i]]
    transformed_data <- tfunction(datos)
    if (cumple_requisitos(transformed_data, dist_name)) {
      dist_fit<-ajustar_distribucion(transformed_data, dist_name)
    }
    if(class(dist_fit)=="fitdist"){
      break
    }
  }
  return(dist_fit)

  
}

test_estadistico <- function(datos, dist_name, dist_fit) {
  
  return(switch(
    dist_name,
    norm = ks.test(datos, "pnorm", mean = dist_fit$estimate[1], sd = dist_fit$estimate[2]),
    exp = ks.test(datos, "pexp", rate = dist_fit$estimate[1]),
    gamma = ks.test(datos, "pgamma", shape = dist_fit$estimate[1], rate = dist_fit$estimate[2]),
    logis = ks.test(datos, "plogis", location = dist_fit$estimate[1], scale = dist_fit$estimate[2]),
    beta = ks.test(datos, "pbeta", shape1 = dist_fit$estimate[1], shape2 = dist_fit$estimate[2]),
    cauchy = ks.test(datos, "pcauchy", location = dist_fit$estimate[1], scale = dist_fit$estimate[2]),
    t = ks.test(datos, "pt", df = dist_fit$estimate[1], ncp = 0),
    pois = chisq.test(tabulate(datos), p = ppois(unique(datos), lambda = dist_fit$estimate[1])),
    weibull = ks.test(datos, "pweibull", shape = dist_fit$estimate[1], scale = dist_fit$estimate[2]),
    chisq = chisq.test(datos, p = pchisq(datos, df = dist_fit$estimate[1])/sum(pchisq(datos, df = dist_fit$estimate[1]))),
    geom = chisq.test(tabulate(datos), p = pgeom(unique(datos), prob = dist_fit$estimate[1])),
    lnorm = ks.test(datos, "plnorm", meanlog = dist_fit$estimate[1], sdlog = dist_fit$estimate[2]),
    invgamma = ks.test(datos, "pinvgamma", shape = dist_fit$estimate[1], rate = dist_fit$estimate[2]),
    stop("Distribución no reconocida")
  ))
}

permutation_test <- function(x, y, rep=1000) {
  observed_diff <- mean(x) - mean(y)
  combined_data <- c(x, y)
  num_permutations <- rep  # Número de permutaciones (ajusta según tus necesidades)
  perm_diffs <- replicate(num_permutations, {
    shuffled_data <- sample(combined_data)
    mean(shuffled_data[1:length(x)]) - mean(shuffled_data[(length(x) + 1):length(combined_data)])
  })
  p_value <- sum(abs(perm_diffs) >= abs(observed_diff)) / num_permutations
  return(p_value)
}


peardist<-function(x,distance){
  #as.dist(1-cor(t(log(x + min(x[x>0])))))
  as.dist(1-cor(t(x),use = "pairwise.complete.obs",method = "pearson"))
}

braydist <- function(x,distance) {
  # Calculate Bray-Curtis dissimilarity
  n <- ncol(x)
  min_values <- apply(x, 2, function(col) min(col[col > 0]))
  x <- t(log(x + min_values))
  numerator <- apply(x, 1, function(row) sum(abs(row - min(row)))/sum(row + min(row)))
  
  # Convert to distance matrix format
  as.dist(1 - numerator)
}
tsd<-function(data,trim=0.1){
  sd(data[data >= quantile(data, trim / 2) & data <= quantile(data, 1 - trim / 2)])
}
