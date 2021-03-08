#' @export
print.estVAS <- function(x, ...) {
  cat("\n")
  cat("Call:", "\n")
  cat("dX\U209C = (\U03B1 - \U03BAX\U209C)dt + \u03C3 dW\U209C", "\n")
  cat("\n")
  cat("Coefficients:", "\n")
  print(round(x$coefficients, 4))
}

#' @export
print.estCEV <- function(x, ...) {
  cat("\n")
  cat("Call:", "\n")
  cat("dX\U209C = (\U03B1 - \U03BAX\U209C)dt + \u03C3 X\U209C\U02E0 dW\U209C", "\n")
  cat("\n")
  cat("Coefficients:", "\n")
  print(round(x$coefficients, 4))
}
