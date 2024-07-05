# mnTS Wrapper ----------------------------------------------------------------
#' mnTS Wrapper
#'
#' This is a simple function that, by default, prints "Hello world". You can
#' customize the text to print (using the \code{to_print} argument) and add
#' an exclamation point (\code{excited = TRUE}).
#'
#' @param par something
#' @param par.fixed More things
#' @param Y inout matrix of
#' @param X input matrix of
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#' @examples
#' hello_world()
#' hello_world(excited = TRUE)
#' hello_world(to_print = "Hi world")
#'
#' @export

mnTS.ml.wrapper <- function(par, par.fixed, Y, X = NULL, Tsample) {
  retval <- mnTS_ml_cpp_list(par, par.fixed, Y, X, Tsample)
  return(retval$LL)
}


