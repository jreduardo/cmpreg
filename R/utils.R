#-----------------------------------------------------------------------
# Verify if all elements of a atomic vector are equal
#' @param x an atomic vector.
#' @return logical, if \code{TRUE} all elements are equal.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom utils combn
#'
all_identical <- function(x) {
  all(apply(combn(x, 2), 2, function(.x) identical(.x[1], .x[2])))
}
