library(R6)

#' @export
StanObj <- R6Class("StanObj",
                         public = list(
                           data = NULL,
                           stan_program = NULL,
                           initialize = function(data, stan_program) {
                             self$data = data
                             self$stan_program = stan_program
                           }
                           ))
