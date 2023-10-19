#' Create \mjseqn{d_{1t}}
#'
#' \loadmathjax
#' Creates the dummy variable \mjseqn{d_{1t} = 1\{t \leq \tau_2\}}
#'
#' @param t An `integer vector` \mjseqn{1, \dots, T}
#' @param tau_2 An `integer scalar`
#'
#' @return Returns an integer vector
#'
dummyfun_1 = function(t, tau_2){
  as.integer(t <= tau_2)
}


#' Create \mjseqn{d_{2t}}
#'
#' \loadmathjax
#' Creates the dummy variable \mjseqn{d_{2t} = 1\{\tau_2 < t < \tau_3\}}
#'
#' @param t An `integer vector` \mjseqn{1, \dots, T}
#' @param tau_2 An `integer scalar`
#' @param tau_3 An `integer scalar`
#'
#' @return Returns an integer vector
#'
dummyfun_2 = function(t, tau_2, tau_3){
  as.integer(t > tau_2 & t < tau_3)
}


#' Create \mjseqn{d_{3t}}
#'
#' \loadmathjax
#' Creates the dummy variable \mjseqn{d_{3t} = 1\{t \geq \tau_3\}}
#'
#' @param t An `integer vector` \mjseqn{1, \dots, T}
#' @param tau_3 An `integer scalar`
#'
#' @return Returns an integer vector
#'
dummyfun_3 = function(t, tau_3){
  as.integer(t >= tau_3)
}


#' Create \mjseqn{z_{1t}}
#'
#' \loadmathjax
#' Creates the dummy variable \mjseqn{z_{1t} = d_{1t} + d_{2t} \frac{\tau_3 - t}{\tau_3 - \tau_2}}
#'
#' @param t An `integer vector` \mjseqn{1, \dots, T}
#' @param tau_2 An `integer scalar`
#' @param tau_3 An `integer scalar`
#'
#' @return Returns an integer vector
#'
#' @export
#'
make_z1_scalar = function(t, tau_2, tau_3){
  dummyfun_1(t, tau_2) + dummyfun_2(t, tau_2, tau_3) * ((tau_3-t)/(tau_3-tau_2))
}

#' Create \mjseqn{z_{3t}}
#'
#' \loadmathjax
#' Creates the dummy variable \mjseqn{z_{3t} = d_{3t} + d_{2t} \frac{t - \tau_2}{\tau_3 - \tau_2}}
#'
#' @param t An `integer vector` \mjseqn{1, \dots, T}
#' @param tau_2 An `integer scalar`
#' @param tau_3 An `integer scalar`
#'
#' @return Returns an integer vector
#'
make_z3_scalar = function(t, tau_2, tau_3){
  dummyfun_3(t, tau_3) + dummyfun_2(t, tau_2, tau_3) * ((t - tau_2)/(tau_3-tau_2))
}


#' Fit the Delventhal model by OLS
#'
#' @param tau_2 tau_2 An `integer scalar`
#' @param tau_3 tau_2 An `integer scalar`
#' @param y `integer vector` The dependent variable, e.g. the crude birth or death rate
#' @param x_t `integer vector` The regressor \mjseqn{x_t}, in the simplest case a vector of ones
#'
#' @return Returns a list of parameters and model diagnostics:
#' * `alpha` Contains the parameters \mjseqn{\alpha_1} and\mjseqn{\alpha_3}
#' * `sigma` Contains the parameters \mjseqn{\sigma^2_1} and\mjseqn{\\sigma^2_3}
#' * `tau` Contains the parameters \mjseqn{\tau_1} and\mjseqn{\tau_3}
#' * `y` The original dependent variable
#' * `y_hat` The predicted values
#' * `RSS` The residual sum of squares
#' @export
#'
#' @examples
#' T_max = 100
#' x_t = rep(1, T_max)
#' tau_2 = 33
#' tau_3 = 67
#' y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
#'           ncol = 1)
#' DT_fit = fit_DT_model(tau_2, tau_3, y, x_t)

fit_DT_model = function(tau_2, tau_3, y, x_t){


  z1_row = make_z1_scalar(t = 1:length(x_t), tau_2, tau_3) * x_t
  z3_row = make_z3_scalar(t = 1:length(x_t), tau_2, tau_3) * x_t

  Z_prime = matrix(c(z1_row, z3_row),
                   nrow = 2, byrow = T)

  Z = t(Z_prime)

  alphas = solve(Z_prime %*% Z) %*% (Z_prime%*%y)

  e_col = matrix(nrow = nrow(Z))
  for (i in 1:nrow(Z)) {
    e_col[i,] =  y[i,] - (Z[i,] %*% alphas)
  }

  e1_col = e_col * z1_row
  e3_col = e_col * z3_row

  sigma_sq_1 = (1/sum(z1_row)) * (t(e1_col) %*% e1_col)
  sigma_sq_3 = (1/sum(z3_row)) * (t(e3_col) %*% e3_col)

  y_hat = predict_DT_model(alpha = as.vector(alphas),
                           tau = c(tau_2, tau_3),
                           as.vector(y))

  RSS = sum((y - y_hat)^2)

  list("alpha" = c("alpha_1" = alphas[1,1],
                   "alpha_3" = alphas[2,1]),
       "sigma_sq" = c("sigma_sq_1" = sigma_sq_1,
                      "sigma_sq_3" = sigma_sq_3),
       "tau" = c("tau_2" = tau_2,
                 "tau_3" = tau_3),
       "y" = as.vector(y),
       "y_hat" = y_hat,
       "RSS" = RSS)

}

#' Predict \mjseqn{y_t} based on estimated parameters
#'
#' @param alpha `numeric vector` of length 2: \mjseqn{(\alpha_1, \alpha_3)}
#' @param tau `numeric vector` of length 2: \mjseqn{(\tau_1, \tau_3)}
#' @param y `numeric vector` with the same length as the original \mjesqn{y_t}
#'
#' @return `numeric vector`
#' @export
#'
#' @examples
#' T_max = 100
#' x_t = rep(1, T_max)
#' tau_2 = 33
#' tau_3 = 67
#' y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
#'           ncol = 1)
#' DT_fit = fit_DT_model(tau_2, tau_3, y, x_t)
#'
#' predict_DT_model(DT_fit$alpha,
#'                  DT_fit$tau,
#'                  DT_fit$y)

predict_DT_model = function(alpha, tau, y){
  scaling = alpha[1] - alpha[2]
  intercept = alpha[2]

  make_z1_scalar(1:length(y),
                 tau[1],
                 tau[2]) *scaling + intercept
}

#' Plot the fitted values
#'
#' @param DT_fit The list returned from a call to \code{\link{predict_DT_model}}
#'
#' @return Returns a `ggplot`
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' T_max = 100
#' x_t = rep(1, T_max)
#' tau_2 = 33
#' tau_3 = 67
#' y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
#'           ncol = 1)
#' DT_fit = fit_DT_model(tau_2, tau_3, y, x_t)
#'
#' plot_DT_fit(DT_fit)
#'
plot_DT_fit = function(DT_fit){
  plotdata = data.frame(y = DT_fit$y,
                        y_hat = DT_fit$y_hat,
                        t = 1:length(DT_fit$y))

  ggplot2::ggplot(plotdata, ggplot2::aes(x = t)) +
    ggplot2::geom_point(ggplot2::aes(y = y)) +
    ggplot2::geom_line(ggplot2::aes(y = y_hat),
                       color = "red") +
    geom_label(x = DT_fit$tau[1]+0.2*length(DT_fit$y),
               y = DT_fit$alpha[1],
               label = paste0("alpha_1 = ", round(DT_fit$alpha[1],2),
                              "\n",
                              "tau_1 = ", round(DT_fit$tau[1])
                              )) +
    geom_label(x = DT_fit$tau[2]-0.2*length(DT_fit$y),
               y = DT_fit$alpha[2],
               label = paste0("alpha_3 = ", round(DT_fit$alpha[2],2),
                              "\n",
                              "tau_3 = ", round(DT_fit$tau[2])
               ))
}

#' Select the best fit by grid searchign over \mjseqn{(\tau_1, \tau_3)}
#'
#' @param y `integer vector` The dependent variable, e.g. the crude birth or death rate
#' @param tau_3_min `integer scalar` Lower bound for \mjseqn{tau_3}
#'
#' @return Returns the same list as \code{\link{fit_DT_model}}
#' @export
#'
#' @examples
#' T_max = 100
#' x_t = rep(1, T_max)
#' tau_2 = 33
#' tau_3 = 67
#' y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
#' ncol = 1)
#' DT_fit = grid_Search_best_DT_fit(y)


grid_Search_best_DT_fit = function(y, tau_3_min = 0){

  max_T = length(y)

  grid = expand.grid(tau_2 = 1:max_T,
                     tau_3 = 1:max_T)
  grid = grid[grid$tau_3 > grid$tau_2,]
  grid = grid[grid$tau_3 >= tau_3_min,]

  y = matrix(y, ncol = 1)

  x_t = rep(1, max_T)

  grid$DT_fit = mapply(fit_DT_model, grid$tau_2, grid$tau_3, MoreArgs = list(y, x_t),
                       SIMPLIFY = F)

  get_RSS = function(list){
    list$RSS
  }

  grid$RSS = mapply(get_RSS, grid$DT_fit)



  grid[grid$RSS == min(grid$RSS),]$DT_fit[[1]]

  # fits_CBR_GER = expand_grid(tau_2 = 1:max_T,
  #                            tau_3 = 1:max_T) %>%
  #   filter(tau_3 > tau_2) %>%
  #   filter(tau_3 >= tau_3_min) %>%
  #   mutate(DT_fit = map2(tau_2, tau_3,
  #                        fit_DT_model,
  #                        y = matrix(y, ncol = 1),
  #                        x_t = rep(1, max_T),
  #                        .progress = progress
  #   )
  #   )
  #
  # fits_CBR_GER %>%
  #   mutate(RSS = map_dbl(DT_fit,
  #                        "RSS")) %>%
  #   arrange(RSS) %>%
  #   first() %>%
  #   pull(DT_fit) %>%
  #   .[[1]]
}

# fit_DT_to_data_dem = function(data_dem){
#   result = data_dem %>%
#     select(Year:Rb, CBRInt, CDRInt) %>%
#     drop_na() %>%
#     pivot_longer(CBRInt:CDRInt) %>%
#     group_by(Sex, Code, Rb, name) %>%
#     summarise(Year = list(Year),
#               value = list(value),
#               max_T = n(),
#               .groups = "drop") %>%
#     mutate(DT_fit = map2(value, max_T,
#                          grid_Search_best_DT_fit,
#                          .progress = TRUE))
#
#   result %>%
#     mutate(alpha_1 = map_dbl(DT_fit, c("alpha", "alpha_1")),
#            alpha_3 = map_dbl(DT_fit, c("alpha", "alpha_3")),
#            tau_2 = map_int(DT_fit, c("tau", "tau_2")),
#            tau_3 = map_int(DT_fit, c("tau", "tau_3")),
#            DT_fit = map(DT_fit, "y_hat")) %>%
#     unnest(cols = c(Year, value, DT_fit)) %>%
#     pivot_wider(id_cols = c(Sex, Code, Rb, Year),
#                 names_from = "name",
#                 values_from = value:tau_3) %>%
#     rename(CBRInt = value_CBRInt,
#            CDRInt = value_CDRInt) %>%
#     left_join(data_dem, .,
#               by = join_by(Year, Code, Rb, Sex, CBRInt, CDRInt))
# }

# T_max = 100
# x_t = rep(1, T_max)
# tau_2 = 33
# tau_3 = 67
# y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
#            ncol = 1)
#
# DT_fit = fit_DT_model(tau_2, tau_3, y, x_t)
#
# plot_DT_fit(DT_fit)

# data = read_csv("_DEV/Dem_Trans_data.csv")
#
# data_GER = data %>%
#   filter(countryname == "Afghanistan") %>%
#   drop_na(cbirth) %>%
#   mutate() %>%
#   mutate(cbirthFitted2 = grid_Search_best_DT_fit(cbirth)$y_hat,
#          .after = cbirthFitted)%>%
#   mutate(cdeathFitted2 = grid_Search_best_DT_fit(cdeath)$y_hat,
#          .after = cdeathFitted)
#
# data_GER %>%
#   ggplot(aes(x = year)) +
#   geom_point(aes(y = cbirth)) +
#   geom_line(aes(y = cbirthFitted),
#             color = "blue") +
#   geom_line(aes(y = cbirthFitted2),
#             color = "red")
#
# data_GER %>%
#   ggplot(aes(x = year)) +
#   geom_point(aes(y = cdeath)) +
#   geom_line(aes(y = cdeathFitted),
#             color = "blue") +
#   geom_line(aes(y = cdeathFitted2),
#             color = "red")

