% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Delventhal_model.R
\name{plot_DT_fit}
\alias{plot_DT_fit}
\title{Plot the fitted values}
\usage{
plot_DT_fit(DT_fit)
}
\arguments{
\item{DT_fit}{The list returned from a call to \code{\link{predict_DT_model}}}
}
\value{
Returns a \code{ggplot}
}
\description{
Plot the fitted values
}
\examples{
T_max = 100
x_t = rep(1, T_max)
tau_2 = 33
tau_3 = 67
y = matrix(make_z1_scalar(1:T_max, tau_2, tau_3)*10 + rnorm(T_max, sd = 0.2),
          ncol = 1)
DT_fit = fit_DT_model(tau_2, tau_3, y, x_t)

plot_DT_fit(DT_fit)

}
