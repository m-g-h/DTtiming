% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Delventhal_model.R
\name{make_z1_scalar}
\alias{make_z1_scalar}
\title{Create \mjseqn{z_{1t}}}
\usage{
make_z1_scalar(t, tau_2, tau_3)
}
\arguments{
\item{t}{An \verb{integer vector} \mjseqn{1, \dots, T}}

\item{tau_2}{An \verb{integer scalar}}

\item{tau_3}{An \verb{integer scalar}}
}
\value{
Returns an integer vector
}
\description{
\loadmathjax
Creates the dummy variable \mjseqn{z_{1t} = d_{1t} + d_{2t} \frac{\tau_3 - t}{\tau_3 - \tau_2}}
}
