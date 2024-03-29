% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_s.t.ponctuel.tps.R
\name{pred_s.t.ponctuel.tps}
\alias{pred_s.t.ponctuel.tps}
\title{Predictions computation}
\usage{
pred_s.t.ponctuel.tps(newdata, object, s, window, event = 1)
}
\arguments{
\item{newdata}{data frame : collected data for a new individual}

\item{object}{lsjm object : estimation of the model}

\item{s}{numeric : the time to begin prediction}

\item{window}{numeric : the side of the prediction window}

\item{event}{integer (0, 1 or 2) : the event of interest for the prediction}
}
\description{
Predictions computation
}
