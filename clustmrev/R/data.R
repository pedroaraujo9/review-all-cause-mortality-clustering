#' Human Mortality Database period life tables (1960--2010)
#'
#' A subset of period life table data from the Human Mortality Database (HMD),
#' covering 30 countries over the period 1960 to 2010 at ages 0 to 110 in
#' 5-year age groups (except for age 0 and age 1, which are given separately).
#' Countries with incomplete coverage over the full period range and
#' sub-national or special populations were excluded.
#'
#' @format A data frame with 36,720 rows and 14 variables:
#' \describe{
#'   \item{\code{year}}{Integer. Calendar year of the life table.}
#'   \item{\code{age}}{Integer. Age group (0, 1, 5, 10, \ldots, 110).}
#'   \item{\code{mx}}{Numeric. Central death rate \eqn{m_x}.}
#'   \item{\code{qx}}{Numeric. Probability of dying between age \eqn{x} and
#'     \eqn{x+n}, i.e. \eqn{q_x}.}
#'   \item{\code{ax}}{Numeric. Average number of years lived in the interval
#'     by those who die, \eqn{a_x}.}
#'   \item{\code{lx}}{Integer. Number of survivors at exact age \eqn{x} out
#'     of a radix of 100,000.}
#'   \item{\code{dx}}{Integer. Number of deaths between age \eqn{x} and
#'     \eqn{x+n}.}
#'   \item{\code{Lx}}{Integer. Person-years lived between age \eqn{x} and
#'     \eqn{x+n}.}
#'   \item{\code{Tx}}{Integer. Total person-years lived above age \eqn{x}.}
#'   \item{\code{ex}}{Numeric. Life expectancy at age \eqn{x}.}
#'   \item{\code{open_interval}}{Logical. Whether the age group is an open
#'     (terminal) interval.}
#'   \item{\code{country_code}}{Character. ISO 3166-1 alpha-3 country code.}
#'   \item{\code{country}}{Character. Country name.}
#'   \item{\code{data_extract}}{Date. Date on which the data were extracted
#'     from the HMD.}
#' }
#'
#' @details
#' The 30 countries included are: Australia, Austria, Belarus, Belgium,
#' Bulgaria, Canada, Czechia, Denmark, Estonia, Finland, France, Hungary,
#' Ireland, Italy, Japan, Latvia, Lithuania, Netherlands, New Zealand, Norway,
#' Poland, Portugal, Russia, Slovakia, Spain, Sweden, Switzerland, U.K.,
#' U.S.A., and Ukraine.
#'
#' The following populations were excluded: Hong Kong, Iceland, Chile,
#' Croatia, Republic of Korea, Luxembourg, East Germany, West Germany,
#' England and Wales (Total and Civilian), Scotland, Northern Ireland,
#' New Zealand Maori, and New Zealand Non-Maori.
#'
#' @source Human Mortality Database. Max Planck Institute for Demographic
#'   Research (Germany), University of California, Berkeley (USA), and French
#'   Institute for Demographic Studies (France). \url{https://www.mortality.org}.
#'   Data extracted on 2024-08-26.
"hmd_data"
