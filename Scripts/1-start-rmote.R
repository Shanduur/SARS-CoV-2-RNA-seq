if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github(c("yihui/servr", "hafen/rmote"))

rmote::start_rmote()
