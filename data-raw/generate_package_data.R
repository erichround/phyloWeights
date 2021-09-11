# EXTERNAL DATA

yin_2020_data <- read.csv("data-raw/yin_2020_data.csv")
usethis::use_data(
  internal = FALSE, 
  overwrite = TRUE,
  yin_2020_data
)