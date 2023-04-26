## code to prepare `Y_series` dataset goes here
# Y_series = read.csv('par0-100000-0th.csv', header = FALSE)$V1
# or
# S0 = C(0.125,0.1,0.25,0.1,-0.7)
# Y_series = crHeston(v_0=S0[3], n_segment=10, par=S0, N=1000, h=1)

usethis::use_data(Y_series, overwrite = TRUE)
