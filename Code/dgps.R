library(tidyverse)
library(parallel)
library(here)

wd <- here()

set.seed(31703)

n_samples <- 10

index_to_coords <- function(n_rows, n_cols, i) {
  c(i %% n_cols, i %/% n_cols + 1)
}

n_rows <- 20
n_cols <- 20

tracts <- 1:(n_rows * n_cols)

tracts_grid <- matrix(tracts, nrow = n_rows, byrow = T)

tract_adjacency <- map(tracts, 
                       \(x) {
                         mat <- matrix(c(0, 0), ncol = 2)
                         if (x %% n_cols != 1) {
                           mat <- rbind(mat, matrix(c(x, x - 1), ncol = 2))
                         }
                         if (x %% n_cols != 0) {
                           mat <- rbind(mat, matrix(c(x, x + 1), ncol = 2))
                         }
                         if (x  >= n_cols ) {
                           mat <- rbind(mat, matrix(c(x, x - n_cols), ncol = 2))
                         }
                         if (x %/% n_cols < n_rows - 1) {
                           mat <- rbind(mat, matrix(c(x, x + n_cols), ncol = 2))
                         }
                         mat[-1,]
                       }) %>% reduce(rbind)

write_csv(tibble(source = tract_adjacency[,1],
                 neighbor = tract_adjacency[,2]),
          str_c(wd, "/out/tract_adjacency.csv"))

#---- Generating true parameter values ----
dest_fes_l_const <- c(-1, -2, -4, -5, -6)
origin_fes_l_const <- c(0, -1, -2, -1.5)
epsilon <- 8
L <- 5e+5

origin_clusters <- map(tracts, \(t) {
  case_when(
    (t %% n_cols <= (n_rows %/% 4) & t %/% n_rows < (n_rows %/% 3)) ~ 1,
    (t %% n_cols > (n_rows %/% 4) & t %/% n_rows < (n_rows %/% 3)) ~ 2,
    (t %% n_cols <= (n_rows %/% 2) & t %/% n_rows >= (n_rows %/% 3)) ~ 3,
    (t %% n_cols > (n_rows %/% 2) & t %/% n_rows >= (n_rows %/% 3)) ~ 4
  )
}) %>% unlist()

dest_clusters <- map(tracts, \(t){
  case_when((t %% n_cols <= (n_rows %/% 5) & t %/% n_rows < (n_rows %/% 3)) ~ 1,
                      (tracts %% n_cols > (n_rows %/% 5) & 
                         tracts %% n_cols <= (n_rows %/% 2) &
                         tracts %/% n_rows < (n_rows %/% 3)) ~ 2,
                      (t %% n_cols > (n_rows %/% 2) & 
                         t %/% n_rows < (n_rows %/% 3)) ~ 3,
                      (t %% n_cols <= (n_rows %/% 2) &
                         t %/% n_rows >= (n_rows %/% 3)) ~ 4,
                      (t %% n_cols > (n_rows %/% 2) & 
                         t %/% n_rows >= (n_rows %/% 3)) ~ 5)
  }) %>% unlist()



locally_constant <- tibble(n = tracts,
                           origin_fe = 
                             origin_fes_l_const[origin_clusters[tracts]]) %>%
  full_join(tibble(k = tracts,
                   dest_fe = 
                     dest_fes_l_const[dest_clusters[tracts]]),
            by = character()) %>%
  select(n, k, origin_fe, dest_fe) %>%
  mutate(fe_sum = origin_fe + dest_fe,
         delta = map2(n, k, \(orr, dest) {
           0.05 * norm(which(tracts_grid == orr, arr.ind = T) -
              which(tracts_grid == dest, arr.ind = T), type = "o")
         }) %>% unlist() + 1,
         prob_numer = exp(fe_sum) / (delta ^ epsilon)) %>%
  mutate(prob = prob_numer / sum(prob_numer))

for (i in 1:n_samples) {
  locally_constant <- locally_constant %>% 
    mutate(ell_kn = rmultinom(1, L, prob) %>% as.integer())
  write_csv(locally_constant %>% select(n, k, ell_kn, delta, origin_fe, dest_fe),
            str_c(wd, "/out/const_nbhd_", i, ".csv")) 
}

#---- Independently generated primitives ----

independent <- tibble(n = tracts,
                      origin_fe = 
                        rnorm(length(tracts),
                              mean = mean(origin_fes_l_const),
                              sd = sd(origin_fes_l_const))) %>%
  full_join(tibble(k = tracts,
                   dest_fe = 
                     rnorm(length(tracts),
                           mean = mean(dest_fes_l_const),
                           sd = sd(dest_fes_l_const))),
            by = character()) %>%
  select(n, k, origin_fe, dest_fe) %>%
  mutate(fe_sum = origin_fe + dest_fe,
         delta = map2(n, k, \(orr, dest) {
           0.05 * norm(which(tracts_grid == orr, arr.ind = T) -
                         which(tracts_grid == dest, arr.ind = T), type = "o")
         }) %>% unlist() + 1,
         prob_numer = exp(fe_sum) / (delta ^ epsilon)) %>%
  mutate(prob = prob_numer / sum(prob_numer))

for (i in 1:n_samples) {
  independent <- independent %>% 
    mutate(ell_kn = rmultinom(1, L, prob) %>% as.integer())
  write_csv(independent %>% select(n, k, ell_kn, delta, origin_fe, dest_fe),
            str_c(wd, "/out/independent_", i, ".csv"))
}
