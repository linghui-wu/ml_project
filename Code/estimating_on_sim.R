library(tidyverse)
library(here)
library(haven)
library(numbers)
library(Matrix)
library(parallel)
library(genlasso)

wd <- here()

source(str_c(wd, "/pglasso_modified.R"))

ind_data_df <- read_csv(str_c(wd, "/out/independent_1.csv"))
const_data_df <- read_csv(str_c(wd, "/out/const_nbhd_1.csv"))

adjacency <- read_csv(str_c(wd, "/out/tract_adjacency.csv"))

origin_tracts <- ind_data_df %>% 
  transmute(origins = k) %>% 
  distinct() %>%
  as.matrix(ncol = 1)

dest_tracts <- ind_data_df %>% 
  transmute(dest = n) %>% 
  distinct() %>%
  as.matrix(ncol = 1)

origin_adj <- adjacency %>%
  filter(source %in% origin_tracts,
         neighbor %in% origin_tracts,
         neighbor > source)

origin_list <- origin_tracts %>% as.matrix()

dest_adj <- adjacency %>%
  filter(source %in% dest_tracts,
         neighbor %in% dest_tracts,
         neighbor > source)

dest_list <- dest_tracts %>% as.matrix()

big_origin_adj_matrix <- matrix(rep(0, length(origin_list)), nrow = 1)
total_origin_size <- length(origin_adj$source)
for (s in 1:length(origin_adj$source)) {
  if (s %% 100 == 0) {
    print(s)
  }
  next_row <- rep(0, length(origin_list))
  next_row[which(origin_list == (origin_adj$source)[s])] <- 1
  next_row[which(origin_list == (origin_adj$neighbor)[s])] <- -1
  big_origin_adj_matrix <- rbind(big_origin_adj_matrix, next_row)
}

sparse_origin_a_mat <- big_origin_adj_matrix[-1,] %>% as("sparseMatrix")

big_dest_adj_matrix <- matrix(rep(0, length(dest_list)), nrow = 1)
total_dest_size <- length(dest_adj$source)
for (s in 1:length(dest_adj$source)) {
  if (s %% 100 == 0) {
    print(s)
  }
  next_row <- rep(0, length(dest_list))
  next_row[which(dest_list == (dest_adj$source)[s])] <- 1
  next_row[which(dest_list == (dest_adj$neighbor)[s])] <- -1
  big_dest_adj_matrix <- rbind(big_dest_adj_matrix, next_row)
}

sparse_dest_a_mat <- big_dest_adj_matrix[-1,] %>% as("sparseMatrix")

total_contrast_matrix <- bdiag(sparse_origin_a_mat,
                               sparse_dest_a_mat)

origin_adj_indices <- map(split(origin_adj %>% as.matrix() %>% t(),
                                rep(1:length(origin_adj$source), each = 2)),
                          \(x) {
                            matrix(c(which(origin_list == x[1]),
                                     which(origin_list == x[2]))
                                   , ncol = 2)
                          }) %>% reduce(rbind)

dest_adj_indices <- map(split(dest_adj %>% as.matrix() %>% t(),
                              rep(1:length(dest_adj$source), each = 2)),
                        \(x) {
                          matrix(c(which(dest_list == x[1]),
                                   which(dest_list == x[2]))
                                 , ncol = 2)
                        }) %>% reduce(rbind)

origin_graph_mat <- sparseMatrix(i = origin_adj_indices[,1],
                                 j = origin_adj_indices[,2],
                                 dims = c(max(origin_adj_indices),
                                          max(origin_adj_indices)))

dest_graph_mat <- sparseMatrix(i = dest_adj_indices[,1],
                               j = dest_adj_indices[,2],
                               dims = c(max(dest_adj_indices),
                                        max(dest_adj_indices)))

total_adjacency_graph <- graph_from_adjacency_matrix(bdiag(origin_graph_mat,
                                                           dest_graph_mat), 
                                                     mode = "undirected")

fe_mat_indices_list <- mclapply(
  (split((ind_data_df %>% select(k, n) %>% as.matrix() %>% 
            t()),
         rep(1:length(ind_data_df$k), each = 2))), 
  \(x) {
    c(which(origin_list == x[1]),
      which(dest_list == x[2]))
  }, mc.cores = 8) 

fe_mat_indices <- fe_mat_indices_list %>% unlist() %>%
  matrix(ncol = 2, byrow = T)

fe_selection_mat <- 
  sparseMatrix(i = 1:length(ind_data_df$k),
               j = fe_mat_indices[,1],
               dims = c(length(ind_data_df$k), 
                        length(origin_list) + 
                          length(dest_list))) -
  sparseMatrix(i = 1:length(ind_data_df$k),
               j = fe_mat_indices[,2] + length(origin_list))

l_1000 <- pglasso_mod(ind_data_df$ell_kn %>% matrix(ncol = 1),
                       log(ind_data_df$delta) %>% matrix(ncol = 1),
                       fe_mat_indices_list,
                       fe_mat_indices,
                       total_contrast_matrix,
                       sparse_origin_a_mat,
                       sparse_dest_a_mat,
                       length(origin_list),
                       length(dest_list),
                       gverbose = F, verbose = T,
                       conv.eps = 0.001, lambda = 1000,
                       mu.init = rep(0.5, 
                                     length(origin_list) + length(dest_list))
                       )

l_50 <- pglasso_mod(ind_data_df$ell_kn %>% matrix(ncol = 1),
                     log(ind_data_df$delta) %>% matrix(ncol = 1),
                     fe_mat_indices_list,
                     fe_mat_indices,
                     total_contrast_matrix,
                     sparse_origin_a_mat,
                     sparse_dest_a_mat,
                     length(origin_list),
                     length(dest_list),
                     gverbose = F, verbose = F,
                     conv.eps = 0.001, lambda = 50,
                     mu.init = rep(0, 
                                   length(origin_list) + length(dest_list))
)

l_10_const <- pglasso_mod(const_data_df$ell_kn %>% matrix(ncol = 1),
                    log(const_data_df$delta) %>% matrix(ncol = 1),
                    fe_mat_indices_list,
                    fe_mat_indices,
                    total_contrast_matrix,
                    sparse_origin_a_mat,
                    sparse_dest_a_mat,
                    length(origin_list),
                    length(dest_list),
                    gverbose = F, verbose = T,
                    conv.eps = 0.0001, lambda = 10,
                    mu.init = rep(0, 
                                  length(origin_list) + length(dest_list))
)

write_csv(tibble(estimates = c(l_10_const$coefficients[2], l_10_const$mu)),
          str_c(wd, "/out/fused_lasso_simulation_const_nbhd_10.csv"))

write_csv(tibble(estimates = c(l_10$coefficients[2], l_10$mu)),
          str_c(wd, "/out/fused_lasso_simulation_ind_10.csv"))
