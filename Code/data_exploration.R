library(tidyverse)
library(here)
library(haven)
library(numbers)
library(Matrix)
library(parallel)
library(genlasso)

source(str_c(wd, "/pglasso_modified.R"))

wd <- here()

nyc_lodes <- read_dta(str_c(wd, 
                            "/data/nyc2010_lodes_wzero_wdelta.dta")) %>%
  mutate(i = as.numeric(i),
         j = as.numeric(j)) %>% 
  filter(!is.na(delta))

write_csv(nyc_lodes,
          str_c(wd, "/data/nyc2010_lodes.csv"))

adjacency <- read_csv(str_c(wd, "/data/nlist_2010.csv")) %>%
  filter(SOURCE_TRACTID <= nyc_lodes %>% 
           summarise(max = max(i, j)) %>% unlist(),
         SOURCE_TRACTID >= nyc_lodes %>% 
           summarise(min = min(i, j)) %>% unlist(),
         NEIGHBOR_TRACTID <= nyc_lodes %>% 
           summarise(max = max(i, j)) %>% unlist(),
         NEIGHBOR_TRACTID >= nyc_lodes %>% 
           summarise(min = min(i, j)) %>% unlist()
  ) %>%
  transmute(source = SOURCE_TRACTID %% 10^8,
            neighbor = NEIGHBOR_TRACTID %% 10^8)

origin_tracts <- nyc_lodes %>% 
  transmute(origins = i %% 10^8) %>% 
  distinct() %>%
  as.matrix(ncol = 1)

dest_tracts <- nyc_lodes %>% 
  transmute(dest = j %% 10^8) %>% 
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

sparse_origin_a_mat <- big_origin_adj_matrix %>% as("sparseMatrix")

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

sparse_dest_a_mat <- big_dest_adj_matrix %>% as("sparseMatrix")

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
  (split((nyc_lodes %>% select(i, j) %>% as.matrix() %>% 
            t()) - 3.6 * 10 ^ 10,
         rep(1:length(nyc_lodes$i), each = 2))), 
  \(x) {
    c(which(origin_list == x[1]),
      which(dest_list == x[2]))
  }, mc.cores = 8) 

fe_mat_indices <- fe_mat_indices_list %>% unlist() %>%
  matrix(ncol = 2, byrow = T)

fe_selection_mat <- 
  sparseMatrix(i = 1:length(nyc_lodes$i),
               j = fe_mat_indices[,1],
               dims = c(length(nyc_lodes$i), 
                        length(origin_list) + 
                          length(dest_list))) +
  sparseMatrix(i = 1:length(nyc_lodes$i),
               j = fe_mat_indices[,2] + length(origin_list))

total_contrast_matrix <- bdiag(sparse_origin_a_mat,
                               sparse_dest_a_mat)

l_1000 <- pglasso_mod(nyc_lodes$X_ij %>% matrix(ncol = 1),
                           nyc_lodes$log_delta %>% matrix(ncol = 1),
                           fe_mat_indices_list,
                           fe_mat_indices,
                           total_contrast_matrix,
                           sparse_origin_a_mat,
                           sparse_dest_a_mat,
                           length(origin_list),
                           length(dest_list),
                           gverbose = T, verbose = T,
                           conv.eps = 0.0001, lambda = 1000,
                      maxsteps = 10000, mu.init = c(rep(-1, length(origin_list)),
                                                    rep(1, length(dest_list)))
                      )

write_csv(tibble(i = origin_list %>% as.vector() + 36*10^9,
                 fe = l_100_unrestricted$mu[1:length(origin_list)]),
                 str_c(wd, "/out/full_100_i.csv"))

write_csv(tibble(j = dest_list %>% as.vector() + 36*10^9,
                 fe = l_100_unrestricted$mu[(length(origin_list)+1):
                                 (length(origin_list)+length(dest_list))]),
          str_c(wd, "/out/full_100_j.csv"))

l_1000 <- read_csv(str_c(wd, "/out/full_1000_i.csv")) %>% 
  select(fe) %>% as.matrix(ncol = 1) %>% 
  append(read_csv(str_c(wd, "/out/full_1000_j.csv")) %>% 
           select(fe) %>% as.matrix(ncol = 1))
