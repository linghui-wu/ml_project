library(tidyverse)
library(here)
library(patchwork)

wd <- here()
plots_out <- str_c(wd, "/out/plots")

ridge_const <- 
  read_csv(str_c(wd, "/out/fused_ridge_simulation_const_nbhd_optimal.csv"))
ridge_ind <- 
  read_csv(str_c(wd, "/out/fused_ridge_simulation_independent_optimal.csv"))

lasso_const <- 
  read_csv(str_c(wd, "/out/fused_lasso_simulation_const_nbhd_500.csv"))
lasso_ind <-
  read_csv(str_c(wd, "/out/fused_lasso_simulation_ind_10.csv"))

ind_data_df <- read_csv(str_c(wd, "/out/independent_1.csv"))
const_data_df <- read_csv(str_c(wd, "/out/const_nbhd_1.csv"))

fe_df <- full_join(tibble(x = 1:20),
                            tibble(y = 1:20),
                            by = character()) %>%
  mutate(ridge_const_o = ridge_const$estimates[2:401],
         ridge_const_d = ridge_const$estimates[402:801],
         ridge_ind_o = ridge_ind$estimates[2:401],
         ridge_ind_d = ridge_ind$estimates[402:801],
         lasso_const_o = lasso_const$estimates[402:801],
         lasso_const_d = lasso_const$estimates[2:401],
         lasso_ind_o = lasso_ind$estimates[402:801],
         lasso_ind_d = lasso_ind$estimates[2:401],
         true_ind_o = ind_data_df %>% group_by(n) %>% 
           summarise(value = mean(origin_fe)) %>% select(value) %>% unlist(),
         true_ind_d = ind_data_df %>% group_by(k) %>% 
           summarise(value = mean(dest_fe)) %>% select(value) %>% unlist(),
         true_const_o = const_data_df %>% group_by(n) %>% 
           summarise(value = mean(origin_fe)) %>% select(value) %>% unlist(),
         true_const_d = const_data_df %>% group_by(k) %>% 
           summarise(value = mean(dest_fe)) %>% select(value) %>% unlist())



#---- Plotting ----
png(str_c(plots_out, "/ridge_o_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_const_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +
  theme(legend.key.height = unit(100, "pt"),
        legend.key.width = unit(100, "pt"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/ridge_d_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_const_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/ridge_o_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/ridge_d_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_ind_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/lasso_o_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_const_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/lasso_d_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_const_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/lasso_o_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/lasso_d_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_ind_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/lasso_o_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/true_d_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_ind_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/true_o_ind.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Origin FEs")
dev.off()

png(str_c(plots_out, "/true_d_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_const_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +   theme(legend.key.height = unit(100, "pt"),         legend.key.width = unit(100, "pt"),         legend.title = element_text(size=20),         legend.text = element_text(size=20),         legend.position = "bottom") +   labs(fill = "Origin FEs") +
  labs(fill = "Dest. FEs")
dev.off()

png(str_c(plots_out, "/true_o_const.png"),
    width = 1200, height = 900)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_const_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053") +
  theme_void() +  
  theme(legend.key.height = unit(100, "pt"),         
        legend.key.width = unit(100, "pt"),         
        legend.title = element_text(size=20),         
        legend.text = element_text(size=20),         
        legend.position = "bottom") +   
  labs(fill = "Origin FEs")
dev.off()

min_const_o <- min(fe_df$lasso_const_o, fe_df$ridge_const_o, fe_df$true_const_o)
max_const_o <- max(fe_df$lasso_const_o, fe_df$ridge_const_o, fe_df$true_const_o)
min_const_d <- min(fe_df$lasso_const_d, fe_df$ridge_const_d, fe_df$true_const_d)
max_const_d <- max(fe_df$lasso_const_d, fe_df$ridge_const_d, fe_df$true_const_d)
min_ind_o <- min(fe_df$lasso_ind_o, fe_df$ridge_ind_o, fe_df$true_ind_o)
max_ind_o <- max(fe_df$lasso_ind_o, fe_df$ridge_ind_o, fe_df$true_ind_o)
min_ind_d <- min(fe_df$lasso_ind_d, fe_df$ridge_ind_d, fe_df$true_ind_d)
max_ind_d <- max(fe_df$lasso_ind_d, fe_df$ridge_ind_d, fe_df$true_ind_d)


png(str_c(plots_out, "/const_o_total.png"),
    width = 2700, height = 1000)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_const_o)) +
  scale_fill_gradient(low = "#E8FDFF", high = "#000053",
                      limits = c(min_const_o, max_const_o)) +
  theme_void() +  
  labs(title = "Fused Lasso Estimates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_const_o)) +
  scale_fill_gradient(low = "#E8FDFF", high = "#000053",
                      limits = c(min_const_o, max_const_o)) +
  theme_void() +  
  labs(title = "Fused Ridge Estimates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_const_o))  +
  scale_fill_gradient(low = "#E8FDFF", high = "#000053",
                      limits = c(min_const_o, max_const_o)) +
  theme_void() +
  labs(title = "True Values") +
  theme(plot.title = element_text(size = 60, hjust = 0.5),
        legend.position = "right",
        legend.key.width = unit(80, "pt"),
        legend.key.height = unit(100, "pt"),
        legend.title = element_text(size=50),         
        legend.text = element_text(size=50)) +   
  labs(fill = "Origin FEs") +
  plot_layout(guides = "collect")
dev.off()

png(str_c(plots_out, "/const_d_total.png"),
    width = 2700, height = 1000)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_const_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_const_d, max_const_d)) +
  theme_void() +  
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_const_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_const_d, max_const_d)) +
  theme_void() +  
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_const_d))  +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_const_d, max_const_d)) +
  theme_void() +
  theme(plot.title = element_text(size = 60, hjust = 0.5),
        legend.position = "right",
        legend.key.width = unit(80, "pt"),
        legend.key.height = unit(100, "pt"),
        legend.title = element_text(size=50),         
        legend.text = element_text(size=50)) +   
  labs(fill = "Dest. FEs") +
  plot_layout(guides = "collect")
dev.off()

png(str_c(plots_out, "/ind_d_total.png"),
    width = 2700, height = 1000)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_ind_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_d, max_ind_d)) +
  theme_void() +  
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_ind_d)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_d, max_ind_d)) +
  theme_void() +  
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_ind_d))  +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_d, max_ind_d)) +
  theme_void() +
  theme(plot.title = element_text(size = 60, hjust = 0.5),
        legend.position = "right",
        legend.key.width = unit(80, "pt"),
        legend.key.height = unit(100, "pt"),
        legend.title = element_text(size=50),         
        legend.text = element_text(size=50)) +   
  labs(fill = "Dest. FEs") +
  plot_layout(guides = "collect")
dev.off()

png(str_c(plots_out, "/ind_o_total.png"),
    width = 2700, height = 1000)
ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = lasso_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_o, max_ind_o)) +
  theme_void() +  
  labs(title = "Fused Lasso Estimates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = ridge_ind_o)) +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_o, max_ind_o)) +
  theme_void() +  
  labs(title = "Fused Ridge Estimates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 60, hjust = 0.5)) +
  ggplot(fe_df, mapping = aes(x, y)) + 
  geom_tile(aes(fill = true_ind_o))  +
  scale_fill_gradient(low = "#E7E7FF", high = "#000053",
                      limits = c(min_ind_o, max_ind_o)) +
  theme_void() +
  labs(title = "True Values") +
  theme(plot.title = element_text(size = 60, hjust = 0.5),
        legend.position = "right",
        legend.key.width = unit(80, "pt"),
        legend.key.height = unit(100, "pt"),
        legend.title = element_text(size=50),         
        legend.text = element_text(size=50)) +   
  labs(fill = "Origin FEs") +
  plot_layout(guides = "collect")
dev.off()

