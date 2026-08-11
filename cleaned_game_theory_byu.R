# load libraries
library(tidyverse)  
library(MASS) # kde 
library(mgcv)  # gam
library(mvtnorm)  # sampling from multivariate normal
library(lpSolve)  # game theory
library(patchwork)  # for combining plots

# set seed
set.seed(1234)

# load data
pks <- read_csv("pk_dat.csv")


## building keeper save distributions to use later in solving the game
# filter to just saved shots
touched_shots <- pks |>
  filter(keeper_touch == 1) 

# creating density of saved shots
k <- kde2d(touched_shots$x, touched_shots$y, n = 200)
k_df <- as.data.frame(expand.grid(x = k$x, y = k$y))
k_df$density <- as.vector(k$z)

# create zones
touched_shots <- touched_shots |>
  dplyr::mutate(zone = case_when(
    x < -1.33 ~ "left",
    x >= -1.33 & x <= 1.33 ~ "center",
    x > 1.33 ~ "right"
  ))

# calculate zone stats
zone_stats <- touched_shots %>%
  group_by(zone) %>%
  summarise(
    mu_x = mean(x, na.rm = TRUE),
    mu_y = mean(y, na.rm = TRUE),
    sd_x = sd(x, na.rm = TRUE),
    sd_y = sd(y, na.rm = TRUE),
    rho  = cor(x, y, use = "complete.obs")
  )

# create list that has the zone stats in easier format
dives <- list(
  left   = with(filter(zone_stats, zone == "left"),   c(mu_x, mu_y)),
  center = with(filter(zone_stats, zone == "center"), c(mu_x, mu_y)),
  right  = with(filter(zone_stats, zone == "right"),  c(mu_x, mu_y))
)

# Define spreads (standard deviations) and correlation
sigma_x <- setNames(zone_stats$sd_x, zone_stats$zone)
sigma_y <- setNames(zone_stats$sd_y, zone_stats$zone)
rho     <- setNames(zone_stats$rho,  zone_stats$zone)

# define goal dimensions
goal_width <- 8
goal_height <- 2.67
goal_x_min <- -4
goal_x_max <- 4
goal_y_max <- goal_height
post_width <- 0.12

# Create a grid of (x, y) values over the goal area
x_vals <- seq(goal_x_min, goal_x_max, length.out = 100)
y_vals <- seq(0, goal_y_max, length.out = 100)
grid <- expand.grid(x = x_vals, y = y_vals)

# Compute density values for each dive direction
df_dive <- do.call(rbind, lapply(names(dives), function(d) {
  mu <- dives[[d]]
  cov_matrix <- matrix(c(sigma_x[d]^2,
                         rho[d] * sigma_x[d] * sigma_y[d],
                         rho[d] * sigma_x[d] * sigma_y[d],
                         sigma_y[d]^2),
                       nrow = 2)
  
  data.frame(
    x = grid$x,
    y = grid$y,
    dive = d,
    density = dmvnorm(cbind(grid$x, grid$y), mean = mu, sigma = cov_matrix)
  )
}))

## goal posts data for geom_polygon
# left post
left_post <- data.frame(
  x = c(-4, -4.11095, -4.11095, -4),
  y = c(0, 0, 2.748, 2.65),
  group = 1
)
# right post
right_post <- data.frame(
  x = c(4, 4.11095, 4.11095, 4),
  y = c(0, 0, 2.748, 2.65),
  group = 1
)
# crossbar
crossbar <- data.frame(
  x = c(-4, -4.11095, 4.11095, 4),
  y = c(2.65, 2.748, 2.748, 2.65)
)


# plot to check
ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_contour(data = df_dive, aes(x = x, y = y, z = density, 
                                   color = factor(dive, 
                                           levels = c("left", "center", "right"),
                                           labels = c("Left", "Center", "Right"))), 
               breaks = seq(0, max(df_dive$density), 
                            length.out = 15)) +  # Contour lines for each dive
  scale_color_manual(values = c("blue", "red", "green")) +  # Three distinct colors for each dive
  theme_minimal() +
  labs(title = "Goalkeeper Dive Coverage",
       x = "Goal Width (X-axis)",
       y = "Goal Height (Y-axis)",
       color = "Dive") +
  theme(aspect.ratio = goal_height / goal_width) 
# note that the distributions with this data is very asymmetric 
# due to small sample size

# create separate data frames for left and right dives
df_left <- df_dive %>% filter(dive == "left")
df_right <- df_dive %>% filter(dive == "right")

# calculate the probability mass
df_left <- df_left %>%
  mutate(prob_mass = density / sum(density)) 

df_right <- df_right %>%
  mutate(prob_mass = density / sum(density)) 

# calculate the cumulative probability
df_left <- df_left %>%
  arrange(desc(density)) %>%
  mutate(cum_prob = cumsum(prob_mass))

df_right <- df_right %>%
  arrange(desc(density)) %>%
  mutate(cum_prob = cumsum(prob_mass))

# create arbitrary quantiles
df_left <- df_left %>%
  mutate(quantile_band = cut(cum_prob,
                             breaks = seq(0, 1, by = 0.1),
                             labels = paste0("Q", 1:10),
                             include.lowest = TRUE))

df_right <- df_right %>%
  mutate(quantile_band = cut(cum_prob,
                             breaks = seq(0, 1, by = 0.1),
                             labels = paste0("Q", 1:10),
                             include.lowest = TRUE))

df_left <- df_left %>%
  mutate(
    save_prob = case_when(
      quantile_band %in% c("Q1", "Q2") ~ 0.90,
      quantile_band %in% c("Q3", "Q4") ~ 0.80,
      quantile_band %in% c("Q5", "Q6", "Q7") ~ 0.70,
      quantile_band == "Q8" ~ 0.60,
      quantile_band == "Q9" ~ 0.50,
      quantile_band == "Q10" ~ 0.00
    )
  )

df_right <- df_right %>%
  mutate(
    save_prob = case_when(
      quantile_band %in% c("Q1", "Q2") ~ 0.90,
      quantile_band %in% c("Q3", "Q4") ~ 0.80,
      quantile_band %in% c("Q5", "Q6", "Q7") ~ 0.70,
      quantile_band == "Q8" ~ 0.60,
      quantile_band == "Q9" ~ 0.50,
      quantile_band == "Q10" ~ 0.00
    )
  )

# create sequence for x and y to create a 40 x 32 grid (1280 grid rectangles)
x_seq <- seq(-4, 4, length.out = 41)
y_seq <- seq(2.67, 0, length.out = 33)

# Pre-compute grid centers
grid_centers <- list()
cell_id <- 1
for (j in 1:(length(x_seq) - 1)) {
  for (i in 1:(length(y_seq) - 1)) {
    x_center <- (x_seq[j] + x_seq[j+1]) / 2
    y_center <- (y_seq[i] + y_seq[i+1]) / 2
    grid_centers[[cell_id]] <- c(x_center, y_center)
    cell_id <- cell_id + 1
  }
}

## a couple of helper functions
# Check if shot is on frame
scoring_probability <- function(x, y) {
  x >= -4 && y <= 4 && y >= 0 && y <= goal_height
}

# Look up save probability from dive map
find_save_prob <- function(x, y, dive_df) {
  idx <- which.min((dive_df$x - x)^2 + (dive_df$y - y)^2)
  prob <- dive_df$save_prob[idx]
  if (is.na(prob)) return(0)
  return(prob)
}

# Simulate shots with execution error
simulate_made_percentage <- function(x_aim, y_aim, sigma_x, sigma_y, 
                                     n_samples, dive_df) {
  mean_vec <- c(x_aim, y_aim)
  cov_matrix <- matrix(c(sigma_x^2, 0, 0, sigma_y^2), nrow=2)
  samples <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = cov_matrix)
  samples[, 2] <- pmax(samples[, 2], 0)  # constrain z above ground
  
  outcomes <- sapply(1:nrow(samples), function(i) {
    x <- samples[i, 1]
    y <- samples[i, 2]
    
    if (!scoring_probability(x, y)) {
      return(0)
    } else {
      save_prob <- find_save_prob(x, y, dive_df)
      return(1 - save_prob)
    }
  })
  
  mean(outcomes) * 100
}

# create heatmap of the goal
generate_goal_heatmap <- function(dive_df, dive_label, sigma_x, sigma_y, n_samples) {
  made_percentages <- numeric(length(grid_centers))
  for (i in 1:length(grid_centers)) {
    aim_x <- grid_centers[[i]][1]
    aim_y <- grid_centers[[i]][2]
    
    made_percentages[i] <- simulate_made_percentage(aim_x, aim_y, sigma_x, 
                                                    sigma_y, n_samples, dive_df)
  }
  
  grid_data <- data.frame(
    x_center = sapply(grid_centers, function(c) c[1]),
    y_center = sapply(grid_centers, function(c) c[2]),
    made_percentage = made_percentages
  )
  
  ggplot(grid_data, aes(x = x_center, y = y_center, fill = made_percentage)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(
      option = "magma",
      name = "Goal %",
      breaks = c(0, 20, 40, 60, 80, 100),
      limits = c(0, 100)) +
    labs(
      title = paste("Goal Probability (", dive_label, " Dive)", sep = ""),
      x = "Y (Goal Width)",
      y = "Z (Height)"
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = goal_height / goal_width,
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 16)
    )
}

# apply the functions
p_left4 <- generate_goal_heatmap(df_left, "Left Side", sigma_x = .85, 
                                 sigma_y = .8, n_samples = 200)
p_right4 <- generate_goal_heatmap(df_right, "Right Side", sigma_x = .85, 
                                  sigma_y = .8, n_samples = 200)

p_left4
p_right4

# store the data from the heat maps
prob_left <- p_left4[["data"]]
prob_right <- p_right4[["data"]]

# combine the 2 data sets
combined2 <- prob_left %>%
  dplyr::select(x_center, y_center) %>%
  bind_cols(
    prob_left   %>% dplyr::select(left = made_percentage),
    prob_right  %>% dplyr::select(right  = made_percentage)
  ) %>% 
  dplyr::select(left, right)

########################### SOLVING THE GAME ###################################
## allowing any location on the goal to be the aiming location
df <- combined2
P <- as.matrix(df[, c("left","right")])
m <- nrow(P)
n <- ncol(P)

# 2 — Keeper LP
# Minimize v  subject to:  P %*% q <= v , sum(q) = 1, q >= 0

# Objective: q1 q2 q3 v  (we minimize v)
f.obj <- c(0, 0, 1)

# Inequality constraints: For each row i: sum_j P[i,j]*q_j - v <= 0
A_ub <- cbind(P, -rep(1, m))
dir_ub <- rep("<=", m)
b_ub <- rep(0, m)

# Equality: sum(q) = 1
A_eq <- matrix(c(1,1,0), nrow = 1)
dir_eq <- "="
b_eq <- 1

A <- rbind(A_ub, A_eq)
dir <- c(dir_ub, dir_eq)
rhs <- c(b_ub, b_eq)

# All variables have lower bound 0 (q1,q2,q3,v >= 0)
lower <- rep(0, 3)

sol <- lp("min", f.obj, A, dir, rhs)

keeper_mix <- sol$solution[1:2]
value <- sol$solution[3]

cat("\n===== Keeper optimal strategy =====\n")
print(keeper_mix)
cat("\n===== Game value (expected scoring %) =====\n")
print(value)

# 3 — Expected payoff for each row
expected <- P %*% keeper_mix

# 4 — Identify support rows (≈ rows where expected payoff = value)
tol <- 1e-8 + 1e-6 * abs(value)
support_idx <- which(abs(expected - value) <= tol)

if (length(support_idx) == 0) {
  maxval <- max(expected)
  support_idx <- which(expected >= maxval - 1e-8)
}

# 5 — Solve kicker distribution restricted to support
k <- length(support_idx)

if (k == 2){
  P_sub <- P[support_idx, ]  # 2x2 submatrix
  
  denom <- (P_sub[1,1] - P_sub[2,1] - P_sub[1,2] + P_sub[2,2])
  numer <- (P_sub[2,2] - P_sub[2,1])
  
  p1 <- numer / denom
  p1 <- max(0, min(1, p1))  # clamp to [0,1] just in case
  p2 <- 1 - p1
  
  kicker_support <- c(p1, p2)
} else if (k == 1) {
  kicker_support <- 1
} else {
  # Variables: p1..pk
  # Constraints: For each keeper column j:
  #     sum_i p_i * P[i,j] = value
  # Plus sum(p_i) = 1
  
  Aeq <- matrix(0, nrow = 3, ncol = k)
  for (j in 1:2)
    Aeq[j, ] <- P[support_idx, j]
  Aeq[3, ] <- 1
  
  beq <- c(value, value, value, 1)
  
  # Objective does not matter — minimize 0
  f.obj.p <- rep(0, k)
  lower.p <- rep(0, k)
  
  sol_p <- lp("min", f.obj.p, Aeq, rep("=",3), beq)
  
  if (sol_p$status == 0)
    kicker_support <- sol_p$solution
  else
    kicker_support <- rep(1/k, k)
}

# 6 — Map to full 1280 rows
kicker_mix <- rep(0, m)
kicker_mix[support_idx] <- kicker_support

# 7 — Attach to data frame and export
df$kicker_prob <- kicker_mix

support_df <- df[support_idx, ]
support_df$row <- support_idx

cat("\n===== Kicker support probabilities =====\n")
print(support_df[, c("row","left","right","kicker_prob")])


# Goalpost parameters
left_post_x <- -4
right_post_x <- 4
post_width <- 0.2
crossbar_height <- 2.67

# Optimal shot coordinates
optimal_x1 <- prob_left[support_idx[1], 1]
optimal_y1 <- prob_left[support_idx[1], 2]

optimal_x2 <- prob_left[support_idx[2], 1]
optimal_y2 <- prob_left[support_idx[2], 2]


# Probabilities and colors
prob_values <- support_df$kicker_prob[1:2]         
colors <- c("blue", "purple")
labels <- paste0(round(prob_values*100,1), "%")      

# Combine into a single data frame for plotting
optimal_df <- data.frame(
  x = c(optimal_x1, optimal_x2),
  y = c(optimal_y1, optimal_y2),
  label = labels,
  fill_color = colors
)

game_prob <- value

# Keeper probabilities and game value
keeper_probs <- c(keeper_mix[1], keeper_mix[2])
keeper_text <- paste0(
  "Keeper Randomization: Left: ", round(keeper_probs[1]*100,1), "% | ",
  "Right: ", round(keeper_probs[2]*100,1), "%\n ",
  "Game value: " , round(game_prob, 2), "%")


# Plot
ggplot() +
  # Goal frame
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  
  # Optimal points
  geom_point(data = optimal_df, aes(x = x, y = y, color = label),
             shape = 21, size = 5, stroke = 1.5, fill = optimal_df$fill_color) +
  
  # Legend mapping
  scale_color_manual(name = "Kicker Randomization",
                     values = setNames(colors, labels)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  
  # Labels with keeper probs in subtitle, game value in x-axis
  labs(
    title = "Optimal Shot Location",
    subtitle = keeper_text,
    x = "X",
    y = "Y"
  ) +
  
  theme_minimal() +
  theme(aspect.ratio = 3 / 10)


############### solving the game without center as an option ###################
## repeat but removing the center since we don't have keeper coverage there
# filter the df_lefts and right
df_left2 <- df_left |>
  filter(x < -1.33 | x > 1.33)

df_right2 <- df_right |>
  filter(x < -1.33 | x > 1.33)

p_left5 <- generate_goal_heatmap(df_left2, "Left Side", sigma_x = .85, 
                                 sigma_y = .8, n_samples = 200)
p_right5 <- generate_goal_heatmap(df_right2, "Right Side", sigma_x = .85, 
                                  sigma_y = .8, n_samples = 200)

p_left5
p_right5


prob_left <- p_left5[["data"]]
prob_right <- p_right5[["data"]]

prob_left2 <- prob_left |>
  filter(x_center < -1.33 | x_center > 1.33)
prob_right2 <- prob_right |>
  filter(x_center < -1.33 | x_center > 1.33)

combined2.1 <- prob_left2 %>%
  dplyr::select(x_center, y_center) %>%
  bind_cols(
    prob_left2   %>% dplyr::select(left = made_percentage),
    prob_right2  %>% dplyr::select(right = made_percentage)
  ) %>% 
  dplyr::select(left, right)


df <- combined2.1
P <- as.matrix(df[, c("left","right")])
m <- nrow(P)
n <- ncol(P)

# 2 — Keeper LP
# Minimize v  subject to:  P %*% q <= v , sum(q) = 1, q >= 0

# Objective: q1 q2 q3 v  (we minimize v)
f.obj <- c(0, 0, 1)

# Inequality constraints: For each row i: sum_j P[i,j]*q_j - v <= 0
A_ub <- cbind(P, -rep(1, m))
dir_ub <- rep("<=", m)
b_ub <- rep(0, m)

# Equality: sum(q) = 1
A_eq <- matrix(c(1,1,0), nrow = 1)
dir_eq <- "="
b_eq <- 1

A <- rbind(A_ub, A_eq)
dir <- c(dir_ub, dir_eq)
rhs <- c(b_ub, b_eq)

# All variables have lower bound 0 (q1,q2,q3,v >= 0)
lower <- rep(0, 3)

sol <- lp("min", f.obj, A, dir, rhs)

keeper_mix <- sol$solution[1:2]
value <- sol$solution[3]

cat("\n===== Keeper optimal strategy =====\n")
print(keeper_mix)
cat("\n===== Game value (expected scoring %) =====\n")
print(value)

# 3 — Expected payoff for each row
expected <- P %*% keeper_mix

# 4 — Identify support rows (≈ rows where expected payoff = value)
tol <- 1e-8 + 1e-6 * abs(value)
support_idx <- which(abs(expected - value) <= tol)

if (length(support_idx) == 0) {
  maxval <- max(expected)
  support_idx <- which(expected >= maxval - 1e-8)
}

# 5 — Solve kicker distribution analytically (2x2 game)
# Make keeper indifferent between columns (left vs right)
# p * P[1,1] + (1-p) * P[2,1] = p * P[1,2] + (1-p) * P[2,2]
# Rearranging: p * (P[1,1] - P[2,1] - P[1,2] + P[2,2]) = P[2,2] - P[2,1]

k <- length(support_idx)

if (k == 2) {
  P_sub <- P[support_idx, ]  # 2x2 submatrix
  
  denom <- (P_sub[1,1] - P_sub[2,1] - P_sub[1,2] + P_sub[2,2])
  numer <- (P_sub[2,2] - P_sub[2,1])
  
  p1 <- numer / denom
  p1 <- max(0, min(1, p1))  # clamp to [0,1] just in case
  p2 <- 1 - p1
  
  kicker_support <- c(p1, p2)
  
} else if (k == 1) {
  kicker_support <- 1
} else {
  # Fallback to LP if support is somehow larger than 2
  Aeq <- matrix(0, nrow = 3, ncol = k)
  for (j in 1:2)
    Aeq[j, ] <- P[support_idx, j]
  Aeq[3, ] <- 1
  beq <- c(value, value, 1)
  sol_p <- lp("min", rep(0, k), Aeq, rep("=", 3), beq)
  kicker_support <- if (sol_p$status == 0) sol_p$solution else rep(1/k, k)
}

# 6 — Map to full 1280 rows
kicker_mix <- rep(0, m)
kicker_mix[support_idx] <- kicker_support

# 7 — Attach to data frame and export
df$kicker_prob <- kicker_mix

support_df <- df[support_idx, ]
support_df$row <- support_idx

cat("\n===== Kicker support probabilities =====\n")
print(support_df[, c("row","left","right","kicker_prob")])


# Goalpost parameters
left_post_x <- -4
right_post_x <- 4
post_width <- 0.2
crossbar_height <- 2.67

# Optimal shot coordinates
optimal_x1 <- prob_left2[support_idx[1], 1]
optimal_y1 <- prob_left2[support_idx[1], 2]

optimal_x2 <- prob_left2[support_idx[2], 1]
optimal_y2 <- prob_left2[support_idx[2], 2]


# Probabilities and colors
prob_values <- support_df$kicker_prob[1:2]         
colors <- c("#9fc8c8", "#298c8c")
labels <- paste0(round(prob_values*100,1), "%")      

# Combine into a single data frame for plotting
optimal_df <- data.frame(
  x = c(optimal_x1, optimal_x2),
  y = c(optimal_y1, optimal_y2),
  label = labels,
  fill_color = colors
)

game_prob <- value

# Keeper probabilities and game value
keeper_probs <- c(keeper_mix[1], keeper_mix[2])
keeper_text <- paste0(
  "Keeper Randomization: Left: ", round(keeper_probs[1]*100,1), "% | ",
  "Right: ", round(keeper_probs[2]*100,1), "%\n ",
  "Game value: " , round(game_prob, 2), "%")

# Plot
ggplot() +
  # Goal frame
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  # Optimal points
  geom_point(data = optimal_df, aes(x = x, y = y, color = label),
             shape = 21, size = 5, stroke = 1.5, fill = optimal_df$fill_color) +
  
  # Legend mapping
  scale_color_manual(name = "Kicker Randomization",
                     values = setNames(colors, labels)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  
  # Labels with keeper probs in subtitle, game value in x-axis
  labs(
    title = "Optimal Shot Location",
    subtitle = keeper_text,
    x = "X",
    y = "Y"
  ) +
  
  theme_minimal() +
  theme(aspect.ratio = 3 / 10)


############## solving but just allowing for the 4 target locations ############
# step 1: get shooter stats aka error distributions
target_stats <- pks |>
  group_by(target) |>
  summarize(n = n(),
            mean_x = mean(error_x),
            mean_y = mean(error_y),
            sd_x = sd(error_x),
            sd_y = sd(error_y))


# step 2: target specific made_percentage function
simulate_made_percentage_empirical <- function(target_name, dive_df, 
                                               shooter_stats, n_samples) {
  stats <- shooter_stats |> filter(target == target_name)
  
  if (nrow(stats) == 0) {
    warning(paste("No data for target:", target_name))
    return(NA)
  }
  
  # Get the target coordinates
  target_coords <- targets |> filter(target_name == !!target_name)
  target_x <- target_coords$x
  target_y <- target_coords$y
  
  # Sample from error distribution (centered at 0 since these are errors)
  mean_vec   <- c(stats$mean_x, stats$mean_y)
  cov_matrix <- matrix(c(stats$sd_x^2, 0, 
                         0, stats$sd_y^2), nrow = 2)
  
  error_samples <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = cov_matrix)
  
  # Translate back to goal coordinates
  actual_x <- error_samples[, 1] + target_x
  actual_y <- error_samples[, 2] + target_y
  actual_y <- pmax(actual_y, 0)  # keep above ground
  
  outcomes <- sapply(1:n_samples, function(i) {
    x <- actual_x[i]
    y <- actual_y[i]
    
    if (!scoring_probability(x, y)) {
      return(0)
    } else {
      save_prob <- find_save_prob(x, y, dive_df)
      return(1 - save_prob)
    }
  })
  
  mean(outcomes) * 100
}

# step 3: build targets data frame
targets <- pks |>
  group_by(target) |>
  summarize(
    x = mean(target_x),
    y = mean(target_y)
  ) |>
  rename(target_name = target)

# step 4: build the payoff matrix P using empirical distributions
# rows = 4 targets, columns = keeper dives left or right
# Using df_left2 and df_right2 from your collegiate keeper data

set.seed(1234)

P <- matrix(NA, nrow = nrow(targets), ncol = 2,
            dimnames = list(targets$target_name, 
                            c("keeper_left", "keeper_right")))

for (i in 1:nrow(targets)) {
  P[i, 1] <- simulate_made_percentage_empirical(
    target_name  = targets$target_name[i],
    dive_df      = df_left2,
    shooter_stats = target_stats,
    n_samples    = 10000
  )
  P[i, 2] <- simulate_made_percentage_empirical(
    target_name  = targets$target_name[i],
    dive_df      = df_right2,
    shooter_stats = target_stats,
    n_samples    = 10000
  )
}

print(P)

# step 5: build df and run LP solver

df <- data.frame(
  x_center = targets$x,
  y_center = targets$y,
  target   = targets$target_name,
  left     = P[, "keeper_left"],
  right    = P[, "keeper_right"]
)

m <- nrow(P)

f.obj <- c(0, 0, 1)

A_ub <- cbind(P, -rep(1, m))
dir_ub <- rep("<=", m)
b_ub <- rep(0, m)

A_eq <- matrix(c(1, 1, 0), nrow = 1)
dir_eq <- "="
b_eq <- 1

A   <- rbind(A_ub, A_eq)
dir <- c(dir_ub, dir_eq)
rhs <- c(b_ub, b_eq)

sol <- lp("min", f.obj, A, dir, rhs)

keeper_mix <- sol$solution[1:2]
value      <- sol$solution[3]

cat("\n===== Keeper optimal strategy =====\n")
print(keeper_mix)
cat("\n===== Game value =====\n")
print(value)

# step 6: kicker distribution 

expected    <- P %*% keeper_mix
tol         <- 1e-8 + 1e-6 * abs(value)
support_idx <- which(abs(expected - value) <= tol)

if (length(support_idx) == 0) {
  maxval      <- max(expected)
  support_idx <- which(expected >= maxval - 1e-8)
}

k <- length(support_idx)

if (k == 2) {
  P_sub <- P[support_idx, ]
  denom <- (P_sub[1,1] - P_sub[2,1] - P_sub[1,2] + P_sub[2,2])
  numer <- (P_sub[2,2] - P_sub[2,1])
  p1    <- max(0, min(1, numer / denom))
  p2    <- 1 - p1
  kicker_support <- c(p1, p2)
} else if (k == 1) {
  kicker_support <- 1
} else {
  Aeq <- matrix(0, nrow = 3, ncol = k)
  for (j in 1:2)
    Aeq[j, ] <- P[support_idx, j]
  Aeq[3, ] <- 1
  beq    <- c(value, value, 1)
  sol_p  <- lp("min", rep(0, k), Aeq, rep("=", 3), beq)
  kicker_support <- if (sol_p$status == 0) sol_p$solution else rep(1/k, k)
}

kicker_mix <- rep(0, m)
kicker_mix[support_idx] <- kicker_support

df$kicker_prob <- kicker_mix
support_df     <- df[support_idx, ]
support_df$row <- support_idx

cat("\n===== Kicker support probabilities =====\n")
print(support_df[, c("row", "left", "right", "kicker_prob")])

# step 7: plot
optimal_x1 <- targets$x[support_idx[1]]
optimal_y1 <- targets$y[support_idx[1]]

optimal_x2 <- targets$x[support_idx[2]]
optimal_y2 <- targets$y[support_idx[2]]


# Probabilities and colors
prob_values <- support_df$kicker_prob[1:2]         
colors <- c("#9fc8c8", "#298c8c")
labels <- paste0(round(prob_values*100,1), "%")      

# Combine into a single data frame for plotting
optimal_df <- data.frame(
  x = c(optimal_x1, optimal_x2),
  y = c(optimal_y1, optimal_y2),
  label = labels,
  fill_color = colors
)

game_prob <- value

# Keeper probabilities and game value
keeper_probs <- c(keeper_mix[1], keeper_mix[2])
keeper_text <- paste0(
  "Keeper Randomization: Left: ", round(keeper_probs[1]*100,1), "% | ",
  "Right: ", round(keeper_probs[2]*100,1), "%\n ",
  "Game value: " , round(game_prob, 2), "%")


ggplot() +
  # Goal frame
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  # Optimal points
  geom_point(data = optimal_df, aes(x = x, y = y, color = label),
             shape = 21, size = 5, stroke = 1.5, fill = optimal_df$fill_color) +
  
  # Legend mapping
  scale_color_manual(name = "Kicker Randomization",
                     values = setNames(colors, labels)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  
  # Labels with keeper probs in subtitle, game value in x-axis
  labs(
    title = "Optimal Shot Location",
    subtitle = keeper_text,
    x = "X",
    y = "Y"
  ) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 0.89, 1.78, 2.67)) +
  theme_minimal() +
  theme(aspect.ratio = 3 / 10)


################# solving the game for individual players ######################
# step 1: get player stats by targets - error distributions
shooter_stats_individual <- pks |>
  group_by(player_id, target) |>
  summarize(n = n(),
            mean_x = mean(error_x),
            mean_y = mean(error_y),
            sd_x = sd(error_x),
            sd_y = sd(error_y),
            .groups = "drop")
print(shooter_stats_individual)

# step 2: modify the simulation function to work with player id
simulate_made_percentage_player <- function(target_name, player, dive_df,
                                            shooter_stats_ind, targets_df,
                                            n_samples) {
  # Get this player's distribution for this target
  stats <- shooter_stats_ind |> 
    filter(player_id == player, target == target_name)
  
  # If player has no data at this target, fall back to population mean
  if (nrow(stats) == 0) {
    warning(paste("No data for player", player, "at target", target_name, 
                  "- using population mean"))
    stats <- pks |>
      filter(target == target_name) |>
      summarize(
        mean_x = mean(error_x, na.rm = TRUE),
        mean_y = mean(error_y, na.rm = TRUE),
        sd_x   = sd(error_x, na.rm = TRUE),
        sd_y   = sd(error_y, na.rm = TRUE)
      )
  }
  
  # Get target coordinates
  target_coords <- targets_df |> filter(target_name == !!target_name)
  target_x <- target_coords$x
  target_y <- target_coords$y
  
  mean_vec   <- c(stats$mean_x, stats$mean_y)
  cov_matrix <- matrix(c(stats$sd_x^2, 0,
                         0, stats$sd_y^2), nrow = 2)
  
  error_samples <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = cov_matrix)
  
  actual_x <- error_samples[, 1] + target_x
  actual_y <- pmax(error_samples[, 2] + target_y, 0)
  
  outcomes <- sapply(1:n_samples, function(i) {
    x <- actual_x[i]
    y <- actual_y[i]
    if (!scoring_probability(x, y)) {
      return(0)
    } else {
      save_prob <- find_save_prob(x, y, dive_df)
      return(1 - save_prob)
    }
  })
  
  mean(outcomes) * 100
}

# step 3: loop over each player and solve the game for them
set.seed(1234)

player_ids <- unique(pks$player_id)
player_results <- list()

for (pid in player_ids) {
  cat("\n========== Player", pid, "==========\n")
  
  # Build payoff matrix for this player
  P_player <- matrix(NA, nrow = nrow(targets), ncol = 2,
                     dimnames = list(targets$target_name,
                                     c("keeper_left", "keeper_right")))
  
  for (i in 1:nrow(targets)) {
    P_player[i, 1] <- simulate_made_percentage_player(
      target_name      = targets$target_name[i],
      player           = pid,
      dive_df          = df_left2,
      shooter_stats_ind = shooter_stats_individual,
      targets_df       = targets,
      n_samples        = 10000
    )
    P_player[i, 2] <- simulate_made_percentage_player(
      target_name      = targets$target_name[i],
      player           = pid,
      dive_df          = df_right2,
      shooter_stats_ind = shooter_stats_individual,
      targets_df       = targets,
      n_samples        = 10000
    )
  }
  
  cat("Payoff matrix:\n")
  print(P_player)
  
  # LP solver — identical to your existing code
  m     <- nrow(P_player)
  f.obj <- c(0, 0, 1)
  A_ub  <- cbind(P_player, -rep(1, m))
  A_eq  <- matrix(c(1, 1, 0), nrow = 1)
  A     <- rbind(A_ub, A_eq)
  dir   <- c(rep("<=", m), "=")
  rhs   <- c(rep(0, m), 1)
  
  sol        <- lp("min", f.obj, A, dir, rhs)
  keeper_mix <- sol$solution[1:2]
  value      <- sol$solution[3]
  
  cat("Keeper optimal strategy:\n")
  print(keeper_mix)
  cat("Game value:\n")
  print(value)
  
  # Kicker distribution
  expected    <- P_player %*% keeper_mix
  tol         <- 1e-8 + 1e-6 * abs(value)
  support_idx <- which(abs(expected - value) <= tol)
  
  if (length(support_idx) == 0) {
    maxval      <- max(expected)
    support_idx <- which(expected >= maxval - 1e-8)
  }
  
  k <- length(support_idx)
  
  if (k == 2) {
    P_sub <- P_player[support_idx, ]
    denom <- (P_sub[1,1] - P_sub[2,1] - P_sub[1,2] + P_sub[2,2])
    numer <- (P_sub[2,2] - P_sub[2,1])
    p1    <- max(0, min(1, numer / denom))
    kicker_support <- c(p1, 1 - p1)
  } else if (k == 1) {
    kicker_support <- 1
  } else {
    Aeq   <- matrix(0, nrow = 3, ncol = k)
    for (j in 1:2) Aeq[j, ] <- P_player[support_idx, j]
    Aeq[3, ] <- 1
    beq   <- c(value, value, 1)
    sol_p <- lp("min", rep(0, k), Aeq, rep("=", 3), beq)
    kicker_support <- if (sol_p$status == 0) sol_p$solution else rep(1/k, k)
  }
  
  kicker_mix <- rep(0, m)
  kicker_mix[support_idx] <- kicker_support
  
  # Store results for this player
  player_results[[as.character(pid)]] <- list(
    player_id   = pid,
    P           = P_player,
    keeper_mix  = keeper_mix,
    value       = value,
    support_idx = support_idx,
    kicker_mix  = kicker_mix,
    targets     = targets
  )
  
  cat("Kicker support:\n")
  print(data.frame(
    target      = targets$target_name[support_idx],
    left        = P_player[support_idx, 1],
    right       = P_player[support_idx, 2],
    kicker_prob = kicker_support
  ))
}

agg_strategy <- c(0.313, 0, 0.687, 0)  # LH, LL, RH, RL

# For each player, compute keeper's best response to aggregate strategy
for (pid in player_ids) {
  P_player <- player_results[[as.character(pid)]]$P
  #cat("Player", pid, "P matrix", P_player, "\n")
  
  # Expected scoring for each keeper action under aggregate kicker strategy
  expected_per_keeper <- t(P_player) %*% agg_strategy
  #cat("Player", pid, "expected per keeper", expected_per_keeper, "\n")
  
  # Keeper best responds by minimizing - picks lowest expected scoring
  value_agg <- min(expected_per_keeper)
  
  cat("Player", pid, "value under aggregate strategy:", value_agg, "\n")
}

### Plot of goalkeeper coverage
df_dive |>
  filter(dive != "center") |>
  ggplot(aes(x = x, y = y, color = dive)) +
  geom_contour(aes(z = density), breaks = seq(0, max(df_dive$density), 
                                              length.out = 15)) +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_color_manual(values = c("left" = "#FF7F50",
                                "right" = "steelblue"),
                     labels = c("left" = "Left",
                                "right" = "Right"),
                     name = "Keeper Dive") +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))


############################### Table Data #####################################
## speed data
pks |>
  group_by(player_id, target) |>
  summarize(n = n(),
            mean_speed = mean(speed_manual),
            sd_speed = sd(speed_manual),
            .groups = "drop")

pks |> 
  group_by(player_id) |>
  summarize(n = n(),
            mean_speed = mean(speed_manual),
            sd_speed = sd(speed_manual))

pks |>
  group_by(target) |>
  summarize(n = n(),
            mean_speed = mean(speed_manual),
            sd_speed = sd(speed_manual))
# overall mean speed and sd
mean(pks$speed_manual)
sd(pks$speed_manual)


## how close are people to the targets (error)
# each players for each target
player_target_summary <- pks |>
  group_by(target, player_id) |>
  summarize(n = n(),
            mean_x = mean(error_x),
            sd_x = sd(error_x),
            mean_y = mean(error_y),
            sd_y = sd(error_y))

# each player overall
pks |>
  group_by(player_id) |>
  summarize(n = n(),
            mean_x = mean(error_x),
            sd_x = sd(error_x),
            mean_y = mean(error_y),
            sd_y = sd(error_y))

# each target overall
overall_target_summary <- pks |>
  group_by(target) |>
  summarize(n = n(),
            mean_x = mean(error_x),
            sd_x = sd(error_x),
            mean_y = mean(error_y),
            sd_y = sd(error_y))


## shooter summary statistics
# shooter 1
pks |>
  filter(player_id == 1) |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 1,
         result == "G") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 1,
         result == "S") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 1,
         result == "F") |>
  group_by(target) |>
  summarize(count = n())

# shooter 2
pks |>
  filter(player_id == 2) |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 2,
         result == "G") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 2,
         result == "S") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 2,
         result == "F") |>
  group_by(target) |>
  summarize(count = n())

# shooter 3
pks |>
  filter(player_id == 3) |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 3,
         result == "G") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 3,
         result == "S") |>
  group_by(target) |>
  summarize(count = n())
pks |>
  filter(player_id == 3,
         result == "F") |>
  group_by(target) |>
  summarize(count = n())


############################### Data Plots #####################################
make_ellipse <- function(mean_x, mean_y, sd_x, sd_y, n_points = 100){
  #scale <- sqrt(qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n_points)
  data.frame(
    x = mean_x + sd_x * cos(theta),
    y = mean_y + sd_y * sin(theta)
  )
}

# player ellipses
player_ellipses <- player_target_summary |>
  rowwise() |>
  do({
    ell <- make_ellipse(.$mean_x, .$mean_y, .$sd_x, .$sd_y)
    ell$target <- .$target
    ell$player_id <- .$player_id
    ell
  })
player_ellipses <- player_target_summary |>
  group_by(player_id, target) |>
  reframe(
    make_ellipse(mean_x, mean_y, sd_x, sd_y)
  )

# overall ellipses
overall_ellipses <- overall_target_summary |>
  rowwise() |>
  do({
    ell <- make_ellipse(.$mean_x, .$mean_y, .$sd_x, .$sd_y)
    ell$target <- .$target
    ell
  })
overall_ellipses <- overall_target_summary |>
  group_by(target) |>
  reframe(
    make_ellipse(mean_x, mean_y, sd_x, sd_y)
  )

# add variable that has the ellipse on the same coordinate system as the goal
player_ellipses <- player_ellipses |>
  mutate(x_coord = case_when(
    target == "LH" ~ (x - 3.3755),
    target == "LL" ~ (x - 3.443),
    target == "RH" ~ (x + 3.388),
    target == "RL" ~ (x + 3.435)
  ),
  y_coord = case_when(
    target == "LH" ~ (y + 2.035),
    target == "LL" ~ (y + 0.617),
    target == "RH" ~ (y + 2.103),
    target == "RL" ~ (y + 0.611)
  ),
  target = factor(target, levels = c("LH", "RH", "LL", "RL"))
  )

overall_ellipses <- overall_ellipses |>
  mutate(x_coord = case_when(
    target == "LH" ~ (x - 3.3755),
    target == "LL" ~ (x - 3.443),
    target == "RH" ~ (x + 3.388),
    target == "RL" ~ (x + 3.435)
  ),
  y_coord = case_when(
    target == "LH" ~ (y + 2.035),
    target == "LL" ~ (y + 0.617),
    target == "RH" ~ (y + 2.103),
    target == "RL" ~ (y + 0.611)
  ))

player_target_summary <- player_target_summary |>
  mutate(x_coord = case_when(
    target == "LH" ~ (mean_x - 3.3755),
    target == "LL" ~ (mean_x - 3.443),
    target == "RH" ~ (mean_x + 3.388),
    target == "RL" ~ (mean_x + 3.435)
  ),
  y_coord = case_when(
    target == "LH" ~ (mean_y + 2.035),
    target == "LL" ~ (mean_y + 0.617),
    target == "RH" ~ (mean_y + 2.103),
    target == "RL" ~ (mean_y + 0.611)
  ),
  target = factor(target, levels = c("LH", "RH", "LL", "RL"))
  ) 

overall_target_summary <- overall_target_summary |>
  mutate(x_coord = case_when(
    target == "LH" ~ (mean_x - 3.3755),
    target == "LL" ~ (mean_x - 3.443),
    target == "RH" ~ (mean_x + 3.388),
    target == "RL" ~ (mean_x + 3.435)
  ),
  y_coord = case_when(
    target == "LH" ~ (mean_y + 2.035),
    target == "LL" ~ (mean_y + 0.617),
    target == "RH" ~ (mean_y + 2.103),
    target == "RL" ~ (mean_y + 0.611)
  ))

# overall error distribution compared to the actual target location
ggplot() +
  geom_path(data = overall_ellipses, 
            aes(x = x_coord, y = y_coord, group = target),
            color = "black", linewidth = 1.5) +
  geom_point(data = overall_target_summary,
             aes(x = x_coord, y = y_coord),
             color = "black", size = 3.5, shape = 18) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)") +
  coord_fixed() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y),
             color = "purple", shape = "cross", size = 1.5, 
             stroke = 2) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67))

# plot with each players error distributions overlaid 
ggplot() +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y),,
             color = "purple", shape = "cross", size = 1.5,
             stroke = 1.5) +
  geom_path(data = player_ellipses,
            aes(x = x_coord, y = y_coord, color = factor(player_id),
                group = interaction(player_id, target)),
            alpha = 0.6,
            linewidth = 0.65) +
  geom_point(data = player_target_summary,
             aes(x = x_coord, y = y_coord, color = factor(player_id)),
             size = 1.75, alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~ factor(target, levels = c("LH", "RH", "LL", "RL"))) +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Player ID") +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0,1.34, 2.67))

# mutating to add the overall as a "player"
overall_target_summary <- overall_target_summary |>
  mutate(player_id = 4) |>
  dplyr::select(target, player_id, n, mean_x, sd_x, 
                mean_y, sd_y, x_coord, y_coord)
target_summary <- player_target_summary |>
  bind_rows(overall_target_summary)

overall_ellipses <- overall_ellipses |>
  mutate(player_id = 4) |>
  dplyr::select(player_id, target, x, y, x_coord, y_coord)
ellipses <- overall_ellipses |>
  bind_rows(player_ellipses)

# adding the average error distribution onto the same plot as the others
ggplot() +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y),,
             color = "purple", shape = "cross", size = 1.5,
             stroke = 1.25) +
  geom_path(data = ellipses, 
            aes(x = x_coord, y = y_coord, 
                group = interaction(player_id, target), 
                color = factor(player_id),
                linetype = factor(player_id)),
            linewidth = 0.75, alpha = 0.8) +
  geom_point(data = target_summary,
             aes(x = x_coord, y = y_coord, color = factor(player_id)),
             size = 1.75, alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~ factor(target, levels = c("LH", "RH", "LL", "RL"))) +
  labs(x = "X (yds)",
       y = "Y (yds)") +
  coord_fixed() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_linetype_manual(values = c("1" = "solid", 
                                   "2" = "solid", 
                                   "3" = "solid", 
                                   "4" = "dashed"),
                        guide = "none") +
  scale_color_manual(values = c("1" = "#F8766D", "2" = "#00BA38", 
                                "3" = "#619CFF", "4" = "black"),
                     labels = c("1" = "1", "2" = "2", 
                                "3" = "3", "4" = "Average"),
                     name = "Player ID") 


## individual error distributions
# shooter 1
ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y,
                 color = target),
             shape = "cross", size = 1.5,
             stroke = 1.5) +
  geom_path(data = player_ellipses |> filter(player_id == 1),
            aes(x = x_coord, y = y_coord,
                group = interaction(target),
                color = target),
            alpha = 0.6,
            linewidth = 0.65) +
  geom_point(data = player_target_summary |> filter(player_id == 1),
             aes(x = x_coord, y = y_coord, color = target),
             size = 1.75, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target",
       caption = "Each ellipse shows the spread of shot locations for a given target zone, with the center point marking the average shot location.\nThe 'X' marks the intended target. Tighter ellipses indicate more consistent shot placement."
  ) +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_color_manual(values = c("LH" = "#F8766D", "LL" = "#7CAE00",
                                "RH" = "#00BFC4", "RL" = "#C77CFF"),
                     labels = c("LH" = "Left High", "LL" = "Left Low",
                                "RH" = "Right High", "RL" = "Right Low")) +
  theme(plot.caption = element_text(hjust = 0.5,
                                    vjust = 0),
        plot.caption.position = "plot")

# shooter 2
ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y,
                 color = target),
             shape = "cross", size = 1.5,
             stroke = 1.5) +
  geom_path(data = player_ellipses |> filter(player_id == 2),
            aes(x = x_coord, y = y_coord,
                group = interaction(target),
                color = target),
            alpha = 0.6,
            linewidth = 0.65) +
  geom_point(data = player_target_summary |> filter(player_id == 2),
             aes(x = x_coord, y = y_coord, color = target),
             size = 1.75, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target",
       caption = "Each ellipse shows the spread of shot locations for a given target zone, with the center point marking the average shot location.\nThe 'X' marks the intended target. Tighter ellipses indicate more consistent shot placement."
  ) +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_color_manual(values = c("LH" = "#F8766D", "LL" = "#7CAE00",
                                "RH" = "#00BFC4", "RL" = "#C77CFF"),
                     labels = c("LH" = "Left High", "LL" = "Left Low",
                                "RH" = "Right High", "RL" = "Right Low")) +
  theme(plot.caption = element_text(hjust = 0.5,
                                    vjust = 0),
        plot.caption.position = "plot")

# shooter 3
ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y,
                 color = target),
             shape = "cross", size = 1.5,
             stroke = 1.5) +
  geom_path(data = player_ellipses |> filter(player_id == 3),
            aes(x = x_coord, y = y_coord,
                group = interaction(target),
                color = target),
            alpha = 0.6,
            linewidth = 0.65) +
  geom_point(data = player_target_summary |> filter(player_id == 3),
             aes(x = x_coord, y = y_coord, color = target),
             size = 1.75, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target",
       caption = "Each ellipse shows the spread of shot locations for a given target zone, with the center point marking the average shot location.\nThe 'X' marks the intended target. Tighter ellipses indicate more consistent shot placement."
  ) +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_color_manual(values = c("LH" = "#F8766D", "LL" = "#7CAE00",
                                "RH" = "#00BFC4", "RL" = "#C77CFF"),
                     labels = c("LH" = "Left High", "LL" = "Left Low",
                                "RH" = "Right High", "RL" = "Right Low")) +
  theme(plot.caption = element_text(hjust = 0.5,
                                    vjust = 0),
        plot.caption.position = "plot")

## projected shot placement (if not saved/deflected)
# shooter 1
pks |>
  filter(player_id == 1) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(aes(x = projected_x, y = projected_y, 
                 color = target, shape = result),
             size = 2, alpha = 0.8) +
  geom_point(data = pks |> group_by(target) |> summarize(mean_x = mean(target_x),
                                                         mean_y = mean(target_y)),
             aes(x = mean_x, y = mean_y, color = target),
             shape = 4, stroke = 1.5, size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target") +
  coord_fixed() +
  scale_shape_manual(values = c("G" = 19, "S" = 1, "F" = 2),
                     labels = c("G" = "Goal", "S" = "Save",
                                "F" = "On Frame"),
                     name = "Result") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67))

# shooter 2
pks |>
  filter(player_id == 2) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(aes(x = projected_x, y = projected_y, 
                 color = target, shape = result),
             size = 2, alpha = 0.8) +
  geom_point(data = pks |> group_by(target) |> summarize(mean_x = mean(target_x),
                                                         mean_y = mean(target_y)),
             aes(x = mean_x, y = mean_y, color = target),
             shape = 4, stroke = 1.5, size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target") +
  coord_fixed() +
  scale_shape_manual(values = c("G" = 19, "S" = 1, "F" = 2),
                     labels = c("G" = "Goal", "S" = "Save",
                                "F" = "On Frame"),
                     name = "Result") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67))

# shooter 3
pks |>
  filter(player_id == 3) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(aes(x = projected_x, y = projected_y, 
                 color = target, shape = result),
             size = 2, alpha = 0.8) +
  geom_point(data = pks |> group_by(target) |> summarize(mean_x = mean(target_x),
                                                         mean_y = mean(target_y)),
             aes(x = mean_x, y = mean_y, color = target),
             shape = 4, stroke = 1.5, size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(x = "X (yds)",
       y = "Y (yds)",
       color = "Target") +
  coord_fixed() +
  scale_shape_manual(values = c("G" = 19, "S" = 1, "F" = 2),
                     labels = c("G" = "Goal", "S" = "Save",
                                "F" = "On Frame"),
                     name = "Result") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67))

## raw shot data and speed plot (combining 4 plots)
# handling the 2 points that are basically on top of each other when plotted
special_points <- pks |> 
  filter(player_id == 1) |> 
  filter(target == "LH" | target == "LL") |> 
  filter(x == -2.960 | x == -2.924) |>
  mutate(jitter_x = jitter(x, amount = 0.1),
         jitter_y = jitter(y, amount = 0.1)) |>
  dplyr::select(x, y, jitter_x, jitter_y, target, result) |>
  mutate(result = fct_collapse(result, "S/F" = c("S", "F")))

s1 <- pks |>
  filter(player_id == 1) |>
  filter(x != -2.960 & x!= -2.924) |>
  mutate(result = fct_collapse(result, "S/F" = c("S", "F"))) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y, 
                 color = target, alpha = 0.8), 
             shape = 4, size = 2, stroke = 1.5) +
  geom_point(aes(x = x, y = y, color = target,
                 shape = result), size = 2, alpha = 0.8) +
  geom_point(data = special_points, aes(x = jitter_x, y = jitter_y, 
                                        color = target, shape = result),
             size = 2, alpha = 0.8) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0, 1.34, 2.67)) +
  scale_shape_manual(values = c("G" = 19, "S/F" = 1),
                     labels = c( "G" = "Goal", "S/F" = "Save/On Frame"),
                     name = "Result:") +
  labs(x = "X (yds)",
       y = "Y (yds)",
       title = "S1", 
       color = "Target:") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.49)) +
  guides(alpha = "none")

s2 <- pks |>
  filter(player_id == 2) |>
  mutate(result = fct_collapse(result, "S/F" = c("S", "F"))) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y, 
                 color = target, alpha = 0.8), 
             shape = 4, size = 2, stroke = 1.5) +
  geom_point(aes(x = x, y = y, color = target,
                 shape = result), size = 2, alpha = 0.8) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0,1.34, 2.67)) +
  scale_shape_manual(values = c("G" = 19, "S/F" = 1),
                     labels = c( "G" = "Goal", "S/F" = "Save/On Frame"),
                     name = "Result") +
  labs(x = "X (yds)",
       y = "",
       title = "S2") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.49)) +
  guides(alpha = "none",
         color = "none",
         shape = "none")


special_points2 <- pks |>
  filter(player_id == 3) |>
  filter(x < -3) |>
  filter(y > 0.67 & y < 1.34)

special_points2 <- special_points2 |>
  mutate(jitter_x = jitter(x, amount = 0.1),
         jitter_y = jitter(y, amount = 0.1)) |>
  mutate(result = fct_collapse(result, "S/F" = c("S", "F")))

s3 <- pks |>
  filter(player_id == 3) |>
  filter(x != -3.350 & x != -3.353) |>
  mutate(result = fct_collapse(result, "S/F" = c("S", "F"))) |>
  ggplot() +
  geom_polygon(data = left_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = right_post, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_polygon(data = crossbar, aes(x = x, y = y),
               fill = "black", color = "black") +
  geom_point(data = pks |> group_by(target) |>
               summarize(mean_target_x = mean(target_x),
                         mean_target_y = mean(target_y)), 
             aes(x = mean_target_x, y = mean_target_y, 
                 color = target, alpha = 0.8), 
             shape = 4, size = 2, stroke = 1.5) +
  geom_point(aes(x = x, y = y, color = target,
                 shape = result), size = 2, alpha = 0.8) +
  geom_point(data = special_points2, aes(x = jitter_x, y = jitter_y,
                                         color = target, shape = result),
             size = 2, alpha = 0.8) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(0,1.34, 2.67)) +
  scale_shape_manual(values = c("G" = 19, "S/F" = 1),
                     labels = c( "G" = "Goal", "S/F" = "Save/On Frame"),
                     name = "Result") +
  labs(x = "X (yds)",
       y = "Y (yds)",
       title = "S3") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.49)) +
  guides(alpha = "none",
         color = "none",
         shape = "none")

speed <- pks |>
  group_by(player_id) |>
  ggplot(aes(x = speed_mph, color = factor(player_id, 
                                              levels = c(3,2,1),
                                              labels = c("S3", "S2", "S1")))) +
  geom_boxplot() +
  geom_text(
    data = pks |>
      group_by(player_id) |>
      summarize(speed_mph = min(speed_mph)),
    aes(x = speed_mph - 2, y = c(0.25, 0, -0.25),
        label = paste0("S", player_id)),
    hjust = -0.3,
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(
    x = "Speed (MPH)",
    y = NULL,
    color = "Player ID",
    title = "Speed by Player"
  ) +
  scale_color_manual(values = c("steelblue", "#3F704D", "#F08080")) +
  scale_x_continuous(breaks = c(40, 44, 48, 52, 56)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  guides(color = "none")


combined <- (s1 + s2 + s3 + speed) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



### individualized speed data plots for each player by target
# change the targets to be factors and add overall a factor
pks2 <- pks |>
  mutate(target = "Overall")

plot_data <- bind_rows(pks2, pks)

plot_data <- plot_data |>
  mutate(target = factor(target, 
                         levels = c("Overall", "RL", "LL", "RH", "LH")))
plot_data <- plot_data |>
  mutate(target = fct_recode(target,
                             "Overall" = "Overall",
                             "Right Low" = "RL",
                             "Left Low" = "LL",
                             "Right High" = "RH",
                             "Left High" = "LH"))


# shooter 1
plot_data |>
  filter(player_id == 1) |>
  ggplot() +
  geom_boxplot(aes(x = speed_mph, y = target,
                   color = target), linewidth = 0.15) +
  geom_point(aes(x = speed_mph, y = target,
                 color = target), stroke = 1.5) +
  geom_point(data = plot_data |> filter(player_id == 1) |> group_by(target) |>
               summarize(med_speed = median(speed_mph)),
             aes(x = med_speed, y = target), color = "black",
             shape = 2, stroke = 1, size = 3) +
  theme_bw() +
  labs(x = "Speed (MPH)",
       y = "Target",
       title = "Speed Distribution by Target",
       subtitle = "Shooter 1",
       caption = "The black triangle represents the median speed for that target") +
  scale_color_manual(values = c("Left High" = "#F3987F", "Left Low" = "#4DBBD5",
                                "Right High" = "#91D1C2", "Right Low" = "#8491B4",
                                "Overall" = "grey40")) +
  guides(color = "none") +
  scale_x_continuous(breaks = c(41, 44, 47, 50, 53, 56))


# shooter 2
plot_data |>
  filter(player_id == 2) |>
  ggplot() +
  geom_boxplot(aes(x = speed_mph, y = target,
                   color = target), linewidth = 0.15) +
  geom_point(aes(x = speed_mph, y = target,
                 color = target), stroke = 1.5) +
  geom_point(data = plot_data |> filter(player_id == 1) |> group_by(target) |>
               summarize(med_speed = median(speed_mph)),
             aes(x = med_speed, y = target), color = "black",
             shape = 2, stroke = 1, size = 3) +
  theme_bw() +
  labs(x = "Speed (MPH)",
       y = "Target",
       title = "Speed Distribution by Target",
       subtitle = "Shooter 2",
       caption = "The black triangle represents the median speed for that target") +
  scale_color_manual(values = c("Left High" = "#F3987F", "Left Low" = "#4DBBD5",
                                "Right High" = "#91D1C2", "Right Low" = "#8491B4",
                                "Overall" = "grey40")) +
  guides(color = "none") +
  scale_x_continuous(breaks = c(45, 47, 49, 51, 53, 55, 57))

# shooter 3
plot_data |>
  filter(player_id == 3) |>
  ggplot() +
  geom_boxplot(aes(x = speed_mph, y = target,
                   color = target), linewidth = 0.15) +
  geom_point(aes(x = speed_mph, y = target,
                 color = target), stroke = 1.5) +
  geom_point(data = plot_data |> filter(player_id == 1) |> group_by(target) |>
               summarize(med_speed = median(speed_mph)),
             aes(x = med_speed, y = target), color = "black",
             shape = 2, stroke = 1, size = 3) +
  theme_bw() +
  labs(x = "Speed (MPH)",
       y = "Target",
       title = "Speed Distribution by Target",
       subtitle = "Shooter 3",
       caption = "The black triangle represents the median speed for that target") +
  scale_color_manual(values = c("Left High" = "#F3987F", "Left Low" = "#4DBBD5",
                                "Right High" = "#91D1C2", "Right Low" = "#8491B4",
                                "Overall" = "grey40")) +
  guides(color = "none") +
  scale_x_continuous(breaks = c(41, 44, 47, 50, 53, 56))

