

# Death age plots -- Lifespan
grouped_df <- results %>%
  group_by(id) %>%
  summarise(across(everything(), last))

# Filter the data for each type of agent
last_row_df_adults <- filter(grouped_df, type == "adult")
last_row_df_juvenile <- filter(grouped_df, type == "juvenile")
last_row_df_eggs <- filter(grouped_df, type == "eggmass")

# Plot a histogram of the lifespan of adults
ggplot(last_row_df_adults, aes(x = Age/365.0)) +
  geom_histogram(bins = 10) +
  labs(x = "Age", y = "Frequency", title = "Lifespan of Adults")

library(dplyr)
library(ggplot2)

# Length at puberty
adults_df <- filter(results, type == "adult")

# Filter take the row where Age == t_puberty
missing_values_df <- filter(adults_df, is.na(Age) | is.na(t_puberty))

# Remove rows from adults_df that are in missing_values_df
filtered_df <- filter(adults_df, !is.na(Age) & !is.na(t_puberty) & Age == t_puberty)

# Create a new plot with a boxplot of Lw
p1 <- ggplot(filtered_df, aes(x = 1, y = Lw)) +
  geom_boxplot(fill = "blue", alpha = 0.75, width = 0.7) +
  labs(title = "Length and Weight at Puberty", y = "Length (cm)") +
  theme(legend.position = "none")

# Add a boxplot of Ww to the right y-axis
p2 <- ggplot(filtered_df, aes(x = 1, y = Ww)) +
  geom_boxplot(fill = "orange", alpha = 0.75, width = 0.7) +
  labs(y = "Weight (g)") +
  theme(legend.position = "none")

# Ultimate size observed in adults
grouped_df <- adults_df %>% group_by(id)
max_size_df <- grouped_df %>% summarise(max_Lw = max(Lw, na.rm = TRUE))

# Histogram of Adults Max Lw
ggplot(max_size_df, aes(x = max_Lw)) +
  geom_histogram(bins = 25, fill = "blue", alpha = 0.75) +
  labs(x = "Max Lw (cm)", y = "Frequency", title = "Adults Max Length (cm)") +
  theme(legend.position = "none")

# Ww and Lw and R for age classes
adults_df <- adults_df %>% mutate(Age_year = ceiling(Age / 365))
adults_df <- na.omit(adults_df)

# Create the boxplots
p1 <- ggplot(adults_df, aes(x = Age_year, y = Lw)) +
  geom_boxplot(fill = "green", alpha = 0.75, width = 0.7) +
  labs(title = "Length of Age classes", y = "Length (cm)", x = "Age (years)") +
  theme(legend.position = "none")

p2 <- ggplot(adults_df, aes(x = Age_year, y = Ww)) +
  geom_boxplot(fill = "orange", alpha = 0.75, width = 0.7) +
  labs(title = "Weight of Age classes", y = "Weight (g)", x = "Age (years)") +
  theme(legend.position = "none")

p3 <- ggplot(adults_df, aes(x = Age_year, y = R)) +
  geom_boxplot(fill = "blue", alpha = 0.75, width = 0.7) +
  labs(title = "R: reproduction investment", y = "Energy (J)", x = "Age (years)") +
  theme(legend.position = "none")

# Ww and Lw for age classes and steps
grouped_adults <- adults_df %>% group_by(Age_year, time)
mean_df <- grouped_adults %>% summarise(Ww_mean = mean(Ww, na.rm = TRUE), Lw_mean = mean(Lw, na.rm = TRUE))

ggplot(mean_df, aes(x = time, y = Ww_mean, color = factor(Age_year))) +
  geom_line() +
  labs(x = "Step", y = "Mean Ww", title = "Mean Ww by Age_year and Step") +
  theme(legend.position = "topleft")

# length frequencies juveniles
juve_df <- filter(results, type == "juvenile")
grouped_df <- juve_df %>% group_by(id)
max_size_df <- grouped_df %>% summarise(max_Lw = max(Lw, na.rm = TRUE))

# Histogram of Juvenile Max Lw
ggplot(max_size_df, aes(x = max_Lw)) +
  geom_histogram(bins = 25, fill = "blue", alpha = 0.75) +
  labs(x = "Max Lw", y = "Frequency", title = "Histogram of Juvenile Max Lw") +
  theme(legend.position = "none")

library(lattice)
library(latticeExtra)

# Assuming 'filtered_df' is your data frame in R with the same structure as in Julia
# Create a boxplot for 'Lw'
p1 <- bwplot(Lw ~ 1, data = filtered_df, ylab = "Length (cm)", main = "Length and Weight at Puberty", box.fill = TRUE, fill = "blue", alpha = 0.75, horizontal = FALSE)

# Create a boxplot for 'Ww'
p2 <- bwplot(Ww ~ 1, data = filtered_df, ylab = "Weight (g)", box.fill = TRUE, fill = "orange", alpha = 0.75, horizontal = FALSE)

# Combine the two plots with a double y-axis
doubleYScale(p1, p2, style1 = 1, style2 = 1, add.ylab2 = TRUE, text = c("Length (cm)", "Weight (g)"))

# demography post processing Multi-spel agents simulations
results <- read.csv("~/PhD/Multi-SPelAgents/results_baseline_last2y_15_0945_multispel_vectorparams.csv")

library(ggplot2)
library(tidyr)
library(dplyr)

# Convert the type and id columns to factors
results$type <- as.factor(results$type)
results$id <- as.factor(results$id)

# Group by time and type, and count the number of agents for each group
summary_df <- results %>%
  group_by(time, type) %>%
  summarise(count = n())

# Plot a line for each type of agent
ggplot(summary_df, aes(x = time, y = count, color = type)) +
  geom_line() +
  labs(x = "Step", y = "Count", title = "Number of Agents by Type") +
  theme(legend.position = "topleft")

# Group by time and type, and sum the Ww for each group
summary_Ww <- results %>%
  group_by(time, type) %>%
  summarise(total_Ww = sum(Ww))

# Plot a line for each type of agent
ggplot(summary_Ww, aes(x = time, y = total_Ww, color = type)) +
  geom_line() +
  labs(x = "Step", y = "Total Biomass", title = "Ww of Agents by Type") +
  theme(legend.position = "topleft")

# Length/age number of adults for each bins

results <- results %>% 
  group_by(Lw, time) %>% 
  summarise(tot_B = sum(Ww, na.rm = TRUE), count = n())

# Group by length bin, then calculate the mean total biomass and count over the year
summary_df <- results %>% 
  group_by(Lw_bin) %>% 
  summarise(mean_tot_B = mean(tot_B, na.rm = TRUE), mean_count = mean(count, na.rm = TRUE))

# Plot the mean total biomass for each length bin
ggplot(summary_df, aes(x = Lw_bin, y = mean_tot_B, fill = mean_tot_B)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Length Bin (cm)", y = "Mean Total Biomass (g)", title = "Mean Total Biomass by Length over the Year") +
  theme_minimal()

# Plot the mean count for each length bin
ggplot(summary_df, aes(x = Lw_bin, y = mean_count, fill = mean_count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Length Bin (cm)", y = "Mean Count", title = "Mean Count by Length over the Year") +
  theme_minimal()

  # do the same for the ages