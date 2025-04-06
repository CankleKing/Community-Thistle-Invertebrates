

                        # SUPER CLEAN CODE 9000 #

################################################################################
#################### Random generation of samples ##############################
################################################################################

# List of codes: Field_(A/B/C if required)_Time-block 
codes <- c(
  "W_A_1", "W_A_2", "W_A_3", 
  "W_B_1", "W_B_2", "W_B_3",
  "S_1", "S_2", "S_3", 
  "BP_1", "BP_2", "BP_3", 
  "BC_A_1", "BC_A_2", "BC_A_3", 
  "BC_B_1", "BC_B_2", "BC_B_3", 
  "BC_C_1", "BC_C_2", "BC_C_3", 
  "M_1", "M_2", "M_3"
)

# Function to randomly choose one code and exclude it from the list
random_choose_exclude <- function(codes_list) {
  if (length(codes_list) == 0) {
    stop("Empty list")
  }
  
  # Randomly choose one code from the list
  chosen_code <- sample(codes_list, 1)
  
  # Exclude the chosen code from the list
  remaining_codes <- codes_list[codes_list != chosen_code]
  
  # Return the chosen code and the remaining codes
  return(list(chosen_code = chosen_code, remaining_codes = remaining_codes))
}

# Set a seed for reproducibility
set.seed(5318008)

# Main loop to randomly choose codes one by one and exclude them from the list
num_codes_to_choose <- length(codes)
chosen_codes <- c()
for (i in 1:num_codes_to_choose) {
  result <- random_choose_exclude(codes)
  chosen_code <- result$chosen_code
  codes <- result$remaining_codes
  chosen_codes <- c(chosen_codes, chosen_code)
}

# Print the randomly chosen codes
elements <- print(chosen_codes)

# Create a data frame from the list
table_df <- data.frame(Order = 1:length(elements), Field_time = elements)

# Print the resulting table to create the random order of fields sampled 
print(table_df)





################################################################################
############################# Merging dataframes ###############################
################################################################################

# Load in the dataframes via "Import Dataset" button in Environment (top right)
# https://github.com/CankleKing/Community-Thistle-Invertebrates

library("dplyr")
# Remove the absence data
species_list_prez <- subset(species_list, count != 0)

# Merge species_list and focal_attributes data frames based on the "code" column
merged_data <- merge(species_list, focal_attributes, by = "code", 
                     all.x = TRUE)

# Check the structure of the merged data frame
str(merged_data)

# Check the number of unique codes and observations after merging
num_unique_codes <- length(unique(merged_data$code))
print(num_unique_codes)
num_observations <- nrow(merged_data)
print(num_observations)

# Select the columns containing wind speeds T0 to T60
wind_speed_columns <- c("wind_t0", "wind_t15", "wind_t30", "wind_t45", 
                        "wind_t60")

# Calculate the average wind speed for each row
merged_data$average_wind_speed <- rowMeans(merged_data[wind_speed_columns], 
                                           na.rm = TRUE)

# Convert the date.x to a date format
merged_data <- merged_data %>%
  mutate(date = as.Date(date.x, format = "%d/%m/%Y"))

# Delete unwanted columns to make a simplified dataframe
merged_data$Written <- NULL
merged_data$written <- NULL
merged_data$pic_obs_samp <- NULL
merged_data$time_sampled_mins <- NULL
merged_data$time_sampled_secs <- NULL
merged_data$Parasitised <- NULL
merged_data$notes <- NULL
merged_data$Key <- NULL
merged_data$`Key pathway`<- NULL
merged_data$...25 <- NULL
merged_data$date.y <- NULL
merged_data$field.y <- NULL
merged_data$time.y <- NULL
merged_data$focal_patch.y <- NULL
merged_data$wind_t0 <- NULL
merged_data$wind_t15 <- NULL
merged_data$wind_t30 <- NULL
merged_data$wind_t45 <- NULL
merged_data$wind_t60 <- NULL
merged_data$wind_bearing <- NULL
merged_data$weather_obvs <- NULL
merged_data$surrounding_veg_type <- NULL
merged_data$date.x <- NULL


# Calculate leaf attributes
merged_data$no_leaves <- merged_data$leaves_alive + merged_data$leaves_dead


####################### Working out distribution data ##########################
library(geometry)
library(tidyr)  
library(ggplot2)
library(sp)
library(geosphere) 

# Count the number of neighbouring thistles for each focal thistle (code)
dist_df <- distribution_values %>%
  group_by(code) %>%
  summarise(num_neigh = n())

# Find the nearest neighbour for each focal thistle 
nearest_neigh <- distribution_values %>%
  group_by(code) %>%
  summarise(min_distance = min(distance, na.rm = TRUE))

# If nearest neighbour is inf replace as blank (W_B_2.7)
nearest_neigh <- nearest_neigh %>%
  mutate(min_distance = as.numeric(replace(min_distance, 
                                           is.infinite(min_distance), NA)))

# Merge the two data frames based on the 'code' column
distribtion_df <- merge(nearest_neigh, dist_df, by = "code", all.x = TRUE)

# Merge species_list and focal_attributes data frames based on the "code" column
merged_data <- merge(merged_data, distribtion_df, by = "code", all.x = TRUE)




#### Performing Convex Hull Analysis 
# Function to convert polar (bearing and distance) to Cartesian coordinates
polar_to_cartesian <- function(bearing, distance) {
  radians <- bearing * (pi / 180)  # Convert degrees to radians
  x <- distance * cos(radians)
  y <- distance * sin(radians)
  return(data.frame(x, y))
}

# Apply the conversion to the dataset
cartesian_coords <- distribution_values %>%
  rowwise() %>%
  mutate(cartesian = list(polar_to_cartesian(bearing, distance))) %>%
  unnest_wider(cartesian)

#### Function to compute hull distance for each focal thistle
# Get unique codes for focal thistles
unique_codes <- unique(cartesian_coords$code)

# Initialize an empty dataframe to store the signed distances
edge_distance_df <- data.frame(code = character(), 
                                 edge_distance = numeric(), 
                                 stringsAsFactors = FALSE)

# Function to generate evenly spaced points along a line segment
generate_points_on_segment <- function(x1, y1, x2, y2, spacing = 0.01) {
  # Calculate the total length of the segment
  segment_length <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
  # Number of points to generate along the segment
  num_points <- ceiling(segment_length / spacing)
  
  # Generate points along the segment
  x_points <- seq(x1, x2, length.out = num_points)
  y_points <- seq(y1, y2, length.out = num_points)
  
  return(data.frame(x = x_points, y = y_points))
}

# Loop through each unique focal thistle code
for (code in unique_codes) {
  # Filter data for the focal thistle
  focal_data <- cartesian_coords %>% filter(code == !!code)
  
  # Extract surrounding thistles' coordinates
  surrounding_thistles <- focal_data[, c("x", "y")]
  
  if (nrow(surrounding_thistles) >= 3) {# Convex hull requires at least 3 points
    # Calculate the convex hull using chull
    hull_indices <- chull(surrounding_thistles$x, surrounding_thistles$y)
    
    # Create a data frame for the convex hull points
    hull_points <- surrounding_thistles[hull_indices, , drop = FALSE]
    
    # Create a new data frame to hold all points along the convex hull edges
    edge_points <- data.frame(x = numeric(0), y = numeric(0))
    
    # Loop over each edge of the convex hull and generate points along the edges
    for (i in seq_len(nrow(hull_points))) {
      # Get the start and end points of the current edge
      x1 <- hull_points$x[i]
      y1 <- hull_points$y[i]
      x2 <- hull_points$x[ifelse(i == nrow(hull_points), 1, i + 1)]
      y2 <- hull_points$y[ifelse(i == nrow(hull_points), 1, i + 1)]
      
      # Generate points along this edge
      segment_points <- generate_points_on_segment(x1, y1, x2, y2, 
                                                   spacing = 0.01)
      
      # Append these points to the edge_points data frame
      edge_points <- rbind(edge_points, segment_points)
    }
    
    # Calculate distances from (0,0) to all edge points
    distances_to_edges <- sqrt((edge_points$x)^2 + (edge_points$y)^2)
    
    # Minimum distance to the hull edge
    min_distance <- min(distances_to_edges)
    
    # Find the coordinates of the point on the hull edge with the min distance
    min_edge_point <- edge_points[which.min(distances_to_edges), ]
    
    # Check if the focal thistle is within the hull
    is_within_hull <- point.in.polygon(0, 0, hull_points$x, hull_points$y) > 0
    
    # Assign sign based on location
    edge_distance <- if (is_within_hull) {
      -min_distance  # Inside hull
    } else {
      min_distance   # Outside hull
    }
    
    # Add the signed distance to the dataframe
    edge_distance_df <- rbind(edge_distance_df, 
                              data.frame(code = code, 
                                         edge_distance = edge_distance))
    
    # Create a plot
    p <- ggplot() +
      geom_point(data = surrounding_thistles, aes(x = x, y = y), 
                 color = "blue") +
      geom_polygon(data = hull_points, aes(x = x, y = y), 
                   fill = "lightblue", alpha = 0.5) +
      geom_point(aes(x = 0, y = 0), color = "red", size = 3) +  # Focal thistle
      geom_segment(aes(x = 0, y = 0, xend = min_edge_point$x, 
                       yend = min_edge_point$y), 
                   color = ifelse(edge_distance > 0, "red", "blue"), 
                   linetype = "dashed") +  # Line to the closest hull edge
      geom_text(aes(x = (min_edge_point$x)/2, y = (min_edge_point$y)/2, 
                    label = paste0("Dist: ", round(edge_distance, 2))), 
                color = "black", size = 3, vjust = -1) +  # Label the distance
      ggtitle(paste("Convex Hull for", code)) +
      xlab("X Coordinate") +
      ylab("Y Coordinate") +
      theme_minimal() +
      coord_equal()
    
    # Show the plot
    print(p)
    
  } else {
    message(paste("Not enough surrounding thistles for focal thistle", code))
    
    # Assign maximum edge distance when not enough surrounding thistles
    edge_distance_df <- rbind(edge_distance_df, 
                              data.frame(code = code, 
                                         edge_distance = 3))
  }
}

# Merge the signed distance data frame with the merged_data data frame
merged_data <- merge(merged_data, edge_distance_df, by = "code", all.x = TRUE)





################## Creating the species richness column ########################

# Count the number of unique species for each focal thistle (code), 
# excluding NA values
focal_richness <- species_list %>%
  group_by(code) %>%
  summarise(number_of_species = n_distinct(species_name[!is.na(species_name)]), 
            .groups = 'drop')

# Print the result
print(focal_richness)


# Merge total abundance and species counts data frames
merged_data <- merge(merged_data, focal_richness, 
                     by = "code", all.x = TRUE)




######################## Creating the Shannon Index column #####################

# Function to calculate Shannon Index
calculate_shannon <- function(abundance) {
  p <- abundance / sum(abundance)  # Calculate proportion
  shannon_index <- -sum(p * log(p + 1e-10), na.rm = TRUE)  # Add small value to 
  return(shannon_index)                                    # avoid log(0)
  }


# Calculate the Shannon Index for each focal thistle
shannon_data <- species_list %>%
  group_by(code) %>%
  summarise(shannon = calculate_shannon(count))

# Set any negative Shannon Index values to 0.0001 (from small proportion)
# Helps to distinguish from absence
shannon_data$shannon[shannon_data$shannon < 0] <- 0.0001

# Merge the Shannon Index data with the merged_data
merged_data <- merge(merged_data, shannon_data, by = "code", all.x = TRUE)





######################## Creating the Simpson Index column #####################

# Function to calculate Simpson Index
calculate_simpson <- function(abundance) {
  p <- abundance / sum(abundance)  
  simpson_index <- sum(p^2, na.rm = TRUE)  
  simpson_index <- 1 - simpson_index # 1 - D makes interpreation more intuitive
  return(simpson_index)              # Making a high D mean higher diversity
}

# Calculate the Simpson Index for each focal thistle
simpson_data <- species_list %>%
  group_by(code) %>%
  summarise(simpson = calculate_simpson(count))

# Merge the Simpson Index data with the merged_data
merged_data <- merge(merged_data, simpson_data, by = "code", all.x = TRUE)




######################## Creating the Evenness Index column ####################

# Function to calculate Evenness Index (Pielou's Evenness)
calculate_evenness <- function(shannon, richness) {
  evenness_index <- shannon / log(richness)
  return(evenness_index)
}

# Calculate the Evenness Index for each focal thistle
merged_data <- merged_data %>%
  group_by(code) %>%
  mutate(
    evenness = ifelse(number_of_species > 1,
                      calculate_evenness(shannon, number_of_species), 
                      0)  # Set to 0 if there's 1 or 0 species
  ) %>%
  ungroup()




######################## Calculating sampling effort ###########################

# Convert the 'Date_sampled' column to date format
sampled_checklist$Date_sampled <- as.Date(sampled_checklist$Date_sampled)

# Omit any rows that where not sampled hence are not a thistle
sampled_checklist <- sampled_checklist[!is.na(sampled_checklist$Date_sampled), ]

# Create the 'nm_sampled_day' column
sampled_checklist <- sampled_checklist %>%
  group_by(Date_sampled) %>%
  mutate(nm_sampled_day = row_number())

# Merging merged_data and sampled_checklist
merged_data <- merge(merged_data, sampled_checklist[, c("code", 
                                           "nm_sampled_day")], by = "code", 
                     all.x = TRUE)

# Calculate effort 
merged_data$effort_day <- merged_data$count / merged_data$nm_sampled_day


# Create a data frame without data on order absence                             
merged_data_prez <- subset(merged_data, count != 0)




####################### Percentage of order occurrence #########################

# Count the total number of occurrences for each order
total_order_count <- species_list_prez %>%
  group_by(order) %>%
  summarise(order_count = n())

# Calculate the total number of observations in the entire dataset
total_observations <- nrow(species_list_prez)

# Calculate the percentage of each order
order_percentage <- total_order_count %>%
  mutate(percentage = (order_count / total_observations) * 100)

# View the result
print(order_percentage)






###################### Create a DF of only major orders ########################

# Clean Na's
omited_order_data <- merged_data[!is.na(merged_data$order), ]

# Create a list of the dominant orders 
orders_of_interest <- c("Araneae", "Coleoptera", "Collembola", "Diptera", 
                        "Hemiptera", "Hymenoptera")

# Creating a new data frame with only the specified orders
selected_order <- omited_order_data[omited_order_data$order %in% 
                                      orders_of_interest, ]

# Recalculate richness, shannon, simpson, and evenness for each patch
selected_order$number_of_species <- NULL
selected_order$shannon <- NULL
selected_order$simpson <- NULL
selected_order$evenness <- NULL

# Count the number of unique species for each order within each focal thistle 
focal_richness_select <- selected_order %>%
  filter(count != 0) %>%  # Exclude rows with zero counts
  group_by(code, order) %>%
  summarise(number_of_species = n_distinct(species_name), .groups = 'drop') %>%
  complete(code, order, fill = list(number_of_species = 0))  # Fill 0 

# Ensure no duplicate rows when merging by summarizing on 'code' and 'order'
focal_richness_select <- focal_richness_select %>%
  group_by(code, order) %>%
  summarise(number_of_species = sum(number_of_species), .groups = 'drop')  

# Print the result
print(focal_richness_select)

# Merge total abundance and species counts data frames for patches
selected_order <- merge(selected_order, focal_richness_select, 
                        by = c("code", "order"), all.x = TRUE)


# Calculate the Shannon Index for each patch
shannon_data_select <- selected_order %>%
  group_by(code) %>%
  summarise(shannon = ifelse(sum(count) > 0, calculate_shannon(count), 0), 
            .groups = 'drop') %>%
  distinct()  # Keep unique rows only

# Set any negative Shannon Index values to a small positive number if needed
shannon_data_select$shannon[shannon_data_select$shannon < 0] <- 0.0001

# Merge the Shannon Index data with the selected_order data frame
selected_order <- merge(selected_order, shannon_data_select, 
                        by = "code", all.x = TRUE)

# Calculate the Simpson Index for each patch
simpson_data_select <- selected_order %>%
  group_by(code) %>%
  summarise(simpson = ifelse(sum(count) > 0, calculate_simpson(count), 0),
            .groups = 'drop') %>%
  distinct()  # Keep unique rows only

# Merge the Simpson Index data with the selected_order data frame
selected_order <- merge(selected_order, simpson_data_select, 
                        by = "code", all.x = TRUE)

# Calculate the Evenness Index for each patch
selected_order <- selected_order %>%
  group_by(code) %>%
  mutate(evenness = ifelse(number_of_species > 1, 
                           calculate_evenness(shannon, number_of_species), 
                           0)) %>%  # Set to 0 if there's 1 or 0 species
  ungroup()

# View the new data frame
View(selected_order)

# Create a data frame without data on species absence
selected_order_prez <- subset(selected_order, count != 0)







################################################################################
############# Create bar plots of abundance & richness by order ################
################################################################################

############################### Species Richness ###############################
# Calculate species richness (number of unique species) for each order and patch 
species_richness_merged <- merged_data_prez %>%
  group_by(order, code, field.x, f_type) %>%
  summarise(SpeciesRichness = n_distinct(species_name), .groups = "drop")

# Box plot for species richness
species_richness_plot_merged <- species_richness_merged %>%
  ggplot(aes(x = order, y = SpeciesRichness)) +
  geom_boxplot(size = 0.55) +
  theme_bw() +
  labs(x = "Order", y = "Species Richness", 
       title = "Species Richness by Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

# Display the plot
species_richness_plot_merged

################################ Abundance #####################################
# Calculate total abundance for each order and patch (code)
abundance <- merged_data_prez %>%
  group_by(order, code, field.x, f_type) %>%
  summarise(Abundance = sum(count), .groups = "drop")

# Box plot for abundance
abundance_plot <- abundance %>%
  ggplot(aes(x = order, y = log(Abundance))) +
  geom_boxplot(size = 0.55) +
  theme_bw() +
  labs(x = "Order", y = "log Abundance", 
       title = "Abundance by Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

# Display the plot
abundance_plot




############################## Create summary tables ###########################

# Table 1: Abundance stats for each invertebrate order per patch
# Create a complete dataset with unique combinations of order and patch
unique_orders_codes <- unique(merged_data[, c("order", "code")])

# Create a complete dataset with all combinations of order and code
complete_data <- unique_orders_codes %>%
  complete(order, code = unique(merged_data$code), 
           fill = list(count = 0, number_of_species = 0))

# Combine the complete dataset with the original data
combined_data <- merge(complete_data, merged_data, by = c("order", "code"), 
                       all.x = TRUE)

# Fill NAs with zeros for count and number_of_species
combined_data$count[is.na(combined_data$count)] <- 0
combined_data$number_of_species[is.na(combined_data$number_of_species)] <- 0

# Calculate the number of unique species for each patch and order
species_counts <- merged_data %>%
  group_by(code, order) %>%
  summarise(number_of_species = n_distinct(species_name), .groups = 'drop')

# Ensure all combinations of orders and patches are represented
species_counts <- species_counts %>%
  complete(order, code = unique(merged_data$code), 
           fill = list(number_of_species = 0))

# Table 2: Species richness statistics for each invertebrate order per patch
species_richness_stats <- species_counts %>%
  group_by(order) %>%
  summarise(
    Mean_Species_Richness = round(mean(number_of_species, na.rm = TRUE), 3),  
    Median_Species_Richness = round(median(number_of_species, na.rm = TRUE), 3),
    SD_Species_Richness = round(sd(number_of_species, na.rm = TRUE), 3),
    Min_Species_Richness = round(min(number_of_species, na.rm = TRUE), 3),
    Max_Species_Richness = round(max(number_of_species, na.rm = TRUE), 3),
    Total_Species_Richness = n_distinct(
      merged_data$species_name[merged_data$order == unique(order)])  
  )

# Calculate overall total species richness across all orders
total_species_richness_stats <- merged_data %>%
  summarise(
    order = "Any",
    Mean_Species_Richness = round(mean(number_of_species, na.rm = TRUE), 3),  
    Median_Species_Richness = round(median(number_of_species, na.rm = TRUE), 3),
    SD_Species_Richness = round(sd(number_of_species, na.rm = TRUE), 3),
    Min_Species_Richness = round(min(number_of_species, na.rm = TRUE), 3),
    Max_Species_Richness = round(max(number_of_species, na.rm = TRUE), 3),
    Total_Species_Richness = n_distinct(merged_data$species_name)  
  )

# Append the overall statistics row to the species richness stats table
species_richness_stats <- bind_rows(species_richness_stats, 
                                    total_species_richness_stats)

# Print the species richness statistics table
print(species_richness_stats)






################################################################################
########################## Summarizing sampling days ###########################
################################################################################


# 1. Total number of unique thistle patches sampled
total_unique_patches <- sampled_checklist %>%
  summarise(total_patches = n_distinct(code))

# 2. Total sampled for each field
total_by_field <- sampled_checklist %>%
  group_by(Field) %>%
  summarise(total_patches = n_distinct(code))

# 3. Total sampled at each time point
total_by_time_point <- sampled_checklist %>%
  filter(!is.na(Time)) %>% # Ignore rows where 'Time' is NA
  group_by(Time) %>%
  summarise(total_patches = n_distinct(code))

# Output the results
cat("Total number of unique thistle patches sampled:", 
    total_unique_patches$total_patches, "\n")
cat("\nTotal sampled for each field:\n")
print(total_by_field)
cat("\nTotal sampled at each time point:\n")
print(total_by_time_point)






################################################################################
############################## Checking the data ###############################
################################################################################
# Save these as numeric
merged_data$number_of_species <- as.numeric(merged_data$number_of_species)
merged_data$nm_sampled_day <- as.numeric(merged_data$nm_sampled_day)
merged_data$edge_distance <- as.numeric(merged_data$edge_distance)

# Checking for skewdness in response variables
hist(merged_data$count) # +ve skew
hist(merged_data$number_of_species)
hist(merged_data$shannon)
hist(merged_data$simpson)
hist(merged_data$evenness) 



# Checking in explanatory variables 
hist(merged_data$height)
hist(merged_data$width) # +ve skew
hist(merged_data$leaves_alive) # +ve skew
hist(merged_data$leaves_dead) # +ve skew
hist(merged_data$no_leaves) # +ve skew
hist(merged_data$flower_bud) # +ve skew
hist(merged_data$flower_bloom) # +ve skew
hist(merged_data$flower_dead) # +ve skew
hist(merged_data$flower_seed) # +ve skew
hist(merged_data$humidity_base)
hist(merged_data$humidity_top)
hist(merged_data$micro_temp)
hist(merged_data$min_surrounding_veg_height) 
hist(merged_data$max_surrounding_veg_height) 
hist(merged_data$average_wind_speed) # +ve skew
hist(merged_data$min_distance) # +ve skew
hist(merged_data$num_neigh) # +ve skew
hist(merged_data$edge_distance)

# Perform log or log +1  transformation on selected variables in merged_data
merged_data$count_log1p <- log1p(merged_data$count) # response

merged_data$width_log <- log(merged_data$width) # natural log transformation
merged_data$leaves_alive_log <- log(merged_data$leaves_alive)
merged_data$no_leaves_log <- log(merged_data$no_leaves) 

merged_data$min_distance_log <- log(merged_data$min_distance)
merged_data$num_neigh_log <- log(merged_data$num_neigh)


merged_data$leaves_dead_log1p <- log1p(merged_data$leaves_dead) # Log + 1
merged_data$flower_bud_log1p <- log1p(merged_data$flower_bud)
merged_data$flower_bloom_log1p <- log1p(merged_data$flower_bloom)
merged_data$flower_dead_log1p <- log1p(merged_data$flower_dead)
merged_data$flower_seed_log1p <- log1p(merged_data$flower_seed)
merged_data$average_wind_speed_log1p <- log1p(merged_data$average_wind_speed)


#  Do the same for selected_order
selected_order$count_log1p <- log1p(selected_order$count) 

selected_order$width_log <- log(selected_order$width) 
selected_order$leaves_alive_log <- log(selected_order$leaves_alive)
selected_order$no_leaves_log <- log(selected_order$no_leaves) 

selected_order$min_distance_log <- log(selected_order$min_distance)
selected_order$num_neigh_log <- log(selected_order$num_neigh)

selected_order$leaves_dead_log1p <- log1p(selected_order$leaves_dead) 
selected_order$flower_bud_log1p <- log1p(selected_order$flower_bud)
selected_order$flower_bloom_log1p <- log1p(selected_order$flower_bloom)
selected_order$flower_dead_log1p <- log1p(selected_order$flower_dead)
selected_order$flower_seed_log1p <- log1p(selected_order$flower_seed)
selected_order$average_wind_speed_log1p <- log1p(
  selected_order$average_wind_speed)


merged_data_prez <- merged_data %>%
  filter(!is.na(species_name) & species_name != "")

selected_order_prez <- selected_order %>%
  filter(!is.na(species_name) & species_name != "")



############################# Correlation matrix ###############################

library(corrplot)

# Select specific columns from merged_data and create a new dataframe
corr_columns <- c("height", "width_log",  "leaves_dead_log1p", 
                       "leaves_alive_log", "no_leaves_log",
                       "max_surrounding_veg_height", 
                       "min_surrounding_veg_height", 
                       "micro_temp", "humidity_base","humidity_top",
                       "average_wind_speed_log1p", "min_distance_log", 
                       "num_neigh_log", "edge_distance")  

# Create a dataframe of variables for the correlation matrix 
corr_data <- merged_data_prez[corr_columns]

# Check that all values are as numeric 
corr_data <- corr_data[sapply(corr_data, is.numeric)]

# Create the correlation_matrix
correlation_matrix <- cor(corr_data, use = "pairwise.complete.obs")

# Create the correlation plot 
corrplot(correlation_matrix, method = "color", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", 
         number.cex = 0.8, 
         addgrid.col = "gray", 
         cl.pos = "n", 
         tl.cex = 1, 
         col = colorRampPalette(c("blue", "white", "red"))(100)) 




################################################################################
####################### Explanatory variable overview ##########################
################################################################################

# Select the relevant columns from merged_data
explanatory_vars <- subset(merged_data, 
                               select = c(code, field.x, time.x,
                                          height, width, leaves_alive, 
                                          leaves_dead, flower_bud, flower_bloom, 
                                          flower_dead, flower_seed, 
                                          humidity_base, humidity_top, 
                                          micro_temp, 
                                          min_surrounding_veg_height, 
                                          max_surrounding_veg_height, 
                                          average_wind_speed, min_distance, 
                                          num_neigh, edge_distance))
# Remove duplicate rows
explanatory_vars <- unique(explanatory_vars)

# Calculate summary statistics for all explanatory variables
summary_stats <- data.frame(
  Variable = c("Height", "Width", "Leaves Alive", "Leaves Dead", 
               "Flower Bud", "Flower Bloom", "Flower Dead", "Flower Seed",
               "Min Distance", "Num Neighbours", "Edge Distance",
               "Humidity Base", "Humidity Top", "Micro Temp", 
               "Min Surrounding Veg Height", "Max Surrounding Veg Height", 
               "Average Wind Speed"),
  
  Mean = c(
    mean(explanatory_vars$height, na.rm = TRUE),
    mean(explanatory_vars$width, na.rm = TRUE),
    mean(explanatory_vars$leaves_alive, na.rm = TRUE),
    mean(explanatory_vars$leaves_dead, na.rm = TRUE),
    mean(explanatory_vars$flower_bud, na.rm = TRUE),
    mean(explanatory_vars$flower_bloom, na.rm = TRUE),
    mean(explanatory_vars$flower_dead, na.rm = TRUE),
    mean(explanatory_vars$flower_seed, na.rm = TRUE),
    mean(explanatory_vars$min_distance, na.rm = TRUE),
    mean(explanatory_vars$num_neigh, na.rm = TRUE),
    mean(explanatory_vars$edge_distance, na.rm = TRUE),
    mean(explanatory_vars$humidity_base, na.rm = TRUE),
    mean(explanatory_vars$humidity_top, na.rm = TRUE),
    mean(explanatory_vars$micro_temp, na.rm = TRUE),
    mean(explanatory_vars$min_surrounding_veg_height, na.rm = TRUE),
    mean(explanatory_vars$max_surrounding_veg_height, na.rm = TRUE),
    mean(explanatory_vars$average_wind_speed, na.rm = TRUE)
  ),
  
  Median = c(
    median(explanatory_vars$height, na.rm = TRUE),
    median(explanatory_vars$width, na.rm = TRUE),
    median(explanatory_vars$leaves_alive, na.rm = TRUE),
    median(explanatory_vars$leaves_dead, na.rm = TRUE),
    median(explanatory_vars$flower_bud, na.rm = TRUE),
    median(explanatory_vars$flower_bloom, na.rm = TRUE),
    median(explanatory_vars$flower_dead, na.rm = TRUE),
    median(explanatory_vars$flower_seed, na.rm = TRUE),
    median(explanatory_vars$min_distance, na.rm = TRUE),
    median(explanatory_vars$num_neigh, na.rm = TRUE),
    median(explanatory_vars$edge_distance, na.rm = TRUE),
    median(explanatory_vars$humidity_base, na.rm = TRUE),
    median(explanatory_vars$humidity_top, na.rm = TRUE),
    median(explanatory_vars$micro_temp, na.rm = TRUE),
    median(explanatory_vars$min_surrounding_veg_height, na.rm = TRUE),
    median(explanatory_vars$max_surrounding_veg_height, na.rm = TRUE),
    median(explanatory_vars$average_wind_speed, na.rm = TRUE)
  ),
  
  SD = c(
    sd(explanatory_vars$height, na.rm = TRUE),
    sd(explanatory_vars$width, na.rm = TRUE),
    sd(explanatory_vars$leaves_alive, na.rm = TRUE),
    sd(explanatory_vars$leaves_dead, na.rm = TRUE),
    sd(explanatory_vars$flower_bud, na.rm = TRUE),
    sd(explanatory_vars$flower_bloom, na.rm = TRUE),
    sd(explanatory_vars$flower_dead, na.rm = TRUE),
    sd(explanatory_vars$flower_seed, na.rm = TRUE),
    sd(explanatory_vars$min_distance, na.rm = TRUE),
    sd(explanatory_vars$num_neigh, na.rm = TRUE),
    sd(explanatory_vars$edge_distance, na.rm = TRUE),
    sd(explanatory_vars$humidity_base, na.rm = TRUE),
    sd(explanatory_vars$humidity_top, na.rm = TRUE),
    sd(explanatory_vars$micro_temp, na.rm = TRUE),
    sd(explanatory_vars$min_surrounding_veg_height, na.rm = TRUE),
    sd(explanatory_vars$max_surrounding_veg_height, na.rm = TRUE),
    sd(explanatory_vars$average_wind_speed, na.rm = TRUE)
  ),
  
  Min = c(
    min(explanatory_vars$height, na.rm = TRUE),
    min(explanatory_vars$width, na.rm = TRUE),
    min(explanatory_vars$leaves_alive, na.rm = TRUE),
    min(explanatory_vars$leaves_dead, na.rm = TRUE),
    min(explanatory_vars$flower_bud, na.rm = TRUE),
    min(explanatory_vars$flower_bloom, na.rm = TRUE),
    min(explanatory_vars$flower_dead, na.rm = TRUE),
    min(explanatory_vars$flower_seed, na.rm = TRUE),
    min(explanatory_vars$min_distance, na.rm = TRUE),
    min(explanatory_vars$num_neigh, na.rm = TRUE),
    min(explanatory_vars$edge_distance, na.rm = TRUE),
    min(explanatory_vars$humidity_base, na.rm = TRUE),
    min(explanatory_vars$humidity_top, na.rm = TRUE),
    min(explanatory_vars$micro_temp, na.rm = TRUE),
    min(explanatory_vars$min_surrounding_veg_height, na.rm = TRUE),
    min(explanatory_vars$max_surrounding_veg_height, na.rm = TRUE),
    min(explanatory_vars$average_wind_speed, na.rm = TRUE)
  ),
  
  Max = c(
    max(explanatory_vars$height, na.rm = TRUE),
    max(explanatory_vars$width, na.rm = TRUE),
    max(explanatory_vars$leaves_alive, na.rm = TRUE),
    max(explanatory_vars$leaves_dead, na.rm = TRUE),
    max(explanatory_vars$flower_bud, na.rm = TRUE),
    max(explanatory_vars$flower_bloom, na.rm = TRUE),
    max(explanatory_vars$flower_dead, na.rm = TRUE),
    max(explanatory_vars$flower_seed, na.rm = TRUE),
    max(explanatory_vars$min_distance, na.rm = TRUE),
    max(explanatory_vars$num_neigh, na.rm = TRUE),
    max(explanatory_vars$edge_distance, na.rm = TRUE),
    max(explanatory_vars$humidity_base, na.rm = TRUE),
    max(explanatory_vars$humidity_top, na.rm = TRUE),
    max(explanatory_vars$micro_temp, na.rm = TRUE),
    max(explanatory_vars$min_surrounding_veg_height, na.rm = TRUE),
    max(explanatory_vars$max_surrounding_veg_height, na.rm = TRUE),
    max(explanatory_vars$average_wind_speed, na.rm = TRUE)
  )
)

# Print the summary statistics table
print(summary_stats)


##################### Create the same table but for each field #################

# Select the relevant columns from merged_data
explanatory_vars <- subset(merged_data, 
                           select = c(code, field.x, time.x,
                                      height, width, leaves_alive, 
                                      leaves_dead, flower_bud, flower_bloom, 
                                      flower_dead, flower_seed, 
                                      humidity_base, humidity_top, 
                                      micro_temp, 
                                      min_surrounding_veg_height, 
                                      max_surrounding_veg_height, 
                                      average_wind_speed, min_distance, 
                                      num_neigh, edge_distance))

# Remove duplicate rows
explanatory_vars <- unique(explanatory_vars)

# Split the data frame by 'field.x'
explanatory_vars_list <- split(explanatory_vars, explanatory_vars$field.x)

# Assign names to the data frames in the list
for (field_name in names(explanatory_vars_list)) {
  assign(paste0("explanatory_vars_", field_name), 
         explanatory_vars_list[[field_name]])
}

# Example to check if the data frames were created
print(ls(pattern = "explanatory_vars_"))

# List to store summary statistics dataframes
summary_stats_list <- list()

# Function to calculate summary statistics
calculate_summary_stats <- function(data) {
  summary_stats <- data.frame(
    Variable = c("Height", "Width", "Leaves Alive", "Leaves Dead", 
                 "Flower Bud", "Flower Bloom", "Flower Dead", "Flower Seed",
                 "Min Distance", "Num Neighbours", "Edge Distance",
                 "Humidity Base", "Humidity Top", "Micro Temp", 
                 "Min Surrounding Veg Height", "Max Surrounding Veg Height", 
                 "Average Wind Speed"),
    
    Mean = c(
      mean(data$height, na.rm = TRUE),
      mean(data$width, na.rm = TRUE),
      mean(data$leaves_alive, na.rm = TRUE),
      mean(data$leaves_dead, na.rm = TRUE),
      mean(data$flower_bud, na.rm = TRUE),
      mean(data$flower_bloom, na.rm = TRUE),
      mean(data$flower_dead, na.rm = TRUE),
      mean(data$flower_seed, na.rm = TRUE),
      mean(data$min_distance, na.rm = TRUE),
      mean(data$num_neigh, na.rm = TRUE),
      mean(data$edge_distance, na.rm = TRUE),
      mean(data$humidity_base, na.rm = TRUE),
      mean(data$humidity_top, na.rm = TRUE),
      mean(data$micro_temp, na.rm = TRUE),
      mean(data$min_surrounding_veg_height, na.rm = TRUE),
      mean(data$max_surrounding_veg_height, na.rm = TRUE),
      mean(data$average_wind_speed, na.rm = TRUE)
    ),
    
    Median = c(
      median(data$height, na.rm = TRUE),
      median(data$width, na.rm = TRUE),
      median(data$leaves_alive, na.rm = TRUE),
      median(data$leaves_dead, na.rm = TRUE),
      median(data$flower_bud, na.rm = TRUE),
      median(data$flower_bloom, na.rm = TRUE),
      median(data$flower_dead, na.rm = TRUE),
      median(data$flower_seed, na.rm = TRUE),
      median(data$min_distance, na.rm = TRUE),
      median(data$num_neigh, na.rm = TRUE),
      median(data$edge_distance, na.rm = TRUE),
      median(data$humidity_base, na.rm = TRUE),
      median(data$humidity_top, na.rm = TRUE),
      median(data$micro_temp, na.rm = TRUE),
      median(data$min_surrounding_veg_height, na.rm = TRUE),
      median(data$max_surrounding_veg_height, na.rm = TRUE),
      median(data$average_wind_speed, na.rm = TRUE)
    ),
    
    SD = c(
      sd(data$height, na.rm = TRUE),
      sd(data$width, na.rm = TRUE),
      sd(data$leaves_alive, na.rm = TRUE),
      sd(data$leaves_dead, na.rm = TRUE),
      sd(data$flower_bud, na.rm = TRUE),
      sd(data$flower_bloom, na.rm = TRUE),
      sd(data$flower_dead, na.rm = TRUE),
      sd(data$flower_seed, na.rm = TRUE),
      sd(data$min_distance, na.rm = TRUE),
      sd(data$num_neigh, na.rm = TRUE),
      sd(data$edge_distance, na.rm = TRUE),
      sd(data$humidity_base, na.rm = TRUE),
      sd(data$humidity_top, na.rm = TRUE),
      sd(data$micro_temp, na.rm = TRUE),
      sd(data$min_surrounding_veg_height, na.rm = TRUE),
      sd(data$max_surrounding_veg_height, na.rm = TRUE),
      sd(data$average_wind_speed, na.rm = TRUE)
    ),
    
    Min = c(
      min(data$height, na.rm = TRUE),
      min(data$width, na.rm = TRUE),
      min(data$leaves_alive, na.rm = TRUE),
      min(data$leaves_dead, na.rm = TRUE),
      min(data$flower_bud, na.rm = TRUE),
      min(data$flower_bloom, na.rm = TRUE),
      min(data$flower_dead, na.rm = TRUE),
      min(data$flower_seed, na.rm = TRUE),
      min(data$min_distance, na.rm = TRUE),
      min(data$num_neigh, na.rm = TRUE),
      min(data$edge_distance, na.rm = TRUE),
      min(data$humidity_base, na.rm = TRUE),
      min(data$humidity_top, na.rm = TRUE),
      min(data$micro_temp, na.rm = TRUE),
      min(data$min_surrounding_veg_height, na.rm = TRUE),
      min(data$max_surrounding_veg_height, na.rm = TRUE),
      min(data$average_wind_speed, na.rm = TRUE)
    ),
    
    Max = c(
      max(data$height, na.rm = TRUE),
      max(data$width, na.rm = TRUE),
      max(data$leaves_alive, na.rm = TRUE),
      max(data$leaves_dead, na.rm = TRUE),
      max(data$flower_bud, na.rm = TRUE),
      max(data$flower_bloom, na.rm = TRUE),
      max(data$flower_dead, na.rm = TRUE),
      max(data$flower_seed, na.rm = TRUE),
      max(data$min_distance, na.rm = TRUE),
      max(data$num_neigh, na.rm = TRUE),
      max(data$edge_distance, na.rm = TRUE),
      max(data$humidity_base, na.rm = TRUE),
      max(data$humidity_top, na.rm = TRUE),
      max(data$micro_temp, na.rm = TRUE),
      max(data$min_surrounding_veg_height, na.rm = TRUE),
      max(data$max_surrounding_veg_height, na.rm = TRUE),
      max(data$average_wind_speed, na.rm = TRUE)
    )
  )
  
  return(summary_stats)
}

# Loop through each subsetted dataframe 
for (field_name in names(explanatory_vars_list)) {
  summary_stats_list[[field_name]] <- 
    calculate_summary_stats(explanatory_vars_list[[field_name]])
}

# Check the first summary stats dataframe
summary_stats_list[["field1"]]






################################################################################
############# Abundance & richness vs. explanatory variable plots ##############
################################################################################

################################### Abundance ##################################

# Make dataframes
all_glm_df <- subset(merged_data_prez, select = c(code, date, field.x, time.x, 
                                                order, height, width_log,
                                                leaves_dead_log1p, 
                                                leaves_alive_log, no_leaves_log,
                                                flower_bud_log1p, 
                                                flower_bloom_log1p, 
                                                flower_seed_log1p, 
                                                flower_dead_log1p, f_type,
                                                grazed,
                                                max_surrounding_veg_height, 
                                                min_surrounding_veg_height, 
                                                micro_temp, humidity_base, 
                                                humidity_top, 
                                                average_wind_speed_log1p, 
                                                min_distance_log, 
                                                num_neigh_log,
                                                edge_distance,
                                                number_of_species, count, 
                                                shannon, simpson, evenness,
                                                count_log1p))

# Presence and absence dattaframe in selected orders dataframe
select_glm_abs <- subset(selected_order, select = c(code, date, field.x, time.x, 
                                          order, height, width_log,
                                          leaves_dead_log1p, 
                                          leaves_alive_log, no_leaves_log,
                                          flower_bud_log1p, 
                                          flower_bloom_log1p, 
                                          flower_seed_log1p, 
                                          flower_dead_log1p, f_type,
                                          grazed,
                                          max_surrounding_veg_height, 
                                          min_surrounding_veg_height, 
                                          micro_temp, humidity_base, 
                                          humidity_top, 
                                          average_wind_speed_log1p, 
                                          min_distance_log, 
                                          num_neigh_log,
                                          edge_distance,
                                          number_of_species, count, 
                                          shannon, simpson, evenness,
                                          count_log1p))


# Loop through column names and create scatter plots for COUNT with the 
# log-transformed explantatory variables for all orders
for (col in names(all_glm_df)) {
  if (col != "count_log1p") { # Exclude the 'count' column itself
    plot_title_count <- paste("Abundance Vs.", col)
    
    # Create count scatter plots
    p_count_all <- ggplot(all_glm_df, 
                          aes_string(x = col, y = "count_log1p")) +
      geom_point() +
      labs(x = col, y = "Abundance") +
      ggtitle(plot_title_count) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_count_all)
  }
}


# ... by selected orders
for (col in names(select_glm_abs)) {
  if (col != "count_log1p") {
    plot_title_count.2 <- paste("Abundance by order Vs.", col)
    
    # Create count scatter plots
    p_count_ord <- ggplot(select_glm_abs, 
                        aes_string(x = col, y = "count_log1p")) +
      facet_wrap(~order) +
      geom_point() +
      labs(x = col, y = "Abundance") +
      ggtitle(plot_title_count.2) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_count_ord)
  }
}

############################## Species richness ################################
# For all species 
for (col_rich in names(all_glm_df)) {
  if (col_rich != "number_of_species") { 
    plot_title_rich <- paste("Species richness Vs.", col_rich)
    
    # Create count scatter plots
    p_rich_all <- ggplot(all_glm_df, 
                          aes_string(x = col_rich, y = "number_of_species")) +
      geom_point() +
      labs(x = col_rich, y = "Abundance") +
      ggtitle(plot_title_rich) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_rich_all)
  }
}


## By order
for (col_rich in names(select_glm_abs)) {
  if (col_rich != "number_of_species") { 
    plot_title_rich.2 <- paste("Species richness by order Vs.", col_rich)
    

    p_rich_ord <- ggplot(select_glm_abs, 
                          aes_string(x = col_rich, y = "number_of_species")) +
      facet_wrap(~order) +
      geom_point() +
      labs(x = col_rich, y = "Species Richness") +
      ggtitle(plot_title_rich.2) +
      geom_smooth(method = "lm", se = TRUE)

    print(p_rich_ord)
  }
}



################################# Shannons #####################################
# For all species 
for (col_shan in names(all_glm_df)) {
  if (col_shan != "shannon") { 
    plot_title_shan <- paste("Shannons Index Vs.", col_shan)
    
    # Create count scatter plots
    p_shan_all <- ggplot(all_glm_df, 
                         aes_string(x = col_shan, y = "shannon")) +
      geom_point() +
      labs(x = col_shan, y = "Shannons Index") +
      ggtitle(plot_title_shan) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_shan_all)
  }
}


## By order
for (col_shan in names(select_glm_abs)) {
  if (col_shan != "shannon") { 
    plot_title_shan.2 <- paste("Shannons Index by order Vs.", col_shan)
    
    
    p_shan_ord <- ggplot(select_glm_abs, 
                         aes_string(x = col_shan, y = "shannon")) +
      facet_wrap(~order) +
      geom_point() +
      labs(x = col_shan, y = "Shannons Index") +
      ggtitle(plot_title_shan.2) +
      geom_smooth(method = "lm", se = TRUE)
    
    print(p_shan_ord)
  }
}

################################# Simpsons #####################################
# For all species 
for (col_simp in names(all_glm_df)) {
  if (col_simp != "simpson") { 
    plot_title_simp <- paste("Simpsons Index Vs.", col_simp)
    
    # Create count scatter plots
    p_simp_all <- ggplot(all_glm_df, 
                         aes_string(x = col_simp, y = "simpson")) +
      geom_point() +
      labs(x = col_simp, y = "Simpsons Index") +
      ggtitle(plot_title_simp) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_simp_all)
  }
}


## By order
for (col_simp in names(select_glm_abs)) {
  if (col_simp != "simpson") { 
    plot_title_simp.2 <- paste("Simpsons Index by order Vs.", col_simp)
    
    
    p_simp_ord <- ggplot(select_glm_abs, 
                         aes_string(x = col_simp, y = "simpson")) +
      facet_wrap(~order) +
      geom_point() +
      labs(x = col_simp, y = "Simpsons Index") +
      ggtitle(plot_title_simp.2) +
      geom_smooth(method = "lm", se = TRUE)
    
    print(p_simp_ord)
  }
}


################################# Evenness #####################################
# For all species 
for (col_even in names(all_glm_df)) {
  if (col_even != "evenness") { 
    plot_title_even <- paste("Evenness Vs.", col_even)
    
    # Create count scatter plots
    p_even_all <- ggplot(all_glm_df, 
                         aes_string(x = col_even, y = "evenness")) +
      geom_point() +
      labs(x = col_even, y = "Evenness") +
      ggtitle(plot_title_even) +
      geom_smooth(method = "lm", se = TRUE)
    
    # Print the scatter plot
    print(p_even_all)
  }
}


## By order
for (col_even in names(select_glm_abs)) {
  if (col_even != "evenness") { 
    plot_title_even.2 <- paste("Evenness by order Vs.", col_even)
    
    
    p_even_ord <- ggplot(select_glm_abs, 
                         aes_string(x = col_even, y = "evenness")) +
      facet_wrap(~order) +
      geom_point() +
      labs(x = col_even, y = "Evenness") +
      ggtitle(plot_title_even.2) +
      geom_smooth(method = "lm", se = TRUE)
    
    print(p_even_ord)
  }
}



################################################################################
############################ Plotting Diversity Indicies #######################
################################################################################
library(ggplot2)
library(gridExtra)
# Boxplot for Shannon Index 
shannon_plot <- ggplot(merged_data, aes(y = shannon)) +
  geom_boxplot(color = "black", fill = "white", width = 0.3) + 
  labs(title = "Shannon Index", y = "Shannon Index") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())  

# Boxplot for Evenness Index
evenness_plot <- ggplot(merged_data, aes(y = evenness)) +
  geom_boxplot(color = "black", fill = "white", width = 0.3) +  
  labs(title = "Evenness Index", y = "Evenness Index") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())  

# Boxplot for Simpson Index
simpson_plot <- ggplot(merged_data, aes(y = simpson)) +
  geom_boxplot(color = "black", fill = "white", width = 0.3) +  
  labs(title = "Simpson Index", y = "Simpson Index") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())  

# Arrange the three boxplots in the same window
grid.arrange(shannon_plot, evenness_plot, simpson_plot, ncol = 3)

# Replace Inf and NaN values to be able to compute
merged_data <- merged_data %>%
  mutate(
    shannon = ifelse(is.infinite(shannon) | is.nan(shannon), NA, shannon),
    evenness = ifelse(is.infinite(evenness) | is.nan(evenness), NA, evenness),
    simpson = ifelse(is.infinite(simpson) | is.nan(simpson), NA, simpson)
  )

# Descriptive statistics for shannon, simpson, and evenness
descriptive_stats <- merged_data %>%
  summarise(
    Mean = c(mean(shannon, na.rm = TRUE), 
             mean(evenness, na.rm = TRUE), 
             mean(simpson, na.rm = TRUE)),
    Median = c(median(shannon, na.rm = TRUE), 
               median(evenness, na.rm = TRUE), 
               median(simpson, na.rm = TRUE)),
    SD = c(sd(shannon, na.rm = TRUE), 
           sd(evenness, na.rm = TRUE), 
           sd(simpson, na.rm = TRUE)),
    Min = c(min(shannon, na.rm = TRUE), 
            min(evenness, na.rm = TRUE), 
            min(simpson, na.rm = TRUE)),
    Max = c(max(shannon, na.rm = TRUE), 
            max(evenness, na.rm = TRUE), 
            max(simpson, na.rm = TRUE))
  )

# Create a new data frame for readability
descriptive_stats <- merged_data %>%
  summarise(
    Mean = round(c(mean(shannon, na.rm = TRUE), 
                   mean(evenness, na.rm = TRUE), 
                   mean(simpson, na.rm = TRUE)), 3),
    Median = round(c(median(shannon, na.rm = TRUE), 
                     median(evenness, na.rm = TRUE), 
                     median(simpson, na.rm = TRUE)), 3),
    SD = round(c(sd(shannon, na.rm = TRUE), 
                 sd(evenness, na.rm = TRUE), 
                 sd(simpson, na.rm = TRUE)), 3),
    Min = round(c(min(shannon, na.rm = TRUE), 
                  min(evenness, na.rm = TRUE), 
                  min(simpson, na.rm = TRUE)), 3),
    Max = round(c(max(shannon, na.rm = TRUE), 
                  max(evenness, na.rm = TRUE), 
                  max(simpson, na.rm = TRUE)), 3)
  )

# Create a new data frame for readability
descriptive_stats_long <- tibble(
  Index = c("Shannon", "Evenness", "Simpson"),
  Mean = descriptive_stats$Mean,
  Median = descriptive_stats$Median,
  SD = descriptive_stats$SD,
  Min = descriptive_stats$Min,
  Max = descriptive_stats$Max
)

# Print descriptive statistics table
print(descriptive_stats_long)




################################################################################
######################### Univariate regression models #########################
################################################################################
# Fit GLM to individual explanatory variables (not to analyse more to just view)
library(MASS)
################################### Abundance ##################################
# For all species
# Create a list of variable names to put through the formula 
expl_vars_count <- names(all_glm_df)[!(names(all_glm_df) %in% 
                                         c("count", "code"))]  

for (var in expl_vars_count) {
  formula <- as.formula(paste("count ~", var))
  model <- glm(formula, data = all_glm_df, family = poisson(link = "log"))
  print(summary(model))
}

# selected orders
expl_vars_count_select <- names(select_glm_abs)[!(names(select_glm_abs) %in% 
                                         c("count", "code"))]  

for (var in expl_vars_count_select) {
  formula <- as.formula(paste("count ~", var))
  model <- glm.nb(formula, data = select_glm_abs)
  print(summary(model))
}


############################## Species richness ################################
# All species
expl_vars_rich <- names(all_glm_df)[!(names(all_glm_df) %in% 
                                       c("number_of_species", "code"))] 

for (var in expl_vars_rich) {
  formula <- as.formula(paste("number_of_species ~", var))
  model <- glm(formula, data = all_glm_df, family = poisson(link = "log"))
  print(summary(model))
}

# selected orders
expl_vars_rich_select <- names(select_glm_abs)[!(names(select_glm_abs) %in% 
                                             c("number_of_species", "code"))] 

for (var in expl_vars_rich_select) {
  formula <- as.formula(paste("number_of_species ~", var))
  model <- glm.nb(formula, data = select_glm_abs)
  print(summary(model))
}

############################## Shannons Index ##################################
# All species
expl_vars_shan <- names(all_glm_df)[!(names(all_glm_df) %in% 
                                        c("shannon", "code"))] 

for (var in expl_vars_shan) {
  formula <- as.formula(paste("shannon ~", var))
  model <- glm(formula, data = all_glm_df, family = gaussian(link = "identity"))
  print(summary(model))
}

# selected orders
expl_vars_shan_select <- names(select_glm_abs)[!(names(select_glm_abs) %in% 
                                                   c("shannon", "code"))] 

for (var in expl_vars_shan_select) {
  formula <- as.formula(paste("shannon ~", var))
  model <- glm(formula, data = select_glm_abs, 
               family = gaussian(link = "identity"))
  print(summary(model))
}

############################## Simpsons Index ##################################
# All species
expl_vars_simp <- names(all_glm_df)[!(names(all_glm_df) %in% 
                                        c("simpson", "code"))] 

for (var in expl_vars_simp) {
  formula <- as.formula(paste("simpson ~", var))
  model <- glm(formula, data = all_glm_df, family = gaussian(link = "identity"))
  print(summary(model))
}

# selected orders
expl_vars_simp_select <- names(select_glm_abs)[!(names(select_glm_abs) %in% 
                                                   c("simpson", "code"))] 

for (var in expl_vars_simp_select) {
  formula <- as.formula(paste("simpson ~", var))
  model <- glm(formula, data = select_glm_abs, 
               family = gaussian(link = "identity"))
  print(summary(model))
}


################################## Evenness ####################################
# All species
expl_vars_even <- names(all_glm_df)[!(names(all_glm_df) %in% 
                                        c("evenness", "code"))] 

for (var in expl_vars_even) {
  formula <- as.formula(paste("evenness ~", var))
  model <- glm(formula, data = all_glm_df, family = gaussian(link = "identity"))
  print(summary(model))
}

# selected orders
expl_vars_even_select <- names(select_glm_abs)[!(names(select_glm_abs) %in% 
                                                   c("evenness", "code"))] 

for (var in expl_vars_even_select) {
  formula <- as.formula(paste("evenness ~", var))
  model <- glm(formula, data = select_glm_abs, 
               family = gaussian(link = "identity"))
  print(summary(model))
}







################################################################################
################## Multivariate GLM with all variables #########################
################################################################################
library(car)
# Fit a GLM with all the variables to see the significance of each explanatory
# variable when all the others are accounted for 

# Fix Evenness calculation producing an Inf value (W_B_2.3)
all_glm_df$evenness[is.infinite(all_glm_df$evenness)] <- 0


################################### Abundance ##################################
# All species
full_model_count <- glm(count ~ . - code, data = all_glm_df, 
                            family = poisson(link = "log"))
summary(full_model_count)


# Selected orders
full_model_count_select <- glm.nb(count ~ . - code, data = select_glm_abs)
summary(full_model_count_select)


############################## Species richness ################################
# All species
full_model_rich <- glm(number_of_species ~ . - code, data = all_glm_df, 
                           family = poisson(link = "log"))
summary(full_model_rich)

# Selected orders
full_model_rich_select <- glm.nb(number_of_species ~ . - code, 
                                  data = select_glm_abs)
summary(full_model_rich_select)

################################ Shannon Index #################################
# All species
full_model_shan <- glm(shannon ~ . - code, data = all_glm_df, 
                       family = gaussian(link = "identity"))
summary(full_model_rich)

# Selected orders
full_model_shan_select <- glm(shannon ~ . - code,
                              data = select_glm_abs, 
                              family = gaussian(link = "identity"))
summary(full_model_shan_select)


################################ Simpson Index #################################
# All species
full_model_simp <- glm(simpson ~ . - code, data = all_glm_df, 
                       family = gaussian(link = "identity"))
summary(full_model_simp)

# Selected orders
full_model_simp_select <- glm(simpson ~ . - code,
                              data = select_glm_abs, 
                              family = gaussian(link = "identity"))
summary(full_model_simp_select)


################################## Evenness ####################################
# All species
full_model_even <- glm(evenness ~ . - code, data = all_glm_df, 
                       family = gaussian(link = "identity"))
summary(full_model_even)

# Selected orders
full_model_even_select <- glm(evenness ~ . - code,
                              data = select_glm_abs, 
                              family = gaussian(link = "identity"))
summary(full_model_even_select)





################################################################################
############################ GLM's and their graphs ############################
################################################################################
library(MASS)

################################### Abundance ##################################

# Create the model for all species                                                           
simple_model_count <- glm( count ~ 1 + 
                             order +
                             height +                     # Patch quality
                             width_log +
                             leaves_dead_log1p +
                             leaves_alive_log +
                             min_surrounding_veg_height + # Environmental 
                             max_surrounding_veg_height +
                             micro_temp +
                             humidity_base + 
                             average_wind_speed_log1p + 
                             min_distance_log +           # Connectivity 
                             edge_distance +
                             num_neigh_log + 
                             date +                       # Unmeasured variance
                             field.x +
                             time.x , 
                           family = poisson(link = "log"),
                           data = all_glm_df )
                          

summary(simple_model_count)              
vif(simple_model_count)

# Patch quality ONLY model for all species
quality_model_count_all <- glm(count ~ 1 + height +                     
                                 width_log +
                                 leaves_dead_log1p +
                                 leaves_alive_log +
                                 date +                       
                                 field.x +
                                 time.x , 
                               family = poisson(link = "log"), 
                               data = all_glm_df)

summary(quality_model_count_all)
vif(quality_model_count_all)


# Connectivity variables model for all species 
connectivity_model_count_all <- glm(count ~  1 + min_distance_log +           
                                      edge_distance +
                                      num_neigh_log +
                                      date +                       
                                      field.x +
                                      time.x, 
                                    family = poisson(link = "log"), 
                                    data = all_glm_df)

summary(connectivity_model_count_all)
vif(connectivity_model_count_all)


# Environmental-only model
environmental_model_count_all <- glm(count ~ 1 + min_surrounding_veg_height +                     
                                        max_surrounding_veg_height +
                                        micro_temp +
                                        humidity_base +
                                        average_wind_speed_log1p +
                                        date +                       
                                        field.x +
                                        time.x, 
                                      family = poisson(link = "log"), 
                                      data = all_glm_df)

summary(environmental_model_count_all)
vif(environmental_model_count_all)


# Create the model for selected orders                                                           
simple_model_count_select <- glm.nb( count ~ 1 + 
                                       order +
                                       height +                     
                                       width_log +
                                       leaves_dead_log1p +
                                       leaves_alive_log +
                                       min_surrounding_veg_height + 
                                       max_surrounding_veg_height +
                                       micro_temp + 
                                       humidity_base + 
                                       average_wind_speed_log1p + 
                                       min_distance_log +           
                                       edge_distance +
                                       num_neigh_log + 
                                       date +                       
                                       field.x +
                                       time.x, 
                                     data = select_glm_abs )
  

summary(simple_model_count_select)                        
vif(simple_model_count_select)

# Patch quality-only model for selected orders
patch_quality_count_select <- glm.nb(count ~ 1 + height +                     
                                       width_log +
                                       leaves_dead_log1p +
                                       leaves_alive_log +
                                       date +                       
                                       field.x +
                                       time.x, 
                                     data = select_glm_abs)

summary(patch_quality_count_select)
vif(patch_quality_count_select)



# Connectivity-only model for selected orders
connectivity_count_select <- glm.nb(count ~  1 + min_distance_log +           
                                      edge_distance +
                                      num_neigh_log +
                                      date +                       
                                      field.x +
                                      time.x, 
                                    data = select_glm_abs)

summary(connectivity_count_select)
vif(patch_quality_count_select)


# Environmental-only model for selected orders
environmental_count_select <- glm.nb(count ~ 1 + min_surrounding_veg_height +                     
                                       max_surrounding_veg_height +
                                       micro_temp + 
                                       humidity_base +
                                       average_wind_speed_log1p +
                                       date +                       
                                       field.x + 
                                       time.x, 
                                     data = select_glm_abs)

summary(environmental_count_select)
vif(environmental_count_select)


# Interaction between quality and connectivity (EXAMPLE HIGH COLLINEARITY)
int_count_select <- glm.nb(count ~ 1 + (min_distance_log +           
                                          edge_distance +
                                          num_neigh_log) * 
                             (height +                     
                                width_log +
                                leaves_dead_log1p +
                                leaves_alive_log) +
                             date +                       
                             field.x +
                             time.x, 
                           data = select_glm_abs)

summary(int_count_select)
vif(int_count_select)

# Create a model with order as an interaction                                                          
order_model_count <- glm.nb( count ~ 1 +
                             height : order +                     
                             width_log : order +
                             leaves_alive_log : order +
                             leaves_dead_log1p : order + 
                             min_surrounding_veg_height : order +                 
                             max_surrounding_veg_height : order +
                             humidity_base : order +
                             average_wind_speed_log1p : order + 
                             min_distance_log : order +           
                             edge_distance : order +
                             num_neigh_log : order +
                             field.x +
                             date +
                             time.x ,
                           data = select_glm_abs )


summary(order_model_count)                        
vif(order_model_count)

# Create a model with order as an interaction                                                          
order_model_count <- glm.nb( count ~ 1 +
                             height : order +                     
                             width_log : order +
                             leaves_alive_log : order +
                             leaves_dead_log1p : order + 
                             min_surrounding_veg_height : order +                 
                             max_surrounding_veg_height : order +
                             humidity_base : order +
                             average_wind_speed_log1p : order + 
                             min_distance_log : order +           
                             edge_distance : order +
                             num_neigh_log : order +
                             field.x +
                             date +
                             time.x ,
                           data = select_glm_abs )

# Patch quality-only model for selected orders with an interaction with order
patch_quality_count_select_ord <- glm.nb(glm.nb(count ~ 1 + height : order + 
                                                  width_log : order + 
                                                  leaves_dead_log1p : order + 
                                                  leaves_alive_log : order + 
                                                  date + 
                                                  field.x +
                                                  time.x, 
                                                  data = select_glm_abs))


summary(patch_quality_count_select_ord)
vif(patch_quality_count_select_ord)



# Connectivity-only model for selected orders with an interaction with order
connectivity_count_select_ord <- glm.nb(count ~ 1 + min_distance_log :order +           
                                      edge_distance : order +
                                      num_neigh_log : order +
                                      date +                       
                                      field.x + 
                                      time.x , 
                                    data = select_glm_abs)

summary(connectivity_count_select_ord)
vif(connectivity_count_select_ord)


# Environmental-only model for selected orders with an interaction with order
environmental_count_select_ord <- glm.nb(count ~ 1 + 
                                       min_surrounding_veg_height : order +                     
                                       max_surrounding_veg_height : order +
                                       micro_temp : order + 
                                       humidity_base : order +
                                       average_wind_speed_log1p : order +
                                       date +                       
                                       field.x + 
                                       time.x, 
                                     data = select_glm_abs)

summary(environmental_count_select_ord)
vif(environmental_count_select_ord)



################################### Richness ###################################

# Create the model for all species                                                           
simple_model_rich <- glm( number_of_species ~ 1 + 
                             order +
                             height +                     
                             width_log +
                             leaves_dead_log1p +
                             leaves_alive_log +
                             min_surrounding_veg_height + 
                             max_surrounding_veg_height +
                             micro_temp +
                             humidity_base + 
                             average_wind_speed_log1p + 
                             min_distance_log +           
                             edge_distance +
                             num_neigh_log + 
                             date +                       
                             field.x +
                             time.x, 
                           family = poisson(link = "log"),
                           data = all_glm_df )


summary(simple_model_rich) 
vif(simple_model_rich)

# Patch quality ONLY model for all species
quality_model_rich_all <- glm(number_of_species ~ 1 + height +                     
                                 width_log +
                                 leaves_dead_log1p +
                                 leaves_alive_log +
                                 date +                       
                                 field.x +
                                 time.x, 
                               family = poisson(link = "log"), 
                               data = all_glm_df)

summary(quality_model_rich_all)
vif(quality_model_rich_all)


# Connectivity variables model for all species 
connectivity_model_rich_all <- glm(number_of_species ~ 1 + min_distance_log +           
                                      edge_distance +
                                      num_neigh_log +
                                      date +                       
                                      field.x +
                                      time.x, 
                                    family = poisson(link = "log"), 
                                    data = all_glm_df)

summary(connectivity_model_rich_all)
vif(connectivity_model_rich_all)


# Environmental-only model
environmental_model_rich_all <- glm(number_of_species ~ 1 + 
                                       min_surrounding_veg_height +                     
                                       max_surrounding_veg_height +
                                       micro_temp + 
                                       humidity_base + 
                                       average_wind_speed_log1p +
                                       date +                       
                                       field.x +
                                       time.x, 
                                     family = poisson(link = "log"), 
                                     data = all_glm_df)

summary(environmental_model_rich_all)
vif(environmental_model_rich_all)


# Create the model for selected orders                                                           
simple_model_rich_select <- glm.nb(number_of_species ~ 1 +  
                                       order +
                                       height +                     
                                       width_log +
                                       leaves_dead_log1p +
                                       leaves_alive_log +
                                       min_surrounding_veg_height + 
                                       max_surrounding_veg_height +
                                       micro_temp + 
                                       humidity_base +
                                       average_wind_speed_log1p + 
                                       min_distance_log +           
                                       edge_distance +
                                       num_neigh_log + 
                                       date +                       
                                       field.x +
                                       time.x, 
                                     data = select_glm_abs )


summary(simple_model_rich_select)                        
vif(simple_model_rich_select)


# Patch quality-only model for selected orders
patch_quality_rich_select <- glm.nb(count ~ 1 + height +                     
                                       width_log +
                                       leaves_dead_log1p +
                                       leaves_alive_log +
                                       date +                       
                                       field.x + 
                                       time.x, 
                                     data = select_glm_abs)

summary(patch_quality_rich_select)
vif(patch_quality_rich_select)                        


# Connectivity-only model for selected orders
connectivity_rich_select <- glm.nb(count ~ 1 + min_distance_log +           
                                      edge_distance +
                                      num_neigh_log +
                                      date +                       
                                      field.x +
                                      time.x, 
                                    data = select_glm_abs)

summary(connectivity_rich_select)
vif(connectivity_rich_select)                        


# Environmental-only model for selected orders
environmental_rich_select <- glm.nb(count ~ 1 + min_surrounding_veg_height +                     
                                       max_surrounding_veg_height +
                                       micro_temp + 
                                       humidity_base +
                                       average_wind_speed_log1p +
                                       date +                       
                                       field.x +
                                       time.x, 
                                     data = select_glm_abs)

summary(environmental_rich_select)
vif(environmental_rich_select)                        



# Create a model with order as an interaction                                                          
order_model_rich <- glm.nb( number_of_species ~ 1 +
                               height :order  +                     
                                  width_log : order +
                                  leaves_dead_log1p : order +
                                  leaves_alive_log : order +
                                  min_surrounding_veg_height : order + 
                                  max_surrounding_veg_height : order+
                                  micro_temp : order + 
                                  humidity_base : order +
                                  average_wind_speed_log1p : order + 
                                  min_distance_log :order  +           
                                  edge_distance : order +
                                  num_neigh_log : order + 
                                  date +                       
                                  field.x +
                                  time.x , 
                             data = select_glm_abs )


summary(order_model_rich)                        
vif(order_model_rich)                        


# Patch quality-only model for selected orders with an interaction with order
patch_quality_rich_select_ord <- glm.nb(number_of_species ~ 1 + height : order +                     
                                                    width_log : order +
                                                    leaves_dead_log1p : order +
                                                    leaves_alive_log : order +
                                                    date +                       
                                                    field.x +
                                                    time.x, 
                                         data = select_glm_abs)

summary(patch_quality_rich_select_ord)
vif(patch_quality_rich_select_ord)                        



# Connectivity-only model for selected orders with an interaction with order
connectivity_rich_select_ord <- glm.nb(number_of_species ~ 1 + 
                                                 min_distance_log : order +           
                                                  edge_distance : order +
                                                  num_neigh_log : order +
                                                  date : order +                       
                                                  field.x : order + 
                                                  time.x, 
                                       data = select_glm_abs)

summary(connectivity_rich_select_ord)
vif(connectivity_rich_select_ord)                        


# Environmental-only model for selected orders with an interaction with order
environmental_rich_select_ord <- glm.nb(number_of_species ~ 1 + 
                                          min_surrounding_veg_height : order +                     
                                          max_surrounding_veg_height : order +
                                          micro_temp : order + 
                                          humidity_base : order +
                                          average_wind_speed_log1p : order +
                                          date +                       
                                          field.x + 
                                          time.x, 
                                         data = select_glm_abs)

summary(environmental_rich_select_ord)
vif(environmental_rich_select_ord)                        


################################### Shannon ####################################

# Create the model for all species                                                           
simple_model_shan <- glm( shannon ~ 1 + 
                            order +
                            height +                     
                            width_log +
                            leaves_dead_log1p +
                            leaves_alive_log +
                            min_surrounding_veg_height + 
                            max_surrounding_veg_height +
                            micro_temp + 
                            humidity_base + 
                            average_wind_speed_log1p + 
                            min_distance_log +           
                            edge_distance +
                            num_neigh_log + 
                            date +                       
                            field.x + 
                            time.x, 
                          family = gaussian(link = "identity"),
                          data = all_glm_df )


summary(simple_model_shan)                        
vif(simple_model_shan)                        

# Patch quality ONLY model for all species
quality_model_shan_all <- glm(shannon ~ 1 + height +                     
                                width_log +
                                leaves_dead_log1p +
                                leaves_alive_log +
                                date +                       
                                field.x + 
                                time.x , 
                              family = gaussian(link = "identity"), 
                              data = all_glm_df)

summary(quality_model_shan_all)
vif(quality_model_shan_all)                        


# Connectivity variables model for all species 
connectivity_model_shan_all <- glm(shannon ~ 1 + min_distance_log +           
                                     edge_distance +
                                     num_neigh_log +
                                     date +                       
                                     field.x + 
                                     time.x, 
                                   family = gaussian(link = "identity"), 
                                   data = all_glm_df)

summary(connectivity_model_shan_all)
vif(connectivity_model_shan_all)                        


# Environmental-only model for all species
environmental_model_shan_all <- glm(shannon ~ 1 + 
                                      min_surrounding_veg_height +                     
                                      max_surrounding_veg_height +
                                      micro_temp + 
                                      humidity_base + 
                                      average_wind_speed_log1p +
                                      date +                       
                                      field.x + 
                                      time.x , 
                                    family = gaussian(link = "identity"), 
                                    data = all_glm_df)

summary(environmental_model_shan_all)
vif(environmental_model_shan_all)                        


# Create the model for selected orders                                                           
simple_model_shan_select <- glm( shannon ~ 1 + 
                                      order +
                                      height +                     
                                      width_log +
                                      leaves_dead_log1p +
                                      leaves_alive_log +
                                      min_surrounding_veg_height + 
                                      max_surrounding_veg_height +
                                      micro_temp + 
                                      humidity_base + 
                                      average_wind_speed_log1p + 
                                      min_distance_log +           
                                      edge_distance +
                                      num_neigh_log + 
                                      date +                       
                                      field.x + 
                                      time.x , 
                                    family = gaussian(link = "identity"),
                                    data = select_glm_abs )


summary(simple_model_shan_select)                        
vif(simple_model_shan_select)                        

# Patch quality-only model for selected orders
patch_quality_shan_select <- glm(shannon ~ 1 + height +                     
                                      width_log +
                                      leaves_dead_log1p +
                                      leaves_alive_log +
                                      date +                       
                                      field.x +
                                      time.x , 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs)

summary(patch_quality_shan_select)
vif(patch_quality_shan_select)                        



# Connectivity-only model for selected orders
connectivity_shan_select <- glm(shannon ~ 1 + min_distance_log +           
                                     edge_distance +
                                     num_neigh_log +
                                     date +                       
                                     field.x + 
                                     time.x , 
                                   family = gaussian(link = "identity"),
                                   data = select_glm_abs)

summary(connectivity_shan_select)
vif(connectivity_shan_select)                        


# Environmental-only model for selected orders
environmental_shan_select <- glm(shannon ~ 1 + min_surrounding_veg_height +                     
                                      max_surrounding_veg_height +
                                      micro_temp + 
                                      humidity_base +
                                      average_wind_speed_log1p +
                                      date +                       
                                      field.x + 
                                      time.x ,
                                    family = gaussian(link = "identity"),
                                    data = select_glm_abs)

summary(environmental_shan_select)
vif(environmental_shan_select)                        


# Create a model with order as an interaction                                                          
order_model_shan <- glm( shannon ~ 1 +
                              height : order +                     
                                 width_log : order +
                                 leaves_dead_log1p : order +
                                 leaves_alive_log : order +
                                 min_surrounding_veg_height : order + 
                                 max_surrounding_veg_height : order +
                                 micro_temp : order + 
                                 humidity_base : order + 
                                 average_wind_speed_log1p : order + 
                                 min_distance_log : order +           
                                 edge_distance : order +
                                 num_neigh_log : order + 
                                 date +                       
                                 field.x +
                                 time.x , 
                            family = gaussian(link = "identity") ,
                            data = select_glm_abs )


summary(order_model_shan)     
vif(order_model_shan)                        


# Patch quality-only model for selected orders with an interaction with order
patch_quality_shan_select_ord <- glm(shannon ~ 1 + 
                                       height:order + 
                                       width_log:order + 
                                       leaves_dead_log1p:order + 
                                       leaves_alive_log:order + 
                                       date + 
                                       field.x +
                                       time.x , 
                                     family = gaussian(link = "identity"), 
                                     data = select_glm_abs)


summary(patch_quality_shan_select_ord)
vif(patch_quality_shan_select_ord)                        



# Connectivity-only model for selected orders with an interaction with order
connectivity_shan_select_ord <- glm(shannon ~ 1 + 
                                        min_distance_log : order +           
                                           edge_distance : order +
                                           num_neigh_log : order +
                                           date +                       
                                           field.x +
                                           time.x, 
                                      family = gaussian(link = "identity") ,
                                      data = select_glm_abs)

summary(connectivity_shan_select_ord)
vif(connectivity_shan_select_ord)                        


# Environmental-only model for selected orders with an interaction with order
environmental_shan_select_ord <- glm(shannon ~ 1 + 
                                          min_surrounding_veg_height : order +                     
                                            max_surrounding_veg_height : order +
                                             micro_temp : order + 
                                             humidity_base : order +
                                             average_wind_speed_log1p : order +
                                             date +                      
                                             field.x +
                                             time.x ,
                                        family = gaussian(link = "identity") ,
                                        data = select_glm_abs)

summary(environmental_shan_select_ord)
vif(environmental_shan_select_ord)                        


################################### Simpson ####################################

# Create the model for all species                                                           
simple_model_simp <- glm( simpson ~ 1 + 
                            order +
                            height +                     
                            width_log +
                            leaves_dead_log1p +
                            leaves_alive_log +
                            min_surrounding_veg_height + 
                            max_surrounding_veg_height +
                            micro_temp + 
                            humidity_base +
                            average_wind_speed_log1p + 
                            min_distance_log +           
                            edge_distance +
                            num_neigh_log + 
                            date +                       
                            field.x + 
                            time.x, 
                          family = gaussian(link = "identity"),
                          data = all_glm_df )


summary(simple_model_simp)                        
vif(simple_model_simp)                        


# Patch quality ONLY model for all species
quality_model_simp_all <- glm(simpson ~ 1 + height +                     
                                width_log +
                                leaves_dead_log1p +
                                leaves_alive_log +
                                date +                       
                                field.x + 
                                time.x, 
                              family = gaussian(link = "identity"), 
                              data = all_glm_df)

summary(quality_model_simp_all)
vif(quality_model_simp_all)                        


# Connectivity variables model for all species 
connectivity_model_simp_all <- glm(simpson ~ 1 + min_distance_log +           
                                     edge_distance +
                                     num_neigh_log +
                                     date +                       
                                     field.x + 
                                     time.x, 
                                   family = gaussian(link = "identity"), 
                                   data = all_glm_df)

summary(connectivity_model_simp_all)
vif(connectivity_model_simp_all)                        


# Environmental-only model for all species
environmental_model_simp_all <- glm(simpson ~ 1 + 
                                      min_surrounding_veg_height +                     
                                      max_surrounding_veg_height +
                                      micro_temp +
                                      humidity_base +
                                      average_wind_speed_log1p +
                                      date +                       
                                      field.x + 
                                      time.x, 
                                    family = gaussian(link = "identity"), 
                                    data = all_glm_df)

summary(environmental_model_simp_all)
vif(environmental_model_simp_all)                        


# Create the model for selected orders                                                           
simple_model_simp_select <- glm( simpson ~ 1 + 
                                   order +
                                   height +                     
                                   width_log +
                                   leaves_dead_log1p +
                                   leaves_alive_log +
                                   min_surrounding_veg_height + 
                                   max_surrounding_veg_height +
                                   micro_temp +
                                   humidity_base +
                                   average_wind_speed_log1p + 
                                   min_distance_log +           
                                   edge_distance +
                                   num_neigh_log + 
                                   date +                       
                                   field.x +
                                   time.x, 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs )


summary(simple_model_simp_select)                        
vif(simple_model_simp_select)                        


# Patch quality-only model for selected orders
patch_quality_simp_select <- glm(simpson ~ 1 + height +                     
                                   width_log +
                                   leaves_dead_log1p +
                                   leaves_alive_log +
                                   date +                       
                                   field.x +
                                   time.x, 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs)

summary(patch_quality_simp_select)
vif(patch_quality_simp_select)                        



# Connectivity-only model for selected orders
connectivity_simp_select <- glm(simpson ~ 1 + min_distance_log +           
                                  edge_distance +
                                  num_neigh_log +
                                  date +                       
                                  field.x + 
                                  time.x, 
                                family = gaussian(link = "identity"),
                                data = select_glm_abs)

summary(connectivity_simp_select)
vif(connectivity_simp_select)                        


# Environmental-only model for selected orders
environmental_simp_select <- glm(simpson ~ 1 + min_surrounding_veg_height +                     
                                   max_surrounding_veg_height +
                                   micro_temp + 
                                   humidity_base + 
                                   average_wind_speed_log1p +
                                   date +                       
                                   field.x + 
                                   time.x,
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs)

summary(environmental_simp_select)
vif(environmental_simp_select)                        


# Create a model with order as an interaction                                                          
order_model_simp <- glm( simpson ~ 1 +
                           height : order +                     
                              width_log : order +
                              leaves_dead_log1p : order +
                              leaves_alive_log : order +
                              min_surrounding_veg_height : order + 
                              max_surrounding_veg_height : order +
                              micro_temp : order + 
                              humidity_base : order +
                              average_wind_speed_log1p : order + 
                              min_distance_log : order +           
                              edge_distance : order +
                              num_neigh_log : order + 
                              date +                      
                              field.x + 
                              time.x , 
                         family = gaussian(link = "identity") ,
                         data = select_glm_abs )


summary(order_model_simp)    
vif(order_model_simp)                        


# Patch quality-only model for selected orders with an interaction with order
patch_quality_simp_select_ord <- glm(simpson ~ 1 + 
                                       height : order + 
                                       width_log : order +
                                       leaves_dead_log1p : order +
                                       leaves_alive_log : order +
                                       date +                       
                                       field.x +
                                       time.x,
                                     family = gaussian(link = "identity") ,
                                     data = select_glm_abs)

summary(patch_quality_simp_select_ord)
vif(patch_quality_simp_select_ord)                        



# Connectivity-only model for selected orders with an interaction with order
connectivity_simp_select_ord <- glm(simpson ~ 1 + 
                                            min_distance_log : order +           
                                               edge_distance : order +
                                               num_neigh_log : order +
                                               date +                       
                                               field.x + 
                                               time.x, 
                                          family = gaussian(link = "identity") ,
                                          data = select_glm_abs)

summary(connectivity_simp_select_ord)
vif(connectivity_simp_select_ord)                        


# Environmental-only model for selected orders with an interaction with order
environmental_simp_select_ord <- glm(simpson ~ 1 + 
                                          min_surrounding_veg_height : order +                     
                                            max_surrounding_veg_height : order +
                                             micro_temp : order +
                                             humidity_base : order + 
                                             average_wind_speed_log1p : order +
                                             date +                       
                                             field.x + 
                                             time.x,
                                        family = gaussian(link = "identity") ,
                                        data = select_glm_abs)

summary(environmental_simp_select_ord)
vif(environmental_simp_select_ord)                        


################################## Evenness ####################################

# Create the model for all species                                                           
simple_model_even <- glm( evenness ~ 1 + 
                            order +
                            height +                     
                            width_log +
                            leaves_dead_log1p +
                            leaves_alive_log +
                            min_surrounding_veg_height + 
                            max_surrounding_veg_height +
                            micro_temp + 
                            humidity_base +
                            average_wind_speed_log1p + 
                            min_distance_log +           
                            edge_distance +
                            num_neigh_log + 
                            date +                       
                            field.x + 
                            time.x, 
                          family = gaussian(link = "identity"),
                          data = all_glm_df )


summary(simple_model_even)                        
vif(simple_model_even)                        

# Patch quality ONLY model for all species
quality_model_even_all <- glm(evenness ~ 1 + height +                     
                                width_log +
                                leaves_dead_log1p +
                                leaves_alive_log +
                                date +                       
                                field.x + 
                                time.x, 
                              family = gaussian(link = "identity"), 
                              data = all_glm_df)

summary(quality_model_even_all)
vif(quality_model_even_all)                        


# Connectivity variables model for all species 
connectivity_model_even_all <- glm(evenness ~ 1 + min_distance_log +           
                                     edge_distance +
                                     num_neigh_log +
                                     date +                       
                                     field.x +
                                     time.x, 
                                   family = gaussian(link = "identity"), 
                                   data = all_glm_df)

summary(connectivity_model_even_all)
vif(connectivity_model_even_all)                        


# Environmental-only model for all species
environmental_model_even_all <- glm(evenness ~ 1 + 
                                      min_surrounding_veg_height +                     
                                      max_surrounding_veg_height +
                                      micro_temp + 
                                      humidity_base + 
                                      average_wind_speed_log1p +
                                      date +                       
                                      field.x + 
                                      time.x, 
                                    family = gaussian(link = "identity"), 
                                    data = all_glm_df)

summary(environmental_model_even_all)
vif(environmental_model_even_all)                        


# Create the model for selected orders                                                           
simple_model_even_select <- glm( evenness ~ 1 + 
                                   order +
                                   height +                     
                                   width_log +
                                   leaves_dead_log1p +
                                   leaves_alive_log +
                                   min_surrounding_veg_height + 
                                   max_surrounding_veg_height +
                                   micro_temp + 
                                   humidity_base + 
                                   average_wind_speed_log1p + 
                                   min_distance_log +           
                                   edge_distance +
                                   num_neigh_log + 
                                   date +                       
                                   field.x +
                                   time.x, 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs )


summary(simple_model_even_select)  
vif(simple_model_even_select)                        


# Patch quality-only model for selected orders
patch_quality_even_select <- glm(evenness ~ 1 + height +                     
                                   width_log +
                                   leaves_dead_log1p +
                                   leaves_alive_log +
                                   date +                       
                                   field.x +
                                   time.x, 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs)

summary(patch_quality_even_select)
vif(patch_quality_even_select)                        



# Connectivity-only model for selected orders
connectivity_even_select <- glm(evenness ~ 1 + min_distance_log +           
                                  edge_distance +
                                  num_neigh_log +
                                  date +                       
                                  field.x + 
                                  time.x, 
                                family = gaussian(link = "identity"),
                                data = select_glm_abs)

summary(connectivity_even_select)
vif(connectivity_even_select)                        



# Environmental-only model for selected orders
environmental_even_select <- glm(evenness ~ 1 + min_surrounding_veg_height +                     
                                   max_surrounding_veg_height +
                                   micro_temp + 
                                   humidity_base +
                                   average_wind_speed_log1p +
                                   date +                       
                                   field.x + 
                                   time.x,
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs)

summary(environmental_even_select)
vif(environmental_even_select)                        


# Create a model with order as an interaction                                                          
order_model_even <- glm( evenness ~ 1 +
                           height : order +                     
                              width_log : order +
                              leaves_dead_log1p : order +
                              leaves_alive_log : order +
                              min_surrounding_veg_height : order + 
                              max_surrounding_veg_height : order +
                              micro_temp : order + 
                              humidity_base : order +
                              average_wind_speed_log1p : order + 
                              min_distance_log : order +           
                              edge_distance : order +
                              num_neigh_log : order + 
                              date +                       
                              field.x +
                              time.x , 
                         family = gaussian(link = "identity") ,
                         data = select_glm_abs )


summary(order_model_simp)    
vif(order_model_simp)                        


# Patch quality-only model for selected orders with an interaction with order
patch_quality_even_select_ord <- glm(evenness ~ 1 + 
                                       height : order + 
                                       width_log : order +
                                       leaves_dead_log1p : order +
                                       leaves_alive_log : order +
                                       date +                       
                                       field.x + 
                                       time.x , 
                                     family = gaussian(link = "identity") ,
                                     data = select_glm_abs)

summary(patch_quality_even_select_ord)
vif(patch_quality_even_select_ord)                        



# Connectivity-only model for selected orders with an interaction with order
connectivity_even_select_ord <- glm(evenness ~ 1 + 
                                            min_distance_log : order +           
                                               edge_distance : order +
                                               num_neigh_log : order +
                                               date +                       
                                               field.x + 
                                               time.x, 
                                          family = gaussian(link = "identity") ,
                                          data = select_glm_abs)

summary(connectivity_even_select_ord)
vif(connectivity_even_select_ord)                        


# Environmental-only model for selected orders with an interaction with order
environmental_even_select_ord <- glm(evenness ~ 1 + 
                                          min_surrounding_veg_height : order +                     
                                            max_surrounding_veg_height : order +
                                             micro_temp : order + 
                                             humidity_base : order +
                                             average_wind_speed_log1p : order +
                                             date +                       
                                             field.x +
                                             time.x,
                                        family = gaussian(link = "identity") ,
                                        data = select_glm_abs)

summary(environmental_even_select_ord)
vif(environmental_even_select_ord)                        


##################################### AIC Testing ##############################

# AIC comparisons for all specified models
response_count_AIC <- AIC(
  full_model_count,
  full_model_count_select,
  simple_model_count,
  quality_model_count_all,
  connectivity_model_count_all,
  environmental_model_count_all,
  simple_model_count_select,
  patch_quality_count_select,
  connectivity_count_select,
  environmental_count_select,
  order_model_count,
  patch_quality_count_select_ord,
  connectivity_count_select_ord,
  environmental_count_select_ord
)

print(response_count_AIC)

# Species Richness
response_rich_AIC <- AIC(
  full_model_rich,
  full_model_rich_select,
  simple_model_rich,
  quality_model_rich_all,
  connectivity_model_rich_all,
  environmental_model_rich_all,
  simple_model_rich_select,
  patch_quality_rich_select,
  connectivity_rich_select,
  environmental_rich_select,
  order_model_rich,
  patch_quality_rich_select_ord,
  connectivity_rich_select_ord,
  environmental_rich_select_ord
)

# Print the AIC values
print(response_rich_AIC)

#  Shannon Index
response_shan_AIC <- AIC(
  full_model_shan,
  full_model_shan_select,
  simple_model_shan,
  quality_model_shan_all,
  connectivity_model_shan_all,
  environmental_model_shan_all,
  simple_model_shan_select,
  patch_quality_shan_select,
  connectivity_shan_select,
  environmental_shan_select,
  order_model_shan,
  patch_quality_shan_select_ord,
  connectivity_shan_select_ord,
  environmental_shan_select_ord
)

# Print the AIC values
print(response_shan_AIC)


# Simpson Index
response_simp_AIC <- AIC(
  full_model_simp,
  full_model_simp_select,
  simple_model_simp,
  quality_model_simp_all,
  connectivity_model_simp_all,
  environmental_model_simp_all,
  simple_model_simp_select,
  patch_quality_simp_select,
  connectivity_simp_select,
  environmental_simp_select,
  order_model_simp,
  patch_quality_simp_select_ord,
  connectivity_simp_select_ord,
  environmental_simp_select_ord
)

# Print the AIC values
print(response_simp_AIC)


# AIC comparisons for models related to Evenness Index
response_even_AIC <- AIC(
  full_model_even,
  full_model_even_select,
  simple_model_even,
  quality_model_even_all,
  connectivity_model_even_all,
  environmental_model_even_all,
  simple_model_even_select,
  patch_quality_even_select,
  connectivity_even_select,
  environmental_even_select,
  order_model_even,
  patch_quality_even_select_ord,
  connectivity_even_select_ord,
  environmental_even_select_ord
)


######################## Creating GLM summary tables ###########################
# Adjust as neccessary
# Extract the estimates, standard errors, z-values, and p-values
# Adjust if using a model that generates t-values 
estimates <- coef(simple_model_simp)
std_errors <- summary(simple_model_simp)$coefficients[, "Std. Error"]
t_values <- summary(simple_model_simp)$coefficients[, "t value"]
p_values <- summary(simple_model_simp)$coefficients[, "Pr(>|t|)"]

# Create a custom data frame with the new Estimate column
custom_table <- data.frame(
  Estimate = estimates,
  Std_Error = std_errors,
  t_value = t_values,
  Pr_z = p_values
)







################################# Plot GLMs ####################################
library(ggeffects)
library(ggplot2)
library(patchwork)  

# Create a model for plotting (without as it cuases issues)
order_model_count_plot <- glm.nb( count ~ 1 +
                               height : order +                     
                               width_log : order +
                               leaves_alive_log : order +
                               leaves_dead_log1p : order + 
                               min_surrounding_veg_height : order +                 
                               max_surrounding_veg_height : order +
                               humidity_base : order +
                               average_wind_speed_log1p : order + 
                               min_distance_log : order +           
                               edge_distance : order +
                               num_neigh_log : order +
                               field.x +
                               time.x ,
                             data = select_glm_abs )


# Create the predicted data to plot from
p1_abun <- ggpredict( order_model_count_plot, terms = c("height", "order") )
p2_abun <- ggpredict( order_model_count_plot, terms = c("width_log", "order") )
p3_abun <- ggpredict( order_model_count_plot, terms = c("num_neigh_log",
                                                        "order") )
p4_abun <- ggpredict( order_model_count_plot, terms = c("edge_distance", 
                                                        "order") )

# Define a function to create the plot with different titles
create_plot <- function(data, title) {
  ggplot(data = data, aes(x = x, y = predicted, col = group, fill = group)) +
    geom_line() +
    xlab("Predictor variables") +
    ylab("Count") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
    facet_wrap(~group) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
    ) +
    ggtitle(title)
}

# Create separate plots with different titles
g1_abun <- create_plot(p1_abun, "height")
g2_abun <- create_plot(p2_abun, "width_log")
g3_abun <- create_plot(p3_abun, "num_neigh_log")
g4_abun <- create_plot(p4_abun, "edge_distance")

# Display the plots
print(g1_abun)
print(g2_abun)
print(g3_abun)
print(g4_abun)


# Create new models for plotting
simple_model_shan_plot <- glm( shannon ~ 1 + 
                            order +
                            height +                     
                            width_log +
                            leaves_dead_log1p +
                            leaves_alive_log +
                            min_surrounding_veg_height + 
                            max_surrounding_veg_height +
                            micro_temp + 
                            humidity_base + 
                            average_wind_speed_log1p + 
                            min_distance_log +           
                            edge_distance +
                            num_neigh_log + 
                            field.x + 
                            time.x, 
                          family = gaussian(link = "identity"),
                          data = all_glm_df )

simple_model_even_plot <- glm( evenness ~ 1 + 
                                   order +
                                   height +                     
                                   width_log +
                                   leaves_dead_log1p +
                                   leaves_alive_log +
                                   min_surrounding_veg_height + 
                                   max_surrounding_veg_height +
                                   micro_temp + 
                                   humidity_base + 
                                   average_wind_speed_log1p + 
                                   min_distance_log +           
                                   edge_distance +
                                   num_neigh_log + 
                                   field.x +
                                   time.x, 
                                 family = gaussian(link = "identity"),
                                 data = select_glm_abs )

simple_model_simp_plot <- glm( simpson ~ 1 + 
                            order +
                            height +                     
                            width_log +
                            leaves_dead_log1p +
                            leaves_alive_log +
                            min_surrounding_veg_height + 
                            max_surrounding_veg_height +
                            micro_temp + 
                            humidity_base +
                            average_wind_speed_log1p + 
                            min_distance_log +           
                            edge_distance +
                            num_neigh_log + 
                            field.x + 
                            time.x, 
                          family = gaussian(link = "identity"),
                          data = all_glm_df )


# Create predicted data for num_neigh_log
p1_shan <- ggpredict(simple_model_shan_plot, terms = c("num_neigh_log"))
p2_even <- ggpredict(simple_model_even_plot, terms = c("num_neigh_log"))
p3_simp <- ggpredict(simple_model_simp_plot, terms = c("num_neigh_log"))

# Create predicted data for edge_distance
p4_shan <- ggpredict(simple_model_shan_plot, terms = c("edge_distance"))
p5_even <- ggpredict(simple_model_even_plot, terms = c("edge_distance"))
p6_simp <- ggpredict(simple_model_simp_plot, terms = c("edge_distance"))

# Function to create a plot for each response variable
create_plot <- function(data, y_label, title) {
  ggplot(data = data, aes(x = x, y = predicted)) +
    geom_line() +
    xlab("Predictor") +
    ylab(y_label) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
    ) +
    ggtitle(title)
}

# Plot for num_neigh_log
plot_shan_num <- create_plot(p1_shan, "Shannon Index", 
                             "Shannon vs num_neigh_log")
plot_even_num <- create_plot(p2_even, "Evenness Index", 
                             "Evenness vs num_neigh_log")
plot_simp_num <- create_plot(p3_simp, "Simpson Index",
                             "Simpson vs num_neigh_log")

# Combine plots for num_neigh_log in a row
num_neigh_plots <- plot_shan_num + plot_even_num + plot_simp_num

# Plot for edge_distance
plot_shan_edge <- create_plot(p4_shan, "Shannon Index", 
                              "Shannon vs edge_distance")
plot_even_edge <- create_plot(p5_even, "Evenness Index", 
                              "Evenness vs edge_distance")
plot_simp_edge <- create_plot(p6_simp, "Simpson Index", 
                              "Simpson vs edge_distance")

# Combine plots for edge_distance in a row
edge_distance_plots <- plot_shan_edge + plot_even_edge + plot_simp_edge

# Display the plots
num_neigh_plots / edge_distance_plots




################################################################################
################################### CCA ########################################
################################################################################
library(vegan)
library(tibble)
library(dplyr)
library(tidyr)
######################## Create the functional trait dataframe #################
traits_df <- species_list_prez[, c("species_name", "feeding", "niche", 
                                   "current_dispersal", "voltinism", 
                                   "size", "stage")]

# Remove duplicate species names
traits_df <- traits_df %>%
       distinct(species_name, .keep_all = TRUE)

# Ensure the species column is set as row names
trait_df <- traits_df %>% column_to_rownames(var = "species_name")


# Remove rows with any NAs
cleaned_data.T <- traits_df %>% filter(complete.cases(.))

# Remove duplicated rows and sum the count for duplicated species
trial_trait_dat <- cleaned_data.T %>%
   group_by(species_name) %>%
   summarise(
    feeding = first(feeding),
    niche = first(niche),
    current_dispersal = first(current_dispersal),
    voltinism = first(voltinism),
    size = first(size),
    stage = first(stage),
    .groups = 'drop'
 )

# Make sure that rows are given by the species name
trial_trait_dat <- as.data.frame(trial_trait_dat)
rownames(trial_trait_dat) <- trial_trait_dat$species_name


##### Create the functional trait dataframe to visulaise order trait space #####
# Name the columns 
trait_columns_ord <- c("feeding", "niche", "current_dispersal", "voltinism",
                   "size", "stage", "order")
trait_columns <- c("feeding", "niche", "current_dispersal", "voltinism",
                       "size", "stage")

# Create a new dataset with the selected columns and species 
traits_df_ord <- species_list_prez[, c("species_name", "feeding", "niche", 
                                       "current_dispersal", "voltinism", 
                                       "size", "stage", "order")]

# Remove duplicate species names
traits_df_ord <- traits_df_ord %>%
  distinct(species_name, .keep_all = TRUE)


# Ensure the species column is set as row names
trait_df_ord <- traits_df_ord %>% column_to_rownames(var = "species_name")



# Remove rows with any NAs
cleaned_data.T_ord <- traits_df_ord %>% filter(complete.cases(.))

# Remove duplicated rows and sum the count for duplicated species
trait_dat_ord <- cleaned_data.T_ord %>%
  group_by(species_name) %>%
  summarise(
    feeding = first(feeding),
    niche = first(niche),
    current_dispersal = first(current_dispersal),
    voltinism = first(voltinism),
    size = first(size),
    stage = first(stage),
    order = first(order),
    .groups = 'drop'
  )

# Make sure that rows are given by the species name
trait_dat_ord <- as.data.frame(trait_dat_ord)
rownames(trait_dat_ord) <- trait_dat_ord$species_name




###################### Create Environmental data matrix ########################

# Merge focal_attributes and distribution data frames based on the "code" column
env_df <- merge(focal_attributes, distribtion_df, by = "code", all.x = TRUE)
env_df <- merge(env_df, edge_distance_df, by = "code", all.x = TRUE)

# Calculate the average wind speed for each row
env_df$average_wind_speed <- rowMeans(env_df[wind_speed_columns], na.rm = TRUE)
                                                                                
# Delete unwanted columns to make a simplified dataframe
env_df$written <- NULL
env_df$flower_bud <- NULL
env_df$flower_bloom <- NULL
env_df$flower_dead <- NULL
env_df$flower_seed<- NULL
env_df$surrounding_veg_type <- NULL
env_df$f_type <- NULL
env_df$focal_patch <- NULL
env_df$wind_t0 <- NULL
env_df$wind_t15 <- NULL
env_df$wind_t30 <- NULL
env_df$wind_t45 <- NULL
env_df$wind_t60 <- NULL
env_df$wind_bearing <- NULL
env_df$weather_obvs <- NULL
env_df$grazed <- NULL


# Make sure that rows are given by the code
rownames(env_df) <- env_df$code
env_df$code <- NULL


# Identify rows with any NA values in env_df (W_B_2.7)
rows_with_na <- which(rowSums(is.na(env_df)) > 0)

# Display the rows that will be removed
print(paste("Rows to be removed:", rows_with_na))
 
# CCA cannot have any NA's so remove the rows that are missing 
env_df<- na.omit(env_df)





######################## Create species data matrix ############################

# Summarise the data to sum the counts of each species per site 
species_data <- merged_data %>%
  group_by(code, species_name) %>%
  summarise(count = sum(count, na.rm = TRUE)) %>%
  ungroup()

# Remove the rows for WB_2.7 as there were no neighbours and CCA won't work
species_data <- species_data[species_data$code != "W_B_2.7",]

# Create the species data matrix into a wide format where rows and sites and 
# columns are species
species_matrix <- species_data %>%
  pivot_wider(names_from = species_name, values_from = count, 
              values_fill = list(count = 0))

# Convert the matrix to a dataframe and add row names
species_matrix <- as.data.frame(species_matrix)
rownames(species_matrix) <- species_matrix$code
species_matrix$code <- NULL

# Extract column names from trait dataframes
columns_to_keep <- row.names(trait_dat_ord)

# Remove species that we don't have traits for 
species_matrix <- species_matrix[, colnames(species_matrix) %in% 
                                          columns_to_keep]

# Identify rows where the sum of species counts is zero as CCA won't work 
zero_sum_codes <- rownames(species_matrix)[rowSums(species_matrix) == 0]

# Print the codes that will be removed
if (length(zero_sum_codes) > 0) {
  print(paste("Codes with zero species counts that will be removed:", 
              zero_sum_codes))
}

# Remove the rows where the sum of species counts is zero
species_matrix <- species_matrix[rowSums(species_matrix) > 0, ]

# Remove the same "code" values from env_df
env_df <- env_df[!rownames(env_df) %in% zero_sum_codes, ]

# Reduce impact of zeros on elucidation distances with Hellinger transformation
species_hell <- decostand(species_matrix, method = "hellinger")






######################### Combine effects analysis #############################
cca_patch_quality <- cca(species_hell ~ width + height + leaves_alive +
                         leaves_dead, data = env_df)

cca_connectivity <- cca(species_hell ~ num_neigh + min_distance + edge_distance, 
                        data = env_df)

cca_environment <- cca(species_hell ~ humidity_base + humidity_top + micro_temp+ 
                      min_surrounding_veg_height + max_surrounding_veg_height +
                      average_wind_speed, data = env_df)

cca_combined_qc <- cca(species_hell ~ width + height + leaves_alive + 
                        leaves_dead +
                        num_neigh + min_distance + edge_distance, 
                        data = env_df)

cca_combined_all <- cca(species_hell ~ width + height + leaves_alive + 
                    leaves_dead + edge_distance + 
                    num_neigh + min_distance + humidity_base + 
                    humidity_top + micro_temp + min_surrounding_veg_height +
                    max_surrounding_veg_height + average_wind_speed, 
                    data = env_df)

# Compare each model
anova(cca_patch_quality, permutations = 9999)
anova(cca_connectivity, permutations = 9999)
anova(cca_environment, permutations = 9999)
anova(cca_combined_qc, permutations = 9999)
anova(cca_combined_all, permutations = 9999)




####################### Visualizing the trait space ############################
# Create the trait space 
trait_scores <- scores(cca_combined_all, display = "species")

# Ensure that the rownames of trial_trait_dat match the rownames of trait_scores
matching_trait_rows <- rownames(trial_trait_dat) %in% rownames(trait_scores)
trial_trait_dat <- trial_trait_dat[matching_trait_rows, ]

# Loop through each trait and create a plot
for (trait in trait_columns) {
  # Create the plot for the current trait
  trait_plot <- ggplot() +
    geom_point(data = as.data.frame(trait_scores), 
               aes(x = CCA1, y = CCA2, 
                   color = factor(trial_trait_dat[[trait]])), size = 3) +
    theme_minimal() +
    scale_color_discrete(name = trait) +
    theme(
      legend.position = c(0.9, 0.8),
      legend.background = element_rect(fill = "white", color = "black"), 
      legend.title = element_text(size = 15), 
      legend.text = element_text(size = 15) 
    )
  
  # Print each plot separately
  print(trait_plot)
}





################## Create the trait space shaded by order ######################
# Create the trait space 
trait_scores <- scores(cca_combined_all, display = "species")

# Ensure that the rownames of trial_trait_dat match the rownames of trait_scores
matching_trait_rows_ord <- rownames(trait_dat_ord) %in% rownames(trait_scores)
trait_dat_ord <- trait_dat_ord[matching_trait_rows_ord, ]

# Convert trait scores to a dataframe and add the corresponding orders
trait_space_df <- as.data.frame(trait_scores)
trait_space_df$Order <- trait_dat_ord$order

# Filter the trait space dataframe to include only the six key Orders
orders_of_interest <- c("Araneae", "Collembola", "Coleoptera", "Diptera", 
                        "Hemiptera", "Hymenoptera")
filtered_trait_space_df <- trait_space_df %>% filter(Order 
                                                     %in% orders_of_interest)

# Function to compute the convex hull for a given Order
compute_convex_hull <- function(data) {
  data[chull(data$CCA1, data$CCA2), ]  # Returns points forming the convex hull
}

# Create an empty list to store convex hulls for each Order
convex_hulls <- list()

# Compute convex hull for each Order
for (order in orders_of_interest) {
  order_data <- trait_space_df %>% filter(Order == order)
  convex_hulls[[order]] <- compute_convex_hull(order_data)
}

# Define a dataframe that combines all convex hulls with an Order column 
hull_data <- do.call(rbind, lapply(names(convex_hulls), function(order) {
  convex_hulls[[order]]$Order <- order
  convex_hulls[[order]]
}))

# Plot the trait space with convex hulls for each Order
trait_space_plot <- ggplot() +
  # Add jittered species points
  geom_jitter(data = trait_space_df, aes(x = CCA1, y = CCA2), 
              color = "black", size = 2, width = 0.1, height = 0.1) +
  
  # Add shaded convex hulls for each Order using ggplot2 default colors
  geom_polygon(data = hull_data, aes(x = CCA1, y = CCA2, fill = Order), 
               alpha = 0.3) +
  
  # Customize plot appearance
  theme_minimal() +
  labs(title = "Functional Trait Space Occupied by Key Orders",
       x = "CCA1", y = "CCA2") +
  
  # Use default ggplot color palette for the Order fill
  scale_fill_discrete(name = "Order") +
  theme(legend.position = "right")  # Adjust legend position as needed

# Print the plot
print(trait_space_plot)



######## Compute the functional trait space areas for each major order #########
library("geometry")  

# Function to compute the convex hull and its area for a given Order
compute_convex_hull_area <- function(data) {
  hull_points <- data[chull(data$CCA1, data$CCA2), ]  # Convex hull points
  hull_area <- geometry::convhulln(as.matrix(hull_points[, c("CCA1", "CCA2")]),
                                   options = "FA")$vol  # Calculate area
  return(hull_area)
}

# Create a list to store convex hull areas for each Order
convex_hull_areas <- list()

# Compute convex hull areas for each Order
for (order in orders_of_interest) {
  order_data <- trait_space_df %>% filter(Order == order)
  convex_hull_areas[[order]] <- compute_convex_hull_area(order_data)
}

# Convert the convex hull areas into a dataframe for easy comparison
convex_hull_areas_df <- data.frame(Order = names(convex_hull_areas),
                                   Area = unlist(convex_hull_areas))

# List the convex hull areas for each Order
convex_hull_areas_df









##################################### Plot the CCA #############################
library(ggplot2)
library(ggrepel)

# List of traits
traits <- c("voltinism", "feeding", "niche", "current_dispersal", "size")

# List of models to plot
models <- list(
  cca_patch_quality = cca_patch_quality,
  cca_connectivity = cca_connectivity,
  cca_environment = cca_environment,
  cca_combined_qc = cca_combined_qc,
  cca_combined_all = cca_combined_all
)


# Make sure species names is as rownames or an ID column
trait_df$Species <- rownames(trait_df)

# Loop through each model
for (model_name in names(models)) {
  
  model <- models[[model_name]]
  
  # Extract species and environmental scores
  species_scores <- scores(model, display = "species") %>% as.data.frame()
  species_scores$Species <- rownames(species_scores)
  
  env_scores <- scores(model, display = "bp")
  
  # Perform ANOVA for each environmental variable
  anova_results <- sapply(rownames(env_scores), function(var) {
    model_temp <- cca(species_hell ~ get(var), data = env_df)
    anova(model_temp, permutations = 999)$`Pr(>F)`[1]
  })
  
  # Create dataframe for environmental variables
  results_df <- data.frame(
    Variable = rownames(env_scores),
    R2 = round(env_scores[, "CCA1"], 2),
    P_value = format.pval(anova_results),
    Label = paste(1:nrow(env_scores))  # Number environmental vectors
  )
  
  # Merge species scores with traits
  species_trait_scores <- left_join(species_scores, trait_df, by = "Species")
  
  # Now loop through each trait
  for (trait in traits) {
    
    # Plot
    model_plot <- ggplot() +
      # Species points colored by trait
      geom_point(data = species_trait_scores, 
                 aes(x = CCA1, y = CCA2, color = .data[[trait]]),
                 size = 2) +
      # Environmental arrows
      geom_segment(aes(x = 0, y = 0, 
                       xend = env_scores[, "CCA1"], 
                       yend = env_scores[, "CCA2"]),
                   arrow = arrow(length = unit(0.3, "cm")), 
                   color = "red") +
      # Text labels for environmental arrows
      geom_text_repel(aes(x = env_scores[, "CCA1"], 
                          y = env_scores[, "CCA2"], 
                          label = results_df$Label),
                      color = "black",
                      size = 4,
                      nudge_x = 0.2,
                      nudge_y = 0.2,
                      fontface = "bold",
                      box.padding = 0.5) +
      # General theme
      theme_minimal() +
      ggtitle(paste("Biplot for", model_name, "-", trait)) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = trait))  # Title for color legend
    
    # Save the plot
    filename <- paste0("biplot_", model_name, "_", trait, ".png")
    ggsave(filename, model_plot, width = 8, height = 6, dpi = 300)
    
    # Print to screen if you want to check
    print(model_plot)
  }
  
  # Save the results table if you want
  assign(paste0(model_name, "_results"), results_df)
}



########## Computing the amount of variance explained by the axis ##############

# Get the eigenvalues from the CCA model
eigenvalues_qual <- cca_patch_quality$CCA$eig
eigenvalues_conn <- cca_connectivity$CCA$eig
eigenvalues_env <- cca_environment$CCA$eig
eigenvalues_comb_qc <- cca_combined_qc$CCA$eig
eigenvalues_comb_all <- cca_combined_all$CCA$eig


# Calculate the proportion of variance explained by each axis
variance_explained_qual <- eigenvalues_qual / sum(eigenvalues_qual)
variance_explained_conn <- eigenvalues_conn / sum(eigenvalues_conn)
variance_explained_env <- eigenvalues_env / sum(eigenvalues_env)
variance_explained_comb_qc <- eigenvalues_comb_qc / sum(eigenvalues_comb_qc)
variance_explained_comb_all <- eigenvalues_comb_all / sum(eigenvalues_comb_all)


# Display the proportion of variance explained by the first two axes
variance_explained_CCA1_qual <- variance_explained_qual[1] * 100
variance_explained_CCA2_qual <- variance_explained_qual[2] * 100

variance_explained_CCA1_conn <- variance_explained_conn[1] * 100
variance_explained_CCA2_conn <- variance_explained_conn[2] * 100

variance_explained_CCA1_env <- variance_explained_env[1] * 100
variance_explained_CCA2_env <- variance_explained_env[2] * 100

variance_explained_CCA1_comb_qc <- variance_explained_comb_qc[1] * 100
variance_explained_CCA2_comb_qc <- variance_explained_comb_qc[2] * 100

variance_explained_CCA1_comb_all <- variance_explained_comb_all[1] * 100
variance_explained_CCA2_comb_all <- variance_explained_comb_all[2] * 100

# Print variance explained by CCA1 and CCA2
cat("Variance explained by CCA1_qual:", variance_explained_CCA1_qual, "%\n")
cat("Variance explained by CCA2_qual:", variance_explained_CCA2_qual, "%\n")

cat("Variance explained by CCA1_conn:", variance_explained_CCA1_conn, "%\n")
cat("Variance explained by CCA2_conn:", variance_explained_CCA2_conn, "%\n")

cat("Variance explained by CCA1_env:", variance_explained_CCA1_env, "%\n")
cat("Variance explained by CCA2_env:", variance_explained_CCA2_env, "%\n")

cat("Variance explained by CCA1_comb_qc:", variance_explained_CCA1_comb_qc, 
    "%\n")
cat("Variance explained by CCA2_comb_qc:", variance_explained_CCA2_comb_qc,
    "%\n")

cat("Variance explained by CCA1_comb_all:", variance_explained_CCA1_comb_all,
    "%\n")
cat("Variance explained by CCA2_comb_all:", variance_explained_CCA2_comb_all,
    "%\n")



######### Create plots of trait centroids and environmental gradients ##########
library(ggrepel)
# Extract species scores from CCA
species_scores <- scores(cca_combined_all, display = "species") %>% 
  as.data.frame()

# Make sure species names are consistent
species_scores$Species <- rownames(species_scores)

# Align trait data with species scores
trait_df$Species <- rownames(trait_df)

# Join traits to species scores
species_traits <- left_join(species_scores, trait_df, by = "Species")

# Define traits to calculate centroids 
traits <- c("feeding", "niche", "current_dispersal", "voltinism", "size", 
            "stage")

# Calculate centroids (average CCA1 and CCA2 for each trait category)
centroids_list <- list()

for (trait in traits) {
  centroids <- species_traits %>%
    group_by(!!sym(trait)) %>%
    summarize(
      CCA1 = mean(CCA1, na.rm = TRUE),
      CCA2 = mean(CCA2, na.rm = TRUE)
    ) %>%
    mutate(Trait = trait) %>%
    rename(Category = !!sym(trait))
  
  centroids_list[[trait]] <- centroids
}

# Combine all centroids into one dataframe
trait_centroids <- bind_rows(centroids_list)

# Plot trait centroids
ggplot(trait_centroids, aes(x = CCA1, y = CCA2, color = Trait)) +
  geom_point(size = 4) +
  geom_text(aes(label = Category), vjust = -0.8, size = 3) +
  geom_text_repel(data = trait_centroids, aes(x = CCA1, y = CCA2, 
                                              label = Category), 
                  size = 3, max.overlaps = 20, box.padding = 0.5, 
                  show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Trait Centroids in CCA Space",
       x = "CCA1 Axis",
       y = "CCA2 Axis",
       color = "Trait") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# fit environmental variables 
envfit_result <- envfit(cca_combined_all, env_df, permutations = 999)  

# Extract arrows (vectors)
env_vectors <- as.data.frame(scores(envfit_result, "vectors"))
env_vectors$Variable <- rownames(env_vectors)

# Scale arrows so they fit nicely
arrow_multiplier <- 1  # You can adjust this value if arrows look too short/long
env_vectors <- env_vectors %>%
  mutate(CCA1 = CCA1 * arrow_multiplier,
         CCA2 = CCA2 * arrow_multiplier)

# Define a custom color palette for each trait
trait_colors <- c(
  "voltinism" = "red",
  "current_dispersal" = "blue",
  "feeding" = "darkgreen",
  "size" = "purple",
  "stage" = "orange",
  "niche" = "brown"
)

# Create a list of unique traits
trait_list <- unique(trait_centroids$Trait)

# Loop through each trait and plot separately
for (trait in trait_list) {
  # Filter for the specific trait
  trait_subset <- trait_centroids %>% filter(Trait == trait)
  
  # Make a new plot for each trait
  p <- ggplot() +
    # Centroids
    geom_point(data = trait_subset, aes(x = CCA1, y = CCA2), 
               color = trait_colors[trait], size = 4) +
    geom_text(data = trait_subset, aes(x = CCA1, y = CCA2, label = Category), 
              vjust = -0.8, size = 3, color = trait_colors[trait]) +
    # Environmental vectors
    geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                 arrow = arrow(length = unit(0.25, "cm")), color = "black") +
    geom_text_repel(data = env_vectors, aes(x = CCA1, y = CCA2, label = Variable),
                    size = 3, color = "black", max.overlaps = 20, 
                    box.padding = 0.5, show.legend = FALSE) +
    theme_minimal() +
    labs(title = paste("Trait Centroids + Environmental Vectors:", trait),
         x = "CCA1 Axis",
         y = "CCA2 Axis") +
    theme(
      legend.position = "none",   # No legend needed since color is fixed
      plot.title = element_text(hjust = 0.5)
    )
  
  # Print the plot to a new window
  print(p)
}
