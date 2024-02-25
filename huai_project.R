#loading packages
library(pacman)
p_load(dplyr, tidyverse, haven, rdrobust, ggplot2, magrittr, rddtools, broom, rddensity)

#loading data 
huai_data = read_dta('huairiver.dta')
huai_data <- huai_data %>%
  mutate(dist_bin = cut(dist_huai, breaks=quantile(dist_huai, probs = seq(0, 1, by = .05), na.rm =
                                                     TRUE)))

#creating binned scatterplot
regression_left <- lm(pm10 ~ dist_huai, data = subset(huai_data, dist_huai < 0))
regression_right <- lm(pm10 ~ dist_huai, data = subset(huai_data, dist_huai >= 0))

#binned scatterplot PM10
bin_plot1 = huai_data %>%
  group_by(dist_bin) %>%
  summarise(dist_huai = mean(dist_huai), pm10 = mean(pm10)) %>%
  ggplot(aes(x = dist_huai, y = pm10)) +
  geom_point(aes(color = ifelse(dist_huai >= 0, "North", "South")), size = 2, alpha = 1) +
  scale_color_manual(values = c("North" = "orange", "South" = "blue")) +
  geom_smooth(data = filter(huai_data, dist_huai <= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_smooth(data = filter(huai_data, dist_huai >= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "purple") +
  labs(
    title = "PM10 as a Function of Degree North of the Huai River",
    color = "",
    x = "Degrees North of the Huai River",
    y = "PM10 Levels (ug/m3)"
  )
bin_plot1

#binned scatterplot for temperature
bin_plot2 = huai_data %>%
  group_by(dist_bin) %>%
  summarise(dist_huai = mean(dist_huai), temp = mean(temp)) %>%
  ggplot(aes(x = dist_huai, y = temp)) +
  geom_point(aes(color = ifelse(dist_huai >= 0, "North", "South")), size = 2, alpha = 1) +
  scale_color_manual(values = c("North" = "orange", "South" = "blue")) +
  geom_smooth(data = filter(huai_data, dist_huai <= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_smooth(data = filter(huai_data, dist_huai >= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "purple") +
  labs(
    title = "Temperature as a Function of Degree North of the Huai River",
    color = "",
    x = "Degrees North of the Huai River",
    y = "Temperature"
  )
bin_plot2

#binned scatterplot for weather
bin_plot3 = huai_data %>%
  group_by(dist_bin) %>%
  summarise(dist_huai = mean(dist_huai), prcp = mean(prcp)) %>%
  ggplot(aes(x = dist_huai, y = prcp)) +
  geom_point(aes(color = ifelse(dist_huai >= 0, "North", "South")), size = 2, alpha = 1) +
  scale_color_manual(values = c("North" = "orange", "South" = "blue")) +
  geom_smooth(data = filter(huai_data, dist_huai <= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_smooth(data = filter(huai_data, dist_huai >= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "purple") +
  labs(
    title = "Precipitation as a Function of Degree North of the Huai River",
    color = "",
    x = "Degrees North of the Huai River",
    y = "Precipitation"
  )
bin_plot3


#binned scatterplot for wind speed
bin_plot4 = huai_data %>%
  group_by(dist_bin) %>%
  summarise(dist_huai = mean(dist_huai), wspd = mean(wspd)) %>%
  ggplot(aes(x = dist_huai, y = wspd)) +
  geom_point(aes(color = ifelse(dist_huai >= 0, "North", "South")), size = 2, alpha = 1) +
  scale_color_manual(values = c("North" = "orange", "South" = "blue")) +
  geom_smooth(data = filter(huai_data, dist_huai <= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_smooth(data = filter(huai_data, dist_huai >= 0), method = "lm", formula = y ~ poly(x), se = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "purple") +
  labs(
    title = "Wind Speed as a Function of Degree North of the Huai River",
    color = "",
    x = "Degrees North of the Huai River",
    y = "Wind Speed"
  )
bin_plot4


#calculating discontinuity effect
#calc by hand 
cutoff_point = 0

regr_hand = rdd_data(y = huai_data$pm10, 
         x = huai_data$dist_huai, 
         cutpoint = 0) %>% 
  rdd_reg_lm(slope = "separate") 
summary(regr_hand)
#calc using rdrobust
# Estimate treatment effect using rdrobust
rd_result = rdrobust(huai_data$pm10, huai_data$dist_huai, 
                      c = cutoff_point, 
                      kernel = "uniform", 
                      bwselect = "mserd", 
                      m = 2,
                      vce = 'hc1')

# Print the treatment effect
summary(rd_result)

#robustness checks
#Symmetric bandwidth
rd_symmetric = rdrobust(huai_data$pm10, huai_data$dist_huai,
                        covs = cbind(huai_data$north_huai * huai_data$dist_huai),
                        c = cutoff_point, 
                        kernel = "Uniform", 
                        m = 2)
summary(rd_symmetric)

#triangular kernel
rd_kernel = rdrobust(huai_data$pm10, huai_data$dist_huai,
                     covs = cbind(huai_data$north_huai * huai_data$dist_huai),
                     c = cutoff_point,
                     kernel = "triangular",
                     bwselect = "mserd", 
                     m = 2,
                     h = 30)
summary(rd_kernel)

# Poly functional form
rd_poly = rdrobust(huai_data$pm10, huai_data$dist_huai, 
                     covs = cbind(huai_data$north_huai * huai_data$dist_huai),
                     c = cutoff_point, 
                     kernel = "Uniform", 
                     bwselect = "mserd", 
                     m = 1,
                     p = 2,
                     h = 30)
summary(rd_poly)

#smallest difference-in-group means estimator
rd_small = rdrobust(y = huai_data$pm10, x = huai_data$dist_huai, 
                    covs = cbind(huai_data$north_huai * huai_data$dist_huai),
                    c = cutoff_point, 
                    h = 1,
                    p = 0,
                    kernel = "Uniform", 
                    bwselect = "msesum", 
                    m = 2)
summary(rd_small)


#make table with regressions
# Extract summary information for each RD model
#Create columns for each test
col1=c(rd_result$Estimate[1], 
       rd_result$se[3], 
       rd_result$ci[3,1], 
       rd_result$ci[3,2], 
       rd_result$bws[1], 
       rd_result$N_h[1])
col2 = c(rd_symmetric$Estimate[1], 
        rd_symmetric$se[3], 
        rd_symmetric$ci[3,1], 
        rd_symmetric$ci[3,2], 
        rd_symmetric$bws[1], 
        rd_symmetric$N_h[1])
col3 = c(rd_kernel$Estimate[1], 
        rd_kernel$se[3], 
        rd_kernel$ci[3,1], 
        rd_kernel$ci[3,2], 
        rd_kernel$bws[1], 
        rd_kernel$N_h[1])
col4 = c(rd_poly$Estimate[1], 
        rd_poly$se[3], 
        rd_poly$ci[3,1], 
        rd_poly$ci[3,2], 
        rd_poly$bws[1], 
        rd_poly$N_h[1])
col5 = c(rd_small$Estimate[1], 
       rd_small$se[3], 
       rd_small$ci[3,1], 
       rd_small$ci[3,2], 
       rd_small$bws[1], 
       rd_small$N_h[1])
#Bind all columns together to create the table
result_table <- cbind(
  col1, col2, col3, col4, col5
)

# Optionally, you can add column names
colnames(result_table) <- c(
  "Regression Discontinuity", "Regression Discontinuity Symmetric", "Regression Discontinuity Kernel", "Regression Discontinuity Linear From", "Regression Discontinuity Small Group Means"
)
# Print or view the resulting table
print(result_table)

#testing covariate smoothness
#testing dist_huai and temp
rd_wind = rdrobust(y = huai_data$temp,
                   x = huai_data$dist_huai,
                   c = cutoff_point, 
                   kernel = "Uniform", 
                   bwselect = "mserd", 
                   p = 1)

summary(rd_wind)
#test dist_huai and precipitation
rd_prcp = rdrobust(y = huai_data$prcp,
              x = huai_data$dist_huai,
              c = cutoff_point,
              kernel = "Uniform",
              bwselect = "mserd",
              p = 1)
summary(rd_prcp)

#test dist_huai and wind
rd_wind = rdrobust(y = huai_data$wspd,
                   x = huai_data$dist_huai,
                   c = cutoff_point,
                   kernel = "Uniform",
                   bwselect = 'mserd',
                   p = 1)
summary(rd_wind)


#manipulation test
rdplotdensity(rdd = rddensity(huai_data$dist_huai, c = 0, p = 1),
              X = huai_data$dist_huai,
              type = "both")

#large confidence intervals indicate potential noise in the data
#overlap in confidence intervals express lack of manipulation between North and South observations
#shows smoothness between groups and minimal concern about manipulation between groups

#9b placebo tests
rd_list = list()

j = 1

for (i in-5:5) {
  out = rdrobust(
    y = huai_data$pm10,
    x = huai_data$dist_huai,
    c = i,
    kernel = "uniform"
  )
  
  ci_df = data.frame(
    "dist" = i,
    "est" = out$Estimate[1],
    "ci_low" = out$ci[1,1],
    "ci_high" = out$ci[1,2]
  )
  rd_list[[j]] = ci_df
  
  j = j + 1
}

ci_df = bind_rows(rd_list)

ggplot(data = ci_df,
       aes(x = dist,
           y = est)) +
  labs(x = "Estimate", y = "Distance from the Huai River") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
               color = "blue",  # Set error bar color
               width = 0.5) + 
  geom_hline(yintercept = 0)
  theme_classic()






