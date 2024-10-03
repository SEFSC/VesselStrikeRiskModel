#Avoid_tot_strike_Sim_Slow_All

library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)
library(stringr)

setwd("~/Desktop/git_RiskModel/RiskModel/All_Code_2024")
o.folder <- "./model_outputs"


# Function to process files
process_files <- function(file_paths, sim_type, class) {
  summaries <- lapply(file_paths, function(file_path) {
    file <- read.csv(file_path, header = TRUE)
    year <- str_sub(file_path, -8, -5)
    
    
    # Pivot data
    tot.strike.pivot <- file %>%
      pivot_longer(cols = c(1:12), names_to = "month", values_to = "rate") %>%
      mutate(month = factor(month, levels = month.abb),
             sim = sim_type)
    
    # Summary statistics
    tot.strike.summary <- tot.strike.pivot %>%
      group_by(boot) %>%
      summarize(tot.mort = sum(rate), tot.deaths = sum(rate * 350)) %>%
      mutate(year = year)
    
    return(tot.strike.summary)
  })
  
  return(bind_rows(summaries))
}


classes <- c("M", "L", "XL")  

average_mortalities_by_class <- NULL

for (c in 1:length(classes)) {
  input.folder <- paste0(o.folder, "/summary_outputs_", classes[c])
  # Define scenarios and file patterns
  scenarios <- c("Status Quo", "Slow All Simulation")
  patterns <- c("tot_strike_Slow_None*", "Avoid_tot_strike_Sim_Slow_All*")
  
  # Process files for each scenario and class
  All_dfs <- lapply(seq_along(scenarios), function(i) {
    scenario <- scenarios[i]
    pattern <- patterns[i]
    
    file_paths <- list.files(input.folder, pattern = pattern, full.names = T)
    file_paths_no_avoid <- file_paths[grepl("No_Avoid", file_paths)]
    file_paths_avoid <- file_paths[!(grepl("No_Avoid", file_paths))]
    
    #avoidance excluded from results
    dfs_no_avoid <- process_files(file_paths_no_avoid, scenario, classes[c])  
    dfs_no_avoid$avoidance <- "Excluded"
    dfs_no_avoid$scenario <- scenarios[i]
    
    #avoidance included in results
    dfs_avoid <- process_files(file_paths_avoid, scenario, classes[c])  
    dfs_avoid$avoidance <- "Included"
    dfs_avoid$scenario <- scenarios[i]
    
    dfs <- rbind(dfs_no_avoid, dfs_avoid)
    
    return(dfs)
  })
  # Combine all dataframes
  all_dfs <- do.call(rbind, All_dfs)
  all_dfs$class <-  classes[c]
  
  # Combine all classes dataframes
  average_mortalities_by_class <- rbind(average_mortalities_by_class, all_dfs)
}

average_mortalities_by_class %>% 
  group_by(class, scenario, avoidance) %>%
  summarise(mean_tot.deaths = mean(tot.deaths), mean_tot.mort.rate = mean(tot.mort))

output.folder <- "./summary_outputs_"
if (!dir.exists(output.folder)) {dir.create(output.folder)}

saveRDS(average_mortalities_by_class, paste0(output.folder, "/average_mortalities_by_class.rds"))

mean_mortalities <- average_mortalities_by_class %>%
  group_by(class, avoidance, scenario) %>%
  summarise(mean_tot.deaths = mean(tot.deaths), sd_tot.deaths = sd(tot.deaths))

mean_mortalities_table <- mean_mortalities %>%
  group_by(class, avoidance) %>%
  mutate(percent_change=((last(mean_tot.deaths)-mean_tot.deaths)/last(mean_tot.deaths)*100))

mean_mortalities_table$class <- factor(mean_mortalities_table$class, levels = c("M", "L", "XL"))

saveRDS(mean_mortalities_table, "mean_mortalities_table.RDS")
write.csv(mean_mortalities_table, "mean_mortalities_table.csv")


average_mortalities_by_class$tot.deaths.r <- round(average_mortalities_by_class$tot.deaths, digits=0)
average_mortalities_by_class$class <- factor(average_mortalities_by_class$class, levels = c("M", "L", "XL"))
average_mortalities_by_class$avoidance <- factor(average_mortalities_by_class$avoidance, levels = c("Included", "Excluded"))

#Plot by each model scenario (Real world (i.e., status quo), SSZ simulation, and slow-all simulation)

mean_mortalities <- average_mortalities_by_class %>%
  group_by(class, avoidance, scenario) %>%
  summarise(mean_tot.deaths = mean(tot.deaths))

mean_mortalities_M_StatusQuo <- mean_mortalities %>% 
  filter (class=="M" & avoidance =="Included" & scenario == "Status Quo")

mean_mortalities_M_SlowAll <- mean_mortalities %>% 
  filter (class=="M" & avoidance =="Included" & scenario == "Slow All Simulation")


Medium_Bar_Mort <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="M" & average_mortalities_by_class$avoidance=="Included",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", position=position_dodge(),  alpha = 0.55, col="black") + 
  theme_bw() + labs(x = "Average Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  geom_vline(xintercept = (mean_mortalities_M_StatusQuo$mean_tot.deaths), color = "blue", linetype="dashed", linewidth = 1.5) +
  geom_vline(xintercept = (mean_mortalities_M_SlowAll$mean_tot.deaths), color = "red", linetype="dashed", linewidth = 1.5) +
  scale_fill_manual(name= "Scenario", values = c( "red", "blue"), labels=c( "Slow-all simulation", "Real world")) +
  ggtitle("Small/Medium (26-65 ft.)")+
  theme(axis.text=element_text(size=12)) +
  scale_x_continuous(limits=c(-1,12), labels = number_format(accuracy = 1)) 

mean_mortalities_L_StatusQuo <- mean_mortalities %>% 
  filter (class=="L" & avoidance =="Included" & scenario == "Status Quo")

mean_mortalities_L_SlowAll <- mean_mortalities %>% 
  filter (class=="L" & avoidance =="Included" & scenario == "Slow All Simulation")


Large_Bar_Mort <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="L" & average_mortalities_by_class$avoidance=="Included",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", position=position_dodge(),  alpha = 0.55, col="black") + 
  theme_bw() + labs(x = "Average Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  geom_vline(xintercept = (mean_mortalities_L_StatusQuo$mean_tot.deaths), color = "blue", linetype="dashed", linewidth = 1.5) +
  geom_vline(xintercept = (mean_mortalities_L_SlowAll$mean_tot.deaths), color = "red", linetype="dashed", linewidth = 1.5) +
  scale_fill_manual(name= "Scenario", values = c( "red", "blue"), labels=c( "Slow-all simulation", "Real world")) +
  ggtitle("Large (65-350 ft.)")+
  theme(axis.text=element_text(size=12)) +
  scale_x_continuous(labels = number_format(accuracy = 1))


mean_mortalities_XL_StatusQuo <- mean_mortalities %>% 
  filter (class=="XL" & avoidance =="Included" & scenario == "Status Quo")

mean_mortalities_XL_SlowAll <- mean_mortalities %>% 
  filter (class=="XL" & avoidance =="Included" & scenario == "Slow All Simulation")


OGV_Bar_Mort <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="XL" & average_mortalities_by_class$avoidance=="Included",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", position=position_dodge(),  alpha = 0.55, col="black") + 
  theme_bw() + labs(x = "Average Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  geom_vline(xintercept = (mean_mortalities_XL_StatusQuo$mean_tot.deaths), color = "blue", linetype="dashed", linewidth = 1.5) +
  geom_vline(xintercept = (mean_mortalities_XL_SlowAll$mean_tot.deaths), color = "red", linetype="dashed", linewidth = 1.5) +
  scale_fill_manual(name= "Scenario", values = c( "red", "blue"), labels=c( "Slow-all simulation", "Real world")) +
  ggtitle("OGV (>350 ft.)")+
  theme(axis.text=element_text(size=12)) +
  scale_x_continuous(labels = number_format(accuracy = 1))


plot_L_M <- ggarrange(Large_Bar_Mort, Medium_Bar_Mort, common.legend=T, legend="none",ncol=1, labels=c("b", "c"))
plot_all_slowAll_vs_SQ <- ggarrange(OGV_Bar_Mort, plot_L_M, common.legend = T, legend="bottom", labels=c("a", ""), widths=c(1.75,1))

plot_all_slowAll_vs_SQ
ggsave("./plot_all_risk_slowAll_vs_slowNone.png", plot_all_slowAll_vs_SQ, dpi=300, height=5, width=10)

#### Avoidance sensitivity analysis plots: including vs. excluding avoidance ####

Medium_Bar_Mort_Sensitivity <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="M",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", aes(fill = factor(scenario), alpha=avoidance), color="black") + #, color="white", alpha=.7, size=.5) +
  theme_bw() + labs(x = "Averge Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  #scale_fill_manual(name= "Scenario", values = c("lightblue", "blue"), labels=c("avoidance", "no avoidance")) +
  #ggtitle("Real world")+
  theme(axis.text=element_text(size=12))  +
  #ylim(c(0,350)) +
  #scale_x_discrete(breaks=seq(0,45,5)) +
  facet_wrap(~scenario, ncol = 1) +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values =c("red", "blue"), name = "Scenario", labels=c( "Slow-all simulation", "Real world")) +
  scale_alpha_manual(values=c(0.4, 1.0), name="Avoidance", labels=c("Included", "Excluded")) +
  ggtitle("Small/Medium (26-65 ft.)")

Large_Bar_Mort_Sensitivity <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="L",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", aes(fill = factor(scenario), alpha=avoidance), color="black") + #, color="white", alpha=.7, size=.5) +
  theme_bw() + labs(x = "Averge Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  #scale_fill_manual(name= "Scenario", values = c("lightblue", "blue"), labels=c("avoidance", "no avoidance")) +
  #ggtitle("Real world")+
  theme(axis.text=element_text(size=12))  +
  #ylim(c(0,350)) +
  #scale_x_discrete(breaks=seq(0,45,5)) +
  facet_wrap(~scenario, ncol = 1) +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values =c("red", "blue"), name = "Scenario", labels=c( "Slow-all simulation", "Real world")) +
  scale_alpha_manual(values=c(0.4, 1.0), name="Avoidance", labels=c("Included", "Excluded")) +
  ggtitle("Large (65-350 ft.)")

OGV_Bar_Mort_Sensitivity <- ggplot(data=average_mortalities_by_class[average_mortalities_by_class$class=="XL",], (aes(x=tot.deaths.r, fill=scenario))) + 
  geom_bar(stat="count", aes(fill = factor(scenario), alpha=avoidance), color="black") + #, color="white", alpha=.7, size=.5) +
  theme_bw() + labs(x = "Averge Mortalities Per Year", y = "Count") + guides(color = FALSE) +
  #scale_fill_manual(name= "Scenario", values = c("lightblue", "blue"), labels=c("avoidance", "no avoidance")) +
  #ggtitle("Real world")+
  theme(axis.text=element_text(size=12))  +
  #ylim(c(0,350)) +
  #scale_x_discrete(breaks=seq(0,45,5)) +
  facet_wrap(~scenario, ncol = 1) +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values =c("red", "blue"), name = "Scenario", labels=c( "Slow-all simulation", "Real world")) +
  scale_alpha_manual(values=c(0.4, 1.0), name="Avoidance", labels=c("Included", "Excluded")) +
  ggtitle("OGV (>350 ft.)")


sensitivity_plot_M_L <-  ggarrange(Large_Bar_Mort_Sensitivity, Medium_Bar_Mort_Sensitivity, ncol=1, common.legend = T, legend="none", labels=c("b", "c"))
sensitivity_plot_all <- ggarrange(OGV_Bar_Mort_Sensitivity, sensitivity_plot_M_L, ncol=2, common.legend = T, legend="bottom", widths=c(2,1), labels=c("a", ""))
sensitivity_plot_all
ggsave("./avoidance_sensitivity_plot_all.png", sensitivity_plot_all, dpi=300, height=5, width=10)

