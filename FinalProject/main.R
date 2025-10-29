# ---------------------------------------------------------------------------- #
# Title: Demand estimation assignment                                          #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# 0 - Import data and functions ----

  # Own functions
    # Utilities: Useful functions
    source("scripts/S1_utilities.r")
    source("scripts/S2_eda.r")
    source("scripts/S3_logits.r")
    source("scripts/S4_mcosts.r")
    source("scripts/S5_Counterfactuals.r")
    # Import libraries
    libcall()
    
  # Import from CSV
    df <- read_csv("data/problem_set_aggregate_data.csv") %>%
      clean_names() %>%
      select(-matches("^unnamed")) %>%
      select(-x1) %>% 
      rename(x1 = x1_2) %>% 
      mutate(
        product_type = case_when(
          reg_cola  == 1 ~ "Regular Cola",
          diet_cola == 1 ~ "Diet Cola",
          non_cola  == 1 ~ "Non-Cola",
          TRUE ~ "Unlabeled"
        )
      )
    
# 1 - Exploratory Data Analysis (EDA) ----
  # Count of markets and products
    table1f(df)
  # General summary statistics
    table2f(df)
  # Firm level activity
    table3f(df)
  # Price Dispersions evidence
    tablea1f(df)
    figure1f(df)
  # Correlations with X's and W's
    figure2f(df)
  
# 2 - Logit Model ----
    Slogit <- logit_ols(df,xvars = c("x1","x2","x3","reg_cola","diet_cola"))

# 3 - IV-Logit Model ----
    IVlogit <- logit_iv(df,xvars = c("x1","x2","x3","reg_cola","diet_cola"),
                        ivars = c("w1","w2"))
    
# 4 - Marginal Costs ----
    mcdata <- margcost(df,fitmod = IVlogit)
    mc_tbl <- mcdata$mcdf
# 5 - Scenario simulation ----
    
    simul <- scenario_simulation(df,mc_tbl = mc_tbl,fitmodel = IVlogit,drop_product_id = 1,doplot=T)
    
    
# - 6 Welfare change ----    
    welfare_table(df, fitmodel = IVlogit, mc_tbl,product_drops = "All",outtex = "results/table4_welfare.tex")
    
# 7 - IV Nested Logit model ----
    NIVlogit <- nested_logit_iv(df,xvars = c("x1","x2","x3"),ivars=c("w1","w2"),
                                intercept = T)
    
    
# 8 - Own and Cross Elasticities ----    
    elast_density_plot(Logit=IVlogit,Nested=NIVlogit)
    
    # Elasticity matrices by market - comparison between models
      # Market 1
      plot_elasticity_matrix(pairs = IVlogit$elast_all,market_id = 1,method="IVlogit")
      plot_elasticity_matrix(pairs = NIVlogit$elast_all,market_id = 1,method="NIVlogit")
      # Market 350
      plot_elasticity_matrix(pairs = IVlogit$elast_all,market_id = 350,method="IVlogit")
      plot_elasticity_matrix(pairs = NIVlogit$elast_all,market_id = 350,method="NIVlogit")
      # Makert 599
      plot_elasticity_matrix(pairs = IVlogit$elast_all,market_id = 599,method="IVlogit")
      plot_elasticity_matrix(pairs = NIVlogit$elast_all,market_id = 599,method="NIVlogit")
      
      
      