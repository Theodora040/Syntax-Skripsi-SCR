# DASHBOARD ANALISIS TINGKAT KEMISKINAN DI PROVINSI JAWA BARAT
# MENGGUNAKAN METODE SPATIAL CLUSTERING REGRESSION

# ==============================================================================
# 1. LOAD REQUIRED PACKAGES
# ==============================================================================

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(corrplot)
library(glmnet)
library(MASS)
library(robustbase)
library(car)
library(lmtest)
library(sandwich)
library(broom)
library(tidyr)
library(patchwork)
library(scales)
library(ggrepel)
library(kableExtra)
library(performance)
library(reshape2)
library(gridExtra)
library(viridis)
library(tibble)
library(sf)
library(leaflet)
library(classInt)

theme_set(theme_minimal(base_size = 12))

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

format_coef <- function(x, var_name) {
  if (grepl("X8|X8D", var_name)) {
    return(format(round(x, 6), scientific = FALSE, nsmall = 6))
  } else {
    return(format(round(x, 4), scientific = FALSE, nsmall = 4))
  }
}

format_t_test <- function(value, var_name, type) {
  if (type == "estimate") {
    if (grepl("X8|X8D", var_name)) {
      return(format(round(value, 6), scientific = FALSE, nsmall = 6))
    } else {
      return(format(round(value, 4), scientific = FALSE, nsmall = 4))
    }
  } else {
    return(format(round(value, 4), scientific = FALSE, nsmall = 4))
  }
}

show_ridge_coef <- function(model, d_var) {
  tryCatch({
    coef_matrix <- as.matrix(coef(model))
    
    if (ncol(coef_matrix) > 1) {
      coef_values <- coef_matrix[, ncol(coef_matrix)]
    } else {
      coef_values <- coef_matrix[, 1]
    }
    
    df <- data.frame(
      Variable = names(coef_values),
      Coefficient = as.numeric(coef_values),
      stringsAsFactors = FALSE
    )
    
    df$Type <- "X - Regularized"
    
    for (i in 1:nrow(df)) {
      var <- df$Variable[i]
      
      if (var == "(Intercept)") {
        df$Type[i] <- "Intercept (No Penalty)"
      } else if (var == d_var) {
        df$Type[i] <- "D (No Penalty)"
      } else if (grepl(paste0(d_var, "$"), var)) {
        df$Type[i] <- "Interaction (X*D) - Regularized"
      }
    }
    
    df <- df[df$Coefficient != 0 | df$Variable == "(Intercept)", ]
    
    df$Coefficient_Formatted <- sapply(1:nrow(df), function(i) {
      x <- df$Coefficient[i]
      var_name <- df$Variable[i]
      
      if (grepl("X8|X8D", var_name)) {
        return(format(round(x, 6), scientific = FALSE, nsmall = 6))
      } else {
        return(format(round(x, 4), scientific = FALSE, nsmall = 4))
      }
    })
    
    return(data.frame(
      Variable = df$Variable,
      Coefficient = df$Coefficient_Formatted,
      Type = df$Type,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    return(data.frame(
      Variable = "Error displaying coefficients",
      Coefficient = "NA",
      Type = "Error"
    ))
  })
}

perform_significance_tests <- function(model, x_vars_main, d_var, x_vars_interaction) {
  f_stat <- summary(model)$fstatistic
  p_val_f <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  
  f_test <- data.frame(
    Test = "F-Test (Overall Model Significance)",
    Statistic = round(f_stat[1], 4),
    df1 = f_stat[2],
    df2 = f_stat[3],
    p_value = ifelse(p_val_f < 0.0001, 
                     format(p_val_f, scientific = TRUE, digits = 4),
                     format(p_val_f, scientific = FALSE, digits = 4)),
    Significance = ifelse(p_val_f < 0.1, "Signifikan", "Tidak Signifikan"),
    stringsAsFactors = FALSE
  )
  
  t_test <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("Variable") %>%
    rename(
      Estimate = "Estimate",
      StdError = "Std. Error",
      t_value = "t value",
      p_value = "Pr(>|t|)"
    ) %>%
    mutate(
      Significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.10 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    rename(
      `Std. Error` = "StdError",
      `t value` = "t_value",
      `Pr(>|t|)` = "p_value"
    )
  
  wald_tests <- list()
  
  if (d_var %in% names(coef(model))) {
    tryCatch({
      wald_intercept <- linearHypothesis(model, paste0("(Intercept) + ", d_var, " = 0"))
      wald_p_val <- wald_intercept$`Pr(>F)`[2]
      
      wald_tests$Intercept <- data.frame(
        Test = "Wald Test - Intercept",
        Hypothesis = paste0("(Intercept) + ", d_var, " = 0"),
        F_Statistic = round(wald_intercept$F[2], 4),
        p_value = ifelse(wald_p_val < 0.0001,
                         format(wald_p_val, scientific = TRUE, digits = 4),
                         format(wald_p_val, scientific = FALSE, digits = 4)),
        Significance = ifelse(wald_p_val < 0.1, "Signifikan", "Tidak Signifikan"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      wald_tests$Intercept <- NULL
    })
  }
  
  for(i in 1:length(x_vars_main)) {
    var_name <- x_vars_main[i]
    interaction_var <- paste0(var_name, d_var)
    
    if(var_name %in% names(coef(model)) && interaction_var %in% names(coef(model))) {
      tryCatch({
        hypothesis <- paste0(var_name, " + ", interaction_var, " = 0")
        wald_result <- linearHypothesis(model, hypothesis)
        wald_p_val <- wald_result$`Pr(>F)`[2]
        
        wald_tests[[var_name]] <- data.frame(
          Test = paste("Wald Test - Pengaruh (", var_name, ")", sep = ""),
          Hypothesis = paste0(var_name, " + ", interaction_var, " = 0"),
          F_Statistic = round(wald_result$F[2], 4),
          p_value = ifelse(wald_p_val < 0.0001,
                           format(wald_p_val, scientific = TRUE, digits = 4),
                           format(wald_p_val, scientific = FALSE, digits = 4)),
          Significance = ifelse(wald_p_val < 0.1, "Signifikan", "Tidak Signifikan"),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        wald_tests[[var_name]] <- NULL
      })
    }
  }
  
  if (length(wald_tests) > 0) {
    wald_df <- do.call(rbind, wald_tests)
    rownames(wald_df) <- NULL
  } else {
    wald_df <- data.frame()
  }
  
  return(list(
    f_test = f_test,
    t_test = t_test,
    wald_tests = wald_df
  ))
}

dt_scroller <- function(x, height = 400) {
  datatable(x,
            options = list(dom = "t",
                           scrollX = TRUE,
                           scrollY = height,
                           paging = FALSE),
            class = "display compact nowrap",
            rownames = FALSE)
}

# ==============================================================================
# 3. CUSTOM CSS
# ==============================================================================

custom_css <- "
/* ============================================================
   FORCE ALL TEXT BLACK
============================================================ */
*, body, h1, h2, h3, h4, h5, h6, p, label, span,
a, strong, button, input, select, textarea {
  color: #000000 !important;
}

/* Fix hover text */
a:hover {
  color: #000000 !important;
}

/* ============================================================
   SIDEBAR COLOR FIX
============================================================ */
.skin-blue .main-sidebar, 
.skin-blue .left-side,
.skin-blue .wrapper .sidebar {
  background-color: #D9ECFF !important;
}

.skin-blue .sidebar-menu > li > a {
  background-color: transparent !important;
  color: #000000 !important;
}

.skin-blue .sidebar-menu > li.active > a {
  background-color: #ADD8E6 !important;
  color: #000000 !important;
}

.skin-blue .sidebar-menu > li:hover > a {
  background-color: #BEE3FF !important;
  color: #000000 !important;
}

.skin-blue .sidebar-menu i {
  color: #000000 !important;
}

.skin-blue .sidebar-menu > li.active > a:after {
  display: none !important;
}

/* ============================================================
   HEADER & BODY STYLE
============================================================ */
.skin-blue .main-header .logo {
  background-color: #2c3e50 !important;
  font-weight: bold;
  color: #FFFFFF !important;
}

.skin-blue .main-header .navbar {
  background-color: #34495e !important;
}

/* ============================================================
   INPUT & TABLES
============================================================ */
.selectize-input, .dataTables_wrapper, table.dataTable td, table.dataTable th {
  color: #000000 !important;
}

/* ============================================================
   BUTTONS
============================================================ */
.btn, .btn-primary, .btn-success, .btn-warning, .btn-info {
  font-weight: bold !important;
  color: #000000 !important;
}

/* ============================================================
   SCROLLBAR
============================================================ */
::-webkit-scrollbar { width: 10px; }
::-webkit-scrollbar-track { background: #f1f1f1; border-radius: 5px; }
::-webkit-scrollbar-thumb { background: #3498db; border-radius: 5px; }
::-webkit-scrollbar-thumb:hover { background: #2980b9; }

/* Style untuk outlier table */
#outliers_table table.dataTable td {
  text-align: center !important;
}

/* Warna untuk jumlah methods detected */
.badge-outlier-2 {
  background-color: #e74c3c !important;
  color: white !important;
}

.badge-outlier-1 {
  background-color: #f39c12 !important;
  color: white !important;
}

/* Warna untuk jenis variabel di tabel ridge */
.bg-intercept {
  background-color: #f8f9fa !important;
}
.bg-d-nopenalty {
  background-color: #d4edda !important;
}
.bg-x-regularized {
  background-color: #fff3cd !important;
}
.bg-interaction {
  background-color: #d1ecf1 !important;
}

/* Style untuk significance stars */
.sig-star {
  font-weight: bold;
  color: #e74c3c;
}

/* ============================================================
   MAP STYLING
============================================================ */
.leaflet-container {
  font-family: inherit !important;
  font-size: 14px !important;
}

.leaflet-popup-content {
  color: #000000 !important;
}

.leaflet-control-attribution {
  background-color: rgba(255, 255, 255, 0.7) !important;
  color: #000000 !important;
}

/* Legend styling */
.info.legend {
  background-color: white !important;
  padding: 10px !important;
  border-radius: 5px !important;
  box-shadow: 0 1px 5px rgba(0,0,0,0.4) !important;
  color: #000000 !important;
}

/* Popup styling untuk polygon map */
.polygon-popup {
  font-family: Arial, sans-serif !important;
  font-size: 12px !important;
  min-width: 200px !important;
}

.polygon-popup h4 {
  margin-top: 0 !important;
  color: #2c3e50 !important;
  border-bottom: 2px solid #3498db !important;
  padding-bottom: 5px !important;
}

.polygon-popup table {
  width: 100% !important;
  border-collapse: collapse !important;
  margin: 5px 0 !important;
}

.polygon-popup td {
  padding: 3px 5px !important;
  border-bottom: 1px solid #eee !important;
}

.polygon-popup tr:last-child td {
  border-bottom: none !important;
}

/* Warna untuk klaster di popup */
.klaster-1 {
  color: #00BCD4 !important;
  font-weight: bold;
}

.klaster-2 {
  color: #E91E63 !important;
  font-weight: bold;
}
"

# ==============================================================================
# 4. UI COMPONENTS
# ==============================================================================

# Header -----------------------------------------------------------------------
header <- dashboardHeader(
  title = tags$div(
    tags$img(
      src    = "https://img.icons8.com/color/48/000000/statistics.png",
      height = "20px",
      style  = "margin-right: 8px;"
    ),
    "SCR Dashboard",
    style = "display: flex; align-items: center;"
  ),
  titleWidth = 300
)

# Sidebar ----------------------------------------------------------------------
sidebar <- dashboardSidebar(
  width = 300,
  sidebarMenu(
    id = "tabs",
    menuItem("Data Input",          tabName = "data",       icon = icon("database")),
    menuItem("Exploratory Analysis",tabName = "explore",    icon = icon("chart-bar")),
    menuItem("Clustering Analysis", tabName = "clustering", icon = icon("map-marked-alt")),
    menuItem("Model Specification", tabName = "model",      icon = icon("cogs")),
    menuItem("Model Estimation",    tabName = "estimation", icon = icon("calculator")),
    menuItem("Diagnostics",         tabName = "diagnostics",icon = icon("stethoscope")),
    menuItem("Robust Methods",      tabName = "robust",     icon = icon("shield-alt")),
    menuItem("Significance Tests",  tabName = "significance", icon = icon("check-circle")),
    menuItem("Results & Reports",   tabName = "results",    icon = icon("file-alt")),
    menuItem("About",               tabName = "about",      icon = icon("info-circle"))
  ),
  hr(),
  
  # BOX UNTUK REGRESSION DATA
  box(
    title       = "Regression Data Upload",
    width       = 12,
    status      = "primary",
    solidHeader = TRUE,
    collapsible = TRUE,
    
    fileInput(
      "file_regression", "Upload Regression Data (Excel)",
      accept      = c(".xlsx", ".xls"),
      buttonLabel = "Browse...",
      placeholder = "No file selected"
    ),
    
    br(),
    
    actionBttn(
      "sample_regression_data", "Use Sample Regression Data",
      style = "gradient",
      color = "warning",
      size  = "sm",
      block = TRUE
    ),
    
    hr(),
    uiOutput("regression_data_info")
  ),
  
  # BOX UNTUK CLUSTERING DATA (SHAPEFILE)
  box(
    title       = "Clustering Data Upload",
    width       = 12,
    status      = "info",
    solidHeader = TRUE,
    collapsible = TRUE,
    collapsed   = TRUE,
    
    h5("Upload file shapefile (.shp, .shx, .dbf):"),
    
    fileInput(
      "shp_file", "Pilih file shapefile",
      accept = c(".shp", ".dbf", ".shx"),
      multiple = TRUE,
      buttonLabel = "Browse...",
      placeholder = "No files selected"
    ),
    
    tags$div(
      style = "font-size: 11px; color: #666; margin-top: 5px;",
      icon("info-circle"),
      " Pilih SEMUA file sekaligus:",
      tags$ul(
        style = "margin: 5px 0; padding-left: 20px;",
        tags$li("File .shp (wajib)"),
        tags$li("File .shx (wajib)"),
        tags$li("File .dbf (wajib)"),
        tags$li("File .prj (opsional)")
      )
    ),
    
    hr(),
    
    actionBttn(
      "sample_polygon_data", "Use Sample Polygon Data",
      style = "gradient",
      color = "success",
      size  = "sm",
      block = TRUE
    )
  ),
  
  # BOX VARIABLE SELECTION (UNTUK REGRESI SAJA)
  box(
    title       = "Regression Model Specification",
    width       = 12,
    status      = "warning",
    solidHeader = TRUE,
    collapsible = FALSE,
    
    selectInput("y_var", "Dependent Variable (Y):", choices = NULL),
    selectInput("d_var", "Cluster Variable (D):", choices = NULL),
    
    pickerInput(
      "x_vars", "Independent Variables (X):",
      choices = NULL,
      options = list(
        `actions-box`          = TRUE,
        `selected-text-format` = "count > 3",
        `count-selected-text`  = "{0} variables selected"
      ),
      multiple = TRUE
    ),
    
    hr(),
    verbatimTextOutput("regression_summary"),
    hr(),
    
    actionBttn(
      "run_regression", "Run Regression Analysis",
      style = "gradient",
      color = "success",
      size  = "md",
      block = TRUE,
      icon  = icon("play")
    )
  ),
  
  # BOX SETTINGS (UNTUK REGRESI SAJA)
  box(
    title       = "Regression Settings",
    width       = 12,
    status      = "danger",
    solidHeader = TRUE,
    collapsible = TRUE,
    collapsed   = TRUE,
    
    awesomeCheckboxGroup(
      "methods", "Estimation Methods:",
      choices  = list(
        "OLS"              = "ols",
        "Ridge Regression" = "ridge",
        "Robust Regression"= "robust",
        "Robust Ridge"     = "robust_ridge"
      ),
      selected = c("ols", "ridge"),
      status   = "primary"
    ),
    
    hr(),
    
    sliderInput(
      "cv_folds", "CV Folds for Ridge:",
      min = 3, max = 10, value = 5
    ),
    
    selectInput(
      "robust_method", "Robust Method:",
      choices  = c("M-estimation" = "M",
                   "MM-estimation" = "MM",
                   "S-estimation"  = "S"),
      selected = "M"
    )
  )
)

# Body -------------------------------------------------------------------------
body <- dashboardBody(
  tags$head(
    tags$style(HTML(custom_css))
  ),
  
  tabItems(
    # ====================== TAB: DATA INPUT ===================================
    tabItem(
      tabName = "data",
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Regression Data Preview",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            
            DTOutput("regression_data_preview"),
            hr(),
            fluidRow(
              valueBoxOutput("n_obs"),
              valueBoxOutput("n_vars"),
              valueBoxOutput("n_missing")
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Regression Data Structure",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            verbatimTextOutput("regression_structure")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Variable Types",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("var_types_plot")
          )
        )
      )
    ),
    
    # ====================== TAB: EXPLORATORY =================================
    tabItem(
      tabName = "explore",
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Descriptive Statistics",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("desc_stats")
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Correlation Matrix",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("corr_plot", height = "500px")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Variable Distributions",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("dist_plot", height = "500px")
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Scatter Plot Matrix",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed   = TRUE,
            plotOutput("scatter_matrix", height = "600px")
          )
        )
      )
    ),
    
    # ====================== TAB: CLUSTERING ANALYSIS =========================
    tabItem(
      tabName = "clustering",
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Clustering Polygon Map",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            
            fluidRow(
              column(6,
                     selectInput("map_color_var", "Color Variable:",
                                 choices = c("Klaster" = "Klaster",
                                             "Tingkat Kemiskinan" = "TngktKm"),
                                 selected = "Klaster")
              ),
              column(6,
                     numericInput("map_zoom", "Zoom Level:", 
                                  value = 8, min = 6, max = 12, step = 0.5)
              )
            ),
            
            hr(),
            
            leafletOutput("cluster_polygon_map", height = "500px"),
            
            tags$div(
              style = "margin-top: 10px; font-size: 12px; color: #666; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
              icon("info-circle"),
              " Klik pada poligon untuk melihat detail wilayah."
            )
          ),
          
          box(
            title       = "Clustering Summary",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            tableOutput("clustering_summary_poly")
          )
        )
      )
    ),
    
    # ====================== TAB: MODEL SPEC ===================================
    tabItem(
      tabName = "model",
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Regression Model Specification",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            
            h4("Model 1: Pooled Regression (Same Slopes)"),
            uiOutput("model1_formula"),
            hr(),
            
            h4("Model 2: Full Interaction (Different Slopes)"),
            uiOutput("model2_formula"),
            hr(),
            
            h4("Model 3: Ridge Regression"),
            uiOutput("model3_formula"),
            hr(),
            
            h4("Model 4: Robust Regression"),
            uiOutput("model4_formula"),
            hr(),
            
            h4("Model 5: Robust Ridge Regression"),
            uiOutput("model5_formula")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Mathematical Specification",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            withMathJax(),
            uiOutput("math_spec")
          ),
          box(
            title       = "Cluster Distribution (from D variable)",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            plotlyOutput("cluster_dist_plot", height = "300px")
          )
        )
      )
    ),
    
    # ====================== TAB: ESTIMATION ===================================
    tabItem(
      tabName = "estimation",
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Model Comparison",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("model_comparison_plot", height = "400px")
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          box(
            title       = "OLS Results",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.methods.includes('ols')",
              DTOutput("ols_results")
            )
          ),
          box(
            title       = "Ridge Regression Results",
            width       = 12,
            status      = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            conditionalPanel(
              condition = "input.methods.includes('ridge')",
              DTOutput("ridge_results")
            )
          )
        ),
        column(
          width = 6,
          box(
            title       = "Coefficient Comparison",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("coef_comparison", height = "500px")
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Performance Metrics",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("performance_table")
          )
        )
      )
    ),
    
    # ====================== TAB: DIAGNOSTICS ==================================
    tabItem(
      tabName = "diagnostics",
      tabBox(
        width = 12,
        title = "Model Diagnostics",
        
        tabPanel(
          "Multicollinearity",
          fluidRow(
            column(
              width = 6,
              box(
                width  = 12,
                status = "info",
                plotOutput("vif_plot", height = "400px")
              )
            ),
            column(
              width = 6,
              box(
                width  = 12,
                status = "info",
                DTOutput("vif_table")
              )
            )
          )
        ),
        tabPanel(
          "Residual Analysis",
          fluidRow(
            column(
              width = 6,
              box(
                width  = 12,
                status = "warning",
                plotOutput("resid_fitted_plot", height = "400px")
              ),
              box(
                width  = 12,
                status = "warning",
                plotOutput("resid_hist", height = "300px")
              )
            ),
            column(
              width = 6,
              box(
                width  = 12,
                status = "warning",
                plotOutput("qq_plot", height = "400px")
              ),
              box(
                width  = 12,
                status = "warning",
                plotOutput("scale_location_plot", height = "300px")
              )
            )
          )
        ),
        tabPanel(
          "Normality Tests",
          fluidRow(
            column(
              width = 12,
              box(
                width = 12,
                status = "info",
                title = "Kolmogorov-Smirnov Test",
                verbatimTextOutput("normality_test")
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              box(
                width = 12,
                status = "info",
                plotOutput("normality_plot", height = "400px")
              )
            )
          )
        ),
        tabPanel(
          "Heteroskedasticity",
          fluidRow(
            column(
              width = 6,
              box(
                width  = 12,
                status = "danger",
                verbatimTextOutput("bp_test")
              )
            ),
            column(
              width = 6,
              box(
                width  = 12,
                status = "danger",
                plotOutput("scale_location_plot2", height = "400px")
              )
            )
          )
        ),
        tabPanel(
          "Autocorrelation",
          fluidRow(
            column(
              width = 6,
              box(
                width  = 12,
                status = "warning",
                verbatimTextOutput("autocorrelation_test")
              )
            ),
            column(
              width = 6,
              box(
                width  = 12,
                status = "warning",
                plotOutput("acf_plot", height = "400px")
              )
            )
          )
        ),
        tabPanel(
          "Outlier Detection",
          fluidRow(
            column(
              width = 6,
              box(
                width  = 12,
                status = "danger",
                plotOutput("cooks_plot", height = "400px")
              )
            ),
            column(
              width = 6,
              box(
                width  = 12,
                status = "danger",
                plotOutput("leverage_plot", height = "400px")
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              box(
                width  = 12,
                status = "primary",
                title = "Cook's Distance Analysis - Two Cut-off Methods",
                verbatimTextOutput("cooks_summary")
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              box(
                width  = 12,
                status = "danger",
                DTOutput("outliers_table")
              )
            )
          )
        )
      )
    ),
    
    # ====================== TAB: ROBUST METHODS ===============================
    tabItem(
      tabName = "robust",
      
      fluidRow(
        box(
          title = "Optimal Lambda Values",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          fluidRow(
            column(3, valueBoxOutput("lambda_min_box", width = 12)),
            column(3, valueBoxOutput("lambda_1se_box", width = 12)),
            column(3, valueBoxOutput("lambda_robust_min_box", width = 12)),
            column(3, valueBoxOutput("lambda_robust_1se_box", width = 12))
          ),
          fluidRow(
            column(3, valueBoxOutput("optimal_mse_box", width = 12)),
            column(3, valueBoxOutput("optimal_mse_1se_box", width = 12)),
            column(3, valueBoxOutput("optimal_mse_robust_box", width = 12)),
            column(3, valueBoxOutput("optimal_mse_robust_1se_box", width = 12))
          )
        )
      ),
      
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Ridge Regression - CV Plot",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            plotOutput("ridge_cv_plot", height = "400px"),
            hr(),
            h5("Lambda Information - Ridge:"),
            verbatimTextOutput("ridge_lambda_info")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Robust Ridge Regression - CV Plot",
            width       = 12,
            status      = "purple",
            solidHeader = TRUE,
            plotOutput("robust_ridge_cv_plot", height = "400px"),
            hr(),
            h5("Lambda Information - Robust Ridge:"),
            verbatimTextOutput("robust_ridge_lambda_info")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Coefficient Paths (Ridge)",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            plotOutput("ridge_path_plot", height = "400px")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Robust Regression Analysis",
            width       = 12,
            status      = "warning",
            solidHeader = TRUE,
            plotOutput("robust_weights_plot", height = "400px"),
            hr(),
            DTOutput("robust_results")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 6,
          box(
            title       = "Comparison: OLS vs Robust",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            plotOutput("ols_vs_robust_plot", height = "400px")
          )
        ),
        column(
          width = 6,
          box(
            title       = "Robust Ridge Coefficients",
            width       = 12,
            status      = "danger",
            solidHeader = TRUE,
            plotOutput("robust_ridge_plot", height = "400px")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Robust Ridge Regression Results",
            width       = 12,
            status      = "danger",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("robust_ridge_results")
          )
        )
      )
    ),
    
    # ====================== TAB: SIGNIFICANCE TESTS ===========================
    tabItem(
      tabName = "significance",
      h2("Uji Signifikansi Parameter untuk Model Terbaik (OLS dengan Interaksi)"),
      
      fluidRow(
        column(
          width = 6,
          box(
            title = tags$b("Uji F - Overall Model Significance"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            tableOutput("f_test_table"),
            hr(),
            verbatimTextOutput("f_test_interpretation")
          )
        ),
        column(
          width = 6,
          box(
            title = tags$b("Uji t - Signifikansi Koefisien Individual"),
            status = "success",
            solidHeader = TRUE,
            width = 12,
            DTOutput("t_test_table"),
            hr(),
            uiOutput("significance_legend")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 12,
          box(
            title = tags$b("Uji Wald - Perbedaan Antar Klaster"),
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            DTOutput("wald_test_table"),
            hr(),
            verbatimTextOutput("wald_test_summary")
          )
        )
      )
    ),
    
    # ====================== TAB: RESULTS & REPORTS ============================
    tabItem(
      tabName = "results",
      fluidRow(
        column(
          width = 8,
          box(
            title       = "Key Findings & Recommendations",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            uiOutput("key_findings")
          ),
          box(
            title       = "Executive Summary",
            width       = 12,
            status      = "info",
            solidHeader = TRUE,
            uiOutput("executive_summary")
          )
        ),
        column(
          width = 4,
          box(
            title       = "Download Reports",
            width       = 12,
            status      = "success",
            solidHeader = TRUE,
            h4("Export Analysis Results"),
            br(),
            downloadButton(
              "download_report_pdf", "PDF Report",
              class = "btn-primary btn-block",
              style = "margin-bottom: 10px;"
            ),
            downloadButton(
              "download_report_html", "HTML Report",
              class = "btn-info btn-block",
              style = "margin-bottom: 10px;"
            ),
            downloadButton(
              "download_coefficients", "Coefficients (CSV)",
              class = "btn-success btn-block",
              style = "margin-bottom: 10px;"
            ),
            downloadButton(
              "download_diagnostics", "Diagnostics (CSV)",
              class = "btn-warning btn-block",
              style = "margin-bottom: 10px;"
            ),
            downloadButton(
              "download_shapefile_data", "Clustering Data (CSV)",
              class = "btn-primary btn-block",
              style = "margin-bottom: 10px;"
            ),
            hr(),
            h4("Export Plots"),
            br(),
            actionButton(
              "export_plots", "Export All Plots",
              class = "btn-primary btn-block",
              icon  = icon("download")
            ),
            br(),
            textInput(
              "report_title", "Report Title:",
              value = "Spatial Clustering Regression Analysis"
            ),
            textAreaInput(
              "report_notes", "Additional Notes:",
              rows = 4,
              placeholder = "Enter any additional notes..."
            )
          ),
          box(
            title       = "Model Selection",
            width       = 12,
            status      = "warning",
            solidHeader = TRUE,
            h4("Recommended Model:"),
            uiOutput("recommended_model"),
            hr(),
            h4("Model Selection Criteria:"),
            tableOutput("model_selection_criteria")
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          box(
            title       = "Detailed Analysis Report",
            width       = 12,
            status      = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed   = TRUE,
            uiOutput("detailed_report")
          )
        )
      )
    ),
    
    # ====================== TAB: ABOUT ========================================
    tabItem(
      tabName = "about",
      fluidRow(
        column(
          width = 12,
          align = "center",
          box(
            title = "About This Dashboard",
            width = 8,
            status = "primary",
            solidHeader = TRUE,
            tags$div(
              style = "text-align: center;",
              tags$img(
                src = "https://img.icons8.com/color/96/000000/statistics.png",
                style = "margin-bottom: 20px;"
              ),
              h2("Spatial Clustering Regression Dashboard"),
              h4("Professional Analysis Tool for Dummy Variable Regression with Spatial Clustering"),
              hr(),
              p("This dashboard provides comprehensive analysis for regression models with dummy variables as clustering approach, with integrated geographic mapping capabilities."),
              br(),
              h4("Features:"),
              tags$ul(
                style = "text-align: left;",
                tags$li("Multiple estimation methods: OLS, Ridge, Robust, Robust Ridge"),
                tags$li("Comprehensive diagnostic tools with α = 10% significance level"),
                tags$li("Interactive visualizations with geographic mapping"),
                tags$li("Separate data handling for regression analysis and clustering visualization"),
                tags$li("Automated model selection based on performance metrics"),
                tags$li("Exportable reports in multiple formats"),
                tags$li("Professional quality output with significance testing (F, t, Wald)"),
                tags$li("Spatial clustering map visualization with shapefile support")
              ),
              br(),
              h4("Methods Implemented:"),
              tags$ul(
                style = "text-align: left;",
                tags$li("Pooled regression with dummy intercept shift"),
                tags$li("Full interaction model with different slopes"),
                tags$li("Ridge regression for multicollinearity"),
                tags$li("Robust regression for outliers"),
                tags$li("Robust ridge regression for both issues"),
                tags$li("Spatial clustering analysis with polygon maps")
              ),
              br(),
              h4("Developed by:"),
              p("Prof. Dr. I Gede Nyoman Mindra Jaya, M.Si."),
              p("Prof. Yuyun Hidayat, MS., Ph.D."),
              p("Theodora Clara"),
              br(),
              p("Department of Statistics, Universitas Padjadjaran"),
              p("For academic and research purposes"),
              hr(),
              tags$footer(
                style = "font-size: 12px; color: #7f8c8d;",
                "© 2026 Spatial Clustering Regression Dashboard"
              )
            )
          )
        )
      )
    )
  )
)

# ==============================================================================
# 5. SERVER FUNCTION
# ==============================================================================

server <- function(input, output, session) {
  
  values <- reactiveValues(
    regression_data = NULL,
    regression_clean = NULL,
    models = list(),
    diagnostics = list(),
    performance = NULL,
    lambda_min = NULL,
    lambda_1se = NULL,
    lambda_robust_min = NULL,
    lambda_robust_1se = NULL,
    optimal_mse_min = NULL,
    optimal_mse_1se = NULL,
    optimal_mse_robust_min = NULL,
    optimal_mse_robust_1se = NULL,
    ridge_penalty_info = NULL,
    significance_tests = NULL,
    shape_data = NULL
  )
  
  # ======================== REGRESSION DATA HANDLING ==========================
  
  observeEvent(input$file_regression, {
    req(input$file_regression)
    tryCatch({
      data <- read_excel(input$file_regression$datapath)
      values$regression_data <- as.data.frame(data)
      
      for(col in names(values$regression_data)) {
        if(is.character(values$regression_data[[col]])) {
          if(any(grepl(",", values$regression_data[[col]]))) {
            values$regression_data[[col]] <- as.numeric(gsub(",", ".", values$regression_data[[col]]))
          }
        }
      }
      
      exclude_patterns <- c("Wilayah", "Region", "Daerah", "Kota", "Kabupaten", 
                            "Longitude", "Latitude", "lon", "lat", "Coordinate", "Location")
      
      all_vars <- names(values$regression_data)
      filtered_vars <- all_vars[!sapply(all_vars, function(x) any(grepl(paste(exclude_patterns, collapse="|"), x, ignore.case = TRUE)))]
      
      if(length(filtered_vars) == 0) {
        filtered_vars <- all_vars
      }
      
      updateSelectInput(session, "y_var", 
                        choices = filtered_vars,
                        selected = ifelse("Y" %in% filtered_vars, "Y", filtered_vars[1]))
      
      updateSelectInput(session, "d_var", 
                        choices = filtered_vars,
                        selected = ifelse("D" %in% filtered_vars, "D", 
                                          ifelse(length(filtered_vars) > 1, filtered_vars[2], filtered_vars[1])))
      
      x_default <- c()
      for(i in 1:8) {
        if(paste0("X", i) %in% filtered_vars) x_default <- c(x_default, paste0("X", i))
        if(paste0("X", i, "D") %in% filtered_vars) x_default <- c(x_default, paste0("X", i, "D"))
        if(paste0("X", i, "0") %in% filtered_vars) x_default <- c(x_default, paste0("X", i, "0"))
      }
      
      if(length(x_default) == 0 && length(filtered_vars) > 2) {
        x_default <- filtered_vars[3:min(8, length(filtered_vars))]
      }
      
      updatePickerInput(session, "x_vars",
                        choices = filtered_vars,
                        selected = x_default)
      
      showNotification("Regression data loaded successfully!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error loading regression data:", e$message), type = "error", duration = 5)
    })
  })
  
  observeEvent(input$sample_regression_data, {
    set.seed(123)
    n <- 100
    
    sample_data <- data.frame(
      Y = rnorm(n, 50, 10),
      X1 = rnorm(n, 10, 2),
      X2 = rnorm(n, 20, 3),
      X3 = rnorm(n, 30, 4),
      X4 = rnorm(n, 40, 5),
      X5 = rnorm(n, 50, 6),
      X6 = rnorm(n, 60, 7),
      X7 = rnorm(n, 70, 8),
      X8 = rnorm(n, 0.0005, 0.0001),
      D = sample(0:1, n, replace = TRUE)
    )
    
    for(i in 1:8) {
      sample_data[[paste0("X", i, "D")]] <- sample_data[[paste0("X", i)]] * sample_data$D
    }
    
    values$regression_data <- sample_data
    
    updateSelectInput(session, "y_var", 
                      choices = names(values$regression_data), 
                      selected = "Y")
    updateSelectInput(session, "d_var", 
                      choices = names(values$regression_data), 
                      selected = "D")
    
    updatePickerInput(session, "x_vars",
                      choices = names(values$regression_data)[!names(values$regression_data) %in% c("Y", "D")],
                      selected = c(paste0("X", 1:8), paste0("X", 1:8, "D")))
    
    showNotification("Sample regression data loaded successfully!", 
                     type = "message", duration = 3)
  })
  
  # ======================== SHAPEFILE HANDLING =============================
  
  observeEvent(input$shp_file, {
    req(input$shp_file)
    
    tryCatch({
      shp_files <- input$shp_file
      
      shp_exts <- c(".shp", ".dbf", ".shx")
      valid_files <- shp_files[grepl(paste(shp_exts, collapse = "|"), shp_files$name), ]
      
      if(nrow(valid_files) < 3) {
        showNotification(
          "File shapefile tidak lengkap! Butuh minimal: .shp, .shx, .dbf",
          type = "error", duration = 5
        )
        return(NULL)
      }
      
      temp_dir <- tempdir()
      
      for(i in 1:nrow(valid_files)) {
        file.copy(
          valid_files$datapath[i],
          file.path(temp_dir, valid_files$name[i])
        )
      }
      
      shp_name <- valid_files$name[grep("\\.shp$", valid_files$name)]
      
      if(length(shp_name) == 0) {
        stop("File .shp tidak ditemukan")
      }
      
      shp_path <- file.path(temp_dir, shp_name)
      values$shape_data <- st_read(shp_path, quiet = TRUE)
      
      if(is.na(st_crs(values$shape_data))) {
        values$shape_data <- st_set_crs(values$shape_data, 4326)
      }
      
      required_cols <- c("Wlyh_rp", "TngktKm", "Klaster")
      missing_cols <- setdiff(required_cols, names(values$shape_data))
      
      if(length(missing_cols) > 0) {
        showNotification(
          paste("Kolom berikut tidak ditemukan:", paste(missing_cols, collapse = ", ")),
          type = "warning", duration = 5
        )
      }
      
      showNotification("Shapefile loaded successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("Error loading shapefile:", e$message), type = "error", duration = 5)
    })
  })
  
  observeEvent(input$sample_polygon_data, {
    sample_data <- data.frame(
      Wilayah = c("Bogor", "Sukabumi", "Cianjur", "Bandung", "Garut", 
                  "Tasikmalaya", "Ciamis", "Kuningan", "Cirebon", "Majalengka",
                  "Sumedang", "Indramayu", "Subang", "Purwakarta", "Karawang",
                  "Bekasi", "Bandung Barat", "Pangandaran", "Kota Bogor",
                  "Kota Sukabumi", "Kota Bandung", "Kota Cirebon","Kota Bekasi",
                  "Kota Depok", "Kota Cimahi", "Kota Tasikmalaya", "Kota Banjar"),
      Wlyh_rp = c("Kab. Bogor", "Kab. Sukabumi", "Kab. Cianjur", "Kab. Bandung", "Kab. Garut",
                  "Kab. Tasikmalaya", "Kab. Ciamis", "Kab. Kuningan", "Kab. Cirebon", "Kab. Majalengka",
                  "Kab. Sumedang", "Kab. Indramayu", "Kab. Subang", "Kab. Purwakarta", "Kab. Karawang",
                  "Kab. Bekasi", "Kab. Bandung Barat", "Kab. Pangandaran", "Kota Bogor", "Kota Sukabumi",
                  "Kota Bandung", "Kota Cirebon", "Kota Bekasi", "Kota Depok", "Kota Cimahi", "Kota Tasikmalaya", 
                  "Kota Banjar"),
      TngktKm = c(7.05, 6.87, 10.14, 6.19, 9.68, 10.23, 7.39, 11.88, 
                  11.00, 10.82, 9.10, 11.93, 9.49, 8.41, 7.86, 4.80, 
                  10.49, 8.75, 6.53, 7.20, 3.87, 9.02, 4.01, 2.34, 
                  4.39, 11.10, 5.85),
      Klaster = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2)
    )
    
    set.seed(123)
    polygons <- list()
    
    base_lon <- 107.6
    base_lat <- -7.0
    
    for(i in 1:nrow(sample_data)) {
      lon_offset <- runif(1, -0.5, 0.5)
      lat_offset <- runif(1, -0.3, 0.3)
      
      lon_min <- base_lon + lon_offset - 0.05
      lon_max <- base_lon + lon_offset + 0.05
      lat_min <- base_lat + lat_offset - 0.03
      lat_max <- base_lat + lat_offset + 0.03
      
      poly <- st_polygon(list(matrix(
        c(lon_min, lat_min,
          lon_max, lat_min,
          lon_max, lat_max,
          lon_min, lat_max,
          lon_min, lat_min),
        ncol = 2, byrow = TRUE
      )))
      
      polygons[[i]] <- st_sfc(poly, crs = 4326)
    }
    
    geometry <- do.call(c, polygons)
    
    values$shape_data <- st_sf(
      sample_data,
      geometry = geometry
    )
    
    showNotification("Sample polygon data loaded successfully!", 
                     type = "message", duration = 3)
  })
  
  # ======================== REGRESSION DATA OUTPUTS ===========================
  
  output$regression_data_preview <- renderDT({
    req(values$regression_data)
    datatable(
      values$regression_data,
      extensions = c("Buttons", "Scroller"),
      options = list(
        dom      = "Bfrtip",
        buttons  = c("copy", "csv", "excel", "pdf"),
        scrollX  = TRUE,
        scrollY  = 400,
        scroller = TRUE,
        pageLength = 10
      ),
      class = 'display compact stripe hover'
    )
  })
  
  output$regression_structure <- renderPrint({
    req(values$regression_data)
    cat("REGRESSION DATA STRUCTURE\n")
    cat("=========================\n\n")
    cat("Dimensions:", nrow(values$regression_data), "rows ×", ncol(values$regression_data), "columns\n\n")
    cat("Variable Types:\n")
    str(values$regression_data, vec.len = 1)
  })
  
  output$var_types_plot <- renderPlotly({
    req(values$regression_data)
    var_types  <- sapply(values$regression_data, function(x) class(x)[1])
    var_counts <- table(var_types)
    plot_ly(
      x    = names(var_counts),
      y    = as.numeric(var_counts),
      type = "bar",
      marker = list(color = c("#3498db", "#2ecc71", "#e74c3c", "#f39c12")),
      text         = as.numeric(var_counts),
      textposition = "auto"
    ) %>%
      layout(
        title = "Variable Types Distribution",
        xaxis = list(title = "Data Type"),
        yaxis = list(title = "Count"),
        hoverlabel = list(font = list(size = 14))
      )
  })
  
  output$n_obs <- renderValueBox({
    req(values$regression_data)
    valueBox(value = nrow(values$regression_data), subtitle = "Observations",
             icon = icon("database"), color = "blue")
  })
  
  output$n_vars <- renderValueBox({
    req(values$regression_data)
    valueBox(value = ncol(values$regression_data), subtitle = "Variables",
             icon = icon("columns"), color = "green")
  })
  
  output$n_missing <- renderValueBox({
    req(values$regression_data)
    missing_count <- sum(is.na(values$regression_data))
    valueBox(
      value    = missing_count,
      subtitle = "Missing Values",
      icon     = icon("exclamation-triangle"),
      color    = ifelse(missing_count > 0, "red", "green")
    )
  })
  
  output$regression_data_info <- renderUI({
    req(values$regression_data)
    tagList(
      h5("Regression Data Information:"),
      tags$ul(
        tags$li(paste("Rows:", nrow(values$regression_data))),
        tags$li(paste("Columns:", ncol(values$regression_data))),
        tags$li(paste("Numeric:", sum(sapply(values$regression_data, is.numeric)))),
        tags$li(paste("Missing:", sum(is.na(values$regression_data))))
      )
    )
  })
  
  output$regression_summary <- renderPrint({
    req(values$regression_data, input$y_var, input$d_var, input$x_vars)
    cat("Regression Model Specification:\n")
    cat("  Y :", input$y_var, "\n")
    cat("  D :", input$d_var, "\n")
    cat("  X :", paste(input$x_vars, collapse = ", "), "\n")
    cat("Methods :", paste(input$methods, collapse = ", "), "\n")
  })
  
  # ======================== SHAPEFILE OUTPUTS ================================
  
  output$shp_attr_preview <- renderDT({
    req(values$shape_data)
    
    attr_data <- st_drop_geometry(values$shape_data)
    
    main_cols <- c("Wlyh_rp", "Klaster", "TngktKm")
    other_cols <- setdiff(names(attr_data), main_cols)
    display_cols <- c(main_cols, other_cols[1:min(3, length(other_cols))])
    
    datatable(
      attr_data[, display_cols, drop = FALSE],
      options = list(
        dom = "t",
        scrollX = TRUE,
        scrollY = 200,
        paging = FALSE
      ),
      rownames = FALSE
    )
  })
  
  # ======================== SHAPEFILE MAP OUTPUT ========================
  
  output$cluster_polygon_map <- renderLeaflet({
    req(values$shape_data)
    
    tryCatch({
      data <- values$shape_data
      
      # Transformasi CRS ke 4326 jika perlu
      current_crs <- st_crs(data)$epsg
      if(!is.na(current_crs) && current_crs != 4326) {
        data <- st_transform(data, 4326)
      }
      
      # Pastikan Klaster sebagai factor
      if("Klaster" %in% names(data)) {
        data$Klaster <- as.factor(data$Klaster)
      }
      
      # --- Color Palette ---
      if(input$map_color_var == "Klaster") {
        color_pal <- colorFactor(
          palette = c("cyan", "pink"),  
          domain = data$Klaster,
          na.color = "grey"
        )
        fill_color <- color_pal(data$Klaster)
      } else if(input$map_color_var == "TngktKm") {
        color_pal <- colorNumeric(
          palette = c("yellow", "orange", "red"),
          domain = data$TngktKm,
          na.color = "grey"
        )
        fill_color <- color_pal(data$TngktKm)
      }
      
      # --- LABEL HOVER (muncul saat kursor diarahkan) ---
      label_content <- sprintf(
        "<strong>%s</strong><br/>
       Klaster: %s<br/>
       Tingkat Kemiskinan: %.2f",
        data$Wlyh_rp,
        as.character(data$Klaster),
        data$TngktKm
      ) %>% lapply(htmltools::HTML)
      
      # --- Popup CONTENT (saat diklik) ---
      popup_content <- sprintf(
        "<div style='font-family: Arial, sans-serif; min-width: 220px;'>
         <h4 style='margin-top: 0; color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 5px;'>%s</h4>
         <table style='width: 100%%; border-collapse: collapse;'>
           <tr><td style='padding: 5px; font-weight: bold;'>Wilayah:</td><td style='padding: 5px;'>%s</td></tr>
           <tr><td style='padding: 5px; font-weight: bold; background-color: #f8f9fa;'>Klaster:</td>
               <td style='padding: 5px; background-color: #f8f9fa;'><span style='font-weight: bold; color: %s;'>%s</span></td></tr>
           <tr><td style='padding: 5px; font-weight: bold;'>Tingkat Kemiskinan:</td><td style='padding: 5px;'>%.2f</td></tr>
         </table>
       </div>",
        data$Wlyh_rp,                                  
        data$Wlyh_rp,                                   
        ifelse(data$Klaster == 1, "cyan", "pink"), 
        as.character(data$Klaster),                    
        data$TngktKm                                   
      )
      
      # --- Leaflet Map (dengan label hover) ---
      map <- leaflet(data) %>%
        addProviderTiles(providers$CartoDB.Positron) %>%
        addPolygons(
          fillColor = fill_color,
          fillOpacity = 0.7,
          color = "white",
          weight = 1,
          highlightOptions = highlightOptions(
            color = "black",
            weight = 3,
            fillOpacity = 0.9,
            bringToFront = TRUE
          ),
          label = label_content,
          labelOptions = labelOptions(
            style = list(
              "font-weight" = "normal",
              "padding" = "6px 10px",
              "border-radius" = "3px",
              "border" = "1px solid #ccc",
              "background-color" = "white",
              "box-shadow" = "0 1px 5px rgba(0,0,0,0.2)"
            ),
            textsize = "12px",
            direction = "auto",
            opacity = 0.9
          ),
          popup = popup_content,
          popupOptions = popupOptions(
            maxWidth = 300,
            minWidth = 220,
            closeButton = TRUE,
            closeOnClick = TRUE
          )
        ) %>%
        setView(lng = 107.6, lat = -7.0, zoom = input$map_zoom)
      
      # --- Legend ---
      if(input$map_color_var == "Klaster") {
        map <- map %>%
          addLegend(
            position = "bottomright",
            pal = color_pal,
            values = ~Klaster,
            title = "Klaster",
            opacity = 0.8,
            labFormat = labelFormat(prefix = "Klaster ")
          )
      } else if(input$map_color_var == "TngktKm") {
        map <- map %>%
          addLegend(
            position = "bottomright",
            pal = color_pal,
            values = ~TngktKm,
            title = "Tingkat Kemiskinan",
            opacity = 0.8,
            labFormat = labelFormat(digits = 2)
          )
      }
      
      return(map)
      
    }, error = function(e) {
      # Error handling
      leaflet() %>%
        addTiles() %>%
        setView(lng = 107.6, lat = -7.0, zoom = 8) %>%
        addControl(
          html = paste0(
            "<div style='padding: 10px; background-color: white; border-radius: 5px;'>
             <h4 style='color: red;'>⚠️ Error Displaying Map</h4>
             <p><b>Error:</b> ", e$message, "</p>
             <p><b>Data Info:</b> ", 
            if(!is.null(values$shape_data)) {
              paste("Features:", nrow(values$shape_data))
            } else {
              "No data loaded"
            }, "</p>
           </div>"
          ),
          position = "topright"
        )
    })
  })
  
  # ======================== CLUSTERING SUMMARY ========================
  
  output$clustering_summary_poly <- renderTable({
    req(values$shape_data)
    
    data_df <- st_drop_geometry(values$shape_data)
    
    if(!"Klaster" %in% names(data_df)) {
      return(data.frame(Message = "Kolom 'Klaster' tidak ditemukan dalam data"))
    }
    
    if(!"TngktKm" %in% names(data_df)) {
      return(data.frame(Message = "Kolom 'TngktKm' tidak ditemukan dalam data"))
    }
    
    summary_df <- data_df %>%
      mutate(Klaster = as.integer(Klaster)) %>%
      group_by(Klaster) %>%
      summarise(
        `Jumlah Daerah` = n(),
        `Rata-rata Kemiskinan` = round(mean(TngktKm, na.rm = TRUE), 2),
        `Min Kemiskinan` = round(min(TngktKm, na.rm = TRUE), 2),
        `Max Kemiskinan` = round(max(TngktKm, na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      arrange(Klaster)
    
    return(summary_df)
  }, striped = TRUE, hover = TRUE, bordered = TRUE, align = 'c')
  
  # ======================== REGRESSION ANALYSIS ===============================
  
  observeEvent(input$run_regression, {
    req(values$regression_data, input$y_var, input$d_var, input$x_vars)
    
    showNotification("Running regression analysis...", type = "warning", duration = 2)
    
    tryCatch({
      selected_vars <- unique(c(input$y_var, input$d_var, input$x_vars))
      data_clean <- values$regression_data %>%
        dplyr::select(all_of(selected_vars)) %>%
        tidyr::drop_na()
      
      values$regression_clean <- data_clean
      
      if(nrow(data_clean) < 10) {
        stop("Not enough complete observations (minimum 10 required)")
      }
      
      x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
      x_vars_interaction <- input$x_vars[grepl(paste0(input$d_var, "$"), input$x_vars)]
      
      formula_pooled <- as.formula(
        paste(input$y_var, "~",
              paste(x_vars_main, collapse = " + "),
              "+", input$d_var)
      )
      
      formula_interaction <- as.formula(
        paste(input$y_var, "~",
              paste(x_vars_main, collapse = " + "),
              "+", input$d_var,
              "+", paste(x_vars_interaction, collapse = " + "))
      )
      
      values$models <- list()
      
      if ("ols" %in% input$methods) {
        values$models$ols_pooled <- lm(formula_pooled, data = data_clean)
        values$models$ols_interaction <- lm(formula_interaction, data = data_clean)
      }
      
      if ("ridge" %in% input$methods) {
        set.seed(123)
        
        X_main <- as.matrix(data_clean[, x_vars_main, drop = FALSE])
        D_vector <- as.matrix(data_clean[[input$d_var]])
        X_interaction <- as.matrix(data_clean[, x_vars_interaction, drop = FALSE])
        
        X_full <- cbind(X_main, D = D_vector, X_interaction)
        y_vector <- data_clean[[input$y_var]]
        
        n_x_main <- length(x_vars_main)
        n_d <- 1
        n_x_interaction <- length(x_vars_interaction)
        
        penalty_vector <- c(
          rep(1, n_x_main),
          0,
          rep(1, n_x_interaction)
        )
        
        colnames(X_full) <- c(x_vars_main, input$d_var, x_vars_interaction)
        
        tryCatch({
          cv_fit <- cv.glmnet(
            X_full, 
            y_vector, 
            alpha = 0, 
            penalty.factor = penalty_vector,
            nfolds = input$cv_folds,
            intercept = TRUE
          )
          
          ridge_optimal <- glmnet(
            X_full, 
            y_vector, 
            alpha = 0, 
            penalty.factor = penalty_vector,
            lambda = cv_fit$lambda.min,
            intercept = TRUE
          )
          
          ridge_full <- glmnet(
            X_full, 
            y_vector, 
            alpha = 0, 
            penalty.factor = penalty_vector,
            intercept = TRUE
          )
          
          values$models$ridge_cv <- cv_fit
          values$models$ridge_full <- ridge_full
          values$models$ridge <- ridge_optimal
          
          values$lambda_min <- cv_fit$lambda.min
          values$lambda_1se <- cv_fit$lambda.1se
          values$optimal_mse_min <- min(cv_fit$cvm)
          values$optimal_mse_1se <- cv_fit$cvm[which(cv_fit$lambda == cv_fit$lambda.1se)]
          
          values$ridge_penalty_info <- data.frame(
            Variable = c("(Intercept)", colnames(X_full)),
            Penalty = c(0, penalty_vector)
          )
          
        }, error = function(e) {
          showNotification(
            paste("Ridge Regression Error:", e$message),
            type = "error", duration = 5
          )
        })
      }
      
      if ("robust" %in% input$methods) {
        tryCatch({
          values$models$robust <- rlm(
            formula_interaction,
            data   = data_clean,
            method = input$robust_method,
            maxit  = 110
          )
        }, error = function(e) {
          showNotification(
            paste("Robust regression failed:", e$message),
            type = "warning", duration = 5
          )
        })
      }
      
      if ("robust_ridge" %in% input$methods) {
        tryCatch({
          X_main <- as.matrix(data_clean[, x_vars_main, drop = FALSE])
          D_vector <- as.matrix(data_clean[[input$d_var]])
          X_interaction <- as.matrix(data_clean[, x_vars_interaction, drop = FALSE])
          
          X_full <- cbind(X_main, D = D_vector, X_interaction)
          y_vector <- data_clean[[input$y_var]]
          
          penalty_vector <- c(
            rep(1, n_x_main),
            0,
            rep(1, n_x_interaction)
          )
          
          colnames(X_full) <- c(x_vars_main, input$d_var, x_vars_interaction)
          
          robust_fit_for_weights <- tryCatch({
            rlm(
              formula_interaction,  
              data   = data_clean,
              method = input$robust_method,
              maxit  = 110
            )
          }, error = function(e) {
            list(w = rep(1, nrow(data_clean)))
          })
          
          weights <- robust_fit_for_weights$w
          
          set.seed(123)
          cv_robust_ridge <- cv.glmnet(
            X_full, 
            y_vector,
            alpha    = 0,
            weights  = weights,
            penalty.factor = penalty_vector,
            nfolds   = input$cv_folds,
            type.measure = "mse",
            intercept = TRUE
          )
          
          values$lambda_robust_min <- cv_robust_ridge$lambda.min
          values$lambda_robust_1se <- cv_robust_ridge$lambda.1se
          values$optimal_mse_robust_min <- min(cv_robust_ridge$cvm)
          values$optimal_mse_robust_1se <- cv_robust_ridge$cvm[which(cv_robust_ridge$lambda == cv_robust_ridge$lambda.1se)]
          
          values$models$robust_ridge <- glmnet(
            X_full, 
            y_vector,
            alpha   = 0,
            lambda  = values$lambda_robust_min,
            weights = weights,
            penalty.factor = penalty_vector,
            intercept = TRUE
          )
          
          values$models$robust_ridge_cv <- cv_robust_ridge
          
        }, error = function(e) {
          showNotification(
            paste("Robust Ridge Regression Error:", e$message),
            type = "error", duration = 5
          )
        })
      }
      
      calculate_performance()
      
      if (!is.null(values$models$ols_interaction)) {
        calculate_diagnostics()
        
        values$significance_tests <- perform_significance_tests(
          values$models$ols_interaction,
          x_vars_main,
          input$d_var,
          x_vars_interaction
        )
      }
      
      showNotification("Regression analysis completed successfully!",
                       type = "message", duration = 3)
      
      updateTabItems(session, "tabs", "estimation")
      
    }, error = function(e) {
      showNotification(paste("Error in regression analysis:", e$message),
                       type = "error", duration = 5)
    })
  })
  
  # ======================== HELPER FUNCTIONS FOR REGRESSION ========================
  
  calculate_performance <- function() {
    req(values$regression_clean)
    performance_list <- list()
    
    if (!is.null(values$models$ols_pooled)) {
      model <- values$models$ols_pooled
      performance_list$OLS_Pooled <- data.frame(
        R2     = summary(model)$r.squared,
        Adj_R2 = summary(model)$adj.r.squared,
        RMSE   = sqrt(mean(residuals(model)^2)),
        MAE    = mean(abs(residuals(model))),
        AIC    = AIC(model),
        BIC    = BIC(model)
      )
    }
    
    if (!is.null(values$models$ols_interaction)) {
      model <- values$models$ols_interaction
      performance_list$OLS_Interaction <- data.frame(
        R2     = summary(model)$r.squared,
        Adj_R2 = summary(model)$adj.r.squared,
        RMSE   = sqrt(mean(residuals(model)^2)),
        MAE    = mean(abs(residuals(model))),
        AIC    = AIC(model),
        BIC    = BIC(model)
      )
    }
    
    if (!is.null(values$models$ridge)) {
      tryCatch({
        model  <- values$models$ridge
        actual <- values$regression_clean[[input$y_var]]
        
        x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
        x_vars_interaction <- input$x_vars[grepl(paste0(input$d_var, "$"), input$x_vars)]
        
        X_main <- as.matrix(values$regression_clean[, x_vars_main, drop = FALSE])
        X_D <- as.matrix(values$regression_clean[[input$d_var]])
        X_interaction <- as.matrix(values$regression_clean[, x_vars_interaction, drop = FALSE])
        
        X_full <- cbind(1, X_main, D = X_D, X_interaction)
        colnames(X_full)[1] <- "(Intercept)"
        
        pred   <- as.numeric(predict(model, newx = X_full[, -1, drop = FALSE]))
        
        performance_list$Ridge <- data.frame(
          R2     = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2),
          Adj_R2 = NA,
          RMSE   = sqrt(mean((actual - pred)^2)),
          MAE    = mean(abs(actual - pred)),
          AIC    = NA,
          BIC    = NA
        )
      }, error = function(e) {
        performance_list$Ridge <- NULL
      })
    }
    
    if (!is.null(values$models$robust)) {
      tryCatch({
        model <- values$models$robust
        performance_list$Robust <- data.frame(
          R2     = NA,
          Adj_R2 = NA,
          RMSE   = sqrt(mean(residuals(model)^2)),
          MAE    = mean(abs(residuals(model))),
          AIC    = NA,
          BIC    = NA
        )
      }, error = function(e) {
        performance_list$Robust <- NULL
      })
    }
    
    if (!is.null(values$models$robust_ridge)) {
      tryCatch({
        model  <- values$models$robust_ridge
        actual <- values$regression_clean[[input$y_var]]
        
        x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
        x_vars_interaction <- input$x_vars[grepl(paste0(input$d_var, "$"), input$x_vars)]
        
        X_main <- as.matrix(values$regression_clean[, x_vars_main, drop = FALSE])
        X_D <- as.matrix(values$regression_clean[[input$d_var]])
        X_interaction <- as.matrix(values$regression_clean[, x_vars_interaction, drop = FALSE])
        
        X_full <- cbind(1, X_main, D = X_D, X_interaction)
        colnames(X_full)[1] <- "(Intercept)"
        
        pred   <- as.numeric(predict(model, newx = X_full[, -1, drop = FALSE]))
        
        performance_list$Robust_Ridge <- data.frame(
          R2     = 1 - sum((actual - pred)^2) / sum((actual - mean(actual))^2),
          Adj_R2 = NA,
          RMSE   = sqrt(mean((actual - pred)^2)),
          MAE    = mean(abs(actual - pred)),
          AIC    = NA,
          BIC    = NA
        )
      }, error = function(e) {
        performance_list$Robust_Ridge <- NULL
      })
    }
    
    if (length(performance_list) > 0) {
      values$performance <- dplyr::bind_rows(performance_list, .id = "Model") %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4)))
    } else {
      values$performance <- NULL
    }
  }
  
  calculate_diagnostics <- function() {
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    
    if (length(input$x_vars) > 1) {
      values$diagnostics$vif <- tryCatch({
        car::vif(model)
      }, error = function(e) {
        NULL
      })
    } else {
      values$diagnostics$vif <- NULL
    }
    
    values$diagnostics$ks <- tryCatch(ks.test(residuals(model), "pnorm"), error = function(e) NULL)
    
    values$diagnostics$bp <- tryCatch(bptest(model), error = function(e) NULL)
    
    values$diagnostics$dw <- tryCatch(durbinWatsonTest(model), error = function(e) NULL)
    
    cooks_dist <- cooks.distance(model)
    values$diagnostics$cooks <- cooks_dist
    
    n <- nrow(model$model)
    p <- length(coef(model))
    
    df1 <- p + 1
    df2 <- n - (p + 1)
    cut_off_f <- ifelse(df2 > 0, qf(0.5, df1 = df1, df2 = df2), 1.0)
    cut_off_fixed <- 1.0
    
    values$diagnostics$cut_off_f <- cut_off_f
    values$diagnostics$cut_off_fixed <- cut_off_fixed
    
    values$diagnostics$outliers_f <- which(cooks_dist > cut_off_f)
    values$diagnostics$outliers_fixed <- which(cooks_dist > cut_off_fixed)
    
    values$diagnostics$outliers_all <- unique(
      c(values$diagnostics$outliers_f,
        values$diagnostics$outliers_fixed)
    )
  }
  
  # ======================== REGRESSION OUTPUTS ================================
  
  output$cluster_dist_plot <- renderPlotly({
    req(values$regression_data, input$d_var)
    
    if(input$d_var %in% names(values$regression_data)) {
      cluster_data <- values$regression_data %>%
        group_by(.data[[input$d_var]]) %>%
        summarise(Count = n(), .groups = "drop")
      
      plot_ly(
        cluster_data,
        x = ~as.factor(.data[[input$d_var]]),
        y = ~Count,
        type = "bar",
        marker = list(color = c("#3498db", "#e74c3c")),
        text = ~Count,
        textposition = "auto"
      ) %>%
        layout(
          title = "Distribution of Cluster Variable (D)",
          xaxis = list(title = paste("Cluster", input$d_var)),
          yaxis = list(title = "Count"),
          showlegend = FALSE
        )
    }
  })
  
  # ======================== MATHEMATICAL SPECIFICATION ========================
  
  output$math_spec <- renderUI({
    withMathJax(
      tags$div(
        style = "line-height: 1.6; font-size: 14px; padding: 15px; background-color: #f8f9fa; border-radius: 5px;",
        
        # Model 1 - Pooled Regression
        tags$p(style = "font-weight: bold; margin-top: 0;", "Model 1: Pooled Regression (Same Slopes)"),
        tags$div(
          style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #3498db; font-size: 15px; overflow-x: auto; white-space: nowrap; margin-bottom: 5px;",
          tags$p(style = "margin: 0; font-family: 'Courier New', monospace; display: inline-block;",
                 HTML("$$Y_i = \\beta_0 + \\sum_{m=1}^{M} \\beta_m X_{mi} + \\sum_{j=1}^{k-1} \\gamma_j D_{ji} + \\varepsilon_i$$")
          )
        ),
        tags$p(style = "font-style: italic; color: #555; margin-top: 2px; margin-bottom: 15px;", 
               "Asumsi: Slope sama untuk semua klaster, hanya intercept yang berbeda"),
        
        # Model 2 - Full Interaction
        tags$p(style = "font-weight: bold;", "Model 2: Full Interaction (Different Slopes)"),
        tags$div(
          style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #e74c3c; font-size: 15px; overflow-x: auto; white-space: nowrap; margin-bottom: 5px;",
          tags$p(style = "margin: 0; font-family: 'Courier New', monospace; display: inline-block;",
                 HTML("$$Y_i = \\beta_0 + \\sum_{m=1}^{M} \\beta_m X_{mi} + \\sum_{j=1}^{k-1} \\gamma_j D_{ji} + \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj} (X_{mi} D_{ji}) + \\varepsilon_i$$")
          )
        ),
        tags$p(style = "font-style: italic; color: #555; margin-top: 2px; margin-bottom: 15px;", 
               "Asumsi: Slope dan intercept berbeda antar klaster"),
        
        # Model 3 - Ridge Regression (SCROLLABLE)
        tags$p(style = "font-weight: bold;", "Model 3: Ridge Regression"),
        tags$div(
          style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #f39c12; font-size: 15px; overflow-x: auto; white-space: nowrap; margin-bottom: 5px;",
          tags$p(style = "margin: 0; font-family: 'Courier New', monospace; display: inline-block;",
                 HTML("$$\\hat{\\beta}^{\\text{ridge}} = \\arg\\min_{\\beta} \\left\\{ \\sum_{i=1}^{n} \\bigl(y_i - \\beta_0 - \\sum_{m=1}^{M} \\beta_m X_{mi} - \\sum_{j=1}^{k-1} \\gamma_j D_{ji} - \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj} (X_{mi} D_{ji}) \\bigr)^2 + \\lambda \\left( \\sum_{m=1}^{M} \\beta_m^2 + \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj}^2 \\right) \\right\\}$$")
          )
        ),
        
        # Model 4 - Robust Regression (SCROLLABLE)
        tags$p(style = "font-weight: bold;", "Model 4: Robust Regression"),
        tags$div(
          style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #2ecc71; font-size: 15px; overflow-x: auto; white-space: nowrap; margin-bottom: 5px;",
          tags$p(style = "margin: 0; font-family: 'Courier New', monospace; display: inline-block;",
                 HTML("$$\\hat{\\beta}^{\\text{robust}} = \\arg\\min_{\\beta} \\sum_{i=1}^{n} \\rho\\left( \\frac{y_i - \\beta_0 - \\sum_{m=1}^{M} \\beta_m X_{mi} - \\sum_{j=1}^{k-1} \\gamma_j D_{ji} - \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj} (X_{mi} D_{ji})}{\\sigma} \\right)$$")
          )
        ),
        
        # Model 5 - Robust Ridge Regression (SCROLLABLE)
        tags$p(style = "font-weight: bold;", "Model 5: Robust Ridge Regression"),
        tags$div(
          style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; border-left: 4px solid #9b59b6; font-size: 15px; overflow-x: auto; white-space: nowrap; margin-bottom: 5px;",
          tags$p(style = "margin: 0; font-family: 'Courier New', monospace; display: inline-block;",
                 HTML("$$\\hat{\\beta}^{\\text{robust-ridge}} = \\arg\\min_{\\beta} \\left\\{ \\sum_{i=1}^{n} w_i \\bigl(y_i - \\beta_0 - \\sum_{m=1}^{M} \\beta_m X_{mi} - \\sum_{j=1}^{k-1} \\gamma_j D_{ji} - \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj} (X_{mi} D_{ji}) \\bigr)^2 + \\lambda \\left( \\sum_{m=1}^{M} \\beta_m^2 + \\sum_{m=1}^{M} \\sum_{j=1}^{k-1} \\delta_{mj}^2 \\right) \\right\\}$$")
          )
        ),
        tags$p(style = "font-style: italic; color: #555; margin-top: 2px; margin-bottom: 15px;", 
               "Kombinasi robust weights (wᵢ) dan ridge penalty (λ) untuk menangani outlier dan multikolinearitas"),
        
        tags$br(),
        tags$hr(),
        
        # ============ NOTASI ============
        tags$p(style = "font-weight: bold; font-size: 16px; margin-bottom: 10px;", "Notasi:"),
        
        tags$div(
          style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px;",
          
          # Kolom Kiri
          tags$div(
            style = "padding-right: 10px;",
            tags$ul(
              style = "list-style-type: disc; padding-left: 20px; margin-top: 0;",
              tags$li(style = "margin-bottom: 8px;", "\\(Y_i\\) : Variabel dependen pada pengamatan ke-\\(i\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\beta_0\\) : Intercept dari model"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\beta_m\\) : Koefisien regresi untuk variabel independen kuantitatif ke-\\(m\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(X_{mi}\\) : Variabel independen kuantitatif ke-\\(m\\) pada pengamatan ke-\\(i\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\gamma_j\\) : Koefisien regresi untuk variabel dummy ke-\\(j\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(D_{ji}\\) : Variabel dummy ke-\\(j\\) pada pengamatan ke-\\(i\\) (dengan \\(D = 0, 1\\))"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\delta_{mj}\\) : Koefisien regresi untuk interaksi variabel independen ke-\\(m\\) dan dummy ke-\\(j\\)")
            )
          ),
          
          # Kolom Kanan
          tags$div(
            style = "padding-left: 10px;",
            tags$ul(
              style = "list-style-type: disc; padding-left: 20px; margin-top: 0;",
              tags$li(style = "margin-bottom: 8px;", "\\(X_{mi} D_{ji}\\) : Interaksi variabel independen ke-\\(m\\) dengan dummy ke-\\(j\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\lambda\\) : Regularization parameter (optimal dari cross-validation)"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\rho(\\cdot)\\) : Robust loss function"),
              tags$li(style = "margin-bottom: 8px;", "\\(w_i\\) : Robust weights"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\sigma\\) : Scale parameter"),
              tags$li(style = "margin-bottom: 8px;", "\\(\\varepsilon_i\\) : Error term / kekeliruan untuk pengamatan ke-\\(i\\)"),
              tags$li(style = "margin-bottom: 8px;", "\\(k\\) : Jumlah klaster")
            )
          )
        )
      )
    )
  })
  
  output$model_comparison_plot <- renderPlotly({
    req(values$performance)
    perf_long <- values$performance %>%
      dplyr::select(Model, R2, RMSE, MAE) %>%
      tidyr::pivot_longer(cols = c(R2, RMSE, MAE),
                          names_to = "Metric",
                          values_to = "Value")
    p <- ggplot(perf_long, aes(x = Model, y = Value, fill = Metric)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("#3498db", "#e74c3c", "#2ecc71")) +
      labs(title = "Model Performance Comparison",
           x = "Model", y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p)
  })
  
  output$ols_results <- renderDT({
    req(values$models$ols_interaction)
    
    ols_summary <- broom::tidy(values$models$ols_interaction) %>%
      mutate(
        estimate = sapply(1:n(), function(i) format_coef(estimate[i], term[i])),
        std.error = sapply(1:n(), function(i) format_coef(std.error[i], term[i])),
        statistic = sapply(1:n(), function(i) format_coef(statistic[i], term[i])),
        p.value = format.pval(p.value, digits = 4),
        Significance = dplyr::case_when(
          as.numeric(p.value) < 0.01  ~ "***",
          as.numeric(p.value) < 0.05  ~ "**",
          as.numeric(p.value) < 0.1   ~ "*",
          TRUE                        ~ ""
        )
      )
    
    dt_scroller(ols_summary, height = 300) %>%
      formatStyle(
        "p.value",
        backgroundColor = styleInterval(
          c(0.01, 0.05, 0.1),
          c("#27ae60", "#2ecc71", "#f1c40f", "#e74c3c")
        )
      )
  })
  
  output$ridge_results <- renderDT({
    req(values$models$ridge, input$d_var)
    
    ridge_coef <- show_ridge_coef(values$models$ridge, input$d_var)
    
    datatable(
      ridge_coef,
      options = list(
        dom = "t",
        scrollX = TRUE,
        scrollY = 300,
        paging = FALSE
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Type',
        backgroundColor = styleEqual(
          c("Intercept (No Penalty)", "D (No Penalty)", "X - Regularized", 
            "Interaction (X*D) - Regularized", "Error"),
          c('#f8f9fa', '#d4edda', '#fff3cd', '#d1ecf1', '#ffcccc')
        )
      )
  })
  
  output$robust_results <- renderDT({
    req(values$models$robust)
    
    tryCatch({
      robust_summary <- broom::tidy(values$models$robust) %>%
        mutate(
          estimate = sapply(1:n(), function(i) format_coef(estimate[i], term[i])),
          std.error = sapply(1:n(), function(i) format_coef(std.error[i], term[i])),
          statistic = sapply(1:n(), function(i) format_coef(statistic[i], term[i]))
        )
      
      dt_scroller(robust_summary, height = 300)
    }, error = function(e) {
      warning_df <- data.frame(
        Warning = "Perfect multicollinearity detected! Cannot estimate Robust Regression with X, D, and X*D together.",
        Recommendation = "Use Ridge Robust Regression instead."
      )
      dt_scroller(warning_df, height = 100)
    })
  })
  
  output$robust_ridge_results <- renderDT({
    req(values$models$robust_ridge, input$d_var)
    
    ridge_coef <- show_ridge_coef(values$models$robust_ridge, input$d_var)
    
    datatable(
      ridge_coef,
      options = list(
        dom = "t",
        scrollX = TRUE,
        scrollY = 300,
        paging = FALSE
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Type',
        backgroundColor = styleEqual(
          c("Intercept (No Penalty)", "D (No Penalty)", "X - Regularized", 
            "Interaction (X*D) - Regularized", "Error"),
          c('#f8f9fa', '#d4edda', '#fff3cd', '#d1ecf1', '#ffcccc')
        )
      )
  })
  
  output$coef_comparison <- renderPlotly({
    req(values$models$ols_interaction)
    ols_coef <- broom::tidy(values$models$ols_interaction) %>%
      dplyr::filter(term != "(Intercept)") %>%
      dplyr::mutate(Method = "OLS") %>%
      dplyr::select(term, estimate, Method)
    
    if (!is.null(values$models$ridge)) {
      tryCatch({
        coef_mat <- as.matrix(coef(values$models$ridge))
        ridge_coef <- data.frame(
          term     = rownames(coef_mat),
          estimate = as.numeric(coef_mat[, 1]),
          row.names = NULL
        ) %>%
          dplyr::filter(term != "(Intercept)") %>%
          dplyr::mutate(Method = "Ridge")
        all_coefs <- dplyr::bind_rows(ols_coef, ridge_coef)
      }, error = function(e) {
        all_coefs <- ols_coef
      })
    } else {
      all_coefs <- ols_coef
    }
    
    p <- ggplot(all_coefs, aes(x = term, y = estimate, fill = Method)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      coord_flip() +
      scale_fill_manual(values = c("#3498db", "#e74c3c")) +
      labs(title = "Coefficient Comparison",
           x = "Variable", y = "Coefficient Value") +
      theme_minimal() +
      theme(legend.position = "bottom")
    ggplotly(p)
  })
  
  output$performance_table <- renderDT({
    req(values$performance)
    num_cols <- names(values$performance)[sapply(values$performance, is.numeric)]
    datatable(
      values$performance,
      options  = list(pageLength = 10, dom = 't'),
      rownames = FALSE
    ) %>%
      formatRound(columns = num_cols, digits = 4)
  })
  
  # ======================== EXPLORATORY OUTPUTS ===============================
  
  output$desc_stats <- renderDT({
    req(values$regression_data)
    numeric_data <- values$regression_data[, sapply(values$regression_data, is.numeric), drop = FALSE]
    if (ncol(numeric_data) > 0) {
      stats <- data.frame(
        Variable = names(numeric_data),
        Mean     = round(colMeans(numeric_data, na.rm = TRUE), 4),
        SD       = round(apply(numeric_data, 2, sd, na.rm = TRUE), 4),
        Min      = round(apply(numeric_data, 2, min, na.rm = TRUE), 4),
        Q1       = round(apply(numeric_data, 2, quantile, 0.25, na.rm = TRUE), 4),
        Median   = round(apply(numeric_data, 2, median, na.rm = TRUE), 4),
        Q3       = round(apply(numeric_data, 2, quantile, 0.75, na.rm = TRUE), 4),
        Max      = round(apply(numeric_data, 2, max, na.rm = TRUE), 4),
        Missing  = colSums(is.na(numeric_data))
      )
      dt_scroller(stats, height = 400) %>%
        formatStyle(
          'Missing',
          backgroundColor = styleInterval(
            c(0, 1),
            c('white', '#FFE5E5', '#FF9999')
          )
        )
    }
  })
  
  output$corr_plot <- renderPlot({
    req(values$regression_data)
    numeric_data <- values$regression_data[, sapply(values$regression_data, is.numeric), drop = FALSE]
    if (ncol(numeric_data) > 1) {
      cor_matrix <- cor(numeric_data, use = "complete.obs")
      corrplot(
        cor_matrix,
        method       = "color",
        type         = "upper",
        order        = "hclust",
        tl.col       = "black",
        tl.srt       = 45,
        addCoef.col  = "black",
        number.cex   = 0.7,
        col          = colorRampPalette(c("#3498db", "white", "#e74c3c"))(200),
        title        = "Correlation Matrix",
        mar          = c(0, 0, 2, 0)
      )
    }
  })
  
  output$dist_plot <- renderPlotly({
    req(values$regression_data, input$y_var)
    plot_vars <- c(input$y_var, input$x_vars[1:min(4, length(input$x_vars))])
    plot_vars <- plot_vars[plot_vars %in% names(values$regression_data)]
    plots <- list()
    for (var in plot_vars) {
      p <- plot_ly(
        values$regression_data,
        x     = ~get(var),
        type  = "histogram",
        nbinsx = 30,
        marker = list(
          color = "#3498db",
          line  = list(color = "#2980b9", width = 1)
        ),
        name = var
      ) %>%
        layout(
          xaxis = list(title = var),
          yaxis = list(title = "Frequency"),
          showlegend = FALSE
        )
      plots[[var]] <- p
    }
    if (length(plots) > 0) {
      subplot(plots,
              nrows  = ceiling(length(plots) / 2),
              shareX = FALSE,
              shareY = FALSE) %>%
        layout(title = "Variable Distributions")
    }
  })
  
  output$scatter_matrix <- renderPlot({
    req(values$regression_data, input$y_var)
    plot_vars <- c(input$y_var, input$x_vars[1:min(5, length(input$x_vars))])
    plot_vars <- plot_vars[plot_vars %in% names(values$regression_data)]
    if (length(plot_vars) > 1) {
      pairs(
        values$regression_data[, plot_vars],
        pch  = 19,
        col  = "#3498db",
        cex  = 0.8,
        main = "Scatter Plot Matrix"
      )
    }
  })
  
  # ======================== MODEL SPECIFICATION OUTPUTS =======================
  
  output$model1_formula <- renderUI({
    req(input$y_var, input$x_vars, input$d_var)
    
    x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
    n_main <- length(x_vars_main)
    
    x1 <- if(n_main >= 1) x_vars_main[1] else "X_1"
    x2 <- if(n_main >= 2) x_vars_main[2] else "X_2"
    x3 <- if(n_main >= 3) x_vars_main[3] else "X_3"
    x4 <- if(n_main >= 4) x_vars_main[4] else "X_4"
    x5 <- if(n_main >= 5) x_vars_main[5] else "X_5"
    x6 <- if(n_main >= 6) x_vars_main[6] else "X_6"
    x7 <- if(n_main >= 7) x_vars_main[7] else "X_7"
    x8 <- if(n_main >= 8) x_vars_main[8] else "X_8"
    
    d_var_name <- input$d_var
    y_var_name <- input$y_var
    
    equation <- paste0(
      y_var_name, " = \\beta_0 + \\beta_1 ", x1, " + \\beta_2 ", x2,
      if(n_main >= 3) paste0(" + \\beta_3 ", x3),
      if(n_main >= 4) paste0(" + \\beta_4 ", x4),
      if(n_main >= 5) paste0(" + \\beta_5 ", x5),
      if(n_main >= 6) paste0(" + \\beta_6 ", x6),
      if(n_main >= 7) paste0(" + \\beta_7 ", x7),
      if(n_main >= 8) paste0(" + \\beta_8 ", x8),
      " + \\gamma ", d_var_name
    )
    
    withMathJax(
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #3498db; font-size: 16px; overflow-x: auto; white-space: normal; word-wrap: break-word;",
        tags$p(style = "margin: 0; font-family: 'Courier New', monospace;",
               HTML(paste0("\\[", equation, "\\]"))
        )
      )
    )
  })
  
  output$model2_formula <- renderUI({
    req(input$y_var, input$x_vars, input$d_var)
    
    x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
    n_main <- length(x_vars_main)
    
    x_vars_interaction <- input$x_vars[grepl(paste0(input$d_var, "$"), input$x_vars)]
    x_original <- gsub(paste0(input$d_var, "$"), "", x_vars_interaction)
    n_inter <- length(x_original)
    
    x1 <- if(n_main >= 1) x_vars_main[1] else "X_1"
    x2 <- if(n_main >= 2) x_vars_main[2] else "X_2"
    x3 <- if(n_main >= 3) x_vars_main[3] else "X_3"
    x4 <- if(n_main >= 4) x_vars_main[4] else "X_4"
    x5 <- if(n_main >= 5) x_vars_main[5] else "X_5"
    x6 <- if(n_main >= 6) x_vars_main[6] else "X_6"
    x7 <- if(n_main >= 7) x_vars_main[7] else "X_7"
    x8 <- if(n_main >= 8) x_vars_main[8] else "X_8"
    
    d_var_name <- input$d_var
    y_var_name <- input$y_var
    
    equation <- paste0(
      y_var_name, " = \\beta_0 + \\beta_1 ", x1, " + \\beta_2 ", x2,
      if(n_main >= 3) paste0(" + \\beta_3 ", x3),
      if(n_main >= 4) paste0(" + \\beta_4 ", x4),
      if(n_main >= 5) paste0(" + \\beta_5 ", x5),
      if(n_main >= 6) paste0(" + \\beta_6 ", x6),
      if(n_main >= 7) paste0(" + \\beta_7 ", x7),
      if(n_main >= 8) paste0(" + \\beta_8 ", x8),
      " + \\gamma ", d_var_name
    )
    
    # Tambah interaksi 
    if(n_inter >= 1) {
      for(i in 1:n_inter) {
        var_name <- if(i <= n_main) x_vars_main[i] else paste0("X_", i)
        equation <- paste0(equation, " + \\delta_", i, " (", var_name, " ", d_var_name, ")")
      }
    }
    
    withMathJax(
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #e74c3c; font-size: 16px; overflow-x: auto; white-space: normal; word-wrap: break-word;",
        tags$p(style = "margin: 0; font-family: 'Courier New', monospace;",
               HTML(paste0("\\[", equation, "\\]"))
        )
      )
    )
  })
  
  output$model3_formula <- renderUI({
    req(input$y_var, input$x_vars)
    x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
    M <- length(x_vars_main)
    
    withMathJax(
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #f39c12;",
        paste0(
          "$$\\text{Minimize: } \\sum_{i=1}^n (y_i - \\hat{y}_i)^2 + \\lambda \\left( \\sum_{m=1}^{", M, "} \\beta_m^2 + \\sum_{m=1}^{", M, "} \\sum_{j=1}^{2-1} \\delta_{mj}^2 \\right)$$"
        )
      )
    )
  })
  
  output$model4_formula <- renderUI({
    req(input$y_var, input$x_vars, input$d_var)
    
    x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
    n_main <- length(x_vars_main)
    
    x_vars_interaction <- input$x_vars[grepl(paste0(input$d_var, "$"), input$x_vars)]
    x_original <- gsub(paste0(input$d_var, "$"), "", x_vars_interaction)
    n_inter <- length(x_original)
    
    x1 <- if(n_main >= 1) x_vars_main[1] else "X_1"
    x2 <- if(n_main >= 2) x_vars_main[2] else "X_2"
    x3 <- if(n_main >= 3) x_vars_main[3] else "X_3"
    x4 <- if(n_main >= 4) x_vars_main[4] else "X_4"
    x5 <- if(n_main >= 5) x_vars_main[5] else "X_5"
    x6 <- if(n_main >= 6) x_vars_main[6] else "X_6"
    x7 <- if(n_main >= 7) x_vars_main[7] else "X_7"
    x8 <- if(n_main >= 8) x_vars_main[8] else "X_8"
    
    d_var_name <- input$d_var
    y_var_name <- input$y_var
    
    equation <- paste0(
      y_var_name, " = \\beta_0 + \\beta_1 ", x1, " + \\beta_2 ", x2,
      if(n_main >= 3) paste0(" + \\beta_3 ", x3),
      if(n_main >= 4) paste0(" + \\beta_4 ", x4),
      if(n_main >= 5) paste0(" + \\beta_5 ", x5),
      if(n_main >= 6) paste0(" + \\beta_6 ", x6),
      if(n_main >= 7) paste0(" + \\beta_7 ", x7),
      if(n_main >= 8) paste0(" + \\beta_8 ", x8),
      " + \\gamma ", d_var_name
    )
    
    # Tambah interaksi
    if(n_inter >= 1) {
      for(i in 1:n_inter) {
        var_name <- if(i <= n_main) x_vars_main[i] else paste0("X_", i)
        equation <- paste0(equation, " + \\delta_", i, " (", var_name, " ", d_var_name, ")")
      }
    }
    
    withMathJax(
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #2ecc71; font-size: 16px; overflow-x: auto; white-space: normal; word-wrap: break-word;",
        tags$p(style = "margin: 0; font-family: 'Courier New', monospace; font-weight: normal; margin-bottom: 0; padding-bottom: 0;",
               HTML(paste0("\\[", equation, "\\]"))
        ),
        tags$div(style = "margin-top: 0; padding-top: 0; font-size: 15px;",
                 HTML("dengan robust weights \\( w_i \\) dari \\( \\rho(\\cdot) \\)")
        )
      )
    )
  })
  
  output$model5_formula <- renderUI({
    req(input$y_var, input$x_vars)
    x_vars_main <- input$x_vars[!grepl(paste0(input$d_var, "$"), input$x_vars)]
    M <- length(x_vars_main)
    
    withMathJax(
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #9b59b6;",
        paste0(
          "$$\\text{Minimize: } \\sum_{i=1}^n w_i (y_i - \\hat{y}_i)^2 + \\lambda \\left( \\sum_{m=1}^{", M, "} \\beta_m^2 + \\sum_{m=1}^{", M, "} \\sum_{j=1}^{2-1} \\delta_{mj}^2 \\right)$$"
        )
      )
    )
  })
  
  # ======================== DIAGNOSTICS OUTPUTS ================================
  
  output$vif_plot <- renderPlot({
    req(values$diagnostics$vif)
    vif_df <- data.frame(
      Variable = names(values$diagnostics$vif),
      VIF = as.numeric(values$diagnostics$vif)
    )
    ggplot(vif_df, aes(x = reorder(Variable, VIF), y = VIF)) +
      geom_bar(stat = "identity", fill = "#3498db", alpha = 0.8) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 1) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "darkred", size = 1) +
      coord_flip() +
      labs(title = "Variance Inflation Factors (VIF)",
           subtitle = "Red lines: VIF > 5 (moderate), VIF > 10 (severe)",
           x = "Variable", y = "VIF") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  output$vif_table <- renderDT({
    req(values$diagnostics$vif)
    vif_df <- data.frame(
      Variable = names(values$diagnostics$vif),
      VIF = round(as.numeric(values$diagnostics$vif), 4),
      Status = ifelse(values$diagnostics$vif > 10, "Severe",
                      ifelse(values$diagnostics$vif > 5, "Moderate", "Acceptable"))
    )
    dt_scroller(vif_df, height = 300) %>%
      formatStyle(
        'Status',
        backgroundColor = styleEqual(
          c("Acceptable", "Moderate", "Severe"),
          c("#2ecc71", "#f39c12", "#e74c3c")
        )
      )
  })
  
  output$resid_fitted_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    resid_df <- data.frame(
      Fitted = fitted(model),
      Residuals = residuals(model)
    )
    ggplot(resid_df, aes(x = Fitted, y = Residuals)) +
      geom_point(color = "#3498db", alpha = 0.7, size = 2) +
      geom_smooth(method = "loess", color = "red", se = FALSE, size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      labs(title = "Residuals vs Fitted Values",
           subtitle = "Checking homoscedasticity",
           x = "Fitted Values", y = "Residuals") +
      theme_minimal()
  })
  
  output$qq_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    ggplot(data.frame(Residuals = residuals(model)), aes(sample = Residuals)) +
      stat_qq(color = "#3498db", size = 2) +
      stat_qq_line(color = "red", size = 1) +
      labs(title = "Normal Q-Q Plot",
           subtitle = "Checking normality of residuals",
           x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal()
  })
  
  output$resid_hist <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    ggplot(data.frame(Residuals = residuals(model)), aes(x = Residuals)) +
      geom_histogram(fill = "#3498db", alpha = 0.7, bins = 30) +
      geom_density(aes(y = ..count..), color = "red", size = 1) +
      labs(title = "Histogram of Residuals",
           x = "Residuals", y = "Frequency") +
      theme_minimal()
  })
  
  output$scale_location_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    resid_df <- data.frame(
      Fitted = fitted(model),
      SqrtAbsResid = sqrt(abs(residuals(model)))
    )
    ggplot(resid_df, aes(x = Fitted, y = SqrtAbsResid)) +
      geom_point(color = "#3498db", alpha = 0.7, size = 2) +
      geom_smooth(method = "loess", color = "red", se = FALSE, size = 1) +
      labs(title = "Scale-Location Plot",
           subtitle = "Checking homoscedasticity (constant variance)",
           x = "Fitted Values", y = "√|Standardized Residuals|") +
      theme_minimal()
  })
  
  output$normality_test <- renderPrint({
    req(values$diagnostics$ks)
    
    cat("KOLMOGOROV-SMIRNOV TEST FOR NORMALITY\n")
    cat("======================================\n\n")
    cat("Test: Residuals ~ Normal distribution\n")
    cat("Alpha level: 10%\n\n")
    cat("Test Statistic (D):", round(values$diagnostics$ks$statistic, 4), "\n")
    cat("p-value:", format.pval(values$diagnostics$ks$p.value, digits = 4), "\n\n")
    
    if (values$diagnostics$ks$p.value < 0.1) {
      cat("CONCLUSION: Reject H0 at 10% significance level\n")
      cat("            Residuals are NOT normally distributed (p < 0.10)\n")
      cat("            Consider using robust methods\n")
    } else {
      cat("CONCLUSION: Fail to reject H0 at 10% significance level\n")
      cat("            Residuals are normally distributed (p >= 0.10)\n")
      cat("            Normality assumption appears valid\n")
    }
  })
  
  output$normality_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    residuals <- residuals(model)
    
    df <- data.frame(Residuals = residuals)
    ggplot(df, aes(x = Residuals)) +
      geom_histogram(aes(y = ..density..), 
                     fill = "#3498db", alpha = 0.5, bins = 30) +
      geom_density(color = "#e74c3c", size = 1.5) +
      stat_function(fun = dnorm, 
                    args = list(mean = mean(residuals), sd = sd(residuals)),
                    color = "#2ecc71", size = 1.5, linetype = "dashed") +
      labs(title = "Residuals Distribution vs Normal Curve",
           subtitle = "Green dashed line: Theoretical normal distribution",
           x = "Residuals", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  })
  
  output$bp_test <- renderPrint({
    req(values$diagnostics$bp)
    cat("BREUSCH-PAGAN TEST FOR HETEROSKEDASTICITY\n")
    cat("==========================================\n\n")
    cat("Alpha level: 10%\n\n")
    cat("Test Statistic (LM):", round(values$diagnostics$bp$statistic, 4), "\n")
    cat("p-value:", format.pval(values$diagnostics$bp$p.value, digits = 4), "\n")
    cat("Degrees of freedom:", values$diagnostics$bp$parameter, "\n\n")
    if (values$diagnostics$bp$p.value < 0.1) {
      cat("CONCLUSION: Reject H0 at 10% significance level\n")
      cat("            Heteroskedasticity detected!\n")
      cat("RECOMMENDATION: Use robust standard errors or robust regression methods.\n")
    } else {
      cat("CONCLUSION: Fail to reject H0 at 10% significance level\n")
      cat("            No heteroskedasticity detected.\n")
      cat("RECOMMENDATION: Homoscedasticity assumption appears valid.\n")
    }
  })
  
  output$scale_location_plot2 <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    resid_df <- data.frame(
      Fitted = fitted(model),
      SqrtAbsResid = sqrt(abs(residuals(model)))
    )
    ggplot(resid_df, aes(x = Fitted, y = SqrtAbsResid)) +
      geom_point(color = "#3498db", alpha = 0.7, size = 2) +
      geom_smooth(method = "loess", color = "red", se = FALSE, size = 1) +
      labs(title = "Scale-Location Plot for Heteroskedasticity Check",
           subtitle = "Non-constant spread indicates heteroskedasticity",
           x = "Fitted Values", y = "√|Standardized Residuals|") +
      theme_minimal()
  })
  
  output$autocorrelation_test <- renderPrint({
    req(values$diagnostics$dw)
    
    cat("DURBIN-WATSON TEST FOR AUTOCORRELATION\n")
    cat("=======================================\n\n")
    cat("Test: Autocorrelation of order 1\n")
    cat("Alpha level: 10%\n\n")
    cat("DW Statistic:", round(values$diagnostics$dw$dw, 4), "\n")
    cat("p-value:", format.pval(values$diagnostics$dw$p, digits = 4), "\n\n")
    
    dw_value <- values$diagnostics$dw$dw
    interpretation <- ""
    
    if (dw_value < 1.5) {
      interpretation <- "Positive autocorrelation (DW < 1.5)"
    } else if (dw_value > 2.5) {
      interpretation <- "Negative autocorrelation (DW > 2.5)"
    } else {
      interpretation <- "No significant autocorrelation (1.5 < DW < 2.5)"
    }
    
    cat("Interpretation:", interpretation, "\n\n")
    
    if (values$diagnostics$dw$p < 0.1) {
      cat("CONCLUSION: Reject H0 at 10% significance level\n")
      cat("            Significant autocorrelation detected (p < 0.10)\n")
      cat("RECOMMENDATION: Consider using time series methods or adding lagged variables\n")
    } else {
      cat("CONCLUSION: Fail to reject H0 at 10% significance level\n")
      cat("            No significant autocorrelation (p >= 0.10)\n")
      cat("RECOMMENDATION: Autocorrelation assumption appears valid\n")
    }
  })
  
  output$acf_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    acf_obj <- acf(residuals(model), plot = FALSE)
    acf_df <- data.frame(
      Lag = acf_obj$lag,
      ACF = acf_obj$acf
    )
    ggplot(acf_df, aes(x = Lag, y = ACF)) +
      geom_segment(aes(xend = Lag, yend = 0), size = 1.5, color = "#3498db") +
      geom_hline(yintercept = 0, color = "black") +
      geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed",
                 color = "red", alpha = 0.5) +
      geom_hline(yintercept = c(-1.96/sqrt(length(residuals(model))), 
                                1.96/sqrt(length(residuals(model)))),
                 linetype = "dotted", color = "darkgreen") +
      labs(title = "Autocorrelation Function (ACF) Plot",
           subtitle = "Red dashed lines: ±0.2, Green dotted lines: 95% confidence bounds",
           x = "Lag", y = "ACF") +
      theme_minimal()
  })
  
  output$cooks_plot <- renderPlot({
    req(values$diagnostics$cooks, 
        values$diagnostics$cut_off_f,
        values$diagnostics$cut_off_fixed)
    
    cooks_df <- data.frame(
      Observation = seq_along(values$diagnostics$cooks),
      CooksD = values$diagnostics$cooks
    )
    
    is_outlier_f <- cooks_df$CooksD > values$diagnostics$cut_off_f
    is_outlier_fixed <- cooks_df$CooksD > values$diagnostics$cut_off_fixed
    
    ggplot(cooks_df, aes(x = Observation, y = CooksD)) +
      geom_segment(aes(xend = Observation, yend = 0), 
                   color = "#3498db", alpha = 0.6) +
      geom_point(aes(color = is_outlier_f, 
                     shape = is_outlier_f), 
                 size = 3) +
      
      geom_hline(yintercept = values$diagnostics$cut_off_f, 
                 linetype = "solid", color = "red", size = 1.5) +
      geom_hline(yintercept = values$diagnostics$cut_off_fixed, 
                 linetype = "dotdash", color = "darkgreen", size = 1.2, alpha = 0.8) +
      
      annotate("text", x = max(cooks_df$Observation) * 0.1, 
               y = values$diagnostics$cut_off_f * 1.1, 
               label = paste("F(α=0.5) =", round(values$diagnostics$cut_off_f, 4)),
               color = "red", hjust = 0, size = 3.5) +
      annotate("text", x = max(cooks_df$Observation) * 0.1, 
               y = values$diagnostics$cut_off_fixed * 1.1, 
               label = "Fixed = 1.0",
               color = "darkgreen", hjust = 0, size = 3.5) +
      
      scale_color_manual(
        name = "Outlier (F-dist)",
        values = c("FALSE" = "#3498db", "TRUE" = "red"),
        labels = c("FALSE" = "Normal", "TRUE" = "Outlier")
      ) +
      scale_shape_manual(
        name = "Outlier (F-dist)",
        values = c("FALSE" = 16, "TRUE" = 17),
        labels = c("FALSE" = "Normal", "TRUE" = "Outlier")
      ) +
      
      labs(
        title = "Cook's Distance dengan Dua Metode Cut-off",
        subtitle = "Garis merah: F-distribution (α=0.5), Hijau: Fixed threshold 1.0",
        x = "Observation", 
        y = "Cook's Distance"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom",
        legend.box = "horizontal"
      )
  })
  
  output$leverage_plot <- renderPlot({
    req(values$models$ols_interaction)
    model <- values$models$ols_interaction
    leverage <- hatvalues(model)
    student_resid <- rstudent(model)
    plot_df <- data.frame(
      Leverage = leverage,
      StudentResid = student_resid,
      CooksD = cooks.distance(model),
      Influential = cooks.distance(model) > 4 / length(cooks.distance(model))
    )
    ggplot(plot_df, aes(x = Leverage, y = StudentResid, size = CooksD, color = Influential)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "#3498db", "TRUE" = "red")) +
      geom_vline(xintercept = 2 * mean(leverage), linetype = "dashed", color = "red") +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "orange") +
      labs(title = "Leverage vs Studentized Residuals",
           subtitle = "Red line: 2×mean leverage, Orange lines: ±2 SD",
           x = "Leverage (Hat values)", y = "Studentized Residuals") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$cooks_summary <- renderPrint({
    req(values$diagnostics$cooks, 
        values$diagnostics$cut_off_f,
        values$diagnostics$cut_off_fixed)
    
    n <- length(values$diagnostics$cooks)
    p <- ifelse(!is.null(values$models$ols_interaction), 
                length(coef(values$models$ols_interaction)), 0)
    
    cat("COOK'S DISTANCE ANALYSIS - TWO CUT-OFF METHODS\n")
    cat("==============================================\n\n")
    
    cat("Model Information:\n")
    cat("- Observations (n):", n, "\n")
    cat("- Parameters (p):", p, "\n")
    cat("- df1 (p+1):", p+1, "\n")
    cat("- df2 (n-(p+1)):", n-(p+1), "\n\n")
    
    cat("Cut-off Values:\n")
    cat("1. F-distribution (α=0.5):", round(values$diagnostics$cut_off_f, 4),
        "   (F(", p+1, ",", n-(p+1), "))\n")
    cat("2. Fixed threshold:", values$diagnostics$cut_off_fixed, "\n\n")
    
    cat("Outlier Detection Results:\n")
    cat("- Method 1 (F-distribution):", 
        length(values$diagnostics$outliers_f), 
        "outliers\n")
    cat("- Method 2 (Fixed 1.0):", 
        length(values$diagnostics$outliers_fixed), 
        "outliers\n")
    cat("- Total unique outliers:", 
        length(values$diagnostics$outliers_all), 
        "outliers\n\n")
    
    cat("Comparison Summary:\n")
    cat("Most conservative (fewest outliers):", 
        min(c(length(values$diagnostics$outliers_f),
              length(values$diagnostics$outliers_fixed))), 
        "outliers\n")
    cat("Most aggressive (most outliers):", 
        max(c(length(values$diagnostics$outliers_f),
              length(values$diagnostics$outliers_fixed))), 
        "outliers\n\n")
    
    if (length(values$diagnostics$outliers_all) > 0) {
      cat("All unique outliers detected (any method):\n")
      cat(paste(sort(values$diagnostics$outliers_all), collapse = ", "), "\n\n")
      
      cat("Detail by method:\n")
      cat("1. F-distribution outliers:", 
          ifelse(length(values$diagnostics$outliers_f) > 0,
                 paste(values$diagnostics$outliers_f, collapse = ", "),
                 "None"), "\n")
      cat("2. Fixed threshold outliers:", 
          ifelse(length(values$diagnostics$outliers_fixed) > 0,
                 paste(values$diagnostics$outliers_fixed, collapse = ", "),
                 "None"), "\n")
    } else {
      cat("No outliers detected by any method.\n")
    }
  })
  
  output$outliers_table <- renderDT({
    req(values$diagnostics$outliers_all, values$regression_clean)
    
    if (length(values$diagnostics$outliers_all) > 0) {
      idx <- values$diagnostics$outliers_all
      
      outliers_df <- values$regression_clean[idx, , drop = FALSE]
      outliers_df$Observation <- idx
      outliers_df$CooksD <- values$diagnostics$cooks[idx]
      
      outliers_df$Exceeds_F <- outliers_df$CooksD > values$diagnostics$cut_off_f
      outliers_df$Exceeds_Fixed <- outliers_df$CooksD > values$diagnostics$cut_off_fixed
      
      outliers_df$Methods_Detected <- 
        rowSums(outliers_df[, c("Exceeds_F", "Exceeds_Fixed")])
      
      display_cols <- c("Observation",
                        input$y_var,
                        input$x_vars[1:min(3, length(input$x_vars))],
                        "CooksD",
                        "Methods_Detected",
                        "Exceeds_F",
                        "Exceeds_Fixed")
      
      display_cols <- display_cols[display_cols %in% names(outliers_df)]
      
      display_df <- outliers_df[, display_cols, drop = FALSE]
      
      names(display_df) <- gsub("Exceeds_", "> ", names(display_df))
      names(display_df)[names(display_df) == "> F"] <- "> F(α=0.5)"
      names(display_df)[names(display_df) == "> Fixed"] <- "> 1.0"
      names(display_df)[names(display_df) == "Methods_Detected"] <- "# Methods"
      
      dt_scroller(display_df, height = 300) %>%
        formatRound(columns = c("CooksD", input$y_var), digits = 4) %>%
        formatStyle(
          'CooksD',
          background = styleColorBar(range(display_df$CooksD, na.rm = TRUE), '#FFE5E5'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        ) %>%
        formatStyle(
          names(display_df)[grepl("> ", names(display_df))],
          backgroundColor = styleEqual(
            c(TRUE, FALSE),
            c('#FF9999', '#99FF99')
          )
        ) %>%
        formatStyle(
          '# Methods',
          backgroundColor = styleInterval(
            c(1, 2),
            c('white', '#f39c12', '#e74c3c')
          )
        )
    }
  })
  
  # ======================== ROBUST METHODS OUTPUTS ==============================
  
  output$ridge_cv_plot <- renderPlot({
    req(values$models$ridge_cv)
    cv_fit <- values$models$ridge_cv
    plot_df <- data.frame(
      log_lambda = log(cv_fit$lambda),
      MSE = cv_fit$cvm,
      SE = cv_fit$cvsd
    )
    ggplot(plot_df, aes(x = log_lambda, y = MSE)) +
      geom_line(color = "#3498db", size = 1) +
      geom_ribbon(aes(ymin = MSE - SE, ymax = MSE + SE),
                  alpha = 0.2, fill = "#3498db") +
      geom_vline(xintercept = log(cv_fit$lambda.min),
                 linetype = "dashed", color = "red", size = 1) +
      geom_vline(xintercept = log(cv_fit$lambda.1se),
                 linetype = "dashed", color = "orange", size = 1) +
      labs(title = "Cross-Validation for Ridge Regression",
           subtitle = paste("Optimal lambda:",
                            round(cv_fit$lambda.min, 4), "(min),",
                            round(cv_fit$lambda.1se, 4), "(1 SE)"),
           x = "log(Lambda)", y = "Mean Squared Error") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$ridge_lambda_info <- renderPrint({
    req(values$lambda_min, values$lambda_1se)
    
    cat("RIDGE REGRESSION LAMBDA VALUES\n")
    cat("===============================\n\n")
    cat("Optimal lambda (min):", round(values$lambda_min, 4), "\n")
    cat("1 SE lambda:", round(values$lambda_1se, 4), "\n")
    cat("MSE at lambda(min):", round(values$optimal_mse_min, 4), "\n")
    cat("MSE at lambda(1SE):", round(values$optimal_mse_1se, 4), "\n")
    
    cat("\nInterpretation:\n")
    cat("- Lambda(min): Minimizes cross-validation error\n")
    cat("- Lambda(1SE): Most regularized model within 1 SE of min\n")
  })
  
  output$ridge_path_plot <- renderPlot({
    req(values$models$ridge_full, values$lambda_min)
    ridge_model <- values$models$ridge_full
    lambda_seq <- ridge_model$lambda
    coef_matrix <- as.matrix(ridge_model$beta)
    coef_df <- as.data.frame(t(coef_matrix))
    coef_df$Lambda <- lambda_seq
    coef_long <- coef_df %>%
      tidyr::pivot_longer(cols = -Lambda,
                          names_to = "Variable",
                          values_to = "Coefficient")
    ggplot(coef_long, aes(x = log(Lambda), y = Coefficient, color = Variable)) +
      geom_line(size = 1) +
      geom_vline(xintercept = log(values$lambda_min),
                 linetype = "dashed", color = "red") +
      labs(title = "Ridge Regression Coefficient Paths",
           subtitle = "Effect of regularization on coefficients",
           x = "log(Lambda)", y = "Coefficient Value") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 3))
  })
  
  output$robust_weights_plot <- renderPlot({
    req(values$models$robust)
    robust_model <- values$models$robust
    weights_df <- data.frame(
      Observation = seq_along(robust_model$w),
      Weight = robust_model$w,
      Residual = residuals(robust_model),
      Influence = ifelse(robust_model$w < 0.5, "Low Weight", "Normal")
    )
    ggplot(weights_df, aes(x = Observation, y = Weight, color = Influence)) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
      scale_color_manual(values = c("Normal" = "#3498db", "Low Weight" = "red")) +
      labs(title = "Observation Weights from Robust Regression",
           subtitle = "Low weights indicate potential outliers",
           x = "Observation", y = "Weight") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$ols_vs_robust_plot <- renderPlot({
    req(values$models$ols_interaction, values$models$robust)
    
    ols_coef <- broom::tidy(values$models$ols_interaction) %>%
      dplyr::mutate(Method = "OLS") %>%
      dplyr::select(term, estimate, Method)
    
    robust_coef <- broom::tidy(values$models$robust) %>%
      dplyr::mutate(Method = "Robust") %>%
      dplyr::select(term, estimate, Method)
    
    compare_coef <- dplyr::bind_rows(ols_coef, robust_coef) %>%
      dplyr::filter(term != "(Intercept)")
    
    ggplot(compare_coef, aes(x = term, y = estimate, fill = Method)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      coord_flip() +
      scale_fill_manual(values = c("#3498db", "#e74c3c")) +
      labs(title = "OLS vs Robust Regression Coefficients",
           subtitle = "Comparison of coefficient estimates",
           x = "Variable", y = "Coefficient Value") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$robust_ridge_cv_plot <- renderPlot({
    req(values$models$robust_ridge_cv)
    cv_fit <- values$models$robust_ridge_cv
    plot_df <- data.frame(
      log_lambda = log(cv_fit$lambda),
      MSE = cv_fit$cvm,
      SE = cv_fit$cvsd
    )
    ggplot(plot_df, aes(x = log_lambda, y = MSE)) +
      geom_line(color = "purple", size = 1) +
      geom_ribbon(aes(ymin = MSE - SE, ymax = MSE + SE),
                  alpha = 0.2, fill = "purple") +
      geom_vline(xintercept = log(cv_fit$lambda.min),
                 linetype = "dashed", color = "red", size = 1) +
      geom_vline(xintercept = log(cv_fit$lambda.1se),
                 linetype = "dashed", color = "orange", size = 1) +
      labs(title = "Cross-Validation for Robust Ridge Regression",
           subtitle = paste("Optimal lambda:",
                            round(cv_fit$lambda.min, 4), "(min),",
                            round(cv_fit$lambda.1se, 4), "(1 SE)"),
           x = "log(Lambda)", y = "Mean Squared Error") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$robust_ridge_lambda_info <- renderPrint({
    req(values$lambda_robust_min, values$lambda_robust_1se)
    
    cat("ROBUST RIDGE REGRESSION LAMBDA VALUES\n")
    cat("=====================================\n\n")
    cat("Optimal lambda (min):", round(values$lambda_robust_min, 4), "\n")
    cat("1 SE lambda:", round(values$lambda_robust_1se, 4), "\n")
    cat("MSE at lambda(min):", round(values$optimal_mse_robust_min, 4), "\n")
    cat("MSE at lambda(1SE):", round(values$optimal_mse_robust_1se, 4), "\n")
    
    cat("\nInterpretation:\n")
    cat("- Lambda(min): Minimizes weighted cross-validation error\n")
    cat("- Lambda(1SE): Most regularized robust model within 1 SE of min\n")
  })
  
  output$robust_ridge_plot <- renderPlot({
    req(values$models$robust_ridge)
    tryCatch({
      coef_mat <- as.matrix(coef(values$models$robust_ridge))
      ridge_coef <- data.frame(
        term = rownames(coef_mat),
        estimate = as.numeric(coef_mat[, 1]),
        row.names = NULL
      ) %>%
        dplyr::filter(term != "(Intercept)")
      
      if (nrow(ridge_coef) > 0) {
        ggplot(ridge_coef, aes(x = reorder(term, estimate), y = estimate)) +
          geom_bar(stat = "identity", fill = "purple", alpha = 0.8) +
          coord_flip() +
          labs(title = "Robust Ridge Regression Coefficients",
               x = "Variable", y = "Coefficient Value") +
          theme_minimal()
      } else {
        ggplot() +
          labs(title = "Robust Ridge Regression Coefficients",
               subtitle = "No coefficients to display") +
          theme_minimal()
      }
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Error in Robust Ridge Coefficients", 
                 size = 6, color = "red") +
        theme_void()
    })
  })
  
  # ======================== VALUE BOXES FOR ROBUST METHODS =====================
  
  output$lambda_min_box <- renderValueBox({
    req(values$lambda_min)
    valueBox(
      value = round(values$lambda_min, 4),
      subtitle = "Ridge Lambda (min)",
      icon = icon("bullseye"),
      color = "red"
    )
  })
  
  output$lambda_1se_box <- renderValueBox({
    req(values$lambda_1se)
    valueBox(
      value = round(values$lambda_1se, 4),
      subtitle = "Ridge Lambda (1 SE)",
      icon = icon("shield-alt"),
      color = "orange"
    )
  })
  
  output$lambda_robust_min_box <- renderValueBox({
    req(values$lambda_robust_min)
    valueBox(
      value = round(values$lambda_robust_min, 4),
      subtitle = "Robust Ridge Lambda (min)",
      icon = icon("shield-alt"),
      color = "purple"
    )
  })
  
  output$lambda_robust_1se_box <- renderValueBox({
    req(values$lambda_robust_1se)
    valueBox(
      value = round(values$lambda_robust_1se, 4),
      subtitle = "Robust Ridge Lambda (1 SE)",
      icon = icon("balance-scale"),
      color = "maroon"
    )
  })
  
  output$optimal_mse_box <- renderValueBox({
    req(values$optimal_mse_min)
    valueBox(
      value = round(values$optimal_mse_min, 4),
      subtitle = "Ridge MSE (min)",
      icon = icon("chart-line"),
      color = "green"
    )
  })
  
  output$optimal_mse_1se_box <- renderValueBox({
    req(values$optimal_mse_1se)
    valueBox(
      value = round(values$optimal_mse_1se, 4),
      subtitle = "Ridge MSE (1 SE)",
      icon = icon("chart-line"),
      color = "olive"
    )
  })
  
  output$optimal_mse_robust_box <- renderValueBox({
    req(values$optimal_mse_robust_min)
    valueBox(
      value = round(values$optimal_mse_robust_min, 4),
      subtitle = "Robust Ridge MSE (min)",
      icon = icon("chart-line"),
      color = "teal"
    )
  })
  
  output$optimal_mse_robust_1se_box <- renderValueBox({
    req(values$optimal_mse_robust_1se)
    valueBox(
      value = round(values$optimal_mse_robust_1se, 4),
      subtitle = "Robust Ridge MSE (1 SE)",
      icon = icon("chart-line"),
      color = "aqua"
    )
  })
  
  # ======================== SIGNIFICANCE TESTS OUTPUTS =======================
  
  output$f_test_table <- renderTable({
    req(values$significance_tests)
    values$significance_tests$f_test %>%
      mutate(
        across(where(is.numeric), ~round(., 4))
      )
  })
  
  output$f_test_interpretation <- renderPrint({
    req(values$significance_tests)
    f_test <- values$significance_tests$f_test
    
    p_val_num <- as.numeric(gsub("< ", "", f_test$p_value))
    if (is.na(p_val_num)) {
      p_val_num <- 0.00001
    }
    
    cat("INTERPRETASI UJI F:\n")
    cat("===================\n\n")
    cat("Hipotesis:\n")
    cat("  H0: Semua koefisien = 0\n")
    cat("  H1: Minimal ada satu koefisien ≠ 0\n\n")
    
    cat("Statistik Uji:\n")
    cat("  F-statistic:", f_test$Statistic, "\n")
    cat("  df1:", f_test$df1, "\n")
    cat("  df2:", f_test$df2, "\n")
    cat("  p-value:", f_test$p_value, "\n\n")
    
    cat("Kesimpulan:\n")
    if (p_val_num < 0.1) {
      cat("  REJECT H0 (p < 0.1)\n")
      cat("  Model secara keseluruhan signifikan secara statistik\n")
      cat("  Minimal ada satu variabel independen yang berpengaruh signifikan terhadap Y\n")
    } else {
      cat("  FAIL TO REJECT H0 (p >= 0.1)\n")
      cat("  Model tidak signifikan secara statistik\n")
      cat("  Tidak cukup bukti bahwa variabel independen mempengaruhi Y\n")
    }
  })
  
  output$t_test_table <- renderDT({
    req(values$significance_tests)
    t_test <- values$significance_tests$t_test %>%
      mutate(
        Estimate = mapply(function(x, y) format_t_test(x, y, "estimate"), 
                          as.numeric(Estimate), Variable),
        `Std. Error` = format(round(as.numeric(`Std. Error`), 4), 
                              scientific = FALSE, nsmall = 4),
        `t value` = format(round(as.numeric(`t value`), 4), 
                           scientific = FALSE, nsmall = 4),
        `Pr(>|t|)` = format(as.numeric(`Pr(>|t|)`), 
                            scientific = FALSE, digits = 4, nsmall = 4)
      )
    
    t_test$`Pr(>|t|)` <- sapply(as.numeric(t_test$`Pr(>|t|)`), function(p) {
      if (is.na(p)) return("NA")
      if (p < 0.0001) {
        return(format(p, scientific = TRUE, digits = 4))
      } else {
        return(format(p, scientific = FALSE, digits = 4, nsmall = 4))
      }
    })
    
    datatable(
      t_test,
      options = list(
        dom = "t",
        scrollX = TRUE,
        scrollY = 400,
        paging = FALSE
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Pr(>|t|)',
        backgroundColor = styleInterval(
          c(0.01, 0.05, 0.1),
          c("#27ae60", "#2ecc71", "#f1c40f", "#e74c3c")
        )
      ) %>%
      formatStyle(
        'Significance',
        fontWeight = 'bold',
        color = styleEqual(
          c("***", "**", "*", ""),
          c("#e74c3c", "#e74c3c", "#e74c3c", "black")
        )
      )
  })
  
  output$significance_legend <- renderUI({
    tagList(
      h5("Keterangan Signifikansi:"),
      tags$ul(
        tags$li(tags$span(class = "sig-star", "***"), "p-value < 0.01 (sangat signifikan)"),
        tags$li(tags$span(class = "sig-star", "**"), "p-value < 0.05 (signifikan)"),
        tags$li(tags$span(class = "sig-star", "*"), "p-value < 0.10 (signifikan)"),
        tags$li("(tanpa bintang): p-value >= 0.1 (tidak signifikan)")
      ),
      hr(),
      h5("Interpretasi Uji t:"),
      tags$ul(
        tags$li("H0: Koefisien = 0 (variabel tidak berpengaruh)"),
        tags$li("H1: Koefisien ≠ 0 (variabel berpengaruh)"),
        tags$li("Jika p-value < 0.1, tolak H0 (variabel signifikan)")
      )
    )
  })
  
  output$wald_test_table <- renderDT({
    req(values$significance_tests)
    wald_tests <- values$significance_tests$wald_tests
    
    if (!is.null(wald_tests) && nrow(wald_tests) > 0) {
      wald_tests <- wald_tests %>%
        mutate(
          F_Statistic = round(as.numeric(F_Statistic), 4),
          Hypothesis = gsub("\\+", "+", Hypothesis)
        )
      
      datatable(
        wald_tests,
        options = list(
          dom = "t",
          scrollX = TRUE,
          scrollY = 300,
          paging = FALSE
        ),
        rownames = FALSE
      ) %>%
        formatStyle(
          'Significance',
          backgroundColor = styleEqual(
            c("Signifikan", "Tidak Signifikan"),
            c("#d4edda", "#f8d7da")
          )
        )
    }
  })
  
  output$wald_test_summary <- renderPrint({
    req(values$significance_tests)
    wald_tests <- values$significance_tests$wald_tests
    
    if (!is.null(wald_tests) && nrow(wald_tests) > 0) {
      cat("INTERPRETASI UJI WALD:\n")
      cat("======================\n\n")
      cat("Hipotesis:\n")
      cat("  H0: Total koefisien di Klaster 2 = 0\n")
      cat("      (Koefisien Klaster 1 + Koefisien Klaster 2 = 0)\n")
      cat("  H1: Total koefisien di Klaster 2 ≠ 0\n\n")
      
      for(i in 1:nrow(wald_tests)) {
        cat(gsub("Wald Test - Perbedaan Antar Klaster ", "", wald_tests$Test[i]), ":\n")
        cat("  H0:", wald_tests$Hypothesis[i], "\n")
        cat("  H1: Total Koefisien ≠ 0 (variabel berpengaruh di klaster pembanding)\n")
        cat("  F-statistic:", wald_tests$F_Statistic[i], "\n")
        cat("  p-value:", wald_tests$p_value[i], "\n")
        cat("  Kesimpulan:", wald_tests$Significance[i], "\n")
        cat("\n")
      }
      
      signif_count <- sum(wald_tests$Significance == "Signifikan")
      total_count <- nrow(wald_tests)
      
      cat("RINGKASAN:\n")
      cat("  Total uji Wald:", total_count, "\n")
      cat("  Variabel signifikan di Klaster 2:", signif_count, "\n")
      cat("  Variabel tidak signifikan di Klaster 2:", total_count - signif_count, "\n")
      
    } else {
      cat("Uji Wald tidak tersedia atau terjadi error.\n")
    }
  })
  
  # ======================== RESULTS & REPORTS OUTPUTS ==========================
  
  output$key_findings <- renderUI({
    req(values$models$ols_interaction)
    has_multicollinearity <- !is.null(values$diagnostics$vif) &&
      any(values$diagnostics$vif > 5)
    has_normality_issue <- !is.null(values$diagnostics$ks) &&
      values$diagnostics$ks$p.value < 0.1
    has_heteroskedasticity <- !is.null(values$diagnostics$bp) &&
      values$diagnostics$bp$p.value < 0.1
    has_autocorrelation <- !is.null(values$diagnostics$dw) &&
      values$diagnostics$dw$p < 0.1
    has_outliers <- !is.null(values$diagnostics$outliers_all) &&
      length(values$diagnostics$outliers_all) > 0
    
    tagList(
      h4("Diagnostic Assessment (α = 10%):"),
      tags$ul(
        if (has_multicollinearity) {
          tags$li(icon("exclamation-triangle", class = "text-warning"),
                  strong("Multicollinearity detected."),
                  "Some VIF values exceed threshold of 10.")
        } else {
          tags$li(icon("check-circle", class = "text-success"),
                  "No significant multicollinearity detected.")
        },
        if (has_normality_issue) {
          tags$li(icon("exclamation-triangle", class = "text-warning"),
                  strong("Non-normal residuals."),
                  "Kolmogorov-Smirnov test indicates deviation from normality (α=10%).")
        } else {
          tags$li(icon("check-circle", class = "text-success"),
                  "Residuals appear normally distributed (α=10%).")
        },
        if (has_heteroskedasticity) {
          tags$li(icon("exclamation-triangle", class = "text-warning"),
                  strong("Heteroskedasticity detected."),
                  "Breusch-Pagan test indicates non-constant variance (α=10%).")
        } else {
          tags$li(icon("check-circle", class = "text-success"),
                  "Homoscedasticity assumption appears valid (α=10%).")
        },
        if (has_autocorrelation) {
          tags$li(icon("exclamation-triangle", class = "text-warning"),
                  strong("Autocorrelation detected."),
                  "Durbin-Watson test indicates correlated residuals (α=10%).")
        } else {
          tags$li(icon("check-circle", class = "text-success"),
                  "No significant autocorrelation detected (α=10%).")
        },
        if (has_outliers) {
          tags$li(icon("exclamation-triangle", class = "text-warning"),
                  strong(paste(length(values$diagnostics$outliers_all),
                               "potential outliers detected.")),
                  "Cook's distance indicates influential observations.")
        } else {
          tags$li(icon("check-circle", class = "text-success"),
                  "No significant outliers detected.")
        }
      ),
      hr(),
      h4("Recommendations:"),
      tags$ul(
        if (has_multicollinearity) {
          tags$li("Use Ridge Regression to address multicollinearity")
        },
        if (has_outliers) {
          tags$li("Use Robust Regression to address outliers")
        },
        if (has_multicollinearity || has_outliers) {
          tags$li("Use Ridge Robust Regression to address multicollinearity and outliers")
        },
        if (has_autocorrelation) {
          tags$li("Consider adding lagged variables or using time series methods")
        },
        if (!has_multicollinearity && !has_outliers && !has_autocorrelation) {
          tags$li("OLS provides reliable estimates for this dataset")
        }
      )
    )
  })
  
  output$executive_summary <- renderUI({
    req(values$models$ols_interaction, values$performance)
    
    best_row <- values$performance %>%
      dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
      dplyr::slice(1)
    best_model <- as.character(best_row$Model)
    best_R2 <- best_row$R2
    best_RMSE <- best_row$RMSE
    
    tagList(
      h4("Executive Summary:"),
      p("This analysis examines the relationship between the dependent variable ",
        strong(input$y_var), " and ", length(input$x_vars),
        " independent variables, with clustering based on ", strong(input$d_var), "."),
      p("The analysis employed multiple regression techniques to address potential issues of multicollinearity and outliers."),
      p("All classical assumption tests were conducted at α = 10% significance level."),
      p("Based on performance metrics, the ", strong(best_model),
        " model appears to provide the best balance of fit and generalization."),
      hr(),
      h4("Key Results:"),
      tags$ul(
        tags$li(paste("Sample size:", nrow(values$regression_clean), "observations")),
        tags$li(paste("Best model:", best_model)),
        tags$li(paste("Best R²:", round(best_R2, 4))),
        tags$li(paste("Best RMSE:", round(best_RMSE, 4)))
      )
    )
  })
  
  output$recommended_model <- renderUI({
    req(values$performance)
    
    best_row <- values$performance %>%
      dplyr::filter(RMSE == min(RMSE, na.rm = TRUE)) %>%
      dplyr::slice(1)
    best_model <- as.character(best_row$Model)
    
    tagList(
      h3(best_model),
      tags$div(
        style = "padding: 10px; background-color: #d4edda; border-radius: 5px;",
        p("Chosen based on the highest R² and the lowest RMSE, with consideration of diagnostic problems.")
      )
    )
  })
  
  output$model_selection_criteria <- renderTable({
    req(values$performance)
    
    values$performance %>%
      dplyr::mutate(
        Score = rowSums(cbind(
          ifelse(RMSE == min(RMSE, na.rm = TRUE), 1, 0),
          ifelse(R2 == max(R2, na.rm = TRUE), 1, 0),
          ifelse(MAE == min(MAE, na.rm = TRUE), 1, 0)
        ), na.rm = TRUE)
      ) %>%
      dplyr::select(Model, R2, RMSE, MAE, Score) %>%
      dplyr::arrange(dplyr::desc(Score)) %>%
      mutate(across(where(is.numeric), ~round(., 4)))
  })
  
  output$detailed_report <- renderUI({
    req(values$models$ols_interaction)
    tagList(
      h5("1. Data Description"),
      p(paste("The analysis used a dataset with", nrow(values$regression_clean),
              "complete observations and",
              length(c(input$x_vars, input$d_var)) + 1,
              "variables.")),
      h5("2. Model Specifications"),
      p("Four main model specifications were tested:"),
      tags$ul(
        tags$li("Pooled regression with dummy intercept shift"),
        tags$li("Full interaction model with different slopes"),
        tags$li("Ridge, Robust, and Robust Ridge models for stability"),
      ),
      h5("3. Estimation Results"),
      p("All models were estimated using maximum likelihood/least squares methods."),
      h5("4. Diagnostic Tests (α = 10%)"),
      p("Comprehensive diagnostic tests were performed including:"),
      tags$ul(
        tags$li("Multicollinearity (VIF)"),
        tags$li("Normality (Kolmogorov-Smirnov)"),
        tags$li("Heteroskedasticity (Breusch-Pagan)"),
        tags$li("Autocorrelation (Durbin-Watson)"),
        tags$li("Outliers (Cook's distance with 2 cut-off methods)")
      ),
      h5("5. Significance Tests"),
      p("Complete significance testing performed:"),
      tags$ul(
        tags$li("F-test for overall model significance"),
        tags$li("t-test for individual coefficient significance"),
        tags$li("Wald test for differences between clusters")
      ),
      h5("6. Recommendations"),
      p("Based on the analysis, specific recommendations are provided for model selection and further analysis.")
    )
  })
  
  # ======================== DOWNLOAD HANDLERS ================================
  
  output$download_coefficients <- downloadHandler(
    filename = function() paste("coefficients_", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      coef_list <- list()
      if (!is.null(values$models$ols_interaction)) {
        coef_list$OLS_Interaction <- broom::tidy(values$models$ols_interaction) %>%
          mutate(Model = "OLS_Interaction")
      }
      if (!is.null(values$models$ridge)) {
        tryCatch({
          coef_mat <- as.matrix(coef(values$models$ridge))
          ridge_coef <- data.frame(
            term     = rownames(coef_mat),
            estimate = as.numeric(coef_mat[, 1]),
            row.names = NULL
          ) %>%
            mutate(Model = "Ridge")
          coef_list$Ridge <- ridge_coef
        }, error = function(e) {
          coef_list$Ridge <- NULL
        })
      }
      if (!is.null(values$models$robust)) {
        tryCatch({
          coef_list$Robust <- broom::tidy(values$models$robust) %>%
            mutate(Model = "Robust")
        }, error = function(e) {
          coef_list$Robust <- NULL
        })
      }
      if (!is.null(values$models$robust_ridge)) {
        tryCatch({
          coef_mat <- as.matrix(coef(values$models$robust_ridge))
          robust_ridge_coef <- data.frame(
            term     = rownames(coef_mat),
            estimate = as.numeric(coef_mat[, 1]),
            row.names = NULL
          ) %>%
            mutate(Model = "Robust_Ridge")
          coef_list$Robust_Ridge <- robust_ridge_coef
        }, error = function(e) {
          coef_list$Robust_Ridge <- NULL
        })
      }
      all_coefs <- dplyr::bind_rows(coef_list)
      write.csv(all_coefs, file, row.names = FALSE)
    }
  )
  
  output$download_report_pdf <- downloadHandler(
    filename = function() paste0("regression_analysis_", Sys.Date(), ".pdf"),
    content = function(file) {
      temp_dir <- tempdir()
      temp_file <- file.path(temp_dir, "report.pdf")
      
      pdf(temp_file, width = 11, height = 8.5)
      
      plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
           xlab = "", ylab = "", main = "")
      text(5, 7, "REGRESSION CLUSTERING ANALYSIS REPORT", cex = 2, font = 2)
      text(5, 6, paste("Generated:", Sys.Date()), cex = 1.2)
      text(5, 5, paste("Dependent Variable:", input$y_var), cex = 1.2)
      text(5, 4, paste("Cluster Variable:", input$d_var), cex = 1.2)
      text(5, 3, paste("Independent Variables:", paste(input$x_vars, collapse = ", ")), cex = 1)
      
      if (!is.null(values$models$ols_interaction)) {
        par(mfrow = c(2, 2))
        plot(values$models$ols_interaction, which = 1)
        plot(values$models$ols_interaction, which = 2)
        plot(values$models$ols_interaction, which = 3)
        plot(values$models$ols_interaction, which = 5)
        
        par(mfrow = c(1, 1))
        plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
             xlab = "", ylab = "", main = "")
        text(5, 7, "MODEL SUMMARY", cex = 1.5, font = 2)
        
        if (!is.null(values$performance)) {
          text(5, 6, paste("Best Model:", values$performance$Model[which.min(values$performance$RMSE)]), cex = 1)
          text(5, 5, paste("Best R²:", round(max(values$performance$R2, na.rm = TRUE), 4)), cex = 1)
          text(5, 4, paste("Best RMSE:", round(min(values$performance$RMSE, na.rm = TRUE), 4)), cex = 1)
        }
      }
      
      dev.off()
      
      file.copy(temp_file, file)
    },
    contentType = "application/pdf"
  )
  
  output$download_report_html <- downloadHandler(
    filename = function() paste0("regression_analysis_", Sys.Date(), ".html"),
    content = function(file) {
      html_content <- paste0(
        "<!DOCTYPE html>
        <html>
        <head>
          <title>Regression Clustering Analysis Report</title>
          <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
            h2 { color: #34495e; margin-top: 30px; }
            .section { margin-bottom: 30px; padding: 20px; background-color: #f8f9fa; border-radius: 5px; }
            table { border-collapse: collapse; width: 100%; margin: 10px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #3498db; color: white; }
            .info-box { background-color: #e8f4f8; padding: 15px; border-left: 4px solid #3498db; margin: 10px 0; }
          </style>
        </head>
        <body>
          <h1>Regression Clustering Analysis Report</h1>
          <div class='info-box'>
            <strong>Generated:</strong> ", Sys.Date(), "<br>
            <strong>Analysis Type:</strong> Dummy Variable Regression Clustering
          </div>
          
          <div class='section'>
            <h2>Model Specification</h2>
            <p><strong>Dependent Variable (Y):</strong> ", input$y_var, "</p>
            <p><strong>Cluster Variable (D):</strong> ", input$d_var, "</p>
            <p><strong>Independent Variables (X):</strong> ", paste(input$x_vars, collapse = ", "), "</p>
            <p><strong>Estimation Methods:</strong> ", paste(input$methods, collapse = ", "), "</p>
          </div>"
      )
      
      if (!is.null(values$performance)) {
        html_content <- paste0(html_content, "
          <div class='section'>
            <h2>Performance Metrics</h2>
            <table>
              <tr>
                <th>Model</th>
                <th>R²</th>
                <th>Adj. R²</th>
                <th>RMSE</th>
                <th>MAE</th>
                <th>AIC</th>
                <th>BIC</th>
              </tr>")
        
        for(i in 1:nrow(values$performance)) {
          html_content <- paste0(html_content, "
              <tr>
                <td>", values$performance$Model[i], "</td>
                <td>", values$performance$R2[i], "</td>
                <td>", ifelse(is.na(values$performance$Adj_R2[i]), "NA", values$performance$Adj_R2[i]), "</td>
                <td>", values$performance$RMSE[i], "</td>
                <td>", values$performance$MAE[i], "</td>
                <td>", ifelse(is.na(values$performance$AIC[i]), "NA", values$performance$AIC[i]), "</td>
                <td>", ifelse(is.na(values$performance$BIC[i]), "NA", values$performance$BIC[i]), "</td>
              </tr>")
        }
        
        html_content <- paste0(html_content, "
            </table>
          </div>")
      }
      
      if (!is.null(values$significance_tests)) {
        html_content <- paste0(html_content, "
          <div class='section'>
            <h2>Significance Tests</h2>
            <h3>F-Test (Overall Model Significance)</h3>
            <p><strong>F-statistic:</strong> ", values$significance_tests$f_test$Statistic, "</p>
            <p><strong>p-value:</strong> ", values$significance_tests$f_test$p_value, "</p>
            <p><strong>Conclusion:</strong> ", values$significance_tests$f_test$Significance, "</p>
          </div>")
      }
      
      if (!is.null(values$shape_data)) {
        cluster_data <- st_drop_geometry(values$shape_data)
        
        if("Klaster" %in% names(cluster_data)) {
          cluster_summary <- cluster_data %>%
            group_by(Klaster) %>%
            summarise(Count = n(), .groups = "drop")
          
          html_content <- paste0(html_content, "
            <div class='section'>
              <h2>Clustering Information</h2>
              <p><strong>Total Regions:</strong> ", nrow(cluster_data), "</p>
              <table>
                <tr>
                  <th>Klaster</th>
                  <th>Count</th>
                  <th>Percentage</th>
                </tr>")
          
          for(i in 1:nrow(cluster_summary)) {
            percentage <- round(cluster_summary$Count[i] / nrow(cluster_data) * 100, 1)
            html_content <- paste0(html_content, "
                <tr>
                  <td>", cluster_summary$Klaster[i], "</td>
                  <td>", cluster_summary$Count[i], "</td>
                  <td>", percentage, "%</td>
                </tr>")
          }
          
          html_content <- paste0(html_content, "
              </table>
            </div>")
        }
      }
      
      html_content <- paste0(html_content, "
          <div class='info-box'>
            <strong>Report generated by:</strong> Regression Clustering Dashboard<br>
            <strong>Timestamp:</strong> ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
          </div>
        </body>
        </html>")
      
      writeLines(html_content, file)
    }
  )
  
  output$download_diagnostics <- downloadHandler(
    filename = function() paste("diagnostics_", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      diag_list <- list()
      
      if (!is.null(values$performance)) {
        diag_list$Performance <- values$performance
      }
      
      if (!is.null(values$diagnostics$vif)) {
        diag_list$VIF <- data.frame(
          Variable = names(values$diagnostics$vif),
          VIF = round(as.numeric(values$diagnostics$vif), 4),
          Status = ifelse(values$diagnostics$vif > 10, "Severe",
                          ifelse(values$diagnostics$vif > 5, "Moderate", "Acceptable"))
        )
      }
      
      if (!is.null(values$diagnostics$ks)) {
        diag_list$Normality <- data.frame(
          Test = "Kolmogorov-Smirnov",
          Statistic = round(values$diagnostics$ks$statistic, 4),
          p_value = format.pval(values$diagnostics$ks$p.value, digits = 4),
          Conclusion = ifelse(values$diagnostics$ks$p.value < 0.1, 
                              "Not Normal (p < 0.10)", "Normal (p >= 0.10)")
        )
      }
      
      if (!is.null(values$diagnostics$bp)) {
        diag_list$Heteroskedasticity <- data.frame(
          Test = "Breusch-Pagan",
          Statistic = round(values$diagnostics$bp$statistic, 4),
          p_value = format.pval(values$diagnostics$bp$p.value, digits = 4),
          df = values$diagnostics$bp$parameter,
          Conclusion = ifelse(values$diagnostics$bp$p.value < 0.1, 
                              "Heteroskedasticity detected", "Homoscedastic")
        )
      }
      
      if (!is.null(values$diagnostics$dw)) {
        diag_list$Autocorrelation <- data.frame(
          Test = "Durbin-Watson",
          Statistic = round(values$diagnostics$dw$dw, 4),
          p_value = format.pval(values$diagnostics$dw$p, digits = 4),
          Conclusion = ifelse(values$diagnostics$dw$p < 0.1, 
                              "Autocorrelation detected", "No autocorrelation")
        )
      }
      
      if (!is.null(values$diagnostics$outliers_all) &&
          length(values$diagnostics$outliers_all) > 0) {
        
        outliers_info <- data.frame(
          Observation = values$diagnostics$outliers_all,
          CooksD = round(values$diagnostics$cooks[values$diagnostics$outliers_all], 4),
          Exceeds_F_Cutoff = values$diagnostics$cooks[values$diagnostics$outliers_all] > values$diagnostics$cut_off_f,
          Exceeds_Fixed_Cutoff = values$diagnostics$cooks[values$diagnostics$outliers_all] > values$diagnostics$cut_off_fixed,
          Methods_Detected = rowSums(cbind(
            values$diagnostics$cooks[values$diagnostics$outliers_all] > values$diagnostics$cut_off_f,
            values$diagnostics$cooks[values$diagnostics$outliers_all] > values$diagnostics$cut_off_fixed
          ))
        )
        
        if (!is.null(values$regression_clean) && nrow(values$regression_clean) > 0) {
          outlier_data <- values$regression_clean[values$diagnostics$outliers_all, , drop = FALSE]
          outliers_info <- cbind(outliers_info, outlier_data)
        }
        
        diag_list$Outliers <- outliers_info
      }
      
      if (!is.null(values$diagnostics$cut_off_f)) {
        diag_list$Cutoffs <- data.frame(
          Cutoff_Type = c("F-distribution (α=0.5)", "Fixed Threshold"),
          Value = c(round(values$diagnostics$cut_off_f, 4), 
                    values$diagnostics$cut_off_fixed)
        )
      }
      
      if (length(diag_list) == 0) {
        all_diag <- data.frame(
          Message = "No diagnostics available. Run regression analysis first.",
          Timestamp = Sys.time()
        )
      } else {
        all_diag <- dplyr::bind_rows(lapply(names(diag_list), function(name) {
          df <- diag_list[[name]]
          df$Diagnostic_Type <- name
          return(df)
        }))
      }
      
      write.csv(all_diag, file, row.names = FALSE)
    }
  )
  
  output$download_shapefile_data <- downloadHandler(
    filename = function() paste("clustering_data_", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      req(values$shape_data)
      data_to_export <- st_drop_geometry(values$shape_data)
      write.csv(data_to_export, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$export_plots, {
    showModal(modalDialog(
      title = "Export Plots",
      "This feature would export all generated plots as image files.",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
}

# ==============================================================================
# 6. RUN APP
# ==============================================================================

shinyApp(
  ui = dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body,
    title = "SCR Dashboard",
    skin  = "blue"
  ),
  server = server
)