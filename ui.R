library(ddtlcm)
library(shiny)
library(shinyFiles)
library(plotly)
library(ape)
library(ggplot2)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
library(DT)
library(stringr)
library(RColorBrewer)

# For plot_summary_test function
library(reshape2)
library(data.table)
library(ggpubr)
library(ggtree)
library(ggtext)

options(warn = -1)

ui = fluidPage(
  id = "fluidpage",
  # Some custom CSS
  tags$head(
    tags$style(HTML("
        /* Smaller font for preformatted text */
        pre, table.table {
          font-size: smaller;
        }

        body {
          min-height: 2000px;
        }

        .option-group {
          border: 1px solid #ccc;
          border-radius: 6px;
          padding: 0px 5px;
          margin: 5px -10px;
          background-color: #f5f5f5;
        }

        .option-header {
          color: #79d;
          text-transform: uppercase;
          margin-bottom: 5px;
        }
        
        td {
        padding: 3px; /* Adjust the padding to control spacing within cells */
        margin: 5px; /* Adjust the margin to control spacing around cells */
        }
      "))
  ),
  
  # New UI design BEGIN
  
  headerPanel("DDT-LCM: Dirichlet Diffusion Tree-Latent Class Model"),
  fluidRow(
    column(12, 
           HTML("<p style='margin-left: 18px; font-size: 15px; color: gray;'>
                This application complements the underlying package, <a href='https://github.com/limengbinggz/ddtlcm/'>ddtlcm</a>, 
                by providing a font-end interface for users to explore and visualize the results of the 'DDT-LCM' model implemented in the package.
                'DDT-LCM' implements a tree-regularized Bayesian latent class analysis, which enhances classical latent class models 
                <a href='https://arxiv.org/abs/2306.04700'>(Li, Stephenson, Wu (2025) Annals of Applied Statistics. In press)</a>.
                It address the long-standing issues of numerical and inferential instabilities for latent class models 
                under weak class separation particularly for small sample sizes. <br>
                <p><small>Current Version: 0.1.0</small></p>
                </p>")
    ),
  ),
  
  
  fluidRow(
    column(12, align = "left",
           checkboxInput("credits_checkbox", "Show References and Maintainer Contacts", TRUE),
           conditionalPanel("input.credits_checkbox == 1",
                            uiOutput("credits_instr"),
           ),
    )
  ),
  
  
                                    

  
  
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      
      div(class = "option-header", "Mode"),
      radioButtons("mode", "Options",
                   choices = c("Simulate Data", 
                               "Upload Raw Data", 
                               "Upload Posterior Samples")),
      wellPanel(
        conditionalPanel("input.mode == 'Simulate Data'",
                         div(class = "option-header", "Model Parameters"),
                         radioButtons("sim_data_src", "Source",
                                      choices = c("Exemplar Parameters", "Upload Parameters"),
                         ),
                         
                         conditionalPanel("input.sim_data_src == 'Exemplar Parameters'",
                                          wellPanel(
                                            fluidRow(
                                              column(width = 11,
                                                     radioButtons("dataset_name", "Scenario",
                                                                  choices = c("Semi-Synthetic Data (mimicking a real data set; see Li, Stephenson, Wu (2023) Arxiv)", "Synthetic Data"),
                                                     ),
                                              )
                                              
                                            ),
                                          )
                         ),
                         conditionalPanel("input.sim_data_src == 'Upload Parameters'",
                                          fluidRow(
                                            column(width = 11,
                                                   numericInput("N", "Number of Individuals, N", value = 496,
                                                                min = 50, max = 1000, step = 1),
                                            )
                                          ),
                                          fluidRow(
                                            column(width = 8,
                                                   fileInput("tree_phylo_file", "Upload the tree phylo")
                                            ),
                                          ),
                                          fluidRow(
                                            column(width = 8,
                                                   fileInput("class_probability_file", "Upload the class probability list")
                                            ),
                                            column(width = 1,
                                                   div(style = "width: 100%; height: 20px;"),
                                                   checkboxInput("class_probability_header", "Header", TRUE)
                                            ),
                                          ),
                                          fluidRow(
                                            column(width = 8,
                                                   fileInput("sigma_by_group_file", "Upload sigma squared by group list")
                                            ),
                                            column(width = 1,
                                                   div(style = "width: 100%; height: 20px;"),
                                                   checkboxInput("sigma_by_group_header", "Header", TRUE)
                                            ),
                                          ),
                         )
                         
        ),
        conditionalPanel("input.mode == 'Upload Raw Data'",
                         div(class = "option-header", "Upload Data"),
                         fluidRow(
                           column(width = 8,
                                  fileInput("data_matrix_file", "Upload the data matrix")
                           ),
                           column(width = 1,
                                  div(style = "width: 100%; height: 20px;"),
                                  checkboxInput("data_matrix_header", "Header", TRUE)
                           ),
                         ),
        ),
        conditionalPanel("(input.mode == 'Simulate Data' & input.sim_data_src == 'Upload Parameters') 
                           | input.mode == 'Upload Raw Data'",
                         fluidRow(
                           column(width = 8,
                                  fileInput("item_membership_file", "Upload the item membership list")
                           ),
                           # column(width = 1,
                           #        div(style = "width: 100%; height: 20px;"),
                           #        checkboxInput("item_membership_header", "Header", FALSE)
                           # ),
                         ),
        ),
        conditionalPanel("input.mode == 'Upload Posterior Samples'",
                         div(class = "option-header", "Upload Data"),
                         fluidRow(
                           column(width = 8,
                                  fileInput("posterior_sample_file", "Upload posterior samples")
                           ),
                         ),
        ),
        conditionalPanel("input.mode != 'Simulate Data' | input.sim_data_src != 'Exemplar Parameters'",
                         radioButtons("item_name_opt", "Optional: Upload the item name list",
                                      choices = c("Default", "Upload"), inline = TRUE),
                         conditionalPanel("input.item_name_opt == 'Upload'",
                                          fluidRow(
                                            column(width = 8,
                                                   fileInput("item_name_file", "")
                                            ),
                                            # column(width = 1,
                                            #        div(style = "width: 100%; height: 20px;"),
                                            #        checkboxInput("item_name_header", "Header", TRUE)
                                            # ),
                                          ),
                         ),
                         # actionButton("read_button", "Read in"),
        ),
        conditionalPanel("input.mode == 'Simulate Data'",
                         fluidRow(
                           column(width = 11,
                                  numericInput("seed_node_param", 
                                               "Random Seed to Generate Node Parameters", 
                                               value = 1,
                                               min = 1, max = 1000, step = 1),
                           ),
                         ),
                         fluidRow(
                           column(width = 11,
                                  numericInput("seed_bin_obs", 
                                               "Random Seed to Generate Binary Observations", 
                                               value = 1,
                                               min = 1, max = 1000, step = 1),
                           ),
                         )
        ),
        actionButton("read_button", "Read in"),
      ),
      
      conditionalPanel("input.mode == 'Simulate Data' | input.mode == 'Upload Raw Data'",
                       wellPanel(
                         div(class = "option-header", "Model Fit Options"),
                         numericInput("K_fit", "Number of Classes, K", value = 6,
                                      min = 1, max = 500, step = 1),
                         numericInput("niter", "Number of Interations",
                                      value = 10, min = 1, max = 100000, step = 1),
                         numericInput("burnin", "Number of Burn-ins",
                                      value = 4, min = 0, max = 99999, step = 1),
                         numericInput("seed_posterior", 
                                      "Random Seed to Generate Posterior Samples",
                                      value = 999, min = 1, max = 1000, step = 1),
                         actionButton("fit_button", "Fit"),
                       ),
      ),
      wellPanel(
        div(class = "option-header", "Output Options"),
        radioButtons("plt_opt", "Plot options",
                     choices = c("all", "profile", "tree"),
                     # inline = TRUE
        ),
      ),
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Truth",
                 # h3("True Parameters: Tree & Response Probabilities"),

                 conditionalPanel("input.mode == 'Simulate Data' & input.sim_data_src == 'Exemplar Parameters'",
                                  fluidRow(h5("The following shows the true tree and item response probabilities for each class; the
                                             detailed values for a more complete list of simulation truths can be found at the 'Parameters' tab.")
                                  ),
                 ),
                 conditionalPanel("input.mode == 'Simulate Data' & input.sim_data_src == 'Upload Parameters'",
                                  checkboxInput("sim_instr_checkbox", "Show instructions", TRUE),
                                  conditionalPanel("input.sim_instr_checkbox == 1",
                                                   uiOutput("sim_instr"),
                                  ),
                 ),
                 fluidRow(
                   column(width = 9,
                          uiOutput("sim_plotui")
                   ),
                   column(width = 3,
                          div(class = "option-group",
                              radioButtons("sim_plot_download_format", "Select download format",
                                           c("png", "jpeg")),
                              downloadButton("download_sim_plot", "Download Plot")
                          )
                   ),
                 ),
                 # fluidRow(
                 #   div(style = "height: 30px;"),
                 #   fluidRow(
                 #     column(width = 12,
                 #            downloadButton("download_sim_data", "Download Simulated Data"),
                 #            ),
                 #   ),
                 # ),
        ),
        tabPanel("Data",
                 h3("Output Data"),
                 # checkboxInput("sim_par_checkbox", "Show Parameter Files", FALSE),
                 conditionalPanel("input.mode",
                                  # conditionalPanel("input.mode != 'Upload Posterior Samples'",
                                  fluidRow(
                                    div(class = "option-header", "Data Matrix"),
                                    # h4("Data Matrix"),
                                    column(width = 9,
                                           dataTableOutput("response_matrix_csv"),
                                           style = "overflow-y: scroll;overflow-x: scroll;"
                                    ),
                                    column(width = 3,
                                           downloadButton("download_response_matrix", "Download Data")
                                    ),
                                  ),
                 ),
                 # fluidRow(
                 #   div(class = "option-header", "Posterior Samples"),
                 #   column(width = 3,
                 #          downloadButton("download_posterior_datatab", "Download RData"),
                 #   ),
                 #   
                 # ),
                 
        ),
        tabPanel("Analysis",
                 h3("Analysis Result: Estimated Tree & Response Probabilities"),
                 h5("You may view only the tree or the response probabilities by selecting the correct 'Plot options' at the bottom of the left panel."),
                 conditionalPanel("input.mode != 'Simulate Data'",
                                  checkboxInput("an_instr_checkbox", "Show instructions", TRUE),
                                  conditionalPanel("input.an_instr_checkbox == 1",
                                                   uiOutput("an_instr"),
                                  ),
                 ), 
                 
                 fluidRow(
                   column(width = 12,
                          # verbatimTextOutput("print_fit_result")
                          textOutput("print_fit_result")
                   ),
                 ),
                 div(style = "height: 20px;"),
                 fluidRow(
                   column(width = 12, 
                          plotlyOutput("fit_plotly")
                   ),
                 ),
                 conditionalPanel("input.mode != 'Upload Posterior Samples'",
                                  div(style = "height: 200px;"),
                                  fluidRow(
                                    downloadButton("download_posterior", "Download Posterior Samples"),
                                  ),
                 ),
                 
                 
        ),
        tabPanel("Parameters",
                
                 conditionalPanel("input.mode != 'Upload Posterior Samples' ",
                                  conditionalPanel("input.mode == 'Simulate Data'",
                                                   fluidRow(h3("Values of true parameters for display and downloads:")),
                                                   fluidRow(
                                                     # div(style = "height: 50px;"),
                                                     div(class = "option-header", "Tree"),
                                                     radioButtons("tree_format", "Format",
                                                                  choices = c("txt", "csv"), inline = TRUE),
                                                     column(width = 9,
                                                            conditionalPanel("input.tree_format == 'txt'",
                                                                             verbatimTextOutput("tree_txt"),
                                                            ),
                                                            conditionalPanel("input.tree_format == 'csv'",
                                                                             dataTableOutput("tree_csv"),
                                                            ),
                                                            
                                                            style = "overflow-y: scroll;overflow-x: scroll;"
                                                     ),
                                                     column(width = 3,
                                                            downloadButton("download_tree", "Download")
                                                     ),
                                                   ),
                                                   fluidRow(
                                                     div(style = "height: 50px;"),
                                                     div(class = "option-header", "Class Probability"),
                                                     column(width = 9,
                                                            dataTableOutput("class_probability_csv"),
                                                            style = "overflow-y: scroll;overflow-x: scroll;"
                                                     ),
                                                     column(width = 3,
                                                            downloadButton("download_class_probability", "Download")
                                                     ),
                                                   ),
                                                   fluidRow(
                                                     div(style = "height: 50px;"),
                                                     div(class = "option-header", "Sigma by Group"),
                                                     column(width = 9,
                                                            dataTableOutput("sigma_by_group_csv"),
                                                            style = "overflow-y: scroll;overflow-x: scroll;"
                                                     ),
                                                     column(width = 3,
                                                            downloadButton("download_sigma_by_group", "Download")
                                                     ),
                                                   ),
                                                   div(style = "height: 50px;"),
                                  ),
                                  fluidRow(
                                    div(class = "option-header", "Item Membership List"),
                                    # h4("Item Membership List"),
                                    column(width = 9,
                                           dataTableOutput("item_membership_csv"),
                                           style = "overflow-y: scroll;overflow-x: scroll;"
                                    ),
                                    column(width = 3,
                                           downloadButton("download_item_membership", "Download")
                                    ),
                                  ),
                 ),
                 conditionalPanel("(input.mode == 'Simulate Data' & input.sim_data_src == 'Exemplar Parameters') | input.item_name_opt == 'Upload'",
                                  fluidRow(
                                    conditionalPanel("input.mode != 'Upload Posterior Samples'",
                                                     div(style = "height: 50px;"),
                                    ),
                                    div(class = "option-header", "Item Name List"),
                                    column(width = 9,
                                           dataTableOutput("item_name_csv"),
                                           style = "overflow-y: scroll;overflow-x: scroll;"
                                    ),
                                    column(width = 3,
                                           downloadButton("download_item_name", "Download Data")
                                    ),
                                  ),
                 ),
        ),
      ),
    ),
    
  ),
  
  fluidRow(
    column(12, align = "center",
           HTML("<p style='font-size: 15px; color: gray;'>
                DDT-LCM Shiny. @2023 MIT LICENSE.
                </p>")
    )
  ),
  
  # fluidRow(
  #   column(12, align = "center",
  #          HTML("<p style='font-size: 15px; color: gray;'>
  #               DDT-LCM Shiny app<br>
  #               Author: Bolin Wu, Briana Stephenson, and Zhenke Wu<br>
  #               Package Maintainer: Mengbing Li (mengbing@umich.edu)<br>
  #               Shinyapp Maintainers: Bolin Wu (bolinw@umich.edu), Zhenke Wu (zhenkewu@umich.edu)<br>
  #               Package Github Page: <a href='https://github.com/limengbinggz/ddtlcm/'>ddtlcm</a><br>
  #               </p>")
  #   )
  # ),
  # 
  # New UI design END
  
)