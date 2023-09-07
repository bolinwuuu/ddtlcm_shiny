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
library(ggbrace)
# library(phylobase)

options(warn = -1)

server = function(input, output, session) {
  data(data_hchs)
  list2env(setNames(data_hchs, names(data_hchs)), envir = globalenv())
  data(data_synthetic)
  list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
  
  observeEvent(input$mode, {
    if (input$mode == "Simulate Data") {
      showTab(inputId = "fluidpage", target = "Simulation")
    } else {
      hideTab(inputId = "fluidpage", target = "Simulation")
      # showTab(inputId = "fluidpage", target = "Analysis")
    }
  })
  
  get_sim_data <- reactive({
    
    validate(
      need(input$mode == "Simulate Data",
           "Please switch to other tabs.")
    )
    
    input$read_button
    N <- 496
    if (input$sim_data_src == "Exemplar Parameters") {
      N <- 496
    } else {
      N <- isolate(as.integer(input$N))
    }
    # set random seed to generate node parameters given the tree
    seed_parameter <- isolate(input$seed_node_param)
    # set random seed to generate multivariate binary observations from LCM
    seed_response <- isolate(input$seed_bin_obs)
    
    return(simulate_lcm_given_tree(get_tree_phylo(), N,
                                   get_class_probability(),
                                   get_item_membership(),
                                   get_sigma_by_group(),
                                   root_node_location = 0,
                                   seed_parameter = seed_parameter,
                                   seed_response = seed_response))
    
  })
  
  read_uploaded_data <- function(file, header, name) {
    validate(
      need(!is.null(file), 
           paste("ERROR: Please upload", name, "file"))
    )
    ext <- tools::file_ext(file)
    
    validate(
      need(ext == "csv",
           paste("Please upload a csv file for", name))
    )
    
    data <- read.csv(file$datapath, header = header)
    
    # TODO: check format
    
    return(data)
    # return(read.csv(file$datapath, header = header))
  }
  
  get_response_matrix <- reactive({
    validate(
      need(input$mode == "Simulate Data" | input$mode == "Upload Raw Data", 
           "ERROR: response_matrix() is NOT called in Simulate Data or Upload Raw Data modes")
    )
    if (input$mode == "Simulate Data") {
      
      sim_data_res <- get_sim_data()
      
      response_matrix <<- sim_data_res$response_matrix
      return(sim_data_res$response_matrix)
      
      # return(response_matrix)
    } else {
      input$read_button
      
      data_file <- isolate(input$data_matrix_file)
      data_header <- isolate(input$data_matrix_header)
      
      uploaded_data <- read_uploaded_data(data_file,
                                          data_header,
                                          "response matrix")
      return(as.matrix(uploaded_data) + 0.0)
    }
  })
  
  get_tree_phylo <- reactive({
    validate(
      need(input$mode == "Simulate Data" | input$mode == "Upload Raw Data",
           "ERROR: tree_phylo() is NOT called in Simulate Data or Upload Raw Data modes")
    )
    if (input$mode == "Simulate Data" & input$sim_data_src == "Exemplar Parameters") {
      if (input$dataset_name == "Synthetic Data") {
        return(data_synthetic$tree_phylo)
      } else {
        return(data_hchs$tree_phylo)
      }
      
    } else {
      input$read_button
      tree_file <- isolate(input$tree_phylo_file)
      validate(
        need(!is.null(tree_file), 
             paste("ERROR: Please upload tree phylo file"))
      )
      ext <- tools::file_ext(tree_file)
      validate(
        need(ext == "csv" | ext == "txt",
             paste("Please upload a csv or a txt file for tree phylo"))
      )
      
      tree_phylo_datapath <- isolate(tree_file$datapath)
      
      return(read.tree(tree_phylo_datapath))
    }
  })
  
  get_class_probability <- reactive({
    validate(
      need(input$mode == "Simulate Data" | input$mode == "Upload Raw Data",
           "ERROR: class_probability() is NOT called in Simulate Data or Upload Raw Data modes")
    )
    if (input$mode == "Simulate Data" & input$sim_data_src == "Exemplar Parameters") {
      if (input$dataset_name == "Synthetic Data") {
        return(data_synthetic$class_probability)
      } else {
        return(data_hchs$class_probability)
      }
    } else {
      input$read_button
      
      class_prob_file <- isolate(input$class_probability_file)
      class_prob_header <- isolate(input$class_probability_header)
      return(c(as.matrix(read_uploaded_data(class_prob_file,
                                            class_prob_header,
                                            "class probability"))))
    }
  })
  
  get_sigma_by_group <- reactive({
    validate(
      need(input$mode == "Simulate Data" | input$mode == "Upload Raw Data",
           "ERROR: sigma_by_group() is NOT called in Simulate Data or Upload Raw Data modes")
    )
    if (input$mode == "Simulate Data" & input$sim_data_src == "Exemplar Parameters") {
      if (input$dataset_name == "Synthetic Data") {
        return(data_synthetic$Sigma_by_group)
      } else {
        return(data_hchs$Sigma_by_group)
      }
    } else {
      input$read_button
      
      sigma_file <- isolate(input$sigma_by_group_file)
      sigma_header <- isolate(input$sigma_by_group_header)
      return(c(as.matrix(read_uploaded_data(sigma_file,
                                            sigma_header,
                                            "sigma by group"))))
    }
    
  })
  
  get_item_membership <- reactive({
    validate(
      need(input$mode == "Simulate Data" | input$mode == "Upload Raw Data",
           "ERROR: item_membership() is NOT called in Simulate Data or Upload Raw Data modes")
    )
    if (input$mode == "Simulate Data" & input$sim_data_src == "Exemplar Parameters") {
      if (input$dataset_name == "Synthetic Data") {
        return(data_synthetic$item_membership_list)
      } else {
        return(data_hchs$item_membership_list)
      }
    } else {
      input$read_button
      
      item_memb_file <- isolate(input$item_membership_file)
      item_memb_header <- isolate(input$item_membership_header)
      
      uploaded_data <- read_uploaded_data(item_memb_file,
                                          item_memb_header,
                                          "item membership list")
      
      return(lapply(as.list(transpose(uploaded_data)),
                    function(vec) vec[!is.na(vec) & vec != ""]))
    }
  })
  
  get_response_prob <- reactive({
    validate(
      need(input$mode == "Simulate Data", 
           "ERROR: response_prob() is NOT called in Simulate Data mode")
    )
    sim_data_res <- get_sim_data()
    return(sim_data_res$response_prob)
  })
  
  
  get_tree_with_parameter <- reactive({
    validate(
      need(input$mode == "Simulate Data", 
           "ERROR: response_prob is NOT called in Simulate Data mode")
    )
    sim_data_res <- get_sim_data()
    return(sim_data_res$tree_with_parameter)
  })
  
  get_item_name_list <- reactive({
    if (input$mode == "Simulate Data" & input$sim_data_src == "Exemplar Parameters") {
      if (input$dataset_name == "Synthetic Data") {
        return(data_synthetic$item_name_list)
      } else {
        return(data_hchs$item_name_list)
      }
      
    } else {
      # item name list option displayed
      input$read_button
      
      opt <- isolate(input$item_name_opt)
      item_name_file <- isolate(input$item_name_file)
      item_name_header <- isolate(input$item_name_header)
      if (opt == "Upload" & !is.null(item_name_file)) {
        
        item_name <- read_uploaded_data(item_name_file,
                                        item_name_header,
                                        "item name list")
        return(lapply(as.list(item_name),
                      function(vec) vec[!is.na(vec) & vec != ""]))
      } else {
        item_memb <- get_result()$setting$item_membership_list
        item_name <- list()
        for (grp in seq_along(item_memb)) {
          item_name[[paste("Group", grp)]] <- paste0("Item ", item_memb[[grp]])
        }
        return(item_name)
        # return(NULL)
      }
    }
  })
  
  sim_plot <- reactive({
    
    validate(
      need(input$mode == "Simulate Data", 
           "Please switch to other tabs.")
    )
    plot_tree_with_heatmap_no_brace(get_tree_with_parameter(), get_response_prob(), get_item_membership())
  })
  
  validate_inputfiles <- reactive({
    if (input$mode == "Simulate Data" & input$sim_data_src == "Upload Parameters") {
      input$read_button
      
      tree_file <- isolate(input$tree_phylo_file)
      class_prob_file <- isolate(input$class_probability_file)
      sigma_file <- isolate(input$sigma_by_group_file)
      item_memb_file <- isolate(input$item_membership_file)
      validate(
        need(!is.null(tree_file),
             "Please upload tree phylo file."),
        need(!is.null(class_prob_file),
             "Please upload class probability file."),
        need(!is.null(sigma_file),
             "Please upload sigma by group file."),
        need(!is.null(item_memb_file),
             "Please upload item membership file."),
        need(!is.null(tree_file) & !is.null(class_prob_file) & 
               !is.null(sigma_file) & !is.null(item_memb_file),
             "Then click 'Read in' button.")
      )
    } else if (input$mode == "Upload Raw Data") {
      input$read_button
      
      data_file <- isolate(input$data_matrix_file)
      item_memb_file <- isolate(input$item_membership_file)
      validate(
        need(!is.null(data_file),
             "Please upload data matrix file."),
        need(!is.null(item_memb_file),
             "Please upload item membership file."),
        need(!is.null(data_file) & !is.null(item_memb_file),
             "Then click 'Read in' button.")
      )
    } else if (input$mode == "Upload Posterior Samples") {
      input$read_button
      
      post_samp_file <- isolate(input$posterior_sample_file)
      validate(
        need(!is.null(post_samp_file),
             "Please upload posterior sample file.\nThen click 'Read in' button.")
      )
    }
  })
  
  output$sim_plotui <- renderUI({
    plotOutput("sim_plot")
  })
  
  output$sim_plot <- renderPlot({
    validate_inputfiles()
    
    sim_plot()
  })
  
  output$download_sim_plot <- downloadHandler(
    filename = function() {
      paste("sim_data.", input$sim_plot_download_format, sep = "")
      # if(input$sim_plot_download_format == "png") {
      #   "sim_data.png"
      # } else if(input$sim_plot_download_format == "jpeg") {
      #   "sim_data.jpeg"
      # }
    },
    content = function(file) {
      if(input$sim_plot_download_format == "png") {
        # png(file)
        # plot(sim_plot())
        # dev.off()
        ggsave(file, plot = sim_plot(), device = "png")
      } else if(input$sim_plot_download_format == "jpeg") {
        # jpeg(file)
        # plot(sim_plot())
        # dev.off()
        ggsave(file, plot = sim_plot(), device = "jpeg")
      }
    }
    # contentType = function() {
      # paste("image/", input$sim_plot_download_format, sep = "")
      
    # }
  )
  
  
  
  niter <- reactive({
    input$fit_button
    
    nit <- isolate(input$niter)
    
    validate(
      need(is.integer(nit) && as.integer(nit > 0), 
           "Iterations must be a positive number")
    )
    return(nit)
  })
  
  burnin <- reactive({
    input$fit_button
    
    nit <- isolate(input$niter)
    
    validate(
      need(is.integer(input$burnin) && as.integer(input$burnin) >= 0, 
           "Burn-ins must be a non-negative integer"),
      need(as.integer(input$burnin) < as.integer(nit), 
           "Burn-ins must be smaller than Iterations")
    )
    input$burnin
  })
  
  get_result <- reactive({
    
    if (input$mode == "Simulate Data" | input$mode == "Upload Raw Data") {
      input$fit_button
      
      K_f <- isolate(input$K_fit)
      resp_mat <- get_response_matrix()
      item_memb <- get_item_membership()
      seed_post <- isolate(input$seed_posterior)
      nit <- niter()
      
      set.seed(seed_post)
      res <- ddtlcm_fit(K = K_f, data = resp_mat,
                        item_membership_list = item_memb,
                        total_iters = nit)
      
      return(res)
      
    } else {
      input$read_button
      
      posterior_file <- isolate(input$posterior_sample_file)
      
      
      validate(
        need(!is.null(posterior_file),
             "ERROR: Please upload posterior sample file")
      )
      ext <- tools::file_ext(posterior_file)
      req(posterior_file)
      validate(
        need(ext == "RData",
             "ERROR: Invalid file; Please upload a .RData file")
      )
      
      load(posterior_file$datapath, envir = .GlobalEnv)
      return(result)
      # data <- switch(ext,
      #                RData = load(input$posterior_sample_file, envir = .GlobalEnv),
      #                validate("ERROR: Invalid file; Please upload a .RData file")
      # )
      
    }
  })
  
  get_summarized_result <- reactive({
    return(summary(get_result(), burnin(), relabel = T, be_quiet = T))
  })
  
  output$print_fit_result <- renderText({
    validate_inputfiles()
    # printed_result <- capture.output(print(get_result()))
    paste(
      gsub("groups.", "groups, with ",
           str_extract(paste(capture.output(print(get_result())), collapse = ""),
                       "DDT-LCM.*samples drawn.")),
      "Posterior means and 95% credible intervals of class probabilities are displayed next to the class labels.",
      sep = "\n"
    )
  })
  
  # output$fit_plotui <- renderUI({
  #   plotOutput("fit_plot")
  # })
  # 
  # fit_plot <- reactive({
  #   # print("FIT_PLOT CALLED")
  #   # plot(x = get_summarized_result(), item_name_list = item_name_list, 
  #   #      plot_option = input$plt_opt)
  #   plot_summary_test(x = get_summarized_result(), item_name_list = item_name_list, 
  #        plot_option = input$plt_opt)
  # })
  # 
  # output$fit_plot <- renderPlot({
  #   fit_plot()
  # })
  
  plot_summary_test <- function(x, plot_option = c("all", "profile", "tree"),
                                item_name_list = NULL,
                                color_palette){
    print("PLOT CALLED")
    plot_option = match.arg(plot_option, c("all", "profile", "tree"))
    tree_with_parameter <- x$tree_map
    response_prob <- matrix(x$response_probs_summary[,"Mean"], nrow = K)
    item_membership_list <- x$setting$item_membership_list
    class_probability <- x$class_probs_summary[,"Mean"]
    plots <- plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list, 
                                    item_name_list, class_probability, color_palette,
                                    return_separate_plots = T)
    probs_lower <- x$response_probs_summary[, "2.5%"]
    probs_higher <- x$response_probs_summary[, "97.5%"]
    response_prob_dat <- plots[["response_prob_dat"]]
    response_prob_dat[, probs_lower := probs_lower]
    response_prob_dat[, probs_higher := probs_higher]
    if (plot_option == "all"){
      plot_lcm <- plots[['lcm_plot']] +     
        geom_errorbar(data = response_prob_dat,
                      aes(ymin=probs_lower, ymax=probs_higher), linewidth=0.2) 
      return(ggarrange(plots[['tree_plot']], plot_lcm, widths = c(1, 1.5), ncol = 2))
    } else if (plot_option == "profile"){
      plot_lcm <- plots[['lcm_plot']] +     
        geom_errorbar(data = response_prob_dat,
                      aes(ymin=probs_lower, ymax=probs_higher), linewidth=0.2) 
      return(plot_lcm)
      # plotly_lcm <- ggplotly(plot_lcm)
      # return(plotly_lcm)
    } else if (plot_option == "tree"){
      return(plots[['tree_plot']])
    } else {
      stop("Please select an option from 'all', 'profile', or 'tree'.")
    }
  }
  
  get_both_plots <- reactive({
    x <- get_summarized_result()

    G <- x$setting$G
    # color_palette = c("#E69F00", "#56B4E9", "#009E73",
    #                   "#000000", "#0072B2", "#D55E00",
    #                   "#CC79A7", "#F0E442", "#999999")
    color_palette <- c(brewer.pal(G, "Dark2"), brewer.pal(G, "Set1"))[1:7]
    
    tree_with_parameter <- x$tree_map
    K <- nrow(x$tree_Sigma)
    response_prob <- matrix(x$response_probs_summary[,"Mean"], nrow = K)
    item_membership_list <- x$setting$item_membership_list
    class_probability <- x$class_probs_summary[,"Mean"]
    class_probability_lower <- x$class_probs_summary[,"2.5%"]
    class_probability_higher <- x$class_probs_summary[,"97.5%"]
    
    plots <- plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list, 
                                    get_item_name_list(), class_probability, class_probability_lower, class_probability_higher, 
                                    color_palette, return_separate_plots = T)
    # plots <- plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list, 
    #                                 get_item_name_list(), class_probability, color_palette,
    #                                 return_separate_plots = T)
    return(plots)
  })
  
  fit_plotly <- reactive({
    plots <- get_both_plots()
    # x <- get_summarized_result()
    plot_option <- input$plt_opt
    if (plot_option == "tree") {
      tree_plot <- plots[['tree_plot']]
      return(get_tree_plotly(tree_plot))
    } else if (plot_option == "profile") {
      profile_plot <- get_profile_plot(plots, get_summarized_result())
      return(get_profile_plotly(profile_plot))
    } else {
      tree_plotly <- get_tree_plotly(plots[['tree_plot']])
      profile_plotly <- get_profile_plotly(get_profile_plot(plots, 
                                                            get_summarized_result())) |>
        layout(
          xaxis = list(
            tickfont = list(
              size = 6
            )
          )
        )
      plotly_all <- subplot(tree_plotly, profile_plotly,
                            titleY = TRUE, titleX = TRUE, 
                            margin = 0.04
      ) |> layout(
        title = list(text = ""),
        # add titles for each plot
        annotations = list( 
          list(x = 0.24, 
               y = 1.1, 
               text = "Tree", 
               showarrow = F, 
               xref='paper', 
               yref='paper'), 
          list(x = 0.83, 
               y = 1.1, 
               text = "Class Profiles", 
               showarrow = F, 
               xref='paper', 
               yref='paper')) 
      )
      
      # disable tree hover info
      for (i in seq_along(tree_plotly$x$data)) {
        plotly_all$x$data[[i]]$hoverinfo <- "none"
      }
      
      return(plotly_all)
    }
  })
  
  get_tree_plotly <- function(tree_plot) {
    K <- sum(tree_plot$data$isTip)
    tree_plotly <- ggplotly(tree_plot) |>
      layout(title = list(text = "Tree"))
    
    # adjust label positions
    tree_plotly$x$data[[length(tree_plotly$x$data)]]$y <- tree_plotly$x$data[[length(tree_plotly$x$data)]]$y + 0.1
    # disable hover info
    tree_plotly$x$layout$hovermode <- FALSE
    
    # add labels for leaf vertices 
    vertex_text <- rep(" ", times = length(tree_plotly$x$data[[1]]$x))
    leaf_indices <- which(round(tree_plotly$x$data[[1]]$x, digits = 5) == 1)
    leaf_indices <- leaf_indices[order(-tree_plotly$x$data[[1]]$y[leaf_indices])]
    vertex_text[leaf_indices] <- paste(" v", 1:K, sep = "")
    tree_plotly$x$data[[1]]$text <- vertex_text
    tree_plotly$x$data[[1]]$mode <- "text"
    return(tree_plotly)
  }
  
  get_profile_plot <- function(plots, x) {
    probs_lower <- x$response_probs_summary[, "2.5%"]
    probs_higher <- x$response_probs_summary[, "97.5%"]
    response_prob_dat <- plots[["response_prob_dat"]]
    response_prob_dat[, probs_lower := probs_lower]
    response_prob_dat[, probs_higher := probs_higher]
    plot_lcm <- plots[['lcm_plot']] +     
      geom_errorbar(data = response_prob_dat,
                    aes(ymin=probs_lower, ymax=probs_higher), linewidth=0.2)
    return(plot_lcm)
  }
  
  get_profile_plotly <- function(profile_plot) {
    
    K <- profile_plot$plot_env$K
    G <- profile_plot$plot_env$G
    
    profile_plot <- profile_plot +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    profile_plotly <- ggplotly(profile_plot, height = 600) |>
      layout(
        title = list(
          text = "Class Profiles",
          y = 1, x = 0.5, xanchor = 'center', yanchor =  'top'
        ),
        legend = list(
          itemclick = FALSE,
          itemdoubleclick = FALSE,
          groupclick = FALSE
        ),
        xaxis = list(
          # tickangle = 45,
          title = list(
            text = "Items",
            standoff = 20,
            font = list(
              size = 15
            )
          ),
          ticktext = get_tick_text(get_item_name_list(), 
                                   unique(profile_plot$data$color)),
          tickfont = list(
            size = 8
          )
        )
      )
    
    # add y-axis title
    if (K < 3) {
      profile_plotly$x$layout[["yaxis"]]$title <- list(
        text = "Item Response Probability",
        font = list(
          size = 15
        ),
        standoff = 10
      )
    } else {
      profile_plotly$x$layout[[paste("yaxis", ceiling(K / 2), sep = "")]]$title <- list(
        text = "Item Response Probability",
        font = list(
          size = 15
        ),
        standoff = 10
      )
    }
    
    # customize legend
    for (i in seq(1, G * K, K)) {
      profile_plotly$x$data[[i]]$name <- names(get_item_name_list())[ceiling(i / K)]
    }
    for (i in setdiff(1:(G * K), seq(1, G * K, K))) {
      profile_plotly$x$data[[i]]$showlegend <- FALSE
    }
    
    # customize error bars
    for (i in seq(G * K + 1, length(profile_plotly$x$data))) {
      profile_plotly$x$data[[i]]$error_y$thickness <- 0.7
      profile_plotly$x$data[[i]]$error_y$width <- 2
    }
    
    profile_plotly <- customize_errorbar(profile_plotly, G, K, 3, 0.7)
    
    return(profile_plotly)
    
  }
  
  customize_errorbar <- function(profile_plotly, G, K, width, thickness) {
    for (i in seq(G * K + 1, length(profile_plotly$x$data))) {
      profile_plotly$x$data[[i]]$error_y$width <- width
      profile_plotly$x$data[[i]]$error_y$thickness <- thickness
    }
    return(profile_plotly)
  }
  
  
  # fit_plotly <- reactive({
  #   # ggplotly(fit_plot())
  #   
  #   profile_plot <- fit_plot() +
  #     theme(axis.title.x = element_blank())
  #   
  #   profile_plotly <- ggplotly(profile_plot, height = 600) |>
  #     layout(
  #       title = list(
  #         y = 1, x = 0.5, xanchor = 'center', yanchor =  'top'
  #       ),
  #       legend = list(
  #         itemclick = FALSE,
  #         itemdoubleclick = FALSE,
  #         groupclick = FALSE
  #       ),
  #       xaxis = list(
  #         # tickangle = 45,
  #         title = list(
  #           text = "Items",
  #           standoff = 20,
  #           font = list(
  #             size = 15
  #           )
  #         ),
  #         ticktext = get_tick_text(item_name_list),
  #         tickfont = list(
  #           size = 10
  #         )
  #       )
  #     )
  #   
  #   # customize legend
  #   G <- length(item_name_list)
  # 
  #   for (i in seq(1, G * K, K)) {
  #     profile_plotly$x$data[[i]]$name <- names(item_name_list)[ceiling(i / K)]
  #   }
  # 
  #   for (i in setdiff(1:(G * K), seq(1, G * K, K))) {
  #     profile_plotly$x$data[[i]]$showlegend <- FALSE
  #   }
  #   
  #   # customize error bars
  #   for (i in seq(G * K + 1, length(profile_plotly$x$data))) {
  #     profile_plotly$x$data[[i]]$error_y$thickness <- 0.7
  #     profile_plotly$x$data[[i]]$error_y$width <- 3
  #   }
  # 
  #   return(profile_plotly)
  # })
  
  get_tick_text <- function(item_name_list, color_palette) {
    # color_palette <- c("#E69F00", "#56B4E9", "#009E73",
    #                  "#000000", "#0072B2", "#D55E00",
    #                  "#CC79A7", "#F0E442", "#999999")
    tick_text <- c()
    for (grp_idx in 1:length(item_name_list)) {
      for (itm in item_name_list[[grp_idx]]) {
        tick_text <- c(tick_text,
                       paste("<span style='color:", color_palette[grp_idx], "'>", itm, "</span>",
                             sep = ""))
        # tick_text <- c(tick_text,
        #                paste("<span style='color:", color_palette[grp_idx], "; font-weight:900'>", itm, "</span>",
        #                      sep = ""))
      }
    }
    return(tick_text)
  }
  
  output$fit_plotly <- renderPlotly({
    validate_inputfiles()
    fit_plotly()
  })
  
  plot_tree_with_heatmap_no_brace <- function(tree_with_parameter, response_prob, item_membership_list){
    K <- phylobase::nTips(tree_with_parameter)
    J <- length(unlist(item_membership_list))
    if (ncol(tree_with_parameter@data) != J){
      stop("The input list 'item_membership_list' should contain indices for all columns of
         node parameters on the tree.")
    }
    branch = branch.length = NULL # Setting the variables to NULL first to avoid NSE notes in R CMD check
    t1 <- ggtree(tree_with_parameter) +
      geom_tiplab(size = 4) +
      geom_nodelab(size = 4) +
      labs(title = paste0("Tree")) +
      theme(plot.title = element_text(hjust = 0.5, margin=margin(b=0), vjust = -1, size = 10)) +
      coord_cartesian(clip="off") +
      # geom_brace(aes_(x=c(0, 1), y=c(0.4, 0), label = paste("Group", 1)),
      #            inherit.data=F, rotate=180, labelsize=3.5, color = "white") +
      theme(plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
      geom_text(aes(x=branch, label=signif(branch.length, 2)), #size=6, 
                vjust=-.3, color = "firebrick")
    # leaf_data <- expit(as.matrix(tree_with_parameter@data[as.character(1:K), 1:J]))
    dat_response_prob <- data.table(value = c(response_prob))
    item = y = value = group_label = NULL # due to NSE notes in R CMD check
    
    dat_response_prob[, class := rep(1:K, J)]
    dat_response_prob[, item := rep(1:J, each = K)]
    
    # reorder rows of the heatmap to match the tree leaves
    class_order <- data.table(t1$data[t1$data$label %in% paste0("v", 1:K), c("label", "y")])
    class_order <- class_order[order(-y),]
    
    class_order <- as.integer(gsub("v", "", class_order$label))
    
    dat_response_prob$class <- factor(dat_response_prob$class, levels = class_order)
    
    plot_prob <- ggplot(data = dat_response_prob, aes(x = item, y = class, fill = value)) +
      # geom_tile(hjust = 0.5, vjust = 0.5, interpolate = FALSE, na.rm = TRUE) +
      geom_tile() +
      labs(x = "Item", y = "Class", fill = "Item Response Probabilities") +
      scale_fill_gradientn(limits = c(0,1), colours=c("gold", "springgreen4", "blue2") ) + 
      theme(plot.title = element_text(hjust = 0.5, vjust = 10, size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position="bottom",
            plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
      coord_cartesian(clip = "off")
    
    G <- length(item_membership_list)
    dat_brace <- c()
    for (g in 1:G) {
      dat_brace <- c(dat_brace, item_membership_list[[g]][1], item_membership_list[[g]][length(item_membership_list[[g]])])
    }
    dat_brace <- data.table(x = dat_brace)
    dat_brace[, y := rep(c(0.4, 0), G)]
    dat_brace[, group_label := rep(1:G, each = 2)]
    plot_prob + geom_brace(aes(dat_brace$x, dat_brace$y, label = group_label), inherit.data=F)
    for (g in 1:G) {
      suppressWarnings({
        plot_prob <- plot_prob +
          geom_brace(aes_(x=c(item_membership_list[[g]][1], item_membership_list[[g]][length(item_membership_list[[g]])]),
                          y=c(0.4, 0), label = g), 
                     inherit.data=F, rotate=180, labelsize=3.5)
      })
    }
    ggarrange(t1, plot_prob, widths = c(0.2, 0.3), ncol = 2, common.legend = T)
  }
  
  # output$download_fit_plot <- downloadHandler(
  #   filename = function() {
  #     paste("fit_data", input$fit_plot_download_format, sep = "")
  #   },
  #   content = function(file) {
  #     if(input$fit_plot_download_format == ".png") {
  #       png(file)
  #       plot(fit_plot())
  #       dev.off()
  #     } else if(input$fit_plot_download_format == ".jpeg") {
  #       jpeg(file)
  #       plot(fit_plot())
  #       dev.off()
  #     }
  #   }
  # )
  
  output$sim_instr <- renderUI({
    wellPanel(
      h4("Instructions:"),
      p("Shortcut: select 'Exemplar Parameters' and switch to 'Parameters' and 'Data' tabs to view and download file templates!"),
      p("The following 4 files should be prepared:"),
      h5("Tree phylo"),
      p("Prepare a csv or txt file containing a tree in Newick format (parenthetic format). An example of a txt file:"),
      p("(((v5:0.14,(v4:0.12,v2:0.12)u6:0.016)u5:0.051,(v6:0.19,(v3:0.16,v1:0.16)u2:0.024)u3:0.0051)u4:0.81)u1;"),
      h5("Class Probability List"),
      p("Prepare a csv file "),
    )
  })
  
  output$an_instr <- renderText({
    if (input$mode == "Simulate Data") {
      
      txt <- "analysis instructions in simulation mode!\ntest newline"
      return(txt)
    } else if (input$mode == "Upload Raw Data") {
      return("analysis instructions in raw data mode!")
    } else {
      return("analysis instructions in posterior sample mode!")
    }
    
  })
  
  output$response_matrix_csv <- renderDataTable({
    get_response_matrix()
  })
  
  # shinyDirChoose(input, "response_matrix_dir", roots = c(home = '~'), filetypes = c('', 'csv'))
  
  output$download_response_matrix <- downloadHandler(
    filename = function() {
      "data_matrix.csv"
    },
    content = function(file) {
      write.csv(get_response_matrix(), file, row.names = FALSE)
    }
  )
  
  output$download_sim_data <- downloadHandler(
    filename = function() {
      "data_matrix.csv"
    },
    content = function(file) {
      write.csv(get_response_matrix(), file, row.names = FALSE)
    }
  )
  
  output$item_membership_csv <- renderDataTable({
    item_memb <- get_item_membership()
    max_length <- max(sapply(item_memb, length))
    padded_item_memb <- lapply(item_memb, function(x) 
      c(x, rep("", max_length - length(x))))
    item_memb_df <- as.data.frame(do.call(rbind, padded_item_memb))
    col_names <- as.character(item_memb_df[1,])
    item_memb_df <- item_memb_df[-1,]
    datatable(item_memb_df, rownames = FALSE, colnames = col_names)
  })
  
  output$download_item_membership <- downloadHandler(
    filename = function() {
      "item_membership_list.csv"
    },
    content = function(file) {
      file.create(file)
      lapply(get_item_membership(), function(x) 
        write.table(t(x), file, append= T, sep=',', col.names = F, row.names = F))
    }
  )
  
  output$tree_txt <- renderText({
    write.tree(get_tree_phylo())
  })
  
  output$tree_csv <- renderDataTable({
    tree_str <- write.tree(get_tree_phylo())
    split_tree <- strsplit(tree_str, ",")
    tree_df <- data.frame(t(unlist(split_tree)))

    datatable(tree_df, colnames = rep("", length(tree_df)))
  })
  
  output$download_tree <- downloadHandler(
    filename = function() {
      paste("tree_phylo.", input$tree_format, sep = "")
    },
    content = function(file) {
      write.tree(get_tree_phylo(), file, append = FALSE, digits = 10, tree.names = FALSE)
    }
  )
  
  output$class_probability_csv <- renderDataTable({
    datatable(data.frame(get_class_probability()), rownames = FALSE, colnames = "x")
  })
  
  output$download_class_probability <- downloadHandler(
    filename = function() {
      "class_probability.csv"
    },
    content = function(file) {
      write.csv(get_class_probability(), file, row.names = FALSE)
    }
  )
  
  output$sigma_by_group_csv <- renderDataTable({
    datatable(data.frame(get_sigma_by_group()), rownames = FALSE, colnames = "x")
  })
  
  output$download_sigma_by_group <- downloadHandler(
    filename = function() {
      "Sigma_by_group.csv"
    },
    content = function(file) {
      write.csv(get_sigma_by_group(), file, row.names = FALSE)
    }
  )
  
  output$item_name_csv <- renderDataTable({
    datatable(filled_item_name(), rownames = FALSE)
  })
  
  output$download_item_name <- downloadHandler(
    filename = function() {
      "item_name_list.csv"
    },
    content = function(file) {
      write.csv(filled_item_name(), file, row.names = FALSE)
    }
  )
  
  filled_item_name <- reactive({
    item_names <- get_item_name_list()
    max_length <- max(sapply(item_names, function(x) length(unlist(x))))
    return(data.frame(lapply(item_names, function(x) c(unlist(x), rep("", max_length - length(unlist(x)))))))
  })
  
  output$download_posterior <- downloadHandler(
    filename = function() {
      "posterior_samples.RData"
    },
    content = function(file) {
      res <- get_result()
      save(res, file = file)
    }
  )
  
  output$download_posterior_datatab <- downloadHandler(
    filename = function() {
      "posterior_samples.RData"
    },
    content = function(file) {
      res <- get_result()
      save(res, file = file)
    }
  )
  
}