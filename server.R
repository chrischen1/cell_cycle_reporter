library(shiny)
library(rgl)
library(rglwidget)
options(rgl.useNULL = TRUE)

source('~/Dropbox/_ChrisProject/workspace/cell_cycle/app/test/api_cycle.R')

data_path <- "~/Dropbox/_ChrisProject/workspace/cell_cycle/rawdata/"
data_files <- list.files(data_path)
exp_data_list <- list()
for(i in 1:length(data_files)){
  x <- read.delim(paste(data_path,data_files[i],sep = ''), sep = "\t",as.is = T)
  datax <- t(x[,2:5])
  colnames(datax) <- x$Time
  exp_data_list[[gsub('Sum_table_(.+)\\.tsv','\\1',data_files[i])]] <- datax
}

shinyServer(function(input, output, session) {
  k_mat <- reactive({
    k_mat <- matrix(as.numeric(unlist(strsplit(input$k_rate,split = ','))),4,4,byrow = T)
    rownames(k_mat) <- colnames(k_mat) <- c('G1','S','G2/M','Other')
    k_mat
  })
  
  init_num <- reactive({
    init_num <-as.numeric(unlist(strsplit(input$init_count,split = ',')))
  })
  data_sim <- reactive({
    data_sim <- simulation(init_num(),k_mat(),duration = input$duration,variation = input$variation)
  })
  duration <- reactive({
    duration = input$duration
  })
  k_mat_solve <- reactive({
    k_mat_solve <- solve_mat(data_sim())
    rownames(k_mat_solve) <- colnames(k_mat_solve) <- c('G1','S','G2/M','Other')
    round(k_mat_solve,2)
  })
  data_pred <- reactive({
    data_pred <- simulation(init_num(),k_mat_solve(),input$duration)
  })
  # for tab2 
  exp_data <- reactive({
    exp_data <- data_smooth(exp_data_list[[input$well]],input$ls_span)
  })
  k_mats_exp <- reactive({
    k_mats_exp <- get_k_mat(exp_data(),u = input$rbf_span)
  })
  data_pred_exp <- reactive({
    data_pred_exp <- simulation(init_num = exp_data()[,1],k_mat = k_mats_exp(),duration = ncol(exp_data())-1)
  })

  
  colfunc <- colorRampPalette(c("red", "blue"))
  open3d()
  scene1 <- reactive({
    plot3d(data_sim()[1,],data_sim()[2,],data_sim()[3,],xlab = 'Count G1',ylab = 'Count S',zlab = 'Count G2/M',col= colfunc(ncol(data_sim())))
    lines3d(data_pred()[1,], data_pred()[2,], data_pred()[3,], col = 'green')
    scene1 <- scene3d()
  })

  scene2 <- reactive({
    plot3d(exp_data()[1,],exp_data()[2,],exp_data()[3,],col = colfunc(ncol(exp_data())),xlab = 'Count G1',ylab = 'Count S',zlab = 'Count G2/M')
    lines3d(data_pred_exp()[1,],data_pred_exp()[2,],data_pred_exp()[3,],col = 'green')
    scene2 <- scene3d()
  })
  rgl.close()
  
  save <- options(rgl.inShiny = TRUE)
  on.exit(options(save))
  output$k_pred <- renderTable({
    cbind(' '=rownames(k_mat_solve()) ,k_mat_solve())
  })
  output$k_mat2_out <- renderTable({
    cbind(' '=rownames(k_mat()) ,k_mat())
  })
  
  output$thewidget1 <- renderRglwidget(
    rglwidget(scene1())
  )
  
  output$thewidget2 <- renderRglwidget(
    rglwidget(scene2())
  )
  output$scap <- renderPlot({
    par(mfrow=c(2,2), mar=c(2,2,2,2))
    plot(1:ncol(exp_data()),apply(exp_data(),2,sum),xlab = 'Time',ylab = 'Total cell count',ylim = c(0,max(apply(exp_data(),2,sum))+200))
    plot_count(exp_data())
    
    plot_rate(k_mats_exp())
    legend('topright',c('G1-S','S-G2/M','G2/M-G1'),lty=1,col=1:3)
    
    plot_rate(k_mats_exp(),ind_list = list(c(4,1),c(4,2),c(4,3)),main = 'Transtions to "Other" state')
    legend('topright',c('G1-other','S-other','G2/M-other'),lty=1,col=c(3,1,2))
  })
})