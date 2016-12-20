library(shiny)
library(rglwidget)

data_path <- "~/Dropbox/_ChrisProject/workspace/cell_cycle/rawdata/"
data_files <- list.files(data_path)
exp_data_list <- list()
for(i in 1:length(data_files)){
  x <- read.delim(paste(data_path,data_files[i],sep = ''), sep = "\t",as.is = T)
  datax <- t(x[,2:5])
  colnames(datax) <- x$Time
  exp_data_list[[gsub('Sum_table_(.+)\\.tsv','\\1',data_files[i])]] <- datax
}

shinyUI(fluidPage(
  mainPanel(
    tabsetPanel(
      tabPanel("Simulation",
               sidebarPanel(h4("Parameters setup"),
                            textInput('init_count', "G1,S,G2/M,Other(separated by comma)", '120,110,110,20'),
                            textInput("k_rate", "Rates in transition matrix(by row):", '0.84,0,0.2,0.01,0.15,0.84,0,0.01,0,0.15,0.89,0.01,0.01,0.01,0.01,0.97'),
                            numericInput("variation", "Variation for simulation in each step:", 0.6),
                            numericInput("duration", "Number of steps to simulate:", 168),
                            
                            submitButton("Update View")),
               mainPanel(
                 h4("Input parameters for transition rate:"),
                 tableOutput("k_mat2_out"),
                 h4("Predicted parameters for transition rate:"),
                 tableOutput("k_pred"),
                 
                 rglwidgetOutput('thewidget1')
                )
               ),
      tabPanel("Experiment",
               sidebarPanel(h4("Parameters setup"),
                            selectInput("well", "Choose a well:", choices = c(names(exp_data_list))),
                            sliderInput("rbf_span", "Size of window(as span in RBF):", min = 0, max = 50, value = 5, step= 0.1),
                            numericInput("ls_span", "Span of loess Smoothing on data:",0.1),
                             
                            submitButton("Update View")),
               mainPanel(
                 rglwidgetOutput('thewidget2'),
                 plotOutput("scap")
               )
               

              )
    )
  )
))

