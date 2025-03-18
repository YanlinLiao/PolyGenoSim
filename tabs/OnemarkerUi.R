output$Onemarker <- renderUI({
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput('choice_onemarker', 'What would you like to simulate?',
                    choices= c("A. SNP array",
                               "B. Sequence reads"), selected = "A. SNP array")
      ),
      
      
      mainPanel(
        tabsetPanel(type="tab",
                    tabPanel('File',
                             h4('All simulations are based on the dosage simulation results'),
                             img(src="example_table1.png", height = 61, width = 247),
                             helpText('The table needs to have only one marker'),
                             fileInput('dsg_onemarker', 'Dosage table(*.csv)',
                                       accept=c('.csv')),
                             helpText('Please specify the foldername that you want to use to 
                                      save your simulated dataset'),
                             textInput("upload_onemarker_fnme", "Folder name", 'Test_OneMarker'),
                    ),
                    tabPanel('Parameters',
                             helpText('If you would like to know more about each parameter, please
                                      check the corresponding tab'),
                             numericInput('ploidy_onemarker',h3('ploidy'), value = 4),
                             
                             numericInput("od_onemarker", h3("overdispersion"), value = 0.03),
                             numericInput("aleB_onemarker", h3("allelic bias"), value = 0.43),
                             
                             #SNP array
                             uiOutput('TotalIntesi_mean'),
                             uiOutput('TotalIntesi_sd'),
                             uiOutput('Background_propotion'),
                             uiOutput('Background_Y'),
                             
                             #sequence reads
                             uiOutput('rD_mean'),
                             uiOutput('rD_dispersion'),
                             uiOutput('seqE_onemarker'),
                             
                             hr(),
                             plotOutput('show_onemarker')
                             
                    ),
                    tabPanel('Simulate',
                             h4('If you are happy with the parameters setting, please click the button below:'),
                             br(),
                             actionButton("run_onemarker_act","Perform Simulation"),
                             helpText("After pressing the Simulate button, 
                                      please wait until the Download button appears below while we prepare your simulation."),
                             
                             br(),
                             textOutput("onemarker_helptxt"),
                             uiOutput("download_onemarker"),
                             uiOutput('redo_onemarker'),
                             
                             hr(),
                             h4('Description of downloaded files'),
                             helpText('Please refer to previous tabs'),
                             
                             #########
                             tagList(
                               tags$head(
                                 tags$link(rel="stylesheet", type="text/css",href="style.css"),
                                 tags$script(type="text/javascript", src = "busy.js")
                               )
                             ),
                             div(class = "busy",  
                                 p("Calculation in progress..."), 
                                 img(src="ajax-loader.gif")
                             )
                             #########
                    )
        )
        
      )
    )
  )
})