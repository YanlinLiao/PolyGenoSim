output$SNParray <- renderUI({
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput('choice_array', 'How you would like to simulate your data?',
                    choices= c("A. From simulated dosage data",
                               "B. From uploaded dosage data file"), selected = "A. From simulated dosage data")
      ),
      mainPanel(
        tabsetPanel(type="tab",
                    tabPanel('File',
                             h4('All simulations are based on the dosage simulation results'),
                             uiOutput('choice_array_suggestion1'),
                             uiOutput('choice_array_suggestion2'),
                             
                             #From already existed dataset
                             uiOutput("check_exist_array"),
                             uiOutput('refresh_array'),
                             uiOutput("File_array_choice"),
                             uiOutput('Replicates_array'),
                             
                             
                             #From an uploaded file
                             uiOutput('File_upload_array'),
                             uiOutput('filename_text_array'),
                             uiOutput('File_upload_name_array')
                             # actionButton("folder_act_array","Continue")
                    ),
                    tabPanel('Ploidy',
                             selectInput('sit_choice_array', 'Which populations did you choose?',
                                         choices= c("A. F1 population",
                                                    "B. F2 population",
                                                    'C. Backcross population',
                                                    'D. Selfing population',
                                                    'E. Provided pedigree'), selected = "A. F1 population")
                             ,
                             uiOutput("ploidy_sug_array"),
                             uiOutput("ploidy_choi1_array"),
                             uiOutput("ploidy_choi2_array")
                    ),
                    tabPanel('Signal intensities',
                             strong('Total intensities'),
                             sliderInput("Intensity", h3("Mean of total intensity"),
                                         min = 100, max = 2000, value = 1000),
                             sliderInput("bmcv", h3("Intensity standard deviation between marker means"),
                                         min = 0, max = 0.1, value = 0.05),
                             h3('Beta distribution parameters to simulate standard deviation for all markers'),
                             sliderInput("b1", h3("Beta 1"),
                                         min = 1, max = 10, value = 2),
                             sliderInput("b2", h3("Beta 2"),
                                         min = 0, max = 10, value = 10),
                             plotOutput('TotalIntensi_plot'),

                             hr(),
                             strong('Background intensities'),
                             sliderInput("l", h3("Amount of background (signal/background)"),
                                         min = 0.5, max = 2, value = 1),
                             sliderInput("m", h3('Y proportion for background'),
                                         min = 0.3, max = 0.7, value = 0.5),
                             plotOutput('background_plot')
                             ),
                    tabPanel('Allelic bias & Overdispersion',
                             helpText("Allelic bias and overdispersion is generated from beta distribution.
                                      In order to find the proper parameter, you can use the beta distribution below to find your parameters."),
                             helpText("Beta 1 and Beta 2 are parameters used for beta distribution"),
                             # sliderInput("beta1_array", "Beta 1", 1, 10,1),
                             # sliderInput("beta2_array", "Beta 2", 1, 10,1),
                             # plotOutput("beta_dist_array"),
                             strong('Allelic bias'),
                             sliderInput("Allelic_1_array", h3("Beta1"),
                                         min = 1, max = 10, value = 10),
                             sliderInput("Allelic_2_array", h3("Beta2"),
                                         min = 1, max = 10, value = 10),
                             plotOutput('allelicarray_plot'),
                             hr(),
                             strong('Overdispersion'),
                             sliderInput("Ovd_1_array", h3("Beta1"),
                                         min = 1, max = 10, value = 5),
                             sliderInput("Ovd_2_array", h3("Beta2"),
                                         min = 1, max = 10, value = 5),
                             plotOutput('odisarray_plot')
                             ),
                    tabPanel('Simulate',
                             h4('If you are happy with the parameters setting, please click the button below:'),
                             br(),
                             actionButton("run_SNParraySim_act","Simulate SNP array"),
                             helpText("After pressing the Simulate button, 
                                      please wait until the Download button appears below while we prepare your simulation."),
                             
                             br(),
                             textOutput("SNParraySim_helptxt"),
                             uiOutput("download_array"),
                             uiOutput('redo_array'),
                             
                             hr(),
                             h4('Description of downloaded files'),
                             helpText('Your downloaded file are a zip file with the name of Date_Time_Prefix.'),
                             helpText('In the folder, it will include folder for each replicate. Under each replicate,
                                      there will be three *.csv files:'),
                             helpText('-GenotypeCallTemplate: It includes a table with columns: MarkerName, SampleName,Read Depth,
                             X,Y,ratio. The X is the reference allele intensities. The Y is the alternative allele intensities. The ratio is Y/(X+Y)'),
                             helpText('-X:A table include all reference allele intensities, where columns are individuals and rows are
                                      all markers'),
                             helpText('-Y:A table include all alternative allele intensities, where columns are individuals and rows are
                                      all markers'),

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