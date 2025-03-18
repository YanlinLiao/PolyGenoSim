output$SeqSim <- renderUI({
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput('choice_seq', 'How you would like to simulate your data?',
                    choices= c("A. From simulated dosage data",
                               "B. From uploaded dosage data file"), selected = "A. From simulated dosage data")
      ),
      
      
      mainPanel(
        tabsetPanel(type="tab",
                    tabPanel('File',
                             h4('All simulations are based on the dosage simulation results'),
                             uiOutput('choice_seq_suggestion1'),
                             uiOutput('choice_seq_suggestion2'),

                             #From already existed dataset
                             uiOutput("check_exist_seq"),
                             uiOutput("refresh_seq"),
                             uiOutput("File_Seq_choice"),
                             uiOutput('Replicates_seq'),

                             #From an uploaded file
                             uiOutput('File_upload_seq'),
                             uiOutput('filename_text_seq'),
                             uiOutput('File_upload_name')
                             # actionButton("folder_act_seq","Continue")
                             ),
                    tabPanel('Ploidy',
                             selectInput('sit_choice_seq', 'Which populations did you choose?',
                                         choices= c("A. F1 population",
                                                    "B. F2 population",
                                                    'C. Backcross population',
                                                    'D. Selfing population',
                                                    'E. Provided pedigree'), selected = "A. F1 population")
                             ,
                             uiOutput("ploidy_sug_seq"),
                             uiOutput("ploidy_choi1_seq"),
                             uiOutput("ploidy_choi2_seq")
                             ),
                    tabPanel('Read Depth',
                             h4('Ideal read depth'),
                             numericInput('rD', "Read depth",
                                          value = 120),
                             h4('Deviations can be added from individual level, marker level'),
                             h4('You can play the widgets below to check what is the range of your read depth if you 
                                choose different deviations'),
                             sliderInput("marker_effect", h3("Deviations (Marker)"),
                                         min = 0, max = 0.05, value = 0.01),
                             sliderInput("individual_effect", h3("Deviations (Individual)"),
                                         min = 0, max = 0.05, value = 0.01),
                             sliderInput("dispersion", h3("Dispersion"),
                                         min = 0, max = 20, value = 1),
                             plotOutput('rD_deviation'),
                             ),
                    tabPanel('Sequencing error',
                             sliderInput("seqE", h3("Sequencing error"),
                                         min = 0, max = 0.1, value = 0.05),
                             plotOutput('seqE_plot')
                             ),
                    tabPanel('Allelic bias & Overdispersion',
                             helpText("Allelic bias and overdispersion is generated from beta distribution.
                                      In order to find the proper parameter, you can use the beta distribution below to find your parameters."),
                             helpText("Beta 1 and Beta 2 are parameters used for beta distribution"),
                             # sliderInput("beta1", "Beta 1", 1, 10,1),
                             # sliderInput("beta2", "Beta 2", 1, 10,1),
                             # plotOutput("beta_dist"),
                             strong('Allelic bias'),
                             sliderInput("Allelic_1", h3("Beta1"),
                                         min = 1, max = 10, value = 1),
                             sliderInput("Allelic_2", h3("Beta2"),
                                         min = 1, max = 30, value = 30),
                             plotOutput('allelicseq_plot'),
                             hr(),
                             strong('Overdispersion'),
                             sliderInput("Ovd_1", h3("Beta1"),
                                         min = 1, max = 10, value = 10),
                             sliderInput("Ovd_2", h3("Beta2"),
                                         min = 1, max = 10, value = 1),
                             plotOutput('odisseq_plot')

                             ),
                    tabPanel('Simulate',
                             h4('If you are happy with the parameters setting, please click the button below:'),
                             br(),
                             actionButton("run_seqSim_act","Simulate sequence reads"),
                             helpText("After pressing the Simulate button, 
                                      please wait until the Download button appears below while we prepare your simulation."),
                             
                             br(),
                             textOutput("seqsim_helptxt"),
                             uiOutput("download_seq"),
                             uiOutput('redo_seq'),
                             
                             hr(),
                             h4('Description of downloaded files'),
                             helpText('Your downloaded file are a zip file with the name of Date_Time_Prefix.'),
                             helpText('In the folder, it will include folder for each replicate. Under each replicate,
                                      there will be three *.csv files:'),
                             helpText('-GenotypeCallTemplate: It includes a table with columns: MarkerName, SampleName,Read Depth,
                             A,a,ratio. The A is the reference allele count. The a is the alternative allele count. The ratio is a/(A+a)'),
                             helpText('-ReadCounts: A table include all reference allele count, where columns are individuals and rows are
                                      all markers'),
                             helpText('-ReadDepth:A table include all read depth, where columns are individuals and rows are
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