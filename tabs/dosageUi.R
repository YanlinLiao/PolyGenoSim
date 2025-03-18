library(shiny)
library(polymapR)

output$dosage <- renderUI({
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        strong("Please specify your need:"),

        selectInput('situation_choice', 'Which populations would you like to simulate?',
                    choices= c("A. F1 population",
                               "B. F2 population",
                               'C. Backcross population',
                               'D. Selfing population',
                               'E. Provided pedigree'), selected = "A. F1 population"),
        uiOutput("ploidy_suggestion"),
        uiOutput("ploidy_choice1"),
        uiOutput("ploidy_choice2")
        
        # actionButton("start_act","Continue"),
        # textOutput("continue_dosage")
      ),
      
      
      mainPanel(
        tabsetPanel(type="tab",
                    tabPanel("Population",
                             ####File####
                             helpText("*Filename is the prefix of the generated files"),
                             textInput('filename', 'Filename prefix',
                                       value = "Test"),
                             
                             numericInput("replicates","Number of replicates",
                                          value = 1),
                             uiOutput('ped_input1'),
                             uiOutput('ped_input2'),
                             uiOutput('ped_input3'),
                             
                             selectInput("mapfun",label = "Mapping function", choices = c("haldane","kosambi"), selected = "haldane")
                        
                      
                             ),
                    tabPanel('Pairing',
                             ####Pairing####
                             helpText("*Degrees of preferential pairing needs to vary between 0 - 1."),
                             helpText(" When it is 0: random pairing"),
                             uiOutput("inheritance_choice1a"),
                             uiOutput("inheritance_choice2a"),
                             helpText("*Frequency of quadrivalents needs to vary between 0 - 1."),
                             helpText(" When it is 0: bivalent formation"),
                             uiOutput("inheritance_choice1b"),
                             uiOutput("inheritance_choice2b")
                             ),
                    tabPanel('Chromosomes',
                             ####Chromosome####
                             numericInput('chromosome', 'Number of chromosomes',
                                          value = 1),
                             numericInput('chr_len', 'Length of chromosome (cM)',
                                          value = 100),
                             numericInput('center', 'Location of centromere (cM)',
                                          value = 50)
                             ),
                    tabPanel('Markers',
                             strong('Markers Parameters'),
                             helpText("Please scroll down to specify the number of each type of markers"),
                             uiOutput('Mrk_chrm')
                    ),
                    tabPanel('Simulate',
                             h4('If you are happy with the parameters setting, please click the button below:'),
                             br(),
                             actionButton("simulate_act","Simulate genotypes"),
                             helpText("After pressing the Simulate button, 
                                      please wait until the Download button appears below while we prepare your simulation."),
                             
                             br(),
                             textOutput("simulation_helptxt"),
                             uiOutput("download_dsg"),
                             
                             hr(),
                             h4('Description of downloaded files'),
                             uiOutput('pedigreesim_description'),
                             uiOutput('pedigreesim_description_csv'),
                             
                             hr(),
                             strong("SNP array OR Sequence reads"),
                             helpText("If you would like to generate SNP array or Sequence reads, please use your simulated
                                      genotypes to continue in the page 'SNP array Simulation' or 'Sequence reads simulation'.
                                      You can choose one of the link below."),
                             actionLink("link_to_SNParray", "SNP array simulation"),
                             br(),
                             actionLink("link_to_SeqSim", "Sequence reads simulation"),
                             
                            
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