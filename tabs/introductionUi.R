library(shiny)
library(polymapR)

output$introduction <- renderUI({
  dashboardBody(
    fluidRow(
    box(title = "Dosage Simulation", background = "black", solidHeader = TRUE,width = 4,
        br(),
        uiOutput("tab"),
        br(),
        actionLink("link_to_tabpanel_dosage", "Simulate dosage")),
    box(title = "SNP array Simulation", background = "black", solidHeader = TRUE,width = 4,
        br(),
        "SNP array simulation: SNP array signals are generated from the discrete dosages.", 
        br(),
        actionLink("link_to_tabpanel_SNParray", "Simulate SNP array")),
    
    box(title = "Sequence Reads Simulation", background = "black", solidHeader = TRUE,width = 4,
        br(),
        "Sequence reads simulation: counts of sequence reads are generated from the discrete dosages", 
        br(),
        actionLink("link_to_tabpanel_SeqSim", "Simulate sequence reads"))
           ))

})