#some of the variables used below have been generated in the global.R
options(bitmapType='cairo')
shinyUI(navbarPage(description[,"Name"], id="tabset",
                   
                   ## tabs for application ##
                   tabPanel("Introduction", value="introduction",  uiOutput("introduction")),
                   tabPanel("Dosage Simulation", value="dosage", id = "dosage", uiOutput("dosage")),
                   tabPanel("SNP array Simulation", value="SNParray",id = "SNParray",  uiOutput("SNParray")),
                   tabPanel("Sequence Reads Simulation", value="SeqSim", id = "SeqSim", uiOutput("SeqSim")),
                   tabPanel("Simulate one marker", value="Onemarker", id = "Onemarker", uiOutput("Onemarker")),
                   ## footer ##
                   footer=footer,
                   theme = "yeti.min.css",
                   # title = a(href = "https://www.wur.nl/nl/Onderzoek-Resultaten/Onderzoeksinstituten/plant-research/Plant-Breeding.htm",
                   #           div(
                   #             img(src = "WUR_RGB_standard.png",
                   #                 height = "40px")
                   #             )
                   #           )
                   # theme = "yeti.min.css"
                   tags$style(type="text/css", css)

))

