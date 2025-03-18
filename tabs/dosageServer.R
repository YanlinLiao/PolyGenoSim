###########################Input files########################
###Population input###
situation_choice <- reactive({
  substring(as.character(input$situation_choice), 1, 1)
}) #can be A,B,C,D,E

###Ploidy input###
output$ploidy_suggestion <- renderUI({
  if(situation_choice() == "A"){
    helpText("Note: Parent 1's ploidy need to be higher than Parent 2's ploidy")
  }
})
output$ploidy_choice1 <- renderUI({
  if(situation_choice() == "A"){
    textInput("ploidy1", "P1: ploidy level", 4)
  }else{
    textInput("ploidy", "ploidy level",4)
  }
})
output$ploidy_choice2 <- renderUI({
  if(situation_choice() == "A"){
    textInput("ploidy2", "P2: ploidy level", value = 4)
  }
})

###Inheritance input###
output$inheritance_choice1a <- renderUI({
  if(situation_choice() == "A"){
    textInput('pairing1', 'P1: Degrees of preferential pairing',
              value = 0)
  }else{
    textInput('pairing', 'Degrees of preferential pairing',
              value = 0)
  }
})
output$inheritance_choice1b <- renderUI({
  if(situation_choice() == "A"){
    textInput('quadrivalents1', 'P1: Frequency of quadrivalents',
              value = 0)
  }else{
    textInput('quadrivalents', 'Frequency of quadrivalents',
              value = 0)
  }
})
output$inheritance_choice2a <- renderUI({
  if(situation_choice() == "A"){
    textInput('pairing2', 'P2: Degrees of preferential pairing',
              value = 0)
  }
})
output$inheritance_choice2b <- renderUI({
  if(situation_choice() == "A"){
    textInput('quadrivalents2', 'P2: Frequency of quadrivalents',
              value = 0)
  }
})

###Pedfile input###
output$ped_input1 <- renderUI({
  if(situation_choice() == "E"){
    helpText('To simulate population based on your own design (pedigree), 
             please upload your pedigree structure design. It needs to be text file named as *.ped.
             It should include three columns: Name, Parent1, Parent2. Each row is one individual. 
             An example is given below:')
  }
})
output$ped_input2 <- renderUI({
  if(situation_choice() == "E"){
    img(src="Ped_example.png", height = 205, width = 212)
  }
})
output$ped_input3 <- renderUI({
  if(situation_choice() == "E"){
    fileInput('ped', 'Pedfile(*.ped)',
              accept=c('.ped'))
  }else{
    numericInput('popnum', "Population size",
                 value = 200)
  }
})

###########################Define Reactive########################
##Ploidy level
#A. fixed ploidy
ploidy <- reactive({
  as.numeric(input$ploidy)
}) 
#B. cross ploidy
ploidy1 <- reactive({
  as.numeric(input$ploidy1)
}) 
ploidy2 <- reactive({
  as.numeric(input$ploidy2)
})
# output$continue_dosage <- renderText({
#   if(input$start_act != 0){
#   "Please continue"
#   }
# })

##Poptype
poptype <- reactive({
  if(situation_choice() == 'A'){
    ptype <- 'F1'
  }else if (situation_choice() == 'B'){
    ptype <- 'F2'
  }else if (situation_choice() == 'C'){
    ptype <- 'BC'
  }else if (situation_choice() == 'D'){
    ptype <- 'Selfing'
  }else{
    ptype <- NULL
  }
  ptype
})

##File parameters
filename_dsg <- reactive({
  as.character(input$filename)
}) 
popnum <- reactive({
  as.numeric(input$popnum)
}) 
replicates <-reactive({
  as.numeric(input$replicates)
}) 
mapfun <- reactive({
  as.character(input$mapfun)
}) 
ped_file <- reactive({
  inFile <- input$ped
  if(is.null(inFile)) return(NULL)
  dt <- read.table(inFile$datapath,header = TRUE)
  return(dt)
})


##Pairing parameters
pairing <- reactive({
  as.numeric(input$pairing)
}) 
quadrivalents <- reactive({
  as.numeric(input$quadrivalents)
}) 
pairing1 <- reactive({
  as.numeric(input$pairing1)
}) 
quadrivalents1 <- reactive({
  as.numeric(input$quadrivalents1)
}) 
pairing2 <- reactive({
  as.numeric(input$pairing2)
}) 
quadrivalents2 <- reactive({
  as.numeric(input$quadrivalents2)
}) 

##Chromosome parameters
chromosome <- reactive({
  as.numeric(input$chromosome)
}) 
chr_len <- reactive({
  as.numeric(input$chr_len)
}) 
center <- reactive({
  as.numeric(input$center)
}) 

##Marker parameters
#according to two parents, generate different MT
MT_choices <- reactive({
  if(situation_choice() == 'D'){
    ch <- as.character(seq(0,ploidy()))
  }else if (situation_choice() %in% c('B','C','E')){
    choices_comb <- expand.grid(seq(0,ploidy()),seq(0,ploidy()))
    ch <- paste0(choices_comb[,1],'x',choices_comb[,2])
  }else{
    choices_comb <- expand.grid(seq(0,ploidy1()),seq(0,ploidy2()))
    ch <- paste0(choices_comb[,1],'x',choices_comb[,2])
  }
  ch
}) 
output$Mrk_chrm <- renderUI({
  lapply(1:length(MT_choices()), function(i) {
    if(situation_choice() == 'D'){
      textInput(paste0('MT_',MT_choices()[i]), paste0('Dosage',MT_choices()[i]),
                value = 20)
    }else{
      textInput(paste0('MT_',MT_choices()[i]), paste0('',MT_choices()[i]),
                value = 20)
    }
  })
})
Mrk_matrix <- reactive({
  t(sapply(1:length(MT_choices()),function(i){
    c(as.numeric(strsplit(MT_choices()[i],'x')[[1]]),
    as.numeric(input[[paste0('MT_',MT_choices()[i])]]))
}))}) 

#######################Perform Simulation#######################
Perform_simulation <- eventReactive(input$simulate_act,{
  if(!"Processing" %in% dir()){
    dir.create("Processing")
  }
  for(i in 1:replicates()){
    if(situation_choice() == "A"){
       PedSim_input_complex(markermat = Mrk_matrix(),
                            LG_number = chromosome(),
                            filename = paste0(filename_dsg(),'_rep',i),
                            folder = 'Processing',
                            ploidy = ploidy1(),
                            ploidy2 = ploidy2(),
                            prefP_P1 = pairing1(),
                            prefP_P2 = pairing2(),
                            quads_P1 = quadrivalents1(),
                            quads_P2 = quadrivalents2(),
                            len = chr_len(),
                            centro = center(),
                            popnum = popnum(),
                            mapfun = mapfun(),
                            marker.positions = NULL)
    }else if (situation_choice() %in% c('B','C')){ #F2 and Backcross
      PedSim_input(markermat = Mrk_matrix(),
                   LG_number = chromosome(),
                   filename = paste0(filename_dsg(),'_rep',i),
                   ploidy = ploidy(),
                   prefP = pairing(),
                   quads = quadrivalents(),
                   len = chr_len(),
                   centro = center(),
                   popnum = popnum(),
                   mapfun = mapfun(),
                   marker.positions = NULL,
                   folder = 'Processing',
                   poptype = poptype(),
                   genfilename = NULL,
                   mapfilename = NULL,
                   chromfilename = NULL,
                   runPedigreeSim = TRUE)
    }else if(situation_choice() == 'D'){  #Selfing
      PedSim_Selfing(markermat = Mrk_matrix(),
                     LG_number = chromosome(),
                     filename = paste0(filename_dsg(),'_rep',i),
                     ploidy = ploidy(),
                     prefP = pairing(),
                     quads = quadrivalents(),
                     len = chr_len(),
                     centro = center(),
                     popnum = popnum(),
                     mapfun = mapfun(),
                     marker.positions = NULL,
                     folder = 'Processing',
                     genfilename = NULL,
                     mapfilename = NULL,
                     chromfilename = NULL,
                     runPedigreeSim = TRUE)
    }else{ #with uploaded Pedigree
      #write the pedfile to directory first
      write.table(ped_file(), file = paste0('Processing/',filename_dsg(),'_rep',i,'.ped'),quote = FALSE, row.names = FALSE, sep = "\t")
      PedSim_pedprovided(markermat = Mrk_matrix(),
                         LG_number = chromosome(),
                         filename = paste0(filename_dsg(),'_rep',i),
                         ploidy = ploidy(),
                         prefP = pairing(),
                         quads = quadrivalents(),
                         len = chr_len(),
                         centro = center(),
                         mapfun = mapfun(),
                         marker.positions = NULL,
                         folder = 'Processing',
                         genfilename = NULL,
                         mapfilename = NULL,
                         chromfilename = NULL,
                         runPedigreeSim = TRUE)
    }
    #create *.csv file
    dsg <- readDatfile(paste0(getwd(),'/Processing/',filename_dsg(),'_rep',i,'_out_alleledose.dat'))
    write.csv(dsg,file = paste0(getwd(),'/Processing/',filename_dsg(),'_rep',i,'.csv'))
    list.of.files <- list.files(paste0(getwd(),'/Processing'), pattern = paste0(filename_dsg(),'_rep',i), full.names=T)
    dir.create(paste0('Processing/',filename_dsg(),'_rep',i))
    file.copy(list.of.files, paste0('Processing/',filename_dsg(),'_rep',i),overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    unlink(list.of.files,recursive = TRUE)
  }
  folders_chosen <- list.files(paste0(getwd(),'/Processing'), pattern = filename_dsg(), full.names=T)
  
  if(!"Dosages_Simulations" %in% dir()){
    dir.create("Dosages_Simulations")
  }
  dir.create(paste0('Dosages_Simulations/',filename_dsg()))
  file.copy(folders_chosen, paste0('Dosages_Simulations/',filename_dsg()), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
  unlink(folders_chosen,recursive = TRUE)
  a <- "Simulated datasets:"
  a
})

output$simulation_helptxt <- renderText({
  Perform_simulation()
})




###################################################Download###################################################
output$download_dsg <- renderUI({
  if(length(Perform_simulation()) > 0){
    downloadButton('download_genotypes', 'Download (*.tar)')
  }
})

output$download_genotypes <- downloadHandler(
  filename <- function() {
    paste(paste0(filename_dsg()),'tar',sep = '.')
  },
  
  content <- function(file) {
    tar(file, paste0("Dosages_Simulations/",filename_dsg()))
  }
)


#######################Link to other panel#######################
observeEvent(input$link_to_SNParray, {
  updateTabsetPanel(session, "tabset", "SNParray")
})
observeEvent(input$link_to_SeqSim, {
  updateTabsetPanel(session, "tabset", "SeqSim")
})


url_dsg <- a("PedigreeSim", href="https://www.wur.nl/en/show/Software-PedigreeSim.htm")
output$pedigreesim_description <- renderUI({
  tagList("The description of downloaded files can be found through:", url_dsg)
})

output$pedigreesim_description_csv <- renderUI({
  helpText("In the downloaded folder, '*.csv' file can be used for the simulation of SNP array OR sequence reads.")
})