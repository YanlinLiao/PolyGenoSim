###########################Upload information########################
###For SNP array
output$TotalIntesi_mean <- renderUI({
  if(choice_onemarker() == "A"){ 
    numericInput("totalintensi_mean", h3("mean of total intensities"), value = 1924)
  }
})

output$TotalIntesi_sd <- renderUI({
  if(choice_onemarker() == "A"){ 
    sliderInput("totalintensi_sd", h3("standard deviation of total intensities"), 
              min = 0, max = 0.1, value = 0.08)
  }
})

output$Background_propotion <- renderUI({
  if(choice_onemarker() == "A"){ 
    sliderInput("background_p", h3("Amount of background (signal/background)"),min = 0.5, max = 3, value = 2.5)
  }
})

output$Background_Y <- renderUI({
  if(choice_onemarker() == "A"){ 
    sliderInput("background_y",h3("Y proportion for background"),min = 0.1, max = 0.9, value = 0.16)
  }
})


###For sequence reads
output$rD_mean <- renderUI({
  if(choice_onemarker() == "B"){ 
    numericInput("rd_mean", h3("mean of read depth"), value = 60)
  }
})

output$rD_dispersion <- renderUI({
  if(choice_onemarker() == "B"){ 
    sliderInput("rd_dispersion", h3("dispersion of read depth"), min = 0, max = 20,value = 1)
  }
})


output$seqE_onemarker <- renderUI({
  if(choice_onemarker() == "B"){ 
    sliderInput("seqE_one", h3("sequencing error"), min = 0, max = 0.1, value = 0.01)
  }
})


###########################Define Reactive########################
#uploaded file
dsgfile_onemarker <- reactive({
  inFile <- input$dsg_onemarker
  if(is.null(inFile)) return(NULL)
  dt <- read.csv(inFile$datapath)
  rownames(dt) <- dt$marker
  dt <- dt[,-1]
  dt <- as.matrix(dt)
  class(dt) <- "integer"
  return(dt)
})
upload_onemarker_fnme <- reactive({
  as.character(input$upload_onemarker_fnme)
})


###Choice pop input
choice_onemarker <- reactive({
  substring(as.character(input$choice_onemarker), 1, 1)
}) #can be A,B

ploidy_onemarker <- reactive({
  as.numeric(input$ploidy_onemarker)
})

od_onemarker <- reactive({
  as.numeric(input$od_onemarker)
})

aleB_onemarker <- reactive({
  as.numeric(input$aleB_onemarker)
})

###SNP array
totalintensi_mean <- reactive({
  as.numeric(input$totalintensi_mean)
})
totalintensi_sd <- reactive({
  as.numeric(input$totalintensi_sd) * totalintensi_mean()
})
background_p <- reactive({
  as.numeric(input$background_p)
})
background_y <- reactive({
  as.numeric(input$background_y)
})

###sequence reads
rd_mean <- reactive({
  as.numeric(input$rd_mean)
})
rd_dispersion <- reactive({
  as.numeric(input$rd_dispersion)
})
seqE_one <- reactive({
  as.numeric(input$seqE_one)
})


##########################RunSimulation########################
onemarker_foldername <- eventReactive(input$run_onemarker_act,{
  prefix <-  gsub(' ','_',gsub(':','',gsub('-','',gsub(' CET','',Sys.time()))))
  f <- paste0(prefix,'_',upload_onemarker_fnme())
  f
})


Run_onemarker <- eventReactive(input$run_onemarker_act,{
  if(choice_onemarker() == 'A'){   #A: SNP array
    if(!"SNParray_Simulations" %in% dir()){
      dir.create("SNParray_Simulations")
    }
    dir.create(paste0("SNParray_Simulations/",onemarker_foldername()))
    results <- GenoSim_oneMarker(dsg = dsgfile_onemarker(),
                                 choice = "A", # A.SNP array, B.sequence reads
                                 ploidy = ploidy_onemarker(), 
                                 ploidy2 = ploidy_onemarker(),
                                 seed_number = 1,
                                 od =od_onemarker(),
                                 ale_b = aleB_onemarker(), #allelic bias (ale_b = Y/X) [0.2-0.8]
                                 avint = totalintensi_mean(), # for total intensities
                                 mcv = totalintensi_sd(), # for total intensities
                                 l = background_p(), #amount of background (l = sign/backg) [0.5-2]
                                 m = background_y(), #balance between marker (SNP array)
                                 ideal_depth,  #ideal read depth
                                 size,# for depth
                                 seq_e)
    for(sheetname in names(results)[1:3]){
      sheet_dat <- results[[sheetname]]
      wb_path <- paste0("SNParray_Simulations/",onemarker_foldername(),"/",sheetname,".csv")
      write.csv(sheet_dat,wb_path)
    }
    write.csv(dsgfile_onemarker(),paste0("SNParray_Simulations/",onemarker_foldername(),"/dsg.csv"))
    }else{
      if(!"SequenceReads_Simulations" %in% dir()){
        dir.create("SequenceReads_Simulations")
      }
      
      dir.create(paste0("SequenceReads_Simulations/",onemarker_foldername()))
      results <- GenoSim_oneMarker(dsg = dsgfile_onemarker(),
                                   choice = "B", # A.SNP array, B.sequence reads
                                   ploidy = ploidy_onemarker(), 
                                   ploidy2 = ploidy_onemarker(),
                                   seed_number = 1,
                                   od =od_onemarker(),
                                   ale_b = aleB_onemarker(), #allelic bias (ale_b = Y/X) [0.2-0.8]
                                   avint, # for total intensities
                                   mcv, # for total intensities
                                   l, #amount of background (l = sign/backg) [0.5-2]
                                   m, #balance between marker (SNP array)
                                   ideal_depth = rd_mean(),  #ideal read depth
                                   size = rd_dispersion(),# for depth
                                   seq_e = seqE_one())
      for(sheetname in names(results)[1:3]){
        sheet_dat <- results[[sheetname]]
        wb_path <- paste0("SequenceReads_Simulations/",onemarker_foldername(),"/",sheetname,".csv")
        write.csv(sheet_dat,wb_path)
      }
      write.csv(dsgfile_onemarker(),paste0("SequenceReads_Simulations/",onemarker_foldername(),"/dsg.csv"))
    }
  a <- "Simulated datasets:"
  a 
})

output$onemarker_helptxt <- renderText({
  Run_onemarker()
})



###################################################Download###################################################
output$download_onemarker <- renderUI({
  if(length(Run_onemarker()) > 0){
    downloadButton('download_onemarker_y', 'Download (*.tar)')
  }
})

download_foldernme_onemarker <- eventReactive(input$run_onemarker_act,{
  f <- onemarker_foldername()
  f
})

output$download_onemarker_y <- downloadHandler(
  filename <- function() {
    paste(paste0(download_foldernme_onemarker()),'tar',sep = '.')
  },
  
  content <- function(file) {
    if(choice_onemarker() == 'A'){
      tar(file, paste0("SNParray_Simulations/",download_foldernme_onemarker()))
    }else{
      tar(file, paste0("SequenceReads_Simulations/",download_foldernme_onemarker()))
    }
   
  }
)


output$redo_onemarker <- renderUI({
  if(length(Run_onemarker()) > 0){
    helpText('If you would like to perform another simulation, please refresh the page and 
          follow the instructions again.')
  }
})



###################################################Output###################################################
get_instant_results <- reactive({
  if(choice_onemarker() == 'A'){   #A: SNP array
    results <- GenoSim_oneMarker(dsg = dsgfile_onemarker(),
                                 choice = "A", # A.SNP array, B.sequence reads
                                 ploidy = ploidy_onemarker(), 
                                 ploidy2 = ploidy_onemarker(),
                                 seed_number = 1,
                                 od =od_onemarker(),
                                 ale_b = aleB_onemarker(), #allelic bias (ale_b = Y/X) [0.2-0.8]
                                 avint = totalintensi_mean(), # for total intensities
                                 mcv = totalintensi_sd(), # for total intensities
                                 l = background_p(), #amount of background (l = sign/backg) [0.5-2]
                                 m = background_y(), #balance between marker (SNP array)
                                 ideal_depth,  #ideal read depth
                                 size,# for depth
                                 seq_e)
  }else{
    results <- GenoSim_oneMarker(dsg = dsgfile_onemarker(),
                                 choice = "B", # A.SNP array, B.sequence reads
                                 ploidy = ploidy_onemarker(), 
                                 ploidy2 = ploidy_onemarker(),
                                 seed_number = 1,
                                 od =od_onemarker(),
                                 ale_b = aleB_onemarker(), #allelic bias (ale_b = Y/X) [0.2-0.8]
                                 avint, # for total intensities
                                 mcv, # for total intensities
                                 l, #amount of background (l = sign/backg) [0.5-2]
                                 m, #balance between marker (SNP array)
                                 ideal_depth = rd_mean(),  #ideal read depth
                                 size = rd_dispersion(),# for depth
                                 seq_e = seqE_one())
  }
  plot_temp <- results$GenotypeCallTemplate
  plot_temp
})

output$show_onemarker <- renderPlot(
  plot(get_instant_results()[,4],get_instant_results()[,5],
       xlab = 'reference allele', ylab = 'alternative allele')
)


