###########################Define Reactive########################
###Choice pop input
choice_seq <- reactive({
  substring(as.character(input$choice_seq), 1, 1)
}) #can be A,B

output$choice_seq_suggestion1 <- renderUI({
  if(choice_seq() == "A"){ 
    #If they would like to use already simulated files
    helpText("Note: the files you would like to use need to be included in the folder 'Dosage_Simulations' ")
  }else{
    #If they want to upload a dosage table
    helpText("Note: your uploaded file need to be a matrix of dosage with markernames as row 
             and individual name as columns. Here is an example:")
  }
})
output$choice_seq_suggestion2 <- renderUI({
  if(choice_seq() == "B"){
    img(src="example_table.png", height = 121, width = 247)
  }
})
#A: choose the correct folder
output$check_exist_seq <- renderUI({
  if(choice_seq() == "A"){
  if(!"Dosages_Simulations" %in% dir()){
    actionLink("link_to_dosages", "Perform dosage simulation first!")
  }
  }
})

output$refresh_seq <- renderUI({
  if(choice_seq() == "A"){
    helpText('If you do not find your simulated dataset. Please refresh the page (ctrl + R)')
  }
})
observeEvent(input$link_to_dosages, {
  updateTabsetPanel(session, "tabset", "dosage")
})
output$File_Seq_choice <- renderUI({
  if(choice_seq() == "A"){
    if("Dosages_Simulations" %in% dir()){
      files_prefix <- dir("Dosages_Simulations")
      files_prefix <- setdiff(files_prefix,c('lib','PedigreeSim.jar'))
      selectInput("File_seq","Please choose the file prefix that you want to use for sequence reads simulation", 
                  choices = files_prefix,
                  selected = files_prefix[1])
    }
  }
})
output$Replicates_seq <- renderUI({
  if(choice_seq() == "A"){
    files_list <- dir(paste0("Dosages_Simulations/",File_seq()))
    helpText(paste0('There are ',length(files_list) ,' replicates in your chosen folder!'))
  }
})

#B: upload the dosage table
output$File_upload_seq <- renderUI({
  if(choice_seq() == "B"){
    fileInput('dsg_seq', 'Dosage table(*.csv)',
              accept=c('.csv'))
  }
})
output$File_upload_name <- renderUI({
  if(choice_seq() == "B"){
    textInput("upload_seq_fnme", "Folder name", 'Test_seq')
  }
})

###Ploidy section
sit_choice_seq <- reactive({
  substring(as.character(input$sit_choice_seq), 1, 1)
}) #can be A,B,C,D,E
output$ploidy_sug_seq <- renderUI({
  if(sit_choice_seq() == "A"){
    helpText("Note: Parent 1's ploidy need to be higher than Parent 2's ploidy")
  }
})
output$ploidy_choi1_seq <- renderUI({
  if(sit_choice_seq() == "A"){
    textInput("ploidy1_seq", "P1: ploidy level", 4)
  }else{
    textInput("ploidy_seq", "ploidy level","4")
  }
})
output$ploidy_choi2_seq <- renderUI({
  if(sit_choice_seq() == "A"){
    textInput("ploidy2_seq", "P2: ploidy level", value = 4)
  }
})

###########################Define Reactive########################
#uploaded file
dsgfile_seq <- reactive({
  inFile <- input$dsg_seq
  if(is.null(inFile)) return(NULL)
  dt <- read.csv(inFile$datapath)
  rownames(dt) <- dt$marker
  dt <- dt[,-1]
  dt <- as.matrix(dt)
  class(dt) <- "integer"
  return(dt)
})
upload_seq_fnme <- reactive({
  as.character(input$upload_seq_fnme)
})
#infolder file
File_seq <- reactive({
  input$File_seq
})
files_list_seq <- reactive({
  dir(paste0("Dosages_Simulations/",File_seq()))
})

##Ploidy level
#A. fixed ploidy
ploidy_seq <- reactive({
  as.numeric(input$ploidy_seq)
}) 
#B. cross ploidy
ploidy1_seq <- reactive({
  as.numeric(input$ploidy1_seq)
}) 
ploidy2_seq <-  reactive({
  as.numeric(input$ploidy2_seq)
}) 


###read depth
rD <- reactive({
  as.numeric(input$rD)
})
marker_effect <- reactive({
  as.numeric(input$marker_effect)
})
individual_effect <- reactive({
  as.numeric(input$individual_effect)
})
dispersion <- reactive({
  as.numeric(input$dispersion)
})
depth_for_plot <- reactive({
  V <- rnorm(200,0,marker_effect()) # marker level effect
  W <- rnorm(200,0,individual_effect())  # individual effect
  s_V <- quantile(V,probs = 0.95)
  s_W <- quantile(W,probs = 0.95)
  
  bt <- log(rD())      # overall mean
  l <- sapply(V, function(x){
    exp(bt + x + W)
  }) 
  
  dep <- t(sapply(1:length(V), function(i){
    rnbinom(n = nrow(l), # per marker
            mu = l[,i], #get the mean for each individual
            size = dispersion()) # dispersion in overall
  }))
  
  return(as.numeric(dep))
})

###sequencing error
fake_dsg <- reactive({
  if(sit_choice_seq() == "A"){
    ploidy_plot <- (ploidy1_seq() + ploidy2_seq())/2
  }else{
    ploidy_plot <- ploidy_seq()
  }
  dsg <- as.numeric(replicate(50,seq(0,ploidy_plot)))
  p <- dsg/ploidy_plot
  return(p)
})
seqE_for_plot <- reactive({
  p <- fake_dsg() * (1 - seqE()) + (1 - fake_dsg()) * seqE()
  rD <- rbinom(length(p), rD(),0.85)
  reference <- p * rD
  alternative <- (1-p) * rD
  return(list('ref' = reference,
              'alt' = alternative))
})
seqE <- reactive({
  as.numeric(input$seqE)
})


###Overdispersion + allelic bias
#for plot
beta1 <- reactive({
  as.numeric(input$beta1)
})
beta2 <- reactive({
  as.numeric(input$beta2)
})

#A/a unbalancement: allelic bias
Allelic_1 <- reactive({
  as.numeric(input$Allelic_1)
})
Allelic_2 <- reactive({
  as.numeric(input$Allelic_2)
})
allelic_for_plot <- reactive({
  aleB <- rbeta(length(fake_dsg()),Allelic_1(),Allelic_2())
  aleBdisset <- t(replicate(1,aleB))
  pu <- fake_dsg()/ (aleBdisset * (1.0 - fake_dsg()) + fake_dsg()) #A/A+B
  rD <- rbinom(length(pu), rD(),0.85)
  reference <- pu * rD
  alternative <- (1-pu) * rD
  return(list('ref' = reference,
              'alt' = alternative))
  
})

#overdispersion
Ovd_1 <- reactive({
  as.numeric(input$Ovd_1)
})
Ovd_2 <- reactive({
  as.numeric(input$Ovd_2)
})
adjust <- function(mu,rho){
  alpha <- mu * (1-rho)/rho
  beta <- (1 - mu) * (1-rho)/rho
  adjusttt <- rbeta(length(mu),shape1 = alpha,shape2 = beta) #random generation of beta distribution
  return(adjusttt)
}
ovd_for_plot <- reactive({
  overdis <- rbeta(length(fake_dsg()),Ovd_1(),Ovd_2())
  overdisset <- t(replicate(1,overdis/50))
  para <- matrix(adjust(mu = fake_dsg(),rho = overdisset),nrow = 1)
  rD <- rbinom(length(para), rD(),0.85)
  reference <- para * rD
  alternative <- (1-para) * rD
  return(list('ref' = reference,
              'alt' = alternative))
})
output$filename_text_seq <- renderUI({
  if(choice_seq() == 'B'){
    helpText('Please specify the foldername that you want to use to 
             save your simulated dataset')
  }
})


##########################RunSimulation########################
seq_foldername <- eventReactive(input$run_seqSim_act,{
  prefix <-  gsub(' ','_',gsub(':','',gsub('-','',gsub(' CET','',Sys.time()))))
  if(choice_seq() == 'A'){
    f <- paste0(prefix,'_',File_seq())
  }else{
    f <- paste0(prefix,'_',upload_seq_fnme())
  }
  f
})


Run_SeqSim <- eventReactive(input$run_seqSim_act,{
  if(!"SequenceReads_Simulations" %in% dir()){
    dir.create("SequenceReads_Simulations")
  }
  if(choice_seq() == 'A'){   #A: from chosen folder: with replicates
    list.of.files <- list.files(paste0("Dosages_Simulations/",File_seq(),'/'), pattern = File_seq(), full.names=T)
    for(i in list.of.files){
      dir.create(paste0("SequenceReads_Simulations","/",gsub('(.+)(/)(.+)(/)(.+)','\\5',i)))
      file <- readDatfile(paste0(i,"/",gsub('(.+)(/)(.+)(/)(.+)','\\5',i),"_out_alleledose.dat"))
      rownames(file) <- file$marker
      file$marker <- NULL
      file <- as.matrix(file)
      if(is.null(ploidy1_seq())){
        results <- GenoSim(dsg = file,
                           choice = "B", # A.SNP array, B.sequence reads
                           ploidy = ploidy_seq(), 
                           ploidy2 = NULL,
                           seed_number = NULL,
                           od = c(Ovd_1(),Ovd_2()), #overdispersion from beta distribution
                           ale_b = c(Allelic_1(),Allelic_2()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                           avint = NULL, # for total intensities
                           bmcv = NULL, # for total intensities
                           b = NULL, # for total intensities
                           l = NULL, #amount of background (l = sign/backg) [0.5-2]
                           m = NULL, #balance between marker (SNP array)
                           ideal_depth = rD(),  #ideal read depth
                           mrk_effect =  marker_effect(), # for depth
                           individual_effect =  individual_effect(),# for depth
                           dispersion = dispersion(),# for depth
                           seq_e = seqE())
      }else{
        results <- GenoSim(dsg = file,
                           choice = "B", # A.SNP array, B.sequence reads
                           ploidy = ploidy1_seq(), 
                           ploidy2 = ploidy2_seq(),
                           seed_number = NULL,
                           od = c(Ovd_1(),Ovd_2()), #overdispersion from beta distribution
                           ale_b = c(Allelic_1(),Allelic_2()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                           avint = NULL, # for total intensities
                           bmcv = NULL, # for total intensities
                           b = NULL, # for total intensities
                           l = NULL, #amount of background (l = sign/backg) [0.5-2]
                           m = NULL, #balance between marker (SNP array)
                           ideal_depth = rD(),  #ideal read depth
                           mrk_effect =  marker_effect(), # for depth
                           individual_effect =  individual_effect(),# for depth
                           dispersion = dispersion(),# for depth
                           seq_e = seqE())
      }
      for(sheetname in names(results)[1:3]){
        sheet_dat <- results[[sheetname]]
        wb_path <- paste0("SequenceReads_Simulations","/",gsub('(.+)(/)(.+)(/)(.+)','\\5',i),"/",sheetname,".csv")
        write.csv(sheet_dat,wb_path)
      }
    }
    list.of.files_move <- list.files("SequenceReads_Simulations/", pattern = paste0(File_seq(),'_'), full.names=T)
    # dir.create(paste0("SequenceReads_Simulations/",File_seq()))
    # file.copy(list.of.files_move, paste0("SequenceReads_Simulations/",File_seq()),overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    dir.create(paste0("SequenceReads_Simulations/",seq_foldername()))
    file.copy(list.of.files_move, paste0("SequenceReads_Simulations/",seq_foldername()),overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    unlink(list.of.files_move,recursive = TRUE)
  }else{
    #B: from upload file: without replicates
    dir.create(paste0("SequenceReads_Simulations/",seq_foldername()))
    if(is.null(ploidy1_seq())){
      results <- GenoSim(dsg = dsgfile_seq(),
                         choice = "B", # A.SNP array, B.sequence reads
                         ploidy = ploidy_seq(), 
                         ploidy2 = NULL,
                         seed_number = NULL,
                         od = c(Ovd_1(),Ovd_2()), #overdispersion from beta distribution
                         ale_b = c(Allelic_1(),Allelic_2()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                         avint = NULL, # for total intensities
                         bmcv = NULL, # for total intensities
                         b = NULL, # for total intensities
                         l = NULL, #amount of background (l = sign/backg) [0.5-2]
                         m = NULL, #balance between marker (SNP array)
                         ideal_depth = rD(),  #ideal read depth
                         mrk_effect =  marker_effect(), # for depth
                         individual_effect =  individual_effect(),# for depth
                         dispersion = dispersion(),# for depth
                         seq_e = seqE())
    }else{
      results <- GenoSim(dsg = dsgfile_seq(),
                         choice = "B", # A.SNP array, B.sequence reads
                         ploidy = ploidy1_seq(), 
                         ploidy2 = ploidy2_seq(),
                         seed_number = NULL,
                         od = c(Ovd_1(),Ovd_2()), #overdispersion from beta distribution
                         ale_b = c(Allelic_1(),Allelic_2()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                         avint = NULL, # for total intensities
                         bmcv = NULL, # for total intensities
                         b = NULL, # for total intensities
                         l = NULL, #amount of background (l = sign/backg) [0.5-2]
                         m = NULL, #balance between marker (SNP array)
                         ideal_depth = rD(),  #ideal read depth
                         mrk_effect =  marker_effect(), # for depth
                         individual_effect =  individual_effect(),# for depth
                         dispersion = dispersion(),# for depth
                         seq_e = seqE())
    }
   
    for(sheetname in names(results)[1:3]){
      sheet_dat <- results[[sheetname]]
      wb_path <- paste0("SequenceReads_Simulations/",seq_foldername(),"/",sheetname,".csv")
      write.csv(sheet_dat,wb_path)
    }
    write.csv(dsgfile_seq(),paste0("SequenceReads_Simulations/",seq_foldername(),"/dsg.csv"))
    
  }
  a <- "Simulated datasets:"
  a 
})

output$seqsim_helptxt <- renderText({
  Run_SeqSim()
})


##########################Output########################
#Deviations effect
output$rD_deviation <- renderPlot(
  hist(depth_for_plot(),col = "#75AADB", border = "white",
       main = "Effect of deviations on read Depth",
       xlab = 'Read depth',
       ylab = 'Frequency')
)

#SeqE effect
output$seqE_plot <- renderPlot(
  plot(seqE_for_plot()$ref,
       seqE_for_plot()$alt,
       xlab = 'Reference allele count', ylab = 'Alternative allele count',
       xlim = c(0,rD() + 10),ylim = c(0, rD() + 10),
       main = 'Effect of sequencing error'
       )
)

#allelic effect
output$allelicseq_plot <- renderPlot(
  plot(allelic_for_plot()$ref,
       allelic_for_plot()$alt,
       xlab = 'Reference allele count', ylab = 'Alternative allele count',
       xlim = c(0,rD() + 10),ylim = c(0, rD() + 10),
       main = 'Effect of allelic bias'
  )
)

#overdispersion effect
output$odisseq_plot <- renderPlot(
  plot(ovd_for_plot()$ref,
       ovd_for_plot()$alt,
       xlab = 'Reference allele count', ylab = 'Alternative allele count',
       xlim = c(0,rD() + 10),ylim = c(0, rD() + 10),
       main = 'Effect of overdispersion'
  )
)


#beta distribution
output$beta_dist <- renderPlot(
  hist(rbeta(10000,beta1(),beta2()),col = "#75AADB", border = "white",xlab = paste0("Beta(",beta1(),",",beta2(),")"), main = "Histogram of Beta distribution")
)

output$redo_seq <- renderUI({
  if(length(Run_SeqSim()) > 0){
    helpText('If you would like to perform another simulation, please refresh the page and 
          follow the instructions again.')
  }
})


###################################################Download###################################################
output$download_seq <- renderUI({
  if(length(Run_SeqSim()) > 0){
    downloadButton('download_Sequence', 'Download (*.tar)')
  }
})

download_foldernme <- eventReactive(input$run_seqSim_act,{
  f <- seq_foldername()
  f
})

output$download_Sequence <- downloadHandler(
  filename <- function() {
    paste(paste0(download_foldernme()),'tar',sep = '.')
  },
  
  content <- function(file) {
    tar(file, paste0("SequenceReads_Simulations/",download_foldernme()))
  }
)


