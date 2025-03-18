###Choice pop input
choice_array <- reactive({
  substring(as.character(input$choice_array), 1, 1)
}) #can be A,B

output$choice_array_suggestion1 <- renderUI({
  if(choice_array() == "A"){ 
    #If they would like to use already simulated files
    helpText("Note: the files you would like to use need to be included in the folder 'Dosage_Simulations' ")
  }else{
    #If they want to upload a dosage table
    helpText("Note: your uploaded file need to be a matrix of dosage with markernames as row 
             and individual name as columns. Here is an example:")
  }
})
output$choice_array_suggestion2 <- renderUI({
  if(choice_array() == "B"){
    img(src="example_table.png", height = 121, width = 247)
  }
})
#A: choose the correct folder
output$check_exist_array <- renderUI({
  if(choice_array() == "A"){
    if(!"Dosages_Simulations" %in% dir()){
      actionLink("link_to_dosages_array", "Perform dosage simulation first!")
    }
  }
})
observeEvent(input$link_to_dosages_array, {
  updateTabsetPanel(session, "tabset", "dosage")
})
output$File_array_choice <- renderUI({
  if(choice_array() == "A"){
    if("Dosages_Simulations" %in% dir()){
      files_prefix <- dir("Dosages_Simulations")
      files_prefix <- setdiff(files_prefix,c('lib','PedigreeSim.jar'))
      selectInput("File_array","Please choose the file prefix that you want to use for sequence reads simulation", 
                  choices = files_prefix,
                  selected = files_prefix[1])
    }
  }
})
output$Replicates_array <- renderUI({
  if(choice_array() == "A"){
    files_list <- dir(paste0("Dosages_Simulations/",File_array()))
    helpText(paste0('There are ',length(files_list) ,' replicates in your chosen folder!'))
  }
})
output$refresh_array <- renderUI({
  if(choice_array() == "A"){
    helpText('If you do not find your simulated dataset. Please refresh the page (ctrl + R)')
  }
})


#B: upload the dosage table
output$File_upload_array <- renderUI({
  if(choice_array() == "B"){
    fileInput('dsg_array', 'Dosage table(*.csv)',
              accept=c('.csv'))
  }
})
output$File_upload_name_array <- renderUI({
  if(choice_array() == "B"){
    textInput("upload_array_fnme", "Folder name", 'Test_array')
  }
})

###Ploidy section
sit_choice_array <- reactive({
  substring(as.character(input$sit_choice_array), 1, 1)
}) #can be A,B,C,D,E
output$ploidy_sug_array <- renderUI({
  if(sit_choice_array() == "A"){
    helpText("Note: Parent 1's ploidy need to be higher than Parent 2's ploidy")
  }
})
output$ploidy_choi1_array <- renderUI({
  if(sit_choice_array() == "A"){
    textInput("ploidy1_array", "P1: ploidy level", 4)
  }else{
    textInput("ploidy_array", "ploidy level",4)
  }
})
output$ploidy_choi2_array <- renderUI({
  if(sit_choice_array() == "A"){
    textInput("ploidy2_array", "P2: ploidy level", value = 4)
  }
})

###########################Define Reactive########################
#uploaded file
dsgfile_array <- reactive({
  inFile <- input$dsg_array
  # if(is.null(inFile)) return(NULL)
  dt <- read.csv(inFile$datapath)
  rownames(dt) <- dt$marker
  dt <- dt[,-1]
  dt <- as.matrix(dt)
  # remove all markers having a missing value in one of the parents
  # dt<-dt[!is.na(dt[,"P1"]) & !is.na(dt[,"P2"]),]
  # convert to integer
  class(dt) <- "integer"
  return(dt)
})

output$array_table <- renderUI({
  if(choice_array() == "B"){
    renderTable(dsgfile_array()[,1:5])
  }
})

upload_array_fnme <- reactive({
  as.character(input$upload_array_fnme)
})
#infolder file
File_array <- reactive({
  input$File_array
})
files_list_array <- reactive({
  dir(paste0("Dosages_Simulations/",File_array()))
})

##Ploidy level
#A. fixed ploidy
ploidy_array <- reactive({
  as.numeric(input$ploidy_array)
}) 
#B. cross ploidy
ploidy1_array <- reactive({
  as.numeric(input$ploidy1_array)
}) 
ploidy2_array <-  reactive({
  as.numeric(input$ploidy2_array)
}) 

##Total intensities
Intensity <- reactive({
  as.numeric(input$Intensity)
})
bmcv <- reactive({ #variation
  as.numeric(input$bmcv)
})
b1 <- reactive({
 as.numeric(input$b1)
})
b2 <- reactive({
  as.numeric(input$b2)
})
intensity_for_plot <- reactive({
  mav <- rnorm(1, Intensity(), Intensity()*bmcv())
  mcv <- rbeta(1, b1(), b2())
  n <- matrix(rnorm(1*1000, mav, mav*mcv), nrow = 1)
  return(n)
})


##fake dsg
fake_dsg_array <- reactive({
  if(sit_choice_array() == "A"){
    ploidy_plot <- (ploidy1_array() + ploidy2_array())/2
  }else{
    ploidy_plot <- ploidy_array()
  }
  dsg <- as.numeric(replicate(50,seq(0,ploidy_plot)))
  p <- dsg/ploidy_plot
  return(p)
})

##background intensities
l <- reactive({
  as.numeric(input$l)
})
m <- reactive({
  as.numeric(input$m)
})
background_for_plot <- reactive({
  n <- rnorm(length(fake_dsg_array()), Intensity(), Intensity()*0.05)
  b <- n/(l()+1) # background intensity
  s <- n - b             # signal intensity
  Ys <- fake_dsg_array() * s
  Y <- Ys + rbinom(n = length(fake_dsg_array()), size = round(b), prob = m())
  # obatin X
  X <- n - Y
  list('ref' = X,
       'alt' = Y)
})



###Overdispersion + allelic bias
#for plot
beta1_array <- reactive({
  as.numeric(input$beta1_array)
})
beta2_array <- reactive({
  as.numeric(input$beta2_array)
})

#A/a unbalancement: allelic bias
Allelic_1_array <- reactive({
  as.numeric(input$Allelic_1_array)
})
Allelic_2_array <- reactive({
  as.numeric(input$Allelic_2_array)
})
allelic_for_plot_array <- reactive({
  aleB <- rbeta(length(fake_dsg_array()),Allelic_1_array(),Allelic_2_array())
  aleBdisset <- t(replicate(1,aleB))
  pu <- fake_dsg_array()/ (aleBdisset * (1.0 - fake_dsg_array()) + fake_dsg_array()) #A/A+B
  intensity <- rnorm(length(fake_dsg_array()), Intensity(), Intensity()*0.05) 
  reference <- pu * intensity
  alternative <- (1-pu) * intensity
  return(list('ref' = reference,
              'alt' = alternative))
  
})

#overdispersion
Ovd_1_array <- reactive({
  as.numeric(input$Ovd_1_array)
})
Ovd_2_array <- reactive({
  as.numeric(input$Ovd_2_array)
})
adjust <- function(mu,rho){
  alpha <- mu * (1-rho)/rho
  beta <- (1 - mu) * (1-rho)/rho
  adjusttt <- rbeta(length(mu),shape1 = alpha,shape2 = beta) #random generation of beta distribution
  return(adjusttt)
}
ovd_for_plot_array <- reactive({
  overdis <- rbeta(length(fake_dsg_array()),Ovd_1_array(),Ovd_2_array())
  overdisset <- t(replicate(1,overdis/50))
  para <- matrix(adjust(mu = fake_dsg_array(),rho = overdisset),nrow = 1)
  intensity <- rnorm(length(fake_dsg_array()), Intensity(), Intensity()*0.05) 
  reference <- para * intensity
  alternative <- (1-para) * intensity
  return(list('ref' = reference,
              'alt' = alternative))
})
output$filename_text_array <- renderUI({
  if(choice_array() == 'B'){
    helpText('Please specify the foldername that you want to use to 
             save your simulated dataset')
  }
})


##########################RunSimulation########################
array_foldername <- eventReactive(input$run_SNParraySim_act,{
  prefix <-  gsub(' ','_',gsub(':','',gsub('-','',gsub(' CET','',Sys.time()))))
  if(choice_array() == 'A'){
    f <- paste0(prefix,'_',File_array())
  }else{
    f <- paste0(prefix,'_',upload_array_fnme())
  }
  f
})

Run_SNParraySim <- eventReactive(input$run_SNParraySim_act,{
  if(!"SNParray_Simulations" %in% dir()){
    dir.create("SNParray_Simulations")
  }
  if(choice_array() == 'A'){   #A: from chosen folder: with replicates
    list.of.files <- list.files(paste0("Dosages_Simulations/",File_array(),'/'), pattern = File_array(), full.names=T)
    for(i in list.of.files){
      dir.create(paste0("SNParray_Simulations","/", gsub('(.+)(/)(.+)(/)(.+)','\\5',i)))
      file <- readDatfile(paste0(i,"/", gsub('(.+)(/)(.+)(/)(.+)','\\5',i),"_out_alleledose.dat"))
      rownames(file) <- file$marker
      file$marker <- NULL
      file <- as.matrix(file)
      if(is.null(ploidy1_array())){
        results <- GenoSim (dsg = file,
                            choice = "A", # A.SNP array, B.sequence reads
                            ploidy = ploidy_array(), 
                            ploidy2 = NULL,
                            seed_number = NULL,
                            od = c(Ovd_1_array(),Ovd_2_array()), #overdispersion from beta distribution
                            ale_b = c(Allelic_1_array(),Allelic_2_array()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                            avint = Intensity(), # for total intensities
                            bmcv = bmcv(), # for total intensities
                            b = c(b1(),b2()), # for total intensities
                            l = l(), #amount of background (l = sign/backg) [0.5-2]
                            m = m(), #balance between marker (SNP array)
                            ideal_depth = NULL,  #ideal read depth
                            mrk_effect = NULL, # for depth
                            individual_effect = NULL,# for depth
                            dispersion = NULL,# for depth
                            seq_e = NULL)
      }else{
          results <- GenoSim (dsg = file,
                              choice = "A", # A.SNP array, B.sequence reads
                              ploidy = ploidy1_array(), 
                              ploidy2 = ploidy2_array(),
                              seed_number = NULL,
                              od = c(Ovd_1_array(),Ovd_2_array()), #overdispersion from beta distribution
                              ale_b = c(Allelic_1_array(),Allelic_2_array()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                              avint = Intensity(), # for total intensities
                              bmcv = bmcv(), # for total intensities
                              b = c(b1(),b2()), # for total intensities
                              l = l(), #amount of background (l = sign/backg) [0.5-2]
                              m = m(), #balance between marker (SNP array)
                              ideal_depth = NULL,  #ideal read depth
                              mrk_effect = NULL, # for depth
                              individual_effect = NULL,# for depth
                              dispersion = NULL,# for depth
                              seq_e = NULL)
      }
      #write to files
      for(sheetname in names(results)[1:3]){
        sheet_dat <- results[[sheetname]]
        wb_path <- paste0("SNParray_Simulations","/",gsub('(.+)(/)(.+)(/)(.+)','\\5',i),"/",sheetname,".csv")
        write.csv(sheet_dat,wb_path)
      }
    }
    list.of.files_move <- list.files("SNParray_Simulations/", pattern = paste0(File_array(),'_'), full.names=T)
    # dir.create(paste0("SNParray_Simulations/",File_array()))
    # file.copy(list.of.files_move, paste0("SNParray_Simulations/",File_array()),overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    dir.create(paste0("SNParray_Simulations/",array_foldername()))
    file.copy(list.of.files_move, paste0("SNParray_Simulations/",array_foldername()),overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)
    unlink(list.of.files_move,recursive = TRUE)
  }else{
    #B: from upload file: without replicates
    dir.create(paste0("SNParray_Simulations/",array_foldername()))
    if(is.null(ploidy2_array())){
      results <-  GenoSim (dsg = dsgfile_array(),
                           choice = "A", # A.SNP array, B.sequence reads
                           ploidy = ploidy_array(), 
                           ploidy2 = NULL,
                           seed_number = NULL,
                           od = c(Ovd_1_array(),Ovd_2_array()), #overdispersion from beta distribution
                           ale_b = c(Allelic_1_array(),Allelic_2_array()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                           avint = Intensity(), # for total intensities
                           bmcv = bmcv(), # for total intensities
                           b = c(b1(),b2()), # for total intensities
                           l = l(), #amount of background (l = sign/backg) [0.5-2]
                           m = m(), #balance between marker (SNP array)
                           ideal_depth = NULL,  #ideal read depth
                           mrk_effect = NULL, # for depth
                           individual_effect = NULL,# for depth
                           dispersion = NULL,# for depth
                           seq_e = NULL)
    }else{
      results <-  GenoSim (dsg = dsgfile_array(),
                           choice = "A", # A.SNP array, B.sequence reads
                           ploidy = ploidy1_array(), 
                           ploidy2 = ploidy2_array(),
                           seed_number = NULL,
                           od = c(Ovd_1_array(),Ovd_2_array()), #overdispersion from beta distribution
                           ale_b = c(Allelic_1_array(),Allelic_2_array()), #allelic bias (ale_b = Y/X+Y) [0.2-0.8]
                           avint = Intensity(), # for total intensities
                           bmcv = bmcv(), # for total intensities
                           b = c(b1(),b2()), # for total intensities
                           l = l(), #amount of background (l = sign/backg) [0.5-2]
                           m = m(), #balance between marker (SNP array)
                           ideal_depth = NULL,  #ideal read depth
                           mrk_effect = NULL, # for depth
                           individual_effect = NULL,# for depth
                           dispersion = NULL,# for depth
                           seq_e = NULL)
    }
    #write to files
    for(sheetname in names(results)[1:3]){
      sheet_dat <- results[[sheetname]]
      wb_path <- paste0("SNParray_Simulations/",array_foldername(),"/",sheetname,".csv")
      write.csv(sheet_dat,wb_path)
    }
    write.csv(dsgfile_array(),paste0("SNParray_Simulations/",array_foldername(),"/dsg.csv"))
    
  }
  a <- "Simulated datasets:"
  a 
})

output$SNParraySim_helptxt <- renderText({
  Run_SNParraySim()
})

##########################Output########################
#Total intensities
output$TotalIntensi_plot <- renderPlot(
  hist(intensity_for_plot(),col = "#75AADB", border = "white",
       main = "Effect of variations on total intensities",
       xlab = 'Total intensities',
       ylab = 'Frequency',nclass = 15)
)

#Background effect
output$background_plot <- renderPlot(
  plot(background_for_plot()$ref,
       background_for_plot()$alt,
       xlab = 'Reference allele intensity', ylab = 'Alternative allele intensity',
       xlim = c(0,Intensity()+200),ylim = c(0, Intensity()+200),
       main = 'Effect of background intensities'
  )
)

#allelic effect
output$allelicarray_plot <- renderPlot(
  plot(allelic_for_plot_array()$ref,
       allelic_for_plot_array()$alt,
       xlab = 'Reference allele intensity', ylab = 'Alternative allele intensity',
       xlim = c(0,Intensity()+200),ylim = c(0, Intensity()+200),
       main = 'Effect of allelic bias'
  )
)

#overdispersion effect
output$odisarray_plot <- renderPlot(
  plot(ovd_for_plot_array()$ref,
       ovd_for_plot_array()$alt,
       xlab = 'Reference allele intensity', ylab = 'Alternative allele intensity',
       xlim = c(0,Intensity()+200),ylim = c(0, Intensity()+200),
       main = 'Effect of overdispersion'
  )
)

#beta distribution
output$beta_dist_array <- renderPlot(
  hist(rbeta(10000,beta1_array(),beta2_array()),col = "#75AADB", border = "white",
       xlab = paste0("Beta(",beta1_array(),",",beta2_array(),")"), 
       main = "Histogram of Beta distribution")
)

output$redo_array <- renderUI({
  if(length(Run_SNParraySim()) > 0){
    helpText('If you would like to perform another simulation, please refresh the page and 
          follow the instructions again.')
  }
})


###################################################Download###################################################
output$download_array <- renderUI({
  if(length(Run_SNParraySim()) > 0){
    downloadButton('download_SNParray', 'Download (*.tar)')
  }
})

download_foldername <- eventReactive(input$run_SNParraySim_act,{
  f <- array_foldername()
  f
})

output$download_SNParray <- downloadHandler(
  filename <- function() {
    paste(paste0(download_foldername()),'tar',sep = '.')
  },
  
  content <- function(file) {
    tar(file, paste0("SNParray_Simulations/",download_foldername()))
  }
)
