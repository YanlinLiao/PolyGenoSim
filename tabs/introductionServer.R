#Make the jump to the page that user want to use
observeEvent(input$link_to_tabpanel_dosage, {
  updateTabsetPanel(session, "tabset", "dosage")
})
observeEvent(input$link_to_tabpanel_SNParray, {
  updateTabsetPanel(session, "tabset", "SNParray")
})
observeEvent(input$link_to_tabpanel_SeqSim, {
  updateTabsetPanel(session, "tabset", "SeqSim")
})

url <- a("PedigreeSim", href="https://www.wur.nl/en/show/Software-PedigreeSim.htm")
output$tab <- renderUI({
  tagList("Dosage simulation: discrete dosages are generated using ", url)
})