#specify how big the input file can accept
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  # On startup, read the URL parameter
  queryParams <- parseQueryString(isolate(session$clientData$url_search))
  if ('tab' %in% names(queryParams))
    updateTabsetPanel(session, 'tabset', selected = queryParams[['tab']])

  # source different tabs
  tabFiles <- list.files(path="tabs/")
  for(tab in tabFiles) {
    path <- file.path("tabs", tab)
    source(path, local=TRUE)
  }

})



