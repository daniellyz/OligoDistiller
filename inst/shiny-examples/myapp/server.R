options(warn=-1)
options(shiny.maxRequestSize = 3000*1024^2)
options(stringsAsFactors = F)

library(shiny)
library(V8)
library(shinyjs)
library(DT, quietly = TRUE)
library(plotly)

shinyServer(function(input, output,clientData, session){

  observeEvent(input$killButton,{
    refresh()})
  
  global <- reactiveValues(OptionAnnotation = NULL)
  
  observeEvent(input$OptionAnnotation, {
    global$OptionAnnotation <- input$OptionAnnotation
  })
  
  output$ParamAnnotation0a <- renderUI({
    if (input$OptionAnnotation == "TDB") {
      downloadLink('downloadDB', 'Download template for database')
    }
  })
  
  output$ParamAnnotation0b <- renderUI({
    if (input$OptionAnnotation == "TDB") {
      fileInput("mDB", "Upload edited database:", multiple = FALSE)
    }
  })
  
  output$ParamAnnotation1a <- renderUI({
    if (input$OptionAnnotation == "T") {
      textInput("mFormula", "FLP chemical formula:", value = "C189H238O119N66P18S4F8")
    }
  })
  
  output$ParamAnnotation1b <- renderUI({
    if (input$OptionAnnotation == "T") {
      downloadLink('downloadTrans', 'Download template for modification list')
    }
  })
  
  output$ParamAnnotation1c <- renderUI({
    if (input$OptionAnnotation == "T") {
      fileInput("mTrans", "Upload edited modification list:", multiple = FALSE)
    }
  })
  
  output$ParamAnnotation2 <- renderUI({
    if (input$OptionAnnotation == "UB") {
      textInput("mBB1", "Building block:", value = "C21 H26.4 O13.2 N7.3 P2 S0.4 F0.9")
    }
  })
  
  output$ParamAnnotation2m <- renderUI({
    if (input$OptionAnnotation == "UB"){
      numericInput("mSigma2", "Maximum isotopic pattern fit deviation (mSigma):",5, min = 0.5, max = 50)
    }
  })
  
  output$ParamAnnotation3 <- renderUI({
    if (input$OptionAnnotation == "UP"){
      selectInput("FLP_type", "FLP type:", choices = c("DNA", "RNA"))
    } 
  })
  
  output$ParamAnnotation3m <- renderUI({
    if (input$OptionAnnotation == "UP"){
      numericInput("mSigma3", "Maximum isotopic pattern fit deviation (mSigma):",5, min = 0.5, max = 50)
    }
  })
  
  output$downloadTrans <- downloadHandler(
    filename = function() {
      "Transformation_list.txt"
    },
    content = function(file) {
      data = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/Transformation_list.txt", sep = "\t", header=T, stringsAsFactors = F)
      write.table(data, file, quote = F, col.names = T, row.names = F, sep = "\t")
    }
  )
  
  output$downloadDB <- downloadHandler(
    filename = function() {
      "ref_target_demo_AB.txt"
    },
    content = function(file) {
      data = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/ref_target_demo_AB.txt", sep = "\t", header=T, stringsAsFactors = F)
      write.table(data, file, quote = F, col.names = T, row.names = F, sep = "\t")
    }
  )
  
  read_MS <-reactive({
    
    query_spectrum = NULL
    valid = 1
    mms = ""

    if (input$blank_file1==""){
      mms="Please paste your mass peaks!"
      valid = 0
    }
    
    if (input$blank_file1!=""){
      
      input_str = input$blank_file1
      input_str = strsplit(input_str,"\n")[[1]]
      
      input_str = lapply(input_str, function(x) strsplit(x, "\\s+")[[1]])
      input_str = lapply(input_str, function(x) x[x!="\t"])
      input_str = lapply(input_str, function(x) x[x!=""])
      
      if (all(sapply(input_str,length)==1)){ # One column situation
        masslist=as.numeric(unlist(input_str))
        masslist=masslist[!is.na(masslist)]
        intlist = rep(100, length(masslist))}
      
      if (any(sapply(input_str,length)>1)){ # >1 column situation
        masslist = as.numeric(sapply(input_str,function(x) x[1]))
        intlist = as.numeric(sapply(input_str,function(x) x[2]))
        valid_peaks = which(!is.na(masslist) & !is.na(intlist))
        masslist = masslist[valid_peaks, drop = FALSE]
        intlist = intlist[valid_peaks, drop = FALSE]
      }
      
      if (length(masslist)<3){
        mms = "The input spectrum must contain at least 3 fragment!"
        valid = 0
      }
      
      if (!is.null(masslist) && valid==1){ 
        query_spectrum = cbind(masslist, intlist)
        if (!input$centroid){
          query_spectrum = centroid_scan(query_spectrum)
        }
      }
    }
    
    list(query_spectrum = query_spectrum, mms = mms,valid=valid)
  })
  
  output$SpectrumShow<- renderPlotly({
    
    query_spectrum = read_MS()$query_spectrum 
    
    if (!is.null(query_spectrum)){
      x = query_spectrum[,1]
      y = query_spectrum[,2]
      
      min_x  = min(x)*0.9
      max_x  = max(x)*1.1
      max_y = max(y)*1.2
      
      plot_ly() %>%
        add_segments(x=~x, xend = ~ x, y = ~0, yend=~y, type="scatter", mode="lines", name="") %>% 
        layout(xaxis = list(zeroline = FALSE, title = "m/z", range = c(min_x, max_x)), 
               yaxis = list(title = "Intensity", range = c(0, max_y))) 
    }
  })
  
  observeEvent(input$exampleButton1, {
    fileText <- paste(readLines("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/example_demo_a.txt"), collapse = "\n")
    updateTextAreaInput(session, "blank_file1", value = fileText)
  })  
  
  process <- eventReactive(input$goButton,{

    results = NULL
    mms = ""
    scan = read_MS()$query_spectrum 
    
    if (!is.null(scan)){

      is_negative = input$polarity
      msms = input$MSMS
      
      if (is_negative){
        polarity = "Negative"
      } else {polarity = "Positive"}
      
      charge_range = input$charge_range
      mz_range = input$mz_range
      mw_range = input$mw_range
      baseline = input$baseline
      mz_error = input$mz_error
      mw_gap = input$mw_gap
      ntheo = mw_gap + 2 # Number of theoretical fragment set a bit higher than mw gap
      
      scan.deconvoluted = process_scan(scan, polarity = polarity, baseline = baseline, 
                MSMS = msms, min_charge = charge_range[1], max_charge = charge_range[2], 
                min_mz = mz_range[1], max_mz = mz_range[2], min_mw = mw_range[1], max_mw = mw_range[2] , 
                mz_error = mz_error, mw_gap = mw_gap)
      scan.deconvoluted.annotated = NULL
      
      if (input$OptionAnnotation == "TDB"){
        
        mDB = input$mDB
        if (!is.null(mDB)){
          DB = read.csv(mDB$datapath, sep = "\t", header = T)
          scan.deconvoluted.annotated = annotate_scan_targeted(
            scan.deconvoluted, formula_flp = "", transformation_list = NULL, mdb = DB, ntheo = ntheo)
        }}
      
      if (input$OptionAnnotation == "T"){
        
        mTrans = input$mTrans
        mFormula = input$mFormula
        if (!is.null(mTrans) & !is.null(mFormula)){
          transformation_list = read.csv(mTrans$datapath, sep = "\t", header = T)
          scan.deconvoluted.annotated = annotate_scan_targeted(
            scan.deconvoluted, formula_flp = mFormula, transformation_list, mdb = NULL, ntheo = ntheo)
        }}
      
      if (input$OptionAnnotation == "UB"){
        
        bblock = input$mBB1
        mSigma = input$mSigma2
        if (input$MSMS){
          min_overlap = 0.4
          min_relative = 0.01 # Less strict filter for MS2 data
        } else {
          min_overlap = 0.6
          min_relative = 0.1 
        }
        
        if (!is.null(bblock) & !is.null(mSigma)){
          scan.deconvoluted.annotated = annotate_scan_untargeted(scan.deconvoluted, bblock, ntheo, min_relative, min_overlap, mSigma)
       }
      }
      
      if (input$OptionAnnotation == "UP"){
        
        bblock = input$FLP_type
        mSigma = input$mSigma3
        if (input$MSMS){
          min_overlap = 0.4
          min_relative = 0.01 # Less strict filter for MS2 data
        } else {
          min_overlap = 0.6
          min_relative = 0.1 
        }
        
        if (!is.null(bblock) & !is.null(mSigma)){
          scan.deconvoluted.annotated = annotate_scan_untargeted(scan.deconvoluted, bblock, ntheo, min_relative, min_overlap, mSigma)
        }
      }
    }
    
    if (!is.null(scan.deconvoluted)){
      mms = paste0("Analysis finished! Please check next tabpanel(s)!")
    } else {mms = "No oligonucleotide feature detected!"}

    list(results1 = scan.deconvoluted, results2 = scan.deconvoluted.annotated, message = mms)
  })

  observeEvent(input$goButton,{

    withProgress({
      if (input$OptionAnnotation == "N"){
        setProgress(message="Deconvoluting data...")
      } else {
        setProgress(message="Deconvoluting and annotating data...")
      }
      Sys.sleep(1)
      setProgress(message=process()$message)
    })
  })

  output$blank_message1<-renderText({process()$message})

  output$table1 <- renderDataTable({
    
    results_table1 = NULL
    results1 = process()$results1
    
    if (!is.null(results1)){
      results_table1 = datatable(results1,escape=rep(FALSE, ncol(results1)), rownames = F)
    }
    return(results_table1)
  })
  
  output$downloadScan <- downloadHandler(
    filename = function() {"Deconvoluted_scan.csv"},
    content = function(file) {
      results1 = process()$results1
      if (!is.null(results1)){
        write.csv(results1, file, sep=",",col.names = T, row.names= F, quote = FALSE)
      }
  })
  
  output$DisplayDeconvoluted<- renderPlotly({
    
    results1 = process()$results1
    deconvoluted_sp = results1[,c("MW", "Response", "Envelop")]
    
    if (!is.null(deconvoluted_sp)){
      x = deconvoluted_sp[,1]
      y = deconvoluted_sp[,2]
      
      if (length(input$table1_rows_selected) > 0) {
        
        inds <- input$table1_rows_selected
        
        # zoom into the complete envelopes that correspond to the selected rows from the DT table
        deconvoluted_sp_tmp <- deconvoluted_sp[deconvoluted_sp$Envelop %in% unique(deconvoluted_sp$Envelop[inds]), -3]
        
        min_x = min(deconvoluted_sp_tmp[, 1])   
        max_x = max(deconvoluted_sp_tmp[, 1])
        small_broadening <- c(max_x-min_x) * 0.05
        min_x <- min_x - small_broadening
        max_x <- max_x + small_broadening
        
        max_y = max(deconvoluted_sp_tmp[, 2])
        
      } else {
        
        min_x = min(x)   
        max_x = max(x)
        max_y <- max(y)
        
      }
      
      plot_ly() %>%
        add_segments(x=~x, xend = ~ x, y = ~0, yend=~y, type="scatter", mode="lines", name="") %>% 
        layout(xaxis = list(zeroline = FALSE, title = "Molecular weight (Da)", range = c(min_x, max_x)), 
               yaxis = list(title = "Intensity", range = c(0, max_y))) 
    }
  })
  
  output$table2 <- renderDataTable({

    results_table2 = NULL
    results2 = process()$results2$feature

    if (!is.null(results2)){
        results_table2 = datatable(results2,escape=rep(FALSE, ncol(results2)), rownames = F, selection = 'single')
    }
    return(results_table2)
  })

  output$downloadFeature <- downloadHandler(
    filename = function() {"Molecular_feature.csv"},
    content = function(file) {
      results2 = process()$results2$feature
      if (!is.null(results2)){
        write.csv(results2, file, sep=",",col.names = T, row.names= F, quote = FALSE)
      }
    }
  )
  
  output$DisplayRawAndReconstructed <- renderPlotly({
    results <- process()$results2
    if (!is.null(results)) {
      test.scan <- read_MS()$query_spectrum
      polarity <- "Negative"
      is_negative <- input$polarity
      if (!is_negative) {
        polarity <- "Positive"
      }
      ntheo <- input$mw_gap + 2
      if (input$OptionAnnotation == "T") {mode <- "targeted"}
      if (input$OptionAnnotation == "TDB") {mode <- "targeted"}
      if (input$OptionAnnotation == "UB") {mode <- "untargeted"}
      if (input$OptionAnnotation == "UP") {mode <- "untargeted"}
      
      yyy <- reconstruct_scan_annotated(test.scan, results,
                                        polarity = polarity, baseline = input$baseline,
                                        mode = mode, mz_error = input$mz_error, ntheo = ntheo
      )
      original_sp <- yyy$original_scan
      reconstructed_sp <- yyy$reconstructed_scan
      idx <- input$table2_rows_selected
      
      p_raw <- p_reconstructed <- NULL
      
      if (!is.null(original_sp)) {
        x_orig <- original_sp[, 1]
        
        y_orig <- original_sp[, 2]
        min_x <- input$mz_range[1] * 0.9
        max_x <- input$mz_range[2] * 1.1
        max_y <- max(y_orig) * 1.2

        p_raw <- plot_ly(source = "p_raw") %>%
          add_segments(x = ~x_orig, xend = ~x_orig, y = ~0, yend = ~y_orig, type = "scatter", mode = "lines", name = "Original spectrum", line = list(color = c("darkblue"))) %>%
          layout(
            xaxis = list(zeroline = FALSE, title = "m/z", range = c(min_x, max_x)),
            yaxis = list(title = "Intensity", range = c(0, max_y))
          )
        
        # add charge state text annotation for the raw plot
        
        if (length(idx) == 1) {
          VT <- which(original_sp$labels == results$feature$FEATURE[idx])
          if (length(VT) > 0) {
            charge_anno <- list(
              x = x_orig[VT],
              y = y_orig[VT],
              text = paste0("<b>", original_sp[VT, 4], "</b>"),
              xref = "x",
              yref = "y",
              showarrow = TRUE,
              startarrowsize = 0.7,
              font = list(size = 15)
            )
            p_raw <- p_raw %>% layout(annotations = charge_anno)
          }
        }
      }
 
      if (!is.null(reconstructed_sp)) {
        x_recon <- reconstructed_sp[, 1]
        
        y_recon <- reconstructed_sp[, 2]
        
        min_x <- input$mz_range[1] * 0.9
        max_x <- input$mz_range[2] * 1.1
        max_y <- max(y_recon) * 1.2
        
        p_reconstructed <- plot_ly(source = "p_reconstructed") %>%
          add_segments(x = ~x_recon, xend = ~x_recon, y = ~0, yend = ~y_recon, type = "scatter", mode = "lines", name = "Reconstructed spectrum", line = list(color = c("brown"))) %>%
          layout(
            xaxis = list(zeroline = FALSE, title = "m/z", range = c(min_x, max_x)),
            yaxis = list(title = "Intensity", range = c(0, max_y))
          )
        
        # add charge state text annotation for the reconstructed plot
        if (length(idx) == 1) {
          VT <- which(reconstructed_sp$labels == results$feature$FEATURE[idx])
          if (length(VT) > 0) {
            charge_anno <- list(
              x = x_recon[VT],
              y = y_recon[VT],
              text = paste0("<b>", reconstructed_sp[VT, 4], "</b>"),
              xref = "x",
              yref = "y",
              showarrow = TRUE,
              startarrowsize = 0.7,
              font = list(size = 15)
            )
            p_reconstructed <- p_reconstructed %>% layout(annotations = charge_anno)
          }
        }
      }
      
      if (!is.null(p_raw) & !is.null(p_reconstructed)){
        subplot(      
          p_raw,
          p_reconstructed,
          shareX = TRUE, 
          titleY = TRUE,
          nrows = 2
          )
      } else if (!is.null(p_raw)) {
        p_raw %>% layout(title = "Raw spectrum")
      } else p_reconstructed %>% layout(title = "Reconstructed spectrum")
    }
  })
})
