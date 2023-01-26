options(warn=-1)
options(shiny.maxRequestSize = 3000*1024^2)
options(stringsAsFactors = F)
#options(shiny.error = recover) 
library(OligoDistiller)
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
  
  output$ParamAnnotation1a <- renderUI({
    if (input$OptionAnnotation == "C1") {
      textInput("mCPD1", "FLP name", value = "Demo A")
    }
  })
  
  output$ParamAnnotation1b <- renderUI({
    if (input$OptionAnnotation == "C1") {
      textInput("mFormula1", "FLP chemical formula:", value = "C189H238O119N66P18S4F8")
    }
  })
  
  output$ParamAnnotation1c<- renderUI({
    if (input$OptionAnnotation == "C1") {
      downloadLink('downloadTrans1', 'Download template for modification list')
    }
  })
  
  output$ParamAnnotation1d<- renderUI({
    if (input$OptionAnnotation == "C1") {
      fileInput("mTrans1", "Upload edited modification list:", multiple = FALSE)
    }
  })
  
  output$ParamAnnotation1e <- renderUI({
    if (input$OptionAnnotation == "C1"){
      selectInput("FLP_type1", "FLP type:", choices = c("DNA", "RNA"))
    } 
  })
  
  output$ParamAnnotation1f <- renderUI({
    if (input$OptionAnnotation == "C1"){
      numericInput("mSigma1", "Maximum isotopic pattern fit deviation (mSigma):", 5, min = 0.5, max = 10)
    }
  })
  
  output$downloadTrans1 <- downloadHandler(
    filename = function() {
      "Transformation_list_jennifer_shortened.txt"
    },
    content = function(file) {
      data = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/Transformation_list_jennifer_shortened.txt", sep = "\t", header=T, stringsAsFactors = F)
      write.table(data, file, quote = F, col.names = T, row.names = F, sep = "\t")
    }
  )
  
  output$ParamAnnotation2a<- renderUI({
    if (input$OptionAnnotation == "C2") {
      downloadLink('downloadDB2', 'Download template for database')
    }
  })
  
  output$ParamAnnotation2b <- renderUI({
    if (input$OptionAnnotation == "C2") {
      fileInput("mDB2", "Upload edited database:", multiple = FALSE)
    }
  })
  
  output$ParamAnnotation2c <- renderUI({
    if (input$OptionAnnotation == "C2"){
      selectInput("FLP_type2", "FLP type:", choices = c("DNA", "RNA"))
    } 
  })
  
  output$ParamAnnotation2d <- renderUI({
    if (input$OptionAnnotation == "C2"){
      numericInput("mSigma2", "Maximum isotopic pattern fit deviation (mSigma):", 5, min = 0.5, max = 20)
    }
  })
  
  output$downloadDB2 <- downloadHandler(
    filename = function() {
      "ref_target_demo_AB.txt"
    },
    content = function(file) {
      data = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/ref_target_demo_AB.txt", sep = "\t", header=T, stringsAsFactors = F)
      write.table(data, file, quote = F, col.names = T, row.names = F, sep = "\t")
    }
  )
  
  output$ParamAnnotation3a <- renderUI({
    if (input$OptionAnnotation == "T") {
      textInput("mCPD3", "FLP name", value = "Demo A")
    }
  })
  
  output$ParamAnnotation3b <- renderUI({
    if (input$OptionAnnotation == "T") {
      textInput("mFormula3", "FLP chemical formula:", value = "C189H238O119N66P18S4F8")
    }
  })
  
  output$ParamAnnotation3c<- renderUI({
    if (input$OptionAnnotation == "T") {
      downloadLink('downloadTrans3', 'Download template for modification list')
    }
  })
  
  output$ParamAnnotation3d<- renderUI({
    if (input$OptionAnnotation == "T") {
      fileInput("mTrans3", "Upload edited modification list:", multiple = FALSE)
    }
  })
  
  output$ParamAnnotation3e<- renderUI({
    if (input$OptionAnnotation == "T") {
      numericInput("mSigma3", "Maximum isotopic pattern fit deviation (mSigma):", 3, min = 0.5, max = 10)
    }
  })
  
  output$downloadTrans3 <- downloadHandler(
    filename = function() {
      "Transformation_list_jennifer_shortened.txt"
    },
    content = function(file) {
      data = read.csv("Transformation_list_jennifer_shortened.txt", sep = "\t", header=T, stringsAsFactors = F)
      write.table(data, file, quote = F, col.names = T, row.names = F, sep = "\t")
    }
  )
  
  output$ParamAnnotation4a <- renderUI({
    if (input$OptionAnnotation == "U"){
      selectInput("FLP_type4", "FLP type:", choices = c("DNA", "RNA"))
    } 
  })
  
  output$ParamAnnotation4b <- renderUI({
    if (input$OptionAnnotation == "U"){
      numericInput("mSigma4", "Maximum isotopic pattern fit deviation (mSigma):", 5, min = 0.5, max = 10)
    }
  })
  
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
      
      if (length(masslist)<5){
        mms = "The input spectrum must contain at least 3 peaks!"
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
      ntheo = input$ntheo
      max_mmw_ppm = input$ppm_error # More tolerance on monoisotopic mass
      
      if (input$MSMS){
        min_overlap = 0.4 # Less strict filter for MS2 data
      } else {
        min_overlap = 0.6
      }
      
      output_process = process_scan(scan, polarity = polarity, MSMS = msms, baseline = baseline, 
                                    min_charge = charge_range[1], max_charge = charge_range[2], 
                                    min_mz = mz_range[1], max_mz = mz_range[2], min_mw = mw_range[1], max_mw = mw_range[2] , 
                                    mz_error = 0.01, mw_gap = 1.1, mw_window = ntheo + 1)

      scan.deconvoluted = output_process$scan_processed_aggregated
      scan.processed = output_process$scan_processed
      
      scan.deconvoluted.annotated = NULL
      
      if (input$OptionAnnotation == "T"){
        
        mTrans = input$mTrans3
        mFormula = input$mFormula3
        mCPD = input$mCPD3
        mSigma = input$mSigma3
        baseline = input$baseline
        if (!is.null(mTrans) & !is.null(mFormula)){
          transformation_list = read.csv(mTrans$datapath, sep = "\t", header = T)
          scan.deconvoluted.annotated = annotate_scan_targeted(
            scan.deconvoluted, formula_flp = mFormula, cpd_flp = mCPD, transformation_list, mdb = NULL, ntheo = ntheo,
            min_overlap = min_overlap, max_msigma = mSigma, max_mmw_ppm = max_mmw_ppm, baseline = baseline)
        }}
      
      if (input$OptionAnnotation == "U"){
        
        bblock = input$FLP_type4
        mSigma = input$mSigma4
        baseline = input$baseline
        
        if (!is.null(bblock) & !is.null(mSigma)){
          scan.deconvoluted.annotated = annotate_scan_untargeted(scan.deconvoluted, bblock, ntheo, min_overlap, mSigma, baseline)
        }
      }
      
      if (input$OptionAnnotation == "C1"){
        
        mTrans = input$mTrans1
        mCPD = input$mCPD1
        mFormula = input$mFormula1
        bblock = input$FLP_type1
        mSigma = input$mSigma1
        baseline = input$baseline
        
        if (!is.null(mTrans) & !is.null(mCPD) & !is.null(mFormula) & !is.null(bblock) & !is.null(mSigma)){
          
          transformation_list = read.csv(mTrans$datapath, sep = "\t", header = T)
          scan.deconvoluted.annotated =  annotate_scan_mix(scan.deconvoluted, ntheo = ntheo, 
                                                           formula_flp = mFormula, cpd_flp = mCPD,
                                                           transformation_list = transformation_list, mdb = NULL, bblock = bblock, 
                                                           min_overlap = min_overlap, max_msigma = mSigma, max_mmw_ppm = max_mmw_ppm, baseline = baseline)
        }
      }
      
      if (input$OptionAnnotation == "C2"){
        
        mDB = input$mDB2
        bblock = input$FLP_type2
        mSigma = input$mSigma2
        baseline = input$baseline
        
        if (!is.null(mDB) & !is.null(bblock) & !is.null(mSigma)){
          
          mDB = read.csv(mDB$datapath, sep = "\t", header = T)
          scan.deconvoluted.annotated =  annotate_scan_mix(scan.deconvoluted, ntheo = ntheo, 
                                                           formula_flp = "", cpd_flp = "",
                                                           transformation_list = NULL, mdb = mDB, bblock = bblock, 
                                                           min_overlap = min_overlap, max_msigma = mSigma, max_mmw_ppm=  max_mmw_ppm, baseline = baseline)
        }
      }
      
      if (!is.null(scan.deconvoluted)){
        mms = paste0("Analysis finished! Please check next tabpanel(s)!")
      } else {mms = "No oligonucleotide feature detected!"}
      
      list(results0 = scan.processed, results1 = scan.deconvoluted, results2 = scan.deconvoluted.annotated, message = mms)
    }
  })
  
  observeEvent(input$goButton,{
    
    withProgress({
      setProgress(message="Deconvoluting and annotating data...")
      Sys.sleep(1)
      setProgress(message=process()$message)
    })
  })
  
  output$blank_message1<-renderText({process()$message})
  
  output$table1 <- renderDataTable({
    
    results_table1 = NULL
    results1 = process()$results1
    
    if (!is.null(results1)){
      results_table1 = datatable(results1, escape=rep(FALSE, ncol(results1)), rownames = F)
    }
    return(results_table1)
  })
  
  table1_proxy <- dataTableProxy("table1")
  
  observeEvent(input$clearRowSelection, {
    table1_proxy %>% selectRows(NULL)
  })
  
  output$downloadScan <- downloadHandler(
    filename = function() {"Scan_Deconvoluted.csv"},
    content = function(file) {
      results1 = process()$results1
      if (!is.null(results1)){
        write.csv(results1, file, sep=",",col.names = T, row.names= F, quote = FALSE)
      }
    })
  
  output$downloadCharged <- downloadHandler(
    filename = function() {"Original_scan_with_charge.csv"},
    content = function(file) {
      results1 = process()$results0
      if (!is.null(results0)){
        write.csv(results0, file, sep=",",col.names = T, row.names= F, quote = FALSE)
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
  
  output$DisplayCharged<- renderPlotly({
    
    charged_sp = process()$results0
    
    if (!is.null(charged_sp)){
      x = charged_sp[,1]
      y = charged_sp[,2]
      y1 = y*1.05
      min_x  = min(x)*0.9
      max_x  = max(x)*1.1
      max_y = max(y)*1.2
      z = charged_sp[,3]
      z[z==0] = ""
      #t <- list(family = "sans serif", size = 14, color = toRGB("grey50"))
      
      plot_ly() %>%
        add_segments(x=~x, xend = ~ x, y = ~0, yend=~y, type="scatter", mode="lines", name="") %>% 
        add_text(x = ~x, y = ~y1, text = z, textposition = "top right", showlegend = FALSE) %>% 
        layout(xaxis = list(zeroline = FALSE, title = "m/z", range = c(min_x, max_x)), 
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
      ntheo <- input$mw_gap 
      if (input$OptionAnnotation == "T") {mode <- "targeted"}
      if (input$OptionAnnotation == "U") {mode <- "untargeted"}
      if (input$OptionAnnotation == "C1") {mode <- "mixed"}
      if (input$OptionAnnotation == "C2") {mode <- "mixed"}
      
      yyy <- reconstruct_scan_annotated(test.scan, results, polarity = polarity, 
                                        mode = mode, mz_error = input$mz_error, ntheo = ntheo)
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
