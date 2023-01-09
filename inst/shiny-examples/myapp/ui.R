options(warn=-1)
options(shiny.maxRequestSize = 3000*1024^2)
options(stringsAsFactors = F)
#options(shiny.error = recover) 

library(shiny)
library(V8)
library(shinyjs)
library(DT, quietly = TRUE)
library(plotly)

shinyUI(navbarPage("OligoDistiller V.01",

        tabPanel("A) Start a Run",
            shinyjs::useShinyjs(),
           # shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),

            column(6,

                h4("Please paste your MS spectrum into the field below:"), 
                textAreaInput("blank_file1", label = '',width=500,height=200),
                   
                checkboxInput("polarity", h5("Data aquired is in negative ion mode."), TRUE),
                checkboxInput("centroid", h5("MS data is already centroided."), TRUE),
                checkboxInput("MSMS", h5("Data aquired is MS/MS data or very noisy."), FALSE),
                
                br(),
                h4("Please estimate the charge state range:"),
                sliderInput("charge_range", "", 1, 30, value = c(3, 12)),

                br(),
                h4("Please define the m/z range in the input spectra:"),
                sliderInput("mz_range", "", 100, 5000, value = c(500, 2000), step = 100),

                br(),
                h4("Please define the molecular weight range in processed spectra:"),
                sliderInput("mw_range", "", 0, 20000, value = c(4000, 12000), step = 1000),

                br(),
                plotlyOutput("SpectrumShow",height = 400, width = 600)
            ),

            column(6,

                h3("MS Processing parameters:"),
                numericInput("baseline", h4("MS spectral baseline"), 100, 0, 10000, 100),
                numericInput("mz_error", h4("MS measurement error (Da)"), 0.01, 0, 0.1, 0.005),
                numericInput("ntheo", h4("Estimated isotope envelop size (Da)"), 10, 0, 20, 1),
        
                br(),
                h3("Annotation parameters:"),
                radioButtons(
                  "OptionAnnotation",
                  "",
                  c("Combined with a modification list" = "C1", 
                    "Combined with a database" = "C2",
                    "Targeted with a modification list" = "T",
                    "Untargeted with Pointless algorithm" = "U")
                ),
                
                uiOutput('ParamAnnotation1a'),
                uiOutput('ParamAnnotation1b'),
                uiOutput('ParamAnnotation1c'),
                uiOutput('ParamAnnotation1d'),
                uiOutput('ParamAnnotation1e'),
                uiOutput('ParamAnnotation1f'),
                uiOutput('ParamAnnotation2a'),
                uiOutput('ParamAnnotation2b'),
                uiOutput('ParamAnnotation2c'),
                uiOutput('ParamAnnotation2d'),
                uiOutput('ParamAnnotation3a'),
                uiOutput('ParamAnnotation3b'),
                uiOutput('ParamAnnotation3c'),
                uiOutput('ParamAnnotation3d'),
                uiOutput('ParamAnnotation3e'),
                uiOutput('ParamAnnotation4a'),
                uiOutput('ParamAnnotation4b'),
                
                tags$head(
                  tags$style(HTML('#goButton{background-color:lightgreen}'))
                ),
                br(),
                actionButton("goButton", "Process",style='padding:6px; font-size:150%'),

                br(),
                tags$head(
                  tags$style(HTML('#killButton{background-color:orange}'))
                ),
                br(),
                actionButton("killButton", "Clear",style='padding:6px; font-size:150%'),
                br(),
                tags$head(
                  tags$style(HTML('#exampleButton1{background-color:lightblue}'))
                ),
                br(),
                actionButton("exampleButton1", "Load Example Demo Strand A",style='padding:6px; font-size:150%'),
                br(),
                br(),
                br(),
                em('Messages from the server:'),
                br(),
                h3(textOutput("blank_message1"))
            )
        ),

        tabPanel("B) Spectral deconvolution results",

                br(),
                dataTableOutput("table1"),
                br(),
                actionButton("clearRowSelection", "Clear row selection"),
                br(),
                br(),
                downloadButton("downloadScan", "Download Deconvoluted Scan"),
                br(),
                br(),
                downloadButton("downloadCharged", "Download peaks with charge"),
                br(),
                br(),
                plotlyOutput("DisplayDeconvoluted", height = 400, width = 800),
                br(),
                br(),
                plotlyOutput("DisplayCharged", height = 400, width = 800)
        ),
        
        tabPanel("C) Mass feature and annotations",
                 
                 br(),
                 dataTableOutput("table2"),
                 br(),
                 downloadButton("downloadFeature", "Download"),
                 br(),
                 br(),
                 #plotlyOutput("DisplayRaw", height = 400, width = 600),
                 #br(),
                 #plotlyOutput("DisplayReconstructed", height = 400, width = 600),
                 #br()
                 plotlyOutput("DisplayRawAndReconstructed", height = 800, width = 800),
                 br()
        )
))

