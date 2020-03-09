library(shiny)
library(dplyr)
library(shinythemes)
library(DT)
library(shinydashboard)
library(visNetwork)
library(readxl)
library(igraph)
library(ggplot2)
#library(shinyjs)
#library(bcrypt)
#library(plotly)
library(shinyWidgets)

source(file.path("sgtc_functions.R"))
source(file.path("functions.R"))
source(file.path("2019_nov20_aug3_runGO.R"))


fdat = read.delim("fdata.txt",stringsAsFactors = F,check.names = F)
noness = fdat  %>% filter(essential == "noness")
ess = fdat  %>% filter(essential == "ess")

phiphop <- readRDS("jan1_phiphop.RDS")
 
xhiphop <- readRDS("jan1_xhiphop.RDS")

xcoinh <- readRDS("jan5_coinhibition_6081_6081.RDS")



xsig <-
  as.matrix(read.delim(
    "2019dec31_signatures_2.76_score_matrix.txt",
    stringsAsFactors = F,
    header = T,
    check.names = F
  ))



xhiphop = xhiphop[,phiphop$name]

we = which(rownames(xhiphop)%in% ess$sgd_gene)
wn = which(rownames(xhiphop)%in% noness$sgd_gene)


xe = xhiphop[we,phiphop$name]
xne = xhiphop[wn,phiphop$name]

pgold <- phiphop %>% filter(gold == "gold")
 
pGO <-
  read.delim(
    "UBC_Novartis_GOenrich.txt",
    stringsAsFactors = F,
    header = T,
    check.names = F
  )

df = read.delim("df_signatures.txt",check.names = F,stringsAsFactors = F)

hits = read.delim("dec25_hiphits.txt",check.names = F,stringsAsFactors = F)

header <- dashboardHeader(title =  "chemogenomics UBC NIBR",titleWidth = 450)

sidebar <- dashboardSidebar(
  sidebarMenu(
   
    menuItem("select fitness profile", tabName="prof", icon = icon("chart-line")),
    menuItem("gold standards", tabName="welc", icon = icon("vial")),
    
    #menuItem("hip vs hop profiles", tabName="fit", icon = icon("chart-line")),
    
    menuItem("hiphop profiles", tabName="hiphop", icon = icon("bullseye")),
    menuItem("GO enrichment", tabName="net", icon = icon("dna")),
    #menuItem("GO enrichments", tabName="goenrich",icon = icon("dna"))
    menuItem("Signatures", tabName="sig", icon = icon("fish"))
  )
)

body <- dashboardBody(
  theme = shinytheme("cerulean"),
  tabItems(
    tabItem(tabName = "welc",
      
      box(title = "log2 ratio by condition",
          plotOutput("coolplot"),status = "primary", solidHeader = T,width = 11),
      box(title = "density z-scores",plotOutput("densplot"),status = "primary", solidHeader = T,width = 11)
      #br(),br()),
    ),#tab
    tabItem(tabName = "prof",
     
      fluidRow(
        
        box(title = "enter gene of interest", 
            
            textInput(
              inputId = "gene",
              label = "",
              value   = "TOR2"
            ), status = "primary", solidHeader = T, width = 3, height = 150),
        
        box(title = "select pool",
            radioButtons("pool1","zygosity",
                     selected = "HIP",choices = c("HIP","HOP"),inline = T),
            status = "primary", solidHeader = T, width = 3, height = 150)
            
             ),
        
    #fluidRow(
        box(title = "Intensity across conditions",
            h4("click and drag to select compounds of interest; view corresponding 
                    fitness plots in hiphop tab"),
            plotOutput("plot4", brush = "brush4", height = 800),status = "primary", solidHeader = T 
            ,width = "80%"),
        
        
        #column(width = 12,
        
        box(title = "experimental detail",dataTableOutput("info"),status = "primary", solidHeader = T 
             ,width = "80%"
        )),#coln fluid tabitem
    
    tabItem("hiphop",
      #p("HIPHOP profiles, GO enrichments & concordance with UBC chemogenomic profiles"),
      fluidRow(
        #column(width = 6,
          box(title = "select compound",status = "primary", solidHeader = T,
            selectizeInput(
              'cmp','',
              choices = phiphop$name,selected = phiphop$name[183],
              multiple = F,
              options = list(
                placeholder = 'Please type to search by compound or select an option below',
                onInitialize = I('function() { this.setValue(""); }')
              )#list
            )#selectize
          )###box
        ),#coln row
     
      fluidRow(
        box(title = "HIP gene target:",status = "primary", solidHeader = T, width = 6,
          
          dataTableOutput('targethip')
          ),#box
        box(title = "PCID gold standards and UBC only:",status = "primary", solidHeader = T, width = 6,
            
            dataTableOutput('pcid')
        )#box
        
      ),#fluid row
      
      #column(width = 10,
      fluidRow(
        
        box(title = "HIP fitness defect scores",status = "primary", 
            solidHeader = T,plotOutput("efd", width = "100%",height = 600)),
        box(title = "HOP fitness defect scores",status = "primary", 
            solidHeader = T,plotOutput("nfd", width = "100%",height = 600)),
        
        column(width = 12,
        box(title = "GO enrichments and concordance with UBC dataset",status = "primary", 
            solidHeader = T, width = "100%",
               dataTableOutput("GO"))
        )#coln
      ),#fluid
    fluidRow(
        
        box(title="HIP genes",
        dataTableOutput('hiptab'),status = "primary", solidHeader = T,width = 6),
        
        box(title="HOP genes",
            dataTableOutput('hoptab'),status = "primary", solidHeader = T,width = 6),
        
        box(title="download HIP",
            downloadButton('downloadhip', 'HIPdownload'),status = "primary", solidHeader = T,width = 6), 
        box(title="download HOP",
            downloadButton('downloadhop', 'HOPdownload'),status = "primary", solidHeader = T,width = 6)
            ),#fluid row
      
    fluidRow(
        
       box(title="Coinhibitory screens",
          dataTableOutput('coinhib'),status = "primary", solidHeader = T,width = "100%"),
            )#fluid row
      
    ),#tabitem
    tabItem("net",
      fluidRow(
          #column(width = 1, align= "left",
              box(title= "select geneset",
                  prettyRadioButtons("pool", label = "background genome:",
                  c("HIP" = "ess", "HOP" = "non","HIPHOP" = "both"), outline = T, fill = F, 
                  status = "primary", shape = "square",bigger = T,
                  selected = "non", inline = T),solidHeader = T,status = "primary",width= 3,height = 150),
                #),#coln
                
          #column(width  = 4,
                       
                       #   # Input: Select quotes ----
              box(title = "fitness score threshold",
                  sliderInput("scorethresh",label = "input score threshold", min = 0, max = 5,
                        value = 1, step = 0.5),
                  status = "primary", solidHeader = T,width = 3,height = 150),
                #),#coln
                
                
          #column(width = 4,
                       
               box(title = "set FDR threshold",
                        sliderInput("fdr",label = "FDR threshold", min = 0, max = 0.5,
                                       value = 0.1, step = 0.05), 
                   status = "primary", solidHeader = T,width = 3,height = 150)
                #)#coln
                
              ),#row
              
      fluidRow(
         column(width  = 3,
              box(title ="enrichment map details",
                           
                           div(style="font-size:16px",
                               tags$strong("Click node to view leading edge genes"),tags$br(),
                               tags$strong("Click node to view enrichment details"),tags$br(),
                               tags$strong("Right click to download image")),
                  status = "primary", solidHeader = T,width="100%",height = 135)
                ),#coln
                
        column(width = 3,
                       #uiOutput("selector"),
              box(title = "select experiment",
                           selectizeInput(
                             'cmp1','',
                             choices = phiphop$name,selected = phiphop$name[183],multiple = F,
                             options = list(
                               placeholder = 'Please type to search or select an option below',
                               onInitialize = I('function() { this.setValue(""); }')
                             )),status = "primary", solidHeader = T,width="100%", height = 135)
                ),#coln
                
       column(width = 3, align="center",
             box(title ="download enrichment",
                           br(),
                           downloadButton('enrich', 'download'),
                 status = "primary", solidHeader = T,width="100%", height = 135)
                       
                )),#coln row
              
    fluidRow(
                
          box(
                  width = 9, title = "GO enrichment network", status = "primary", solidHeader = TRUE,
                  visNetworkOutput("network_proxy",width = "100%",height = 900)
                ),
                
          box(title = "top contributing genes",width = 3,
                    uiOutput("ui1"),solidHeader = T,status = "primary"),
                
          box(title = "enrichment details",width = 3,
                    uiOutput("ui2"),solidHeader = T,status = "primary",background = "navy")
                
                
              ),#row tab
      ),#tabItem
      
    tabItem("sig",
            fluidRow(
             
              
              column(width  = 4,
                     
                     box(title = "fitness score threshold (minimum = 3)",
                         sliderInput("scorethresh2",label = "input score threshold", min = 3, max = 10,
                                     value = 1, step = 0.5),
                         status = "primary", solidHeader = T,width="100%")
                     ),#coln
              
              
              column(width = 4,
                     
                     box(title = "set FDR threshold",
                         sliderInput("fdr2",label = "FDR cutoff", min = 0, max = 0.5,
                                     value = 0.2, step = 0.05), status = "primary", solidHeader = T,width="100%")
                    )#coln
              
               ),#row
            
           
            
            fluidRow(
             
           
              column(width  = 4,
                     box(title ="enrichment map details",
                         
                         div(style="font-size:16px",
                             tags$strong("Click node to view leading edge genes"),tags$br(),
                             tags$strong("Click node to view enrichment details"),tags$br(),
                             tags$strong("Right click to download image")),status = "primary", solidHeader = T,width="100%")
                   ),#coln
              
              column(width = 4,
                    
                     box(title = "select signature",
                         selectizeInput(
                           'resp','',
                           choices = colnames(xsig),multiple = F,
                           options = list(
                             placeholder = 'Please type to search or select an option below',
                             onInitialize = I('function() { this.setValue(""); }')
                           )),status = "primary", solidHeader = T,width="100%")
                      )
                    ),#coln row
              
            
            fluidRow(
              
              box(
                width = 9, title = "GO enrichment network for response signature", 
                status = "primary", solidHeader = TRUE,
                visNetworkOutput("network_proxy2",width = "100%",height = 800)
                  ),
              
             
              box(title = "enrichment details",width = 3,
                  uiOutput("ui4"),solidHeader = T,status = "primary",background = "navy"),
    
              
              box(title = "top contributing genes",width = 3,
                  uiOutput("ui3"),solidHeader = T,status = "primary"),
              
              
              box(title = "gene response signature", width = 3,
                  uiOutput("ui5"), status = "primary", solidHeader = T
              )
              
              
            ),#row tab
    )#tabItem
    
    
  ))
    

ui <- dashboardPage(header, sidebar, body)
# 02 - Server -------------------------------------------------------------

server <- function(input, output, session) { 
  output$menu <- renderMenu({
    sidebarMenu(
      menuItem("welcome", tabName="welc", icon = icon("vial")),
      
      menuItem("select fitness profile", tabName="prof", icon = icon("chart-line")),
      #menuItem("hip vs hop profiles", tabName="fit", icon = icon("chart-line")),
      menuItem("hiphop profiles", tabName="hiphop", icon = icon("bullseye"))
      
    )
  })
  #####go
  ig <- reactive({ 
    req(input$cmp1)
    req(input$pool)
    if(input$pool == "non") { fdata = xne} else if  (input$pool == "ess") 
    {fdata = xe} else fdata = xhiphop
    
    w = which(colnames(fdata) %in% input$cmp1)
    validate(
      need(length(w) != 0, "Please select a compound")
    )
    req(length(w)!=0)
    df = data.frame(gene = rownames(fdata), score = fdata[,w],index = 0,stringsAsFactors = F)
    w = which(df$score >= input$scorethresh)
    validate(need(length(w)!=0, message = "No scores above threshold"))
    df$index[w]=1
    df = df %>% arrange(desc(score))
    curr_exp = "tst"
    df = df[,c('index','score','gene')]
    w = which(df$index ==1)
   
    
    tst = myrun_go_enrich1(fdrThresh = input$fdr, curr_exp = "tst",score = df, bp_path = "BP_geneNames_current.gmt")
    validate(
      need(!is.null(tst$edgeMat), "No GO enrichment, try another compound or pool")
    )
    #this stop worked 
    #if(is.null(tst$edgeMat))  {stop("No Go enrichment")}
    n = tst$enrichInfo
    e = tst$edgeMat
    w = which(names(n) == "id")
    coln = ncol(n) - 1
    n = n[,c(w,1:coln)]
    
   
    w = which(names(e) == "label")
    let = graph_from_data_frame(e[,-w],vertices = n,directed = F)
    v=set_vertex_attr(let,"label",value=n$formattedLabel)
    v=set_vertex_attr(v,"color.background",value=n$cluster)
    v=set_vertex_attr(v,"label.color",value="black")
    v=set_vertex_attr(v,"label.family",value=1)
    v=set_vertex_attr(v,"label.font",value=3)
    v=set_edge_attr(v,"color",value="black")
    
    vis = toVisNetworkData(v)
    vis$nodes$label = vis$nodes$formattedLabel
    vis$nodes = vis$nodes %>% arrange(label)
    vis$edges$color = "black"
    vis$nodes$color.border="black"
    vis$nodes$font.face ="Courier"
    vis$nodes$shape = "dot"
    vis$nodes$font.size = 27
    vis$nodes$borderWidth=1
    w = which.min(vis$nodes$FDR)
    if(length(w)>0) vis$nodes$font.size[w] = 33
    if(length(w)>0) vis$nodes$borderWidth[w] = 5
    vis
  })
  
  output$network_proxy <- renderVisNetwork({
    validate(
      need(!is.na(ig()$nodes),"lone 298 ig nodes"))
    w =  nrow(ig()$nodes)
    
    validate(
      need(w != 0, "line 302 GO enrichment")
    )
    
    vis = ig()
    g = grep("color",names(vis$nodes))
    names(vis$nodes)[g] = "color.background"
    
    vis$nodes$color.border="black"
    vis$nodes$FDR = signif(vis$nodes$FDR,2)
    visNetwork(vis$nodes, vis$edges) %>% 
      visNodes(shadow=list(enabled=T,size=25),borderWidth=1) %>%
      visOptions(nodesIdSelection = list(enabled = TRUE, 
                                         style = 'width: 150px;color: darkblue;'),
                 selectedBy = list(variable="FDR",
                                   style = 'width: 150px;color: darkblue;')) %>%
      visIgraphLayout()  %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
  })
  
  observe({
    nodes_selection <- input$networkid_selected
    visNetworkProxy("network_proxy_select") %>%
      visSelectNodes(id = nodes_selection)
  })
  
  observeEvent(input$current_node_id, {
    req(input$current_node_id)
    
    vis = ig()
    n = ig()$nodes
    w =  nrow(n)
   
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id)
    n = vis$nodes[w,]
  })
  
  observeEvent(input$id,
               {
                 vis = ig()
                 w = which(vis$nodes$label == input$id)
                 n = vis$nodes[w,]
                 
               })
  
  output$table1 <- DT::renderDataTable({
    req(input$current_node_id)
    vis = ig()
    n = ig()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id)
    term = vis$nodes$formattedLabel[w]
    nam = c("term","nGenes","geneSetFraction","FDR")
    m = match(nam,names(vis$nodes))
    vis$nodes$FDR = signif(vis$nodes$FDR,2)
    names(vis$nodes)[m] = c("GO term","geneSet size","% of geneSet","FDR")
    n = vis$nodes[w,m]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    t = t(n[,2:4])
    datatable(t,width=220,caption = htmltools::tags$caption(term,
                                                            style = 'caption-side: top; 
      text-align: center; color:black;background:white;font-weight:bold;'),
              options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F),
              height = 400,colnames = "") %>%
      formatStyle( target = "row", color = "black",backgroundColor = "white",
                   columns = c(1,2),fontWeight = 'bold')
  })
  
  output$bar <- renderPlot({
    req(input$current_node_id)
    vis = ig()
    n = ig()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id)
    n = vis$nodes[w,]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    s6 = mygenebarplot(n$overlapGenes)
    
    barplot(s6[[1]]$score,names.arg = s6[[1]]$gene,las=1,horiz=T,col="dodgerblue")
  })
  
  hgt = reactive({
    
    req(input$current_node_id)
    vis = ig()
    n = ig()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id)
    n = vis$nodes[w,]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    o = mygenebarplot(n$overlapGenes)
    height = mybarheight(o[[1]])
    height = height*2
  })
  
  output$ui1 = renderUI({
    #req(input$current_node_id),
    if (!is.null(input$current_node_id)){
      plotOutput("bar", height = hgt())
    }
  })
  
  output$ui2 = renderUI({
    req(ig()$nodes)
    if (!is.null(input$current_node_id))  
    {
      hr()
      dataTableOutput("table1")
    }
  })
  
  # output$gotable <- DT::renderDataTable({
  #   
  #   if(input$pool == "non") { fdata = xne} else if  (input$pool == "ess") 
  #   {fdata = xe} else fdata = xhiphop
  #   w = which(colnames(fdata) %in% input$cmp1)
  #   validate(
  #     need(length(w) != 0, "Please select a valid pool")
  #   )
  #   req(length(w)!=0)
  #   df = data.frame(gene = rownames(fdata), score = fdata[,w],index = 0,stringsAsFactors = F)
  #   w = which(df$score >= input$scorethresh)
  #   validate(need(length(w)!=0, message = "No scores above threshold"))
  #   df$index[w]=1
  #   df = df %>% arrange(desc(score))
  #   curr_exp = "tst"
  #   df = df[,c('index','score','gene')]
  #   w = which(df$index ==1)
  #   
  #   tst = myrun_go_enrich1(fdrThresh = input$fdr, curr_exp = "tst",score = df, bp_path = "BP_geneNames_current.gmt")
  #   
  #   req(input$cmp1)
  #   req(cmp1())
  #   req(input$pool)
  #   drug = cmp1()
  #   if(input$pool == "non") { fdata = xne} else if  (input$pool == "ess") 
  #   {fdata = xe} else fdata = xhiphop
  #   w = which(colnames(fdata) %in% input$cmp1)
  #   validate(
  #     need(length(w) != 0, "Please select a valid pool")
  #   )
  #   req(length(w)!=0)
  #   df = data.frame(gene = rownames(fdata), score = fdata[,w],index = 0,stringsAsFactors = F)
  #   w = which(df$score >= input$scorethresh)
  #   validate(need(length(w)!=0, message = "No scores above threshold"))
  #   df$index[w]=1
  #   df = df %>% arrange(desc(score))
  #   curr_exp = "tst"
  #   df = df[,c('index','score','gene')]
  #   w = which(df$index ==1)
  #   tst = myrun_go_enrich1(fdrThresh = input$fdr, curr_exp = "tst",score = df, bp_path = "BP_geneNames_current.gmt")
  #   enrich = tst$enrichInfo
  #   s = sapply(enrich,is.numeric)
  #   ws = which(s == T)
  #   enrich[,ws] = signif(enrich[,ws],digits = 2)
  #   datatable(enrich, rownames = F,options = list(autoWidth = T))
  # })
  
  outenrich = reactive({
    req(input$cmp1)
    #req(cmp1())
    req(input$pool)
    #drug = cmp1()
    
    if(input$pool == "non") { fdata = xne} else if  (input$pool == "ess") 
    {fdata = xe} else fdata = xhiphop
    w = which(colnames(fdata) %in% input$cmp1)
    validate(
      need(length(w) != 0, "Please select a compound")
    )
    req(length(w)!=0)
    df = data.frame(gene = rownames(fdata), score = fdata[,w],index = 0,stringsAsFactors = F)
    w = which(df$score >= input$scorethresh)
    validate(need(length(w)!=0, message = "No scores above threshold"))
    df$index[w]=1
    df = df %>% arrange(desc(score))
    curr_exp = "tst"
    df = df[,c('index','score','gene')]
    w = which(df$index ==1)
    
    tst = myrun_go_enrich1(fdrThresh = input$fdr, curr_exp = "tst",score = df, bp_path = "BP_geneNames_current.gmt")
    enrich = tst$enrichInfo
    enrich
  }
  )  
  
  output$enrich <- downloadHandler(
    filename = paste0("enrich:", Sys.Date(), ".txt"),
    content = function(file) {
      write.table(outenrich(), file, row.names = T,sep="\t")
    }
  )
  ###end go
  ig2 <- reactive({ 
    #req(input$resp)
    
    w = which(colnames(xsig) %in% input$resp)
    validate(
      need(length(w) != 0, "Please select a signature")
    )
    req(length(w)!=0)
    df = data.frame(gene = rownames(xsig), score = xsig[,w],index = 0,stringsAsFactors = F)
    w = which(df$score >= input$scorethresh2)
    validate(need(length(w)!=0, message = "No scores above threshold"))
    df$index[w]=1
    df = df %>% arrange(desc(score))
    curr_exp = "tst"
    df = df[,c('index','score','gene')]
    w = which(df$index ==1)
    
    
    tst = myrun_go_enrich1(fdrThresh = input$fdr2, curr_exp = "tst",score = df, bp_path = "BP_geneNames_current.gmt")
    validate(
      need(!is.null(tst$edgeMat), "No GO enrichment, try another signature")
    )
   
    n = tst$enrichInfo
    e = tst$edgeMat
    w = which(names(n) == "id")
    coln = ncol(n) - 1
    n = n[,c(w,1:coln)]
    
    
    w = which(names(e) == "label")
    let = graph_from_data_frame(e[,-w],vertices = n,directed = F)
    v=set_vertex_attr(let,"label",value=n$formattedLabel)
    v=set_vertex_attr(v,"color.background",value=n$cluster)
    v=set_vertex_attr(v,"label.color",value="black")
    v=set_vertex_attr(v,"label.family",value=1)
    v=set_vertex_attr(v,"label.font",value=3)
    v=set_edge_attr(v,"color",value="black")
    
    vis = toVisNetworkData(v)
    vis$nodes$label = vis$nodes$formattedLabel
    vis$nodes = vis$nodes %>% arrange(label)
    vis$edges$color = "black"
    vis$nodes$color.border="black"
    vis$nodes$font.face ="Courier"
    vis$nodes$shape = "dot"
    vis$nodes$font.size = 27
    vis$nodes$borderWidth=1
    w = which.min(vis$nodes$FDR)
    if(length(w)>0) vis$nodes$font.size[w] = 33
    if(length(w)>0) vis$nodes$borderWidth[w] = 5
    vis
  })
  
  output$network_proxy2 <- renderVisNetwork({
    validate(
      need(!is.na(ig2()$nodes),"lone 298 ig nodes"))
    w =  nrow(ig2()$nodes)
    
    validate(
      need(w != 0, "line 302 GO enrichment")
    )
    
    vis = ig2()
    g = grep("color",names(vis$nodes))
    names(vis$nodes)[g] = "color.background"
    
    vis$nodes$color.border="black"
    vis$nodes$FDR = signif(vis$nodes$FDR,2)
    visNetwork(vis$nodes, vis$edges) %>% 
      visNodes(shadow=list(enabled=T,size=25),borderWidth=1) %>%
      visOptions(nodesIdSelection = list(enabled = TRUE, 
                                         style = 'width: 150px;color: darkblue;'),
                 selectedBy = list(variable="FDR",
                                   style = 'width: 150px;color: darkblue;')) %>%
      visIgraphLayout()  %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id2', nodes.nodes);
                ;}")
  })
  
  observe({
    nodes_selection <- input$networkid_selected
    visNetworkProxy("network_proxy2_select") %>%
      visSelectNodes(id = nodes_selection)
  })
  
  observeEvent(input$current_node_id2, {
    req(input$current_node_id2)
   
    vis = ig2()
    n = ig2()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id2)
    n = vis$nodes[w,]
  })
  
  
  
  observeEvent(input$id2,
               {
                 vis = ig2()
                 w = which(vis$nodes$label == input$id2)
                 n = vis$nodes[w,]
                 
               })
  
  output$table2 <- DT::renderDataTable({
    req(input$current_node_id2)
    vis = ig2()
    n = ig2()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id2)
    term = vis$nodes$formattedLabel[w]
    nam = c("term","nGenes","geneSetFraction","FDR")
    m = match(nam,names(vis$nodes))
    vis$nodes$FDR = signif(vis$nodes$FDR,2)
    names(vis$nodes)[m] = c("GO term","geneSet size","% of geneSet","FDR")
    n = vis$nodes[w,m]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    t = t(n[,2:4])
    datatable(t,width=220,caption = htmltools::tags$caption(term,
                                                            style = 'caption-side: top; 
      text-align: center; color:black;background:white;font-weight:bold;'),
              options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F),
              height = 400,colnames = "") %>%
      formatStyle( target = "row", color = "black",backgroundColor = "white",
                   columns = c(1,2),fontWeight = 'bold')
  })
  
  output$bar2 <- renderPlot({
    req(input$current_node_id2)
    vis = ig2()
    n = ig2()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id2)
    n = vis$nodes[w,]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    s6 = mygenebarplot(n$overlapGenes)
    
    barplot(s6[[1]]$score,names.arg = s6[[1]]$gene,las=1,horiz=T,col="dodgerblue")
  })
  
  output$bar3 <- renderPlot({
    w = which(df$signature %in% input$resp)
    req(input$resp)
    validate(
      need(length(w)!= 0, "Please choose a signature")
    )
    par(mar = c(3,6,2,1))
    
    df = df[w,]
    # w = which(d != 0)
    # d = d[w,,drop = F]
    # d = data.frame(gene = rownames(d),score = d[,1],stringsAsFactors = F)
    # g = grep("NIBR",)
    d = df %>% arrange(desc(score))
    nrows = nrow(d)
    
    if(nrows > 10) d = d[1:10,]
    d = d %>% arrange(score)
    d$score = as.numeric(d$score)
    barplot(d$score,names.arg = d$gene,las=1,horiz=T,col="navy",main = paste(input$resp,d$site[1], sep=":"))
  })
  
  hgt3 = reactive({
    
    w = which(df$signature %in% input$resp)
    req(input$resp)
    validate(
      need(length(w)!= 0, "Please choose a signature")
    )
    d = df[w,]
    d = d %>% arrange(score)
    nrows = nrow(d)
    if(nrows > 10) d = d[1:10,]
    height = mybar(d)
    height = height*2
  })
  
  
  hgt2 = reactive({
    
    req(input$current_node_id2)
    vis = ig2()
    n = ig2()$nodes
    w =  nrow(n)
    
    validate(
      need(w != 0, "No GO enrichment")
    )
    req(n)
    w = which(vis$nodes$id == input$current_node_id2)
    n = vis$nodes[w,]
    validate(need(nrow(n)!=0, message = "please wait for graph to update"))
    o = mygenebarplot(n$overlapGenes)
    height = mybarheight(o[[1]])
    height = height*2
  })
  
  output$ui3 = renderUI({
    #req(input$current_node_id),
    if (!is.null(input$current_node_id2)){
      plotOutput("bar2", height = hgt2())
    }
  })
  
  output$ui5 = renderUI({
    #req(input$current_node_id),
    if (!is.null(input$resp)){
      plotOutput("bar3", height = hgt3())
    }
  })
  
  output$ui4 = renderUI({
    req(ig2()$nodes)
    if (!is.null(input$current_node_id2))  
    {
      hr()
      dataTableOutput("table2")
    }
  })
  
  
  
  ###end sig
  ###
  #### done
  output$coolplot <- renderPlot({
    par(pch=19)
    
    #req(input$pool3)
    req(input$gene)
    gene = toupper(input$gene)
    pool = input$pool1
    if(pool== "HIP") xdat = xe
    else {if(pool == "HOP") xdat = xne }
    
    par(mar = c(9,4,2,1))
    w = which(rownames(xdat) == toupper(input$gene))
    validate(
      need(length(w) == 1 ,message = 
          "please enter a valid gene"))
    myplotpts(xdat[,pgold$name], row = gene, group = pgold$cond, cex = 1.4)
  })
  
 #######################
  ##########################
 
  data_for_plot4 <- reactive({
    
    gene = toupper(input$gene)
    pool = input$pool1
    if(pool== "HIP") xdat = xe
    else {if(pool == "HOP") xdat = xne }
    
    
    #req(input$pool1)
    req(input$gene)
    gene = toupper(input$gene)
    
    w = which(rownames(xdat) == gene)
    validate(
      need(length(w) == 1 ,message = 
          "please enter a valid gene"))
    if(length(w) > 0) {

      mx2 = mymeltdf(xhiphop,row = gene,df = phiphop)
     
      mx2
    }
    
    mx2
  })
  ##########################
  output$plot4 <- renderPlot({
    x4 = data_for_plot4()
    mx2 = x4
      med = median(mx2$fitness_defect)
      #wx2 = which(mx2$gene %in% hits$target)
      w = which(hits$target %in% mx2$gene[1])
      
      wx = which(mx2$name %in% hits$name[w])
      
      mx2$fitness_defect = round(mx2$fitness_defect,2)
      mx2$sig = 0
      mx2$shape = 17
      if(length(wx) > 0) mx2$shape[wx] = 19
      mx2$shape = factor(mx2$shape)
      if(length(wx) > 0) mx2$sig[wx] = 1
      wsig = which(mx2$sig == T)
    #med = median(x4$fitness[w])
        g  = ggplot(mx2,aes(x =  name,y=fitness_defect, col = factor(sig))) + theme_bw() +
             geom_point(aes(size = 4,shape = mx2$shape, col = factor(mx2$sig)))

        g1 = g + theme(legend.position="none") +
          theme(panel.grid.minor =   element_blank()) +
          theme(panel.grid.major = element_blank()) +
          theme(axis.ticks = element_blank(), axis.text.x =   element_blank())

        g1 = g1 + labs(y = "fitness defect score") + labs(x = "drug") +
          geom_hline(yintercept=median(mx2$fitness_defect),color = "black",
            linetype = "dashed",size =1) +
          ggtitle(as.character(mx2$gene[1]))
        #if (length(wsig) > 0) {
        g2 = g1 + geom_text_repel(
          data = subset(mx2, sig == TRUE),
          aes(x =  name,y=fitness_defect, label = name),
          point.padding = 0.25,
          segment.alpha = 0.2, col = "black"
        )
    #g2 = ggplotly(g1, source= "subset") %>% layout(dragmode = "select")
    #g2 = event_register(g2, 'plotly_selected')
   g2
  })
  
  output$targethip <- renderDataTable({
    w = which(hits$name %in% input$cmp)
    datatable(hits[w,c("site","SGD")],rownames = F,options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F),
              height = 400,escape = F)
  })
  
  output$pcid <- renderDataTable({
    w = which(phiphop$name %in% input$cmp)
    datatable(phiphop[w,c("drug","PCID")],rownames = F,options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F,
                                                                    columnDefs = list(list(className = 'dt-left', targets = '_all'))),
              height = 400,escape = F)
  })
  

  
  output$densplot <- renderPlot({
    par(pch=19)
    
    #req(input$pool3)
    req(input$gene)
    gene = toupper(input$gene)
    pool = input$pool1
    if(pool== "HIP") xdat = xe
    else {if(pool == "HOP") xdat = xne }
    
    #par(mar = c(9,4,2,1))
    w = which(rownames(xdat) == toupper(input$gene))
    validate(
      need(length(w) == 1 ,message = 
             "please enter a valid gene"))
    x = xdat[w,]
    z = (x - median(x))/mad(x)
    dens = density(z)
    
    lens = length(z)/512
    plot(dens$x,dens$y,col="dodgerblue",cex = 1.5)
    wz = which(dens$x>= 3.09)
    #wz = w/512
    points(dens$x[wz],dens$y[wz],col="red",cex = 1.5)
    
  })
  ##########################
  brush <- reactive({
    mx2 = data_for_plot4()
    mx2$fitness_defect = round(mx2$fitness_defect,2)
    mx = brushedPoints(mx2, input$brush4, xvar = "name", yvar = "fitness_defect")
    mx3 = mx$name
    
    mx3
  })
  
  observeEvent(brush(),{
    drug = brush()
   
    w = which(phiphop$name %in% drug)
    
    x = phiphop %>% filter(name %in% drug) #%>% select(name) 
    x = phiphop[w,]
    
    x = x$name
    
    if(length(drug) == 0) x = phiphop$name
    
    
    updateSelectizeInput(session,'cmp',label = "select condition",
      choices = x, selected = x[1])#if there weren't shared conditions, this would have to be specified wrt x, e.g. x[1]
    updateSelectizeInput(session,'cmp1',label = "",
                         choices = x, selected = x[1])
  },ignoreInit =T, ignoreNULL = F)
  
 
  
  xdat1 = reactive({
    if(input$pool1 == "HOP") xdat = xne 
    else if(input$pool1 == "HIP") {xdat = xe } 
    else if(input$pool1 == "HIPHOP") {xdat = xhiphop }
    xdat
  })
  

  output$fd <- renderPlot({
    
    par(mar = c(9,3,2,1))
    pdat = phiphop
    xdat = xdat1()
    w1 = which(pdat$name %in% input$cmp)
    validate(
      need(length(w1) == 1 ,message = 
          "please enter a valid condition"))
   
    p10(xdat,w1[1],pch =17)
    
    
  })
  
  output$efd = renderPlot({
    validate(
      need(input$cmp,message = 
          "please select condition"))
    cond1 = input$cmp
    
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    we = which(rownames(xhiphop)%in% ess$sgd_gene)
    
    xess = xhiphop[we,]
    
    wx = which(colnames(xess) %in% cond1)
    
    colnames(xess)[wx] = paste0("HIP|",colnames(xess)[wx])
    
    
    p10(xess,wx,pch = 19)
    
    
  })
  
  output$nfd = renderPlot({
    validate(
      need(input$cmp,message = 
          "please select condition"))
    
    
    cond1 = input$cmp
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    wne = which(rownames(xhiphop)%in% noness$sgd_gene)
    
    xnon = xhiphop[wne,]
    wx = which(colnames(xnon) %in% cond1)
    colnames(xnon)[wx] = paste0("HOP|",colnames(xnon)[wx])
    p10(xnon,wx,pch = 19)
    
    
  })
  output$info <- renderDataTable({
    mx2 = data_for_plot4()
    mx2$fitness_defect = round(mx2$fitness_defect,2)
    mx2 = mx2[,c("screen","gene","fitness_defect","drug","dose")]
    
    brushedPoints(mx2, input$brush4, xvar = "screen", yvar = "fitness_defect")
    
  },rownames = F, escape = F)
  
  out <- reactive({
    mx2 = data_for_plot4()
    
    dx2 = dcast(mx2, formula = exp~gene,value.var = "fitness_defect")
    dx2[,2] = round(dx2[,2],2)
    dx2
  })
  
outfd = reactive({
  pdat = phiphop
  xdat = xdat1()
  w1 = which(pdat$name %in% input$cond)
  validate(
    need(length(w1) == 1 ,message = 
        "please enter a valid condition"))
 
  #if(length(w1) == 0) p10(xdat,185,pch = 17)
  fd = xdat[,w1[1],drop = F]
  colnames(fd) = input$cond
  fd
})

outhip = reactive({
  
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
        "please enter a valid condition"))
  
  #if(length(w1) == 0) p10(xdat,185,pch = 17)
  fd = xe[,w,drop = F]
  colnames(fd) = paste0("HIP:",input$cmp)
  fd
})

output$downloadhip <- downloadHandler(
  filename = paste0("HIP:",input$cmp, Sys.Date(), ".txt"),
  content = function(file) {
    write.table(as.data.frame(outhip()), file, row.names = T,sep="\t",quote=F)
  }
)
outhop = reactive({
  
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
        "please enter a valid condition"))
  
  #if(length(w1) == 0) p10(xdat,185,pch = 17)
  fd = xne[,w,drop = F]
  colnames(fd) = paste0("HOP:",input$cmp)
  fd
})

output$hoptab = renderDataTable({
  
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  
  #if(length(w1) == 0) p10(xdat,185,pch = 17)
  fd = xne[,w,drop = F]
  colnames(fd) = paste0("HOP:",input$cmp)
  fd = fd[order(fd[,1],decreasing = T),,drop=F]
  w = which(fdat$sgd_gene %in% rownames(fd))
  
  f = fdat[w,c("sgd_orf","sgd_gene","descriptor")]
  wf = which(duplicated(f$sgd_gene))
  f = f[-wf,]
  
  wf = which(rownames(fd) %in% f$sgd_gene)
  m = match(rownames(fd)[wf],f$sgd_gene)
  f = f[m,]
  f$sgd_gene = paste0("<a href=https://www.yeastgenome.org/locus/",f$sgd_gene,">",f$sgd_gene,"</a>")
  names(f)=c("ORF","GENE","descriptor")
  f
},rownames = F,escape = F,options=list(pageLength=10))

output$hiptab = renderDataTable({
  
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  
  #if(length(w1) == 0) p10(xdat,185,pch = 17)
  fd = xe[,w,drop = F]
  colnames(fd) = paste0("HIP:",input$cmp)
  fd = fd[order(fd[,1],decreasing = T),,drop=F]
  w = which(fdat$sgd_gene %in% rownames(fd))
  
  f = fdat[w,c("sgd_orf","sgd_gene","descriptor")]
  wf = which(duplicated(f$sgd_gene))
  f = f[-wf,]
  
  wf = which(rownames(fd) %in% f$sgd_gene)
  m = match(rownames(fd)[wf],f$sgd_gene)
  f = f[m,]
  f$sgd_gene = paste0("<a href=https://www.yeastgenome.org/locus/",f$sgd_gene,">",f$sgd_gene,"</a>")
  names(f)=c("ORF","GENE","descriptor")
  f
},escape = F,rownames=F,options=list(pageLength=10))


output$coinhib = renderDataTable({
  xcoinh = xcoinh[,phiphop$name]
  xcoinh = xcoinh[phiphop$name,]
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  
  d = xcoinh[,w,drop = F]
  
  d = d[order(d[,1],decreasing = T),,drop=F]
  print(paste("dim d=",dim(d)))
  df = data.frame(screen = rownames(d)[1:50], coinhibition = d[1:50,1],stringsAsFactors = F)
  
  m = match(df$screen,phiphop$name)
  df$coinhibition = round(df$coinhibition,2)
  df$compound = phiphop$drug[m]
  #df$screen = paste0("<a href='#hiphop'>", df$screen, "</a>")
  df$pcid = phiphop$PCID[m]
  
  df
  
},escape = F,rownames=F,options=list(pageLength=10),selection = "single",server = F)
 
observeEvent(input$coinhib_rows_selected, {
  row <- input$coinhib_rows_selected 
  
  xcoinh = xcoinh[,phiphop$name]
  xcoinh = xcoinh[phiphop$name,]
  w = which(phiphop$name %in% input$cmp)
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  
  d = xcoinh[,w,drop = F]
  
  d = d[order(d[,1],decreasing = T),,drop=F]
  
  df = data.frame(screen = rownames(d)[1:10], coinhibition = d[1:10,1],stringsAsFactors = F)
  
  m = match(df$screen,phiphop$name)
  df$coinhibition = round(df$coinhibition,2)
  df$compound = phiphop$drug[m]
  #df$screen = paste0("<a href='#hiphop'>", df$screen, "</a>")
  df$pcid = phiphop$PCID[m]
  
  df
  print(row)
  #output$text <- renderText({paste("X =", X[row, "x"], "Y =", X[row, "y"])})
  #row = row$coinhibition
  updateSelectizeInput(session,'cmp',label = "",
                       choices = df$screen[row], selected = df$screen[row][1])
})  

output$downloadhop <- downloadHandler(
  filename = paste0("HOP:",input$cmp, Sys.Date(), ".txt"),
  content = function(file) {
    write.table(as.data.frame(outhop()), file, row.names = T,sep="\t",quote=F)
  }
)

output$go_table <- renderDataTable({
  
  z = pGO
 
  w = which(z$condition %in% input$cond)
  
  zw = z[w,]
  
  zw
  
},escape = F,rownames = F,options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F))

output$GO <- renderDataTable({
  
  z = pGO
  
  w = which(tolower(z$condition) %in% tolower(input$cmp))
  
  zw = z[w,]
  
  zw
  
},escape = F,rownames = F,options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F))

}

shinyApp(ui, server)
