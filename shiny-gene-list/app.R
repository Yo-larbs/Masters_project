#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
install.packages("shinythemes")
library(shiny)
library(shinydashboard)
library(bslib)
library(shinythemes)
library(plotly)
#install.packages("htmlwidgets")
#knitr::purl("../scripts/cell_line_classification/cleaning_and_classification_vector.Rmd",output)
#source(output)
#knitr::purl("scripts/cell_line_classification/cleaning_and_classification_virus.Rmd",output)
#source(output)
main_menu=(
  
  card(
    card_header("PCA Dashboard", 
                tooltip(
                  bsicons::bs_icon("question-circle"),
                  "click on  points for cell line specifics",
                  placement = "right"
                ),
                popover( bsicons::bs_icon("gear"),
                         actionButton("reset", label = "reset",icon = icon("rewind")),placement = "right"
                         
                ),
    ),
    layout_sidebar(fillable = TRUE,
                   sidebar = sidebar(
                     sliderInput("bins",label = "number of clusters:",min = 1,max = 50,value = 8),
                     radioButtons("assay",label = "which assay?",choices = c("vector","virus"),selected = "vector"),
                     input_switch("contributions","overlay parameter contributions?",value = F)
                     ),
                   card_body( plotlyOutput("PCA",height = "85vh",width = "100%",fill = T),fill = T
                   )
                   
                                 ),card_footer("designed by Yoofi Larbi",tooltip(
                                   bsicons::bs_icon("Exclamation-circle"),
                                   "Special thanks to Dr.Alessandra Vigilante and Dr. Luis Apolonia",
                                   placement = "right"
                                 )))
  )

  
specific_graph=(
  navset_tab(id = "assay_selection",header = actionButton("back_button",label = "",icon = icon("circle-xmark")),
    nav_panel("Vector_cell",
              card(card_header("Box plot"),
                   card_body(plotlyOutput("Box_plot"),fill = T)
                
              ),
              card(card_header("logistic plot"),
                     card_body(plotlyOutput("logistic_plot"),fill = T)
              ),
              card(card_header("logarithmic plot"),
                     card_body(plotlyOutput("logarithmic_plot"),fill = T)
              )),
    
    nav_panel("Virus_cell", 
              card(card_header("Box plot"),
                card_body(plotlyOutput("Box_plot_virus",fill = T))
    ),
    card(card_header("logistic plot"),
         card_body(plotlyOutput("logistic_plot_virus"))
    ),
    card(card_header("logarithmic plot"),
         card_body(plotlyOutput("logarithmic_plot_virus"))
    )))
  )
  
  
  
  # Application title
  
  # Show a plot of the generated distribution

# Define UI for application that draws a histogram
ui <- page_fillable(uiOutput("dynamic_Ui"))


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  
  
  selected_group <- reactiveVal(NULL)
  
  last_click <- reactiveVal(NULL)
  
  Specific_cell=reactiveVal(NULL)
  
  ui_refresh <- reactiveVal(0)
  
  observeEvent(Specific_cell(), {
    ui_refresh(ui_refresh() + 1)
  })
  
  
  output$dynamic_Ui <- renderUI({
    
    if(is.null(Specific_cell())){
      ui_refresh(0)
      main_menu
      
      
    }
    else{
      specific_graph
      
      
    }
  })
  
  output$cell_name=renderUI({
    headerPanel(Specific_cell())
  })
  # this is a reactive for PCA which specifies which assay would be used depending on the radio button
  PCA_value <- reactiveVal(NULL)
  
  # this is for the overlay of the PCA contributions
  contr_df=reactiveVal(NULL)
  
  
    

#This checks if its virus or vector   
    assay_obs_virus= observe({
      req(input$assay, input$bins)
      if (input$assay=="virus"){
        updateSliderInput(session,"bins",label = "number of clusters:",min = 1,max = 50,value = 10)
        update_switch("contributions",label = "overlay parameter contributions?",value = F )
        set.seed(39)
        clusters <- kmeans(ind_coord_virus[, c("Dim.1", "Dim.2")], centers = input$bins)
        ind_coord_virus$group <- as.factor(clusters$cluster)
        PCA_value(ind_coord_virus)
        contr_df(as.data.frame(res_var_virus$coord))
        selected_group(NULL)
        assay_obs_virus$suspend()
        assay_obs_vector$resume()
    }})
    
    assay_obs_vector= observe({
      req(input$assay, input$bins)
      if (input$assay=="vector"){
        updateSliderInput(session,"bins",label = "number of clusters:",min = 1,max = 50,value = 8)
        update_switch("contributions",label = "overlay parameter contributions?",value = F )
        set.seed(18)
        clusters <- kmeans(ind_coord[, c("Dim.1", "Dim.2")], centers = input$bins)
        selected_group(NULL)
        ind_coord$group <- as.factor(clusters$cluster)
        PCA_value(ind_coord)    
        contr_df(as.data.frame(res_var$coord))
        assay_obs_vector$suspend()
        assay_obs_virus$resume()
      }})
    
    observeEvent(input$bins,{
      if (input$assay=="vector"){
      clusters <- kmeans(ind_coord[, c("Dim.1", "Dim.2")], centers = input$bins)
      selected_group(NULL)
      ind_coord$group <- as.factor(clusters$cluster)
      PCA_value(ind_coord)    
    }
      else{
        clusters <- kmeans(ind_coord_virus[, c("Dim.1", "Dim.2")], centers = input$bins)
      ind_coord_virus$group <- as.factor(clusters$cluster)
      PCA_value(ind_coord_virus)
      }})

    
    
    #PCA plot 
    #PCA plot function
    
    output$PCA <- renderPlotly({
      df=PCA_value()
      dummy=ui_refresh()
        # generate bins based on input$bins from ui.R
      PCA_obj=ggplot(data = df,aes(x=Dim.1,y = Dim.2,label=row.names(df),text= paste("cell_line:",row.names(df),"<br>PC1:",signif(Dim.1,3),"<br>PC2:",signif(Dim.2,3))))+
        geom_point(size = 2,aes(colour = group))+
        theme_classic()+
        xlab("PC1")+
        ylab("PC2")
      
      if (!is.null(selected_group())){
        PCA_obj=PCA_obj+
        geom_text(aes(label=ifelse(group==selected_group(),row.names(df),'')))
      }
        
     
       p=ggplotly(PCA_obj,source = "PCA",tooltip = "text")%>%
         event_register(event = "plotly_legendclick")%>%
         event_register(event = "plotly_click")%>%
        layout(legend = list(itemclick = F,itemdoubleclick = F))
       
       
       
       
       if (input$contributions){
         df=contr_df()
         for (i in 1:nrow(df)){
         p=p%>%layout(annotations=list(
         list(
           x=(df$Dim.1[i]*3.5),
           y=(df$Dim.2[i]*3.5),
           showarrow=TRUE,
           xref = "x",
           yref = "y",
           axref = "x",
           ayref = "y",
           arrowhead = 2,
           arrowsize = 1,
           ax=0,
           ay=0,
           opacity=0.5,
           arrowcolor='blue'
          
         )
         ))
         
         
      p=p%>%layout(annotations=list(
        list( 
      text=row.names(df)[i],
      font = list(size = 15),
      xanchor = "left",  # optional: controls horizontal alignment
      yanchor = "middle",
      xshift = 5,        # extra pixel offset to avoid overlap
      #yshift = 0,
      x =  df$Dim.1[i] *3.5 + 0,
      y =  df$Dim.2[i] *3.5 + 0.15,
      showarrow=F,
      xref = "x",
      yref = "y"
      )))
         }}
    p
    })
    
    output$Box_plot = renderPlotly({
      req(Specific_cell())
        box_plot_vector[[Specific_cell()]]
      
        
    })
    
    output$Box_plot_virus = renderPlotly({
      req(Specific_cell())
      box_plot_virus[[Specific_cell()]]
      
    })
    
    
    output$logarithmic_plot = renderPlotly({
      req(Specific_cell())
        logarithmic_plot_vector[[Specific_cell()]]

    })
    
    output$logarithmic_plot_virus = renderPlotly({
      req(Specific_cell())
      logarithmic_plot_virus[[Specific_cell()]]
      
      
    })
    
    
    output$logistic_plot = renderPlotly({
      req(Specific_cell())
        logistic_plot_vector[[Specific_cell()]]
      
    })
    
    output$logistic_plot_virus = renderPlotly({
      req(Specific_cell())
      logistic_plot_virus[[Specific_cell()]]

    })
    
    
    observeEvent(event_data("plotly_legendclick",source = "PCA",priority = "event"),{
      ed=event_data("plotly_legendclick",source = "PCA")
      print(as.integer(ed[[1]]["legendgroup"]))
      
      if (!is.null(selected_group()) && selected_group()==as.integer(ed[[1]]["legendgroup"])){
        selected_group(NULL)
      }
      else{
      selected_group(as.integer(ed[[1]]["legendgroup"]) )
      }
    })
    
    observeEvent(event_data("plotly_click",source = "PCA",priority = "event"),{
      ed2=event_data("plotly_click",source = "PCA")
      print(((ed2)))#["name"]))
      if (input$assay=="vector"){
        df=ind_coord
        df$cell_line=row.names(df)
        cell_key= inner_join(df,ed2,join_by(Dim.1==x,Dim.2==y))
      }
      else if(input$assay=="virus"){
        df=ind_coord_virus
        df$cell_line=row.names(df)
        cell_key= inner_join(df,ed2,join_by(Dim.1==x,Dim.2==y))
      }
      print(cell_key$cell_line)
      Specific_cell(cell_key$cell_line)
    })
  
    observeEvent(input$reset,{
      selected_group(NULL)
      Specific_cell(NULL)
      ui_refresh(ui_refresh()+1)
      updateSliderInput(session,"bins",label = "number of clusters:",min = 1,max = 50,value = 8)
      updateRadioButtons(session,"assay",label = "which assay?",choices = c("vector","virus"),selected = "vector")
      update_switch("contributions",label = "overlay parameter contributions?",value = F )
      
      
      
    })
    observeEvent(input$back_button,{
      Specific_cell(NULL)
      selected_group(NULL)
      ui_refresh(ui_refresh()+1)
      
      
    })
}
# Run the application 
shinyApp(ui = ui, server = server)
