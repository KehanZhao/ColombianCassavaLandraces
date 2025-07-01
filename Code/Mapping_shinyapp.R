library(shiny)
library(dplyr)
library(plotly)
library(data.table)
setwd("D:/WorkSpace/Rstudio_Workspace/Map_PCA")


EnvPhenotypes <- fread("Colombia_Batch1234_Phenos.txt")

# Cluster <- read.table("Colombia_Kinship_Clusters.txt",header = T)
# colnames(Cluster)[1] <- "Accession"
# EnvPhenotypes <- left_join(EnvPhenotypes,Cluster)
# EnvPhenotypes$Cluster <- as.factor(EnvPhenotypes$Cluster)
Phenotypes <- colnames(EnvPhenotypes)[c(2:9,11:ncol(EnvPhenotypes))]
colnames(EnvPhenotypes)[1] <- "Accession"

EnvPhenotypes$bio1 <- EnvPhenotypes$bio1/10
EnvPhenotypes$bio12 <- EnvPhenotypes$bio12/10

write.csv(EnvPhenotypes,"Map_Shiny_input.csv")
# Define the UI

ui <- fluidPage(
  titlePanel("Colombia Cassava Locations"),
  sidebarLayout(
    sidebarPanel(
      selectInput("color_by", "Color points by:", choices = Phenotypes),
      selectInput("variable_x", "Variable X:", choices = Phenotypes),
      selectInput("variable_y", "Variable Y:", choices = Phenotypes)
    ),
    mainPanel(
      plotlyOutput("map"),
      plotlyOutput("scatterplot")
    )
  )
)

# Define the server
server <- function(input, output) {
  # Read the data
  cassava_data <- reactive({
    read.csv("Map_Shiny_input.csv")
  })
  # Create the map plot
  output$map <- renderPlotly({
    data_to_plot <- cassava_data()
    
    plot_ly(data_to_plot, 
            type = "scattermapbox",
            lat = ~Latitude,
            lon = ~Longitude,
            mode = "markers",
            marker = list(
              size = 8,
              color = data_to_plot[[input$color_by]],
              colorscale = "Viridis",
              showscale = TRUE
            ),
            text = ~Accession,
            hoverinfo = "text") %>%
      layout(mapbox = list(
        style = "open-street-map",
        center = list(lat = 4.5709, lon = -74.2973),
        zoom = 4.5
      ))
  })
  
  # Create the scatterplot
  output$scatterplot <- renderPlotly({
    data_to_plot <- cassava_data()
    
    plot_ly(data_to_plot, 
            x = ~get(input$variable_x),
            y = ~get(input$variable_y),
            mode = "markers",
            marker = list(
              size = 8,
              color = ~get(input$color_by),
              colorscale = "Viridis",
              showscale = TRUE
            ),
            text = ~Accession,
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = input$variable_x),
        yaxis = list(title = input$variable_y)
      )
  })
  
  
}

# Run the application
shinyApp(ui = ui, server = server)


Colombia <- read.csv("Map_Shiny_input.csv")
#Colombia[,grepl("bio",colnames(Colombia))] <- Colombia[,grepl("bio",colnames(Colombia))]/10 

ggplot(Colombia,aes(x=Annual.Mean.Temperature))+
  geom_histogram(bins=20,color = "#000000", fill = "#0099F8")+
  labs(
    title = "Annual Mean Temperature",
    caption = "Source: bioclim",
    x = "Degrees (C)",
    y = "Count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(color = "#000000", size = 20, face = "bold", hjust = 0.5),
    plot.caption = element_text(face = "italic"),
    axis.text = element_text(color = "#000000",size = 14),
    axis.title = element_text(color = "#000000",size = 16)
  )




#####testing


meantemp <- lm(Annual.Mean.Temperature ~ PC1+PC2+PC3+PC4+PC5,data=EnvPhenotypes)
summary(meantemp)

