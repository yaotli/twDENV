library(plotly)

Sys.setenv("plotly_username"="yaoli")
Sys.setenv("plotly_api_key"="ylcqmx9pdd")

#3D plot for temporal-spatial analysis #https://plot.ly/~yaoli/2/ii-ib-ia/
td<-plot_ly(gisp, x = Longitude, y = Latitude, z = Date,
            text = paste("Clarity: ", Group),
            type="scatter3d", mode="markers", color= Group)

plotly_POST(td, filename = "Temporal-spatial distribution of DENV2", sharing = "public")