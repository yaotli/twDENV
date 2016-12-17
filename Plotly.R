library(plotly)

Sys.setenv("plotly_username"="yaoli")
Sys.setenv("plotly_api_key"="ylcqmx9pdd")

#3D plot for temporal-spatial analysis #https://plot.ly/~yaoli/2/ii-ib-ia/
td<-plot_ly(gisp, x = Longitude, y = Latitude, z = Date,
            text = paste("Clarity: ", Group),
            type="scatter3d", mode="markers", color= Group)

plotly_POST(td, filename = "Temporal-spatial distribution of DENV2", sharing = "public")


# ploty
# https://plot.ly/~yaoli/4/ia-ib-ii/

library(plotly)

Sys.setenv("plotly_username"="yaoli")
Sys.setenv("plotly_api_key"="ylcqmx9pdd")

p<-plot_ly(dbx_a, x = ~onset, y = ~V5, z = ~log10(VL), color = ~clade) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Onset'),
                      yaxis = list(title = 'Variation'),
                      zaxis = list(title = 'log10(VL)')))

plotly_POST(p, filename = "VL vs Onset vs vatiation", sharing = "public")