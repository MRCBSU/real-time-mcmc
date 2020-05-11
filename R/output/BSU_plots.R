make.plots <- function(projections, ylab = "", data = NULL, y.format = ".3s") {
  Eng_projection <- sum.all(projections) %>% add.quantiles(NULL, QUANTILES)
  Eng_data <- sum.all.data(data)
  plot.height <- 420 + 150
  date <- ymd(date.data)
  lines <- list(
    list(
      type = "line", 
      y0 = 0, 
      y1 = 1, 
      yref = "paper",
      x0 = ymd(date.data), 
      x1 = ymd(date.data), 
      line = list(color = "red"),
      hoverinfo = "Today"
    ),
    list(
      type = "line", 
      y0 = 0, 
      y1 = 1, 
      yref = "paper",
      x0 = ymd(20200323), 
      x1 = ymd(20200323), 
      line = list(color = "blue"),
      hoverinfo = "Lockdown"
    )
  )
  create.subplot <- function(projections, data) {
    plot <- projections %>%
      pivot_wider(names_from = quantile) %>%
      plot_ly(x = ~date, width = 800, height = plot.height) %>%
      add_ribbons(ymin = ~`0.025`, ymax = ~`0.975`, color = I("lightblue2"), alpha = 0.25) %>%
      add_lines(y = ~`0.5`, color = I("black")) %>%
      layout(shapes = lines, showlegend = FALSE, xaxis = list(title = "Date"),
             hovermode = "x unified")
    if (is.null(data)) return(plot)
    plot %>%
      add_markers(
        data = data,
        x = ~date,
        y = ~True,
        color = I("red"),
        hovertemplate = "%{x}: %{y:.3s}<extra>Actual report</extra>"
      )
  }
  return(create.subplot(Eng_projection, Eng_data))
}

make.plots(infections) %>% layout(yaxis = list(title = "Number of new infections"))
make.plots(noisy_deaths, data = dth.dat) %>% layout(yaxis = list(title = "Number of deaths"))