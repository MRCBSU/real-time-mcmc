library(lubridate)
library(plotly)
if (!exists("sum.all")) source(file.path(Rfile.loc, "results_api.R"))
QUANTILES <- c(0.025, 0.5, 0.975)
title <- function(name) {
  list(
    text = str_replace_all(name, "_", " "),
    xref = "paper",
    yref = "paper",
    yanchor = "top",
    xanchor = "left",
    x = 0.1,
    y = 0.9,
    font = list(size = 20),
    showarrow = FALSE
  )
}

create.base.subplot <- function(data, num.rows, subplot_title) {
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
    ),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20200511), 
        x1 = ymd(20200511), 
        line = list(color = "blue"),
        hoverinfo = "Lockdown"
    )
  )
  plot.height <- num.rows * 420 + 150
  return(
    plot_ly(data, x = ~date, width = 800, height = plot.height) %>%
      layout(shapes = lines, annotations = title(subplot_title), showlegend = FALSE,
             hovermode = "x unified")
  )
}

make.plots <- function(projections, ylab = "", by = NULL, data = NULL,
                       y.format = ".3s", combine.to.England = sum.all,
                       combine.data.to.England = sum.all.data,
                       project.forwards = !(is.null(data) && external), x.label = "") {

  if (is.null(combine.to.England)) {
    Eng_projection <- Eng_data <- NULL
  } else {
    Eng_projection <- combine.to.England(projections) %>% add.quantiles(NULL, QUANTILES)
    Eng_data <- combine.data.to.England(data)
  }
  projections <- get.aggregated.quantiles(projections, by, c(0.025, 0.5, 0.975)) %>%
    rename(by = all_of(by))
  if (!project.forwards) {
    if (!is.null(Eng_projection)) {
      Eng_projection <- filter(Eng_projection, date <= ymd(date.data))
    }
    projections <- filter(projections, date <= ymd(date.data))
  }
  plot.names <- unique(projections$by)
  num.plots <- length(plot.names)
  if (!is.null(Eng_projection)) num.plots <- num.plots + 1
  num.rows <- ceiling(num.plots/2)
  date <- ymd(date.data)
  create.subplot <- function(projections, subplot_title, data) {
    plot <- projections %>%
      pivot_wider(names_from = quantile) %>%
      create.base.subplot(num.rows, subplot_title) %>%
      add_ribbons(ymin = ~`0.025`, ymax = ~`0.975`, color = I("lightblue2"), alpha = 0.25,
                hovertemplate = paste0("%{x}: %{y:", y.format, "}<extra>Upper 95% CrI</extra>")) %>%
      add_lines(y = ~`0.025`, alpha = 0,   # An extra trace just for hover text
                hovertemplate = paste0("%{x}: %{y:", y.format, "}<extra>Lower 95% CrI</extra>")) %>%
      add_lines(y = ~`0.5`, color = I("black"),
                hovertemplate = paste0("%{x}: %{y:", y.format, "}<extra>Median</extra>")) %>%
      layout(xaxis = list(title = x.label))
    if (is.null(data)) {
      return(plot)
    } else {
      if (external) {
        data.hover <- "%{x} <extra>Observed deaths</extra>"
      } else {
        data.hover <- "%{x}: %{y:.3s}<extra>Observed deaths</extra>"
      }
      plot %>%
        add_markers(
          data = data,
          x = ~date,
          y = ~True,
          color = I("red"),
          hovertemplate = data.hover
        )
    }
  }
  plots <- NULL
  if (!is.null(Eng_projection)) plots <- list("Total" = create.subplot(Eng_projection, "Total", Eng_data))
  for(subplot in plot.names) {
    if (is.null(data)) {
      plot.data <- NULL
    } else {
      plot.data <- data %>%
        filter(!!sym(by) == subplot) %>%
        group_by(date) %>%
        summarise(True = sum(value))
    }
    plots[[subplot]] <- projections %>%
      filter(by == subplot) %>%
      create.subplot(subplot, plot.data)
  }
  return(subplot(plots, nrows = num.rows))
}
