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
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20200323), 
        x1 = ymd(20200511), 
		fillcolor = "black",
		opacity = 0.15,
        hoverinfo = "First national lockdown"
    ),
    list(
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20201105), 
        x1 = ymd(20201202), 
		fillcolor = "black",
		opacity = 0.15,
        hoverinfo = "Second national lockdown"
    ),
    list(
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210105), 
        x1 = ymd(20210308), 
		fillcolor = "black",
		opacity = 0.15,
        hoverinfo = "Third national lockdown"
    ),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210308), 
        x1 = ymd(20210308), 
        line = list(color = "green"),
        hoverinfo = "Step one of roadmap"
	),
    list(
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210308), 
        x1 = ymd(20210329), 
		fillcolor = "black",
		opacity = 0.12,
        hoverinfo = "Third national lockdown"
    ),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210329), 
        x1 = ymd(20210329), 
        line = list(color = "green"),
        hoverinfo = "Step one of roadmap"
	),
    list(
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210329), 
        x1 = ymd(20210412), 
		fillcolor = "black",
		opacity = 0.09,
        hoverinfo = "Third national lockdown"
    ),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210412), 
        x1 = ymd(20210412), 
        line = list(color = "green"),
        hoverinfo = "Step two of roadmap"
	),
    list(
        type = "rect", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210412), 
        x1 = ymd(20210517), 
		fillcolor = "black",
		opacity = 0.06,
        hoverinfo = "Third national lockdown"
	),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210517), 
        x1 = ymd(20210517), 
        line = list(color = "green"),
        hoverinfo = "Step three of roadmap"
	),
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = ymd(20210719), 
        x1 = ymd(20210719), 
        line = list(color = "green"),
        hoverinfo = "Step four of roadmap"
	)
  )
  lines <- lines[sapply(lines, function(x) x$x0 %in% data$date)]
  plot.height <- num.rows * 420 + 150
  return(
    plot_ly(data, x = ~date, width = 800, height = plot.height) %>%
      layout(shapes = lines, annotations = title(subplot_title), showlegend = FALSE,
             hovermode = "x unified")
  )
}

make.plots <- function(projections, ylab = "", by = NULL, data = NULL,
                       y.format = ".3s}", combine.to.England = sum.all,
                       combine.data.to.England = sum.all.data,
                       project.forwards = !(is.null(data) && external), x.label = "",
                       denoms = NULL, y.percent = FALSE) {

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
  create.subplot <- function(projections, subplot_title, data, denom) {
    plot_dat <- projections %>%
      pivot_wider(names_from = quantile)
    if (!is.null(denom)) {
      plot_dat <- plot_dat %>%
        mutate_at(vars(starts_with("0")), list(prop = ~./denom))
      y.format <- paste0(y.format, " (%{text:.2%})")
    } else {
      plot_dat <- plot_dat %>%
        mutate_at(vars(starts_with("0")), list(prop = ~1))
    }

    plot <- plot_dat %>%
      create.base.subplot(num.rows, subplot_title) %>%
      add_ribbons(ymin = ~`0.025`, ymax = ~`0.975`, text = ~`0.975_prop`, color = I("lightblue2"), alpha = 0.25,
                hovertemplate = paste0("%{x}: %{y:", y.format, "<extra>Upper 95% CrI</extra>")) %>%
      add_lines(y = ~`0.025`, alpha = 0, text = ~`0.025_prop`,   # An extra trace just for hover text
                hovertemplate = paste0("%{x}: %{y:", y.format, "<extra>Lower 95% CrI</extra>")) %>%
      add_lines(y = ~`0.5`, color = I("black"), text = ~`0.5_prop`,
                hovertemplate = paste0("%{x}: %{y:", y.format, "<extra>Median</extra>")) %>%
      layout(xaxis = list(title = x.label))
    if (y.percent) plot <- plot %>% layout(yaxis = list(tickformat = "%"))
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
  if (!is.null(Eng_projection)) {
    if (is.null(denoms)) {
      Eng_denom <- NULL
    } else {
      Eng_denom <- sum(denoms$denom)
    }
    plots <- list("England" = create.subplot(Eng_projection, "England", Eng_data, Eng_denom))
  }
  for(subplot in plot.names) {
    if (is.null(denoms)) {
      denom <- NULL
    } else {
      denom <- denoms %>%
        filter(name == subplot) %>%
        `$`(denom)
    }
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
      create.subplot(subplot, plot.data, denom)
  }
  return(plotly::subplot(plots, nrows = num.rows))
}
