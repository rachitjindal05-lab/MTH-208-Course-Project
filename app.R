library(shiny)
library(shinydashboard)
library(bslib)
library(quantmod)      
library(TTR)
library(DT)
library(xts)
library(zoo)
library(yfR)
library(PerformanceAnalytics)
library(quadprog)
library(ggplot2)
gold_dark_theme <- bs_theme(
  version = 5,
  bg = "#0d0d0d",
  fg = "#f5d76e",
  primary = "#d4af37",
  secondary = "#2e2e2e",
  base_font    = font_google("Poppins"),
  heading_font = font_google("Poppins"),
  code_font    = font_google("JetBrains Mono"),
  "navbar-bg" = "#1a1a1a",
  "body-bg" = "#0d0d0d",
  "card-bg" = "#151515",
  "card-border-color" = "#2c2c2c",
  "link-color" = "#ffcf40"
)
ui <- dashboardPage(
  dashboardHeader(title = "OHLCV + Portfolio Dashboard (MTH208 Group Project)"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Technical Analysis RSI/MA/EMA", tabName = "stock", icon = icon("chart-line")),
      menuItem("Portfolio Optimize (Markowitz Model)", tabName = "portfolio", icon = icon("chart-pie")),
      menuItem("Option Chain Data", tabName = "options", icon = icon("money-bill-trend-up")),
      menuItem("Option Chain + Black–Scholes", tabName = "options_bsm", icon = icon("money-bill-trend-up"))
    )
  ),
  dashboardBody(
    theme = gold_dark_theme,
    tabItems(
      tabItem(tabName = "stock",
              fluidRow(
                box(width = 3, title = "Controls", solidHeader = TRUE,
                    textInput("ticker", "Ticker", "AAPL"),
                    dateRangeInput("dates", "Date range",
                                   start = Sys.Date() - 365, end = Sys.Date()),
                    selectInput("indicator", "Indicator",
                                c("SMA", "EMA", "Bollinger Bands", "RSI", "MACD")),
                    numericInput("n", "Period", 20, min = 2),
                    actionButton("go", "Load / Refresh", class = "btn btn-warning")
                ),
                box(width = 9,
                    tabsetPanel(
                      tabPanel("Candlestick Chart", plotOutput("chart", height = "600px")),
                      tabPanel("OHLCV Table", DTOutput("table"))
                    ))
              )
      ),
      tabItem(tabName = "portfolio",
              fluidRow(
                box(width = 4, title = "Portfolio Inputs", solidHeader = TRUE,
                    textAreaInput("tickers_multi", "Enter stock tickers (comma separated)",
                                  "AAPL,MSFT,GOOG,AMZN,TSLA"),
                    dateRangeInput("dates_multi", "Date range",
                                   start = Sys.Date() - 365, end = Sys.Date()),
                    actionButton("load_multi", "Load & Optimize", class = "btn btn-warning")
                ),
                box(width = 8, title = "Efficient Frontier", solidHeader = TRUE,
                    plotOutput("frontier", height = "500px"),
                    DTOutput("weights_table")
                )
              )
      ),
      tabItem(tabName = "options",
              fluidRow(
                box(width = 3, title = "Option Chain Controls", solidHeader = TRUE,
                    textInput("opt_ticker", "Ticker", "AAPL"),
                    actionButton("fetch_opt", "Fetch Option Chain", class = "btn btn-warning")
                ),
                box(width = 9, title = "Option Chain Data", solidHeader = TRUE,
                    tabsetPanel(
                      # FIX: output IDs must be unique in the UI
                      tabPanel("Calls", DTOutput("calls_table")),
                      tabPanel("Puts", DTOutput("puts_table"))
                    )
                )
              )
      ),
      tabItem(tabName = "options_bsm",
              fluidRow(
                box(width = 3, title = "Option Chain Inputs", solidHeader = TRUE,
                    textInput("opt_ticker_bsm", "Ticker", "AAPL"),
                    numericInput("r", "Risk-free rate (annual, %)", 5, min = 0, step = 0.1),
                    numericInput("days", "Days to expiry", 30, min = 1),
                    numericInput("sigma_guess", "Fallback volatility (%, used if IV missing)", 20, min = 1, max = 300),
                    actionButton("fetch_opt_bsm", "Fetch & Compute", class = "btn btn-warning")
                ),
                box(width = 9, title = "Option Chain + BSM Comparison", solidHeader = TRUE,
                    tabsetPanel(
                      tabPanel("Calls", DTOutput("calls_table_bsm")),
                      tabPanel("Puts", DTOutput("puts_table_bsm")),
                      tabPanel("IV vs Strike", plotOutput("iv_plot", height = "500px"))
                    )
                )
              )
      ),
      tabItem(tabName = "about",
              fluidRow(
                box(width = 12, title = "About This App", solidHeader = TRUE,
                    HTML('<h3 style="color:#00ff99;">OHLCV + Portfolio + Options Dashboard</h3>
<p style="color:#0080ff;">This dashboard combines stock technical analysis, Markowitz portfolio optimization, and option chain analysis with theoretical Black–Scholes pricing.</p>
<ul style="color:#ff2a00;">
<li>View and analyze historical OHLCV data with indicators like SMA, EMA, RSI, MACD, and Bollinger Bands.</li>
<li>Construct and visualize efficient frontiers using mean-variance optimization.</li>
<li>Compare live option chain data with Black–Scholes theoretical prices and implied volatility plots.</li>
</ul>
<p style="color:#f5d76e;">Developed and designed by <strong style="color:#FFD700;"> MTH208 Course Group Project</strong></p>
<p style="color:#f5d76e;">Powered by R, Shiny, Quantmod, and yfR.</p>')
                )
              )
      )
      
    )
  )
)
server <- function(input, output, session) {
  bsm_price <- function(S, K, T, r, sigma, type = c("call", "put")) {
    type <- match.arg(type)
    d1 <- (log(S/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    if (type == "call") {
      return(S * pnorm(d1) - K * exp(-r*T) * pnorm(d2))
    } else {
      return(K * exp(-r*T) * pnorm(-d2) - S * pnorm(-d1))
    }
  }
  clean_iv <- function(x) {
    
    if (is.null(x)) return(NULL)
    if (is.numeric(x)) return(x)
    suppressWarnings({
      y <- gsub("%", "", as.character(x))
      y <- as.numeric(y)
    })
    y / 100
  }
  
  first_nonempty_chain <- function(ch) {
    if (!is.null(ch$calls) && !is.null(ch$puts)) return(ch)
    if (is.list(ch)) {
      for (i in seq_along(ch)) {
        if (!is.null(ch[[i]]$calls) && !is.null(ch[[i]]$puts)) return(ch[[i]])
      }
    }
    NULL
  }
  fetch_stock <- eventReactive(input$go, {
    req(input$ticker, input$dates)
    df <- try(yf_get(input$ticker,
                     first_date = input$dates[1],
                     last_date  = input$dates[2],
                     do_cache   = FALSE),
              silent = TRUE)
    validate(need(!inherits(df, "try-error") && nrow(df) > 0, "Data unavailable"))
    xt <- xts(df[, c("price_open", "price_high", "price_low",
                     "price_close", "volume")],
              order.by = df$ref_date)
    colnames(xt) <- c("Open", "High", "Low", "Close", "Volume")
    xt <- na.locf(xt); xt <- na.locf(xt, fromLast = TRUE)
    xt
  })
  
  output$chart <- renderPlot({
    x <- fetch_stock(); req(NROW(x) > 0)
    custom_theme <- chart_theme()
    custom_theme$col$bg     <- "#0d0d0d"
    custom_theme$col$grid   <- "#2b2b2b"
    custom_theme$col$up.col <- "#32CD32"
    custom_theme$col$dn.col <- "#FF4500"
    custom_theme$col$border <- "#FFD700"
    
    chart_Series(x, name = input$ticker, type = "candlesticks", theme = custom_theme)
    add_Vo()
    
    ind <- input$indicator
    if (!is.null(ind) && ind != "None") {
      if (ind == "SMA") {
        add_TA(SMA(Cl(x), n = input$n), on = 1, col = "#FFD700")
      } else if (ind == "EMA") {
        add_TA(EMA(Cl(x), n = input$n), on = 1, col = "#FFA500")
      } else if (ind == "Bollinger Bands") {
        bb <- BBands(HLC(x), n = input$n)
        add_TA(bb$mavg, on = 1, col = "#FFD700")
        add_TA(bb$up,   on = 1, col = "#DC143C")
        add_TA(bb$dn,   on = 1, col = "#DC143C")
      } else if (ind == "RSI") {
        add_TA(RSI(Cl(x), n = input$n), on = NA, col = "#1E90FF")
      } else if (ind == "MACD") {
        mc <- MACD(Cl(x))
        add_TA(mc$macd,   on = NA, col = "#FFD700")
        add_TA(mc$signal, on = NA, col = "#FF6347")
      }
    }
  })
  
  output$table <- renderDT({
    x <- fetch_stock(); req(x)
    df <- data.frame(Date = index(x), coredata(x))
    datatable(df, options = list(pageLength = 10))
  })
  
  port_data <- eventReactive(input$load_multi, {
    tickers <- unique(unlist(strsplit(gsub(" ", "", input$tickers_multi), ",")))
    rets <- list()
    for (t in tickers) {
      df <- try(yf_get(t,
                       first_date = input$dates_multi[1],
                       last_date  = input$dates_multi[2],
                       do_cache   = FALSE),
                silent = TRUE)
      if (inherits(df, "try-error") || nrow(df) == 0) next
      xt <- xts(df$price_close, order.by = df$ref_date)
      xt <- na.locf(xt); xt <- na.locf(xt, fromLast = TRUE)
      rets[[t]] <- dailyReturn(xt)
    }
    validate(need(length(rets) >= 2, "Need at least two valid tickers"))
    merge_xts <- na.omit(do.call(merge, rets))
    colnames(merge_xts) <- names(rets)
    merge_xts
  })
  
  output$frontier <- renderPlot({
    R <- port_data(); req(R)
    mu <- colMeans(R); covmat <- cov(R)
    N <- ncol(R)
    target_seq <- seq(min(mu), max(mu), length.out = 50)
    port_sd <- numeric(length(target_seq))
    port_mu <- numeric(length(target_seq))
    for (i in seq_along(target_seq)) {
      target <- target_seq[i]
      Dmat <- 2 * covmat
      dvec <- rep(0, N)
      Amat <- cbind(rep(1, N), mu)
      bvec <- c(1, target)
      sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 2), silent = TRUE)
      if (inherits(sol, "try-error")) next
      w <- sol$solution
      port_sd[i] <- sqrt(t(w) %% covmat %% w)
      port_mu[i] <- sum(mu * w)
    }
    df <- na.omit(data.frame(sd = port_sd, mu = port_mu))
    ggplot(df, aes(x = sd, y = mu)) +
      geom_line(linewidth = 1.3) +
      geom_point(size = 2) +
      labs(title = "Markowitz Efficient Frontier",
           x = "Risk (Std Dev)", y = "Expected Return") +
      
      theme_minimal(base_size = 14) +
      theme(
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text        = element_text(color = "black"),
        axis.title       = element_text(color = "black"),
        plot.title       = element_text(color = "#004080", face = "bold", size = 16),
        panel.grid.major = element_line(color = "#cccccc"),
        panel.grid.minor = element_line(color = "#e6e6e6")
      )
    
  })
  
  output$weights_table <- renderDT({
    R <- port_data(); req(R)
    mu <- colMeans(R); covmat <- cov(R); N <- ncol(R)
    Dmat <- 2 * covmat; dvec <- rep(0, N)
    Amat <- cbind(rep(1, N)); bvec <- 1
    sol <- try(solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
    if (inherits(sol, "try-error")) {
      return(datatable(data.frame(Message = "Optimization failed")))
    }
    wt <- sol$solution / sum(sol$solution)
    datatable(data.frame(Ticker = colnames(R), Weight = round(wt, 3)))
  })
  opt_chain <- eventReactive(input$fetch_opt, {
    req(input$opt_ticker)
    tryCatch({
      ch <- getOptionChain(Symbols = input$opt_ticker)
      ch <- first_nonempty_chain(ch)
      validate(need(!is.null(ch), "No option data available"))
      ch
    }, error = function(e) {
      showNotification(paste("Error fetching data:", e$message), type = "error")
      NULL
    })
  })
  
  output$calls_table <- renderDT({
    d <- opt_chain(); req(d)
    datatable(d$calls, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$puts_table <- renderDT({
    d <- opt_chain(); req(d)
    datatable(d$puts, options = list(pageLength = 10, scrollX = TRUE))
  })
  opt_chain_bsm <- eventReactive(input$fetch_opt_bsm, {
    req(input$opt_ticker_bsm)
    tryCatch({
      ch <- getOptionChain(Symbols = input$opt_ticker_bsm)
      ch <- first_nonempty_chain(ch)
      validate(need(!is.null(ch), "No option data available"))
      S0 <- NA_real_
      q <- try(getQuote(input$opt_ticker_bsm), silent = TRUE)
      if (!inherits(q, "try-error") && nrow(q) > 0) {
        if ("Last" %in% names(q)) S0 <- as.numeric(q$Last)
        if (is.na(S0) && "Last.Trade.Price.Only" %in% names(q)) S0 <- as.numeric(q$Last.Trade.Price.Only)
      }
      if (is.na(S0)) {
        px <- try(yf_get(input$opt_ticker_bsm, first_date = Sys.Date() - 5, last_date = Sys.Date(), do_cache = FALSE), silent = TRUE)
        if (!inherits(px, "try-error") && nrow(px) > 0) {
          S0 <- as.numeric(tail(px$price_close, 1))
        }
      }
      validate(need(!is.na(S0), "Unable to fetch underlying price"))
      
      Texp <- input$days / 365
      r    <- input$r / 100
      sig_fallback <- input$sigma_guess / 100
      
      calls <- ch$calls
      puts  <- ch$puts
      iv_calls <- clean_iv(calls$IV)
      if (is.null(iv_calls)) iv_calls <- clean_iv(calls$ImpVol)
      if (is.null(iv_calls)) iv_calls <- rep(NA_real_, nrow(calls))
      
      iv_puts <- clean_iv(puts$IV)
      if (is.null(iv_puts)) iv_puts <- clean_iv(puts$ImpVol)
      if (is.null(iv_puts)) iv_puts <- rep(NA_real_, nrow(puts))
      
      # Compute theoretical prices
      calls$IV_num <- ifelse(is.na(iv_calls), sig_fallback, iv_calls)
      puts$IV_num  <- ifelse(is.na(iv_puts),  sig_fallback, iv_puts)
      
      calls$BSM_Theoretical <- mapply(function(K, v) bsm_price(S0, K, Texp, r, v, "call"), calls$Strike, calls$IV_num)
      puts$BSM_Theoretical  <- mapply(function(K, v) bsm_price(S0, K, Texp, r, v, "put"),  puts$Strike,  puts$IV_num)
      
      list(calls = calls, puts = puts, S0 = S0)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      NULL
    })
  })
  
  output$calls_table_bsm <- renderDT({
    d <- opt_chain_bsm(); req(d)
    cols <- intersect(c("Strike", "Last", "Bid", "Ask", "IV", "ImpVol", "IV_num", "BSM_Theoretical"), names(d$calls))
    datatable(d$calls[, cols, drop = FALSE], options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$puts_table_bsm <- renderDT({
    d <- opt_chain_bsm(); req(d)
    cols <- intersect(c("Strike", "Last", "Bid", "Ask", "IV", "ImpVol", "IV_num", "BSM_Theoretical"), names(d$puts))
    datatable(d$puts[, cols, drop = FALSE], options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$iv_plot <- renderPlot({
    d <- opt_chain_bsm(); req(d)
    if (!"IV_num" %in% names(d$calls)) return(NULL)
    ggplot(d$calls, aes(x = Strike, y = IV_num)) +
      geom_line(linewidth = 1.2) +
      geom_point() +
      labs(title = paste0("Implied Volatility vs Strike — ", input$opt_ticker_bsm),
           x = "Strike Price", y = "Implied Volatility (decimal)") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text        = element_text(color = "black"),
        axis.title       = element_text(color = "black"),
        plot.title       = element_text(color = "#004080", face = "bold", size = 16),
        panel.grid.major = element_line(color = "#cccccc"),
        panel.grid.minor = element_line(color = "#e6e6e6")
      )
  })
}
shinyApp(ui, server)