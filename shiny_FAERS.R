install.packages("shiny")

library(shiny)
library(bslib)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(gt)
library(ggplot2)

DATA_DIR <- "data_raw"

detect_delim <- function(file) {
  first_line <- readLines(file, n = 1)
  if (grepl("\\|", first_line)) return("|")
  if (grepl("\\$", first_line)) return("$")
  if (grepl("\\!", first_line)) return("!")
  return("|")
}

drug_delim <- detect_delim(file.path(DATA_DIR, "DRUG25Q3.txt"))
reac_delim <- detect_delim(file.path(DATA_DIR, "REAC25Q3.txt"))

drug_raw <- read_delim(
  file.path(DATA_DIR, "DRUG25Q3.txt"),
  delim = drug_delim,
  escape_double = FALSE,
  quote = "",
  trim_ws = TRUE,
  col_types = cols(.default = col_character())
)

reac_raw <- read_delim(
  file.path(DATA_DIR, "REAC25Q3.txt"),
  delim = reac_delim,
  escape_double = FALSE,
  quote = "",
  trim_ws = TRUE,
  col_types = cols(.default = col_character())
)

drug_clean <- drug_raw %>%
  mutate(
    drugname = str_to_upper(str_squish(drugname)),
    role_cod = str_to_upper(str_squish(role_cod))
  )

reac_clean <- reac_raw %>%
  mutate(pt = str_to_upper(str_squish(pt)))

drugs_by_case <- drug_clean %>%
  group_by(primaryid) %>%
  summarise(drugs = list(unique(drugname)), .groups = "drop")

reactions_by_case <- reac_clean %>%
  group_by(primaryid) %>%
  summarise(reactions = list(unique(pt)), .groups = "drop")

case_base <- drugs_by_case %>%
  left_join(reactions_by_case, by = "primaryid")

ror_ci <- function(a, b, c, d) {
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  ror <- (a * d) / (b * c)
  se  <- sqrt(1/a + 1/b + 1/c + 1/d)
  tibble(
    Measure  = "ROR",
    Estimate = ror,
    Lower    = exp(log(ror) - 1.96 * se),
    Upper    = exp(log(ror) + 1.96 * se)
  )
}

prr_ci <- function(a, b, c, d) {
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  risk1 <- a / (a + b)
  risk0 <- c / (c + d)
  prr  <- risk1 / risk0
  se   <- sqrt((1/a - 1/(a+b)) + (1/c - 1/(c+d)))
  tibble(
    Measure  = "PRR",
    Estimate = prr,
    Lower    = exp(log(prr) - 1.96 * se),
    Upper    = exp(log(prr) + 1.96 * se)
  )
}

my_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Inter")
)

ui <- fluidPage(
  theme = my_theme,
  navbarPage("FAERS Signal Detection Dashboard",
             
             tabPanel("Signal Detection",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Select Inputs"),
                          selectInput("drug", "Drug:", choices = sort(unique(drug_clean$drugname))),
                          selectInput("event", "Event (PT):", choices = sort(unique(reac_clean$pt))),
                          actionButton("run", "Run Signal Detection", class = "btn btn-primary w-100"),
                          hr(),
                          downloadButton("download_results", "Download Results", class = "w-100"),
                          width = 3
                        ),
                        
                        mainPanel(
                          h3("2×2 Table"),
                          gt_output("matrix2x2"),
                          
                          h3("Disproportionality Metrics (ROR / PRR)"),
                          gt_output("signal_metrics"),
                          
                          h3("2×2 Cell Counts Plot"),
                          plotOutput("plot_counts", height = "300px"),
                          
                          width = 9
                        )
                      )
             ),
             
             tabPanel("Raw Output",
                      verbatimTextOutput("debug")
             ),
             
             tabPanel("About",
                      h4("FAERS Signal Detection System"),
                      p("This dashboard performs disproportionality analysis using FAERS data."),
                      p("Includes case aggregation, 2×2 table generation, and ROR/PRR metrics.")
             )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    DRUG_PATTERN  <- input$drug
    EVENT_PATTERN <- input$event
    
    flagged <- case_base %>%
      mutate(
        has_drug  = map_lgl(drugs,     ~ any(str_detect(.x, DRUG_PATTERN))),
        has_event = map_lgl(reactions, ~ any(str_detect(.x, EVENT_PATTERN)))
      )
    
    a <- sum(flagged$has_drug & flagged$has_event)
    b <- sum(flagged$has_drug & !flagged$has_event)
    c <- sum(!flagged$has_drug & flagged$has_event)
    d <- sum(!flagged$has_drug & !flagged$has_event)
    
    two_by_two <- tibble(
      Exposure = c("Drug Present", "Drug Absent"),
      `Event Present` = c(a, c),
      `Event Absent`  = c(b, d)
    )
    
    signals <- bind_rows(
      ror_ci(a, b, c, d),
      prr_ci(a, b, c, d)
    )
    
    list(table = two_by_two, signals = signals)
  })
  
  output$matrix2x2 <- render_gt({
    req(results())
    results()$table %>% gt() %>% tab_header(title = "2×2 Contingency Table")
  })
  
  output$signal_metrics <- render_gt({
    req(results())
    results()$signals %>% gt() %>%
      fmt_number(columns = c(Estimate, Lower, Upper), decimals = 3)
  })
  
  output$plot_counts <- renderPlot({
    req(results())
    df <- results()$table %>%
      pivot_longer(cols = -Exposure, names_to = "Outcome", values_to = "Count")
    
    ggplot(df, aes(x = Outcome, y = Count, fill = Exposure)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#1F77B4", "#2CA02C")) +
      theme_minimal(base_size = 14)
  })
  
  output$debug <- renderPrint({
    results()
  })
  
  output$download_results <- downloadHandler(
    filename = function() "faers_signal_results.csv",
    content = function(file) {
      write.csv(results()$signals, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
