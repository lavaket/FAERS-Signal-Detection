
## Faers Drugs Safety Signal Detection

packages <- c("dplyr", "readr", "stringr", "purrr", "tidyr", "ggplot2")
to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(to_install)) install.packages(to_install)


install.packages("knitr")
install.packages("xfun")

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

drug_file <- "DRUG25Q3.txt"
reac_file <- "REAC25Q3.txt"


DRUG_PATTERN  <- "OZEMPIC"
EVENT_PATTERN <- "PANCREATITIS"


# Directory
DATA_DIR <- "data_raw"


file.exists("data_raw/DRUG25Q3.txt")
file.exists("data_raw/REAC25Q3.txt")


getwd()
 
load_faers <- function(filename) {
  read_delim(
    file = paste0(DATA_DIR, "/", filename),
    delim = "$",
    trim_ws = TRUE,
    col_types = cols(.default = col_character())
  )
}


drug_raw <- load_faers(drug_file)
reac_raw <- load_faers(reac_file)

cat("Files Loaded:\n")
print(head(drug_raw))
print(head(reac_raw))

names(drug_raw)


drug_clean <- drug_raw %>%
  mutate(
    drugname = str_to_upper(str_squish(drugname)),
    role_cod = str_to_upper(str_squish(role_cod))
  )

reac_clean <- reac_raw %>%
  mutate(
    pt = str_to_upper(str_squish(pt))
  )

cat("Cleaning complete.\n")


# Collapse to case level lists of drugs + reactions


drugs_by_case <- drug_clean %>%
  group_by(primaryid) %>%
  summarise(drugs = list(unique(drugname)), .groups = "drop")

reactions_by_case <- reac_clean %>%
  group_by(primaryid) %>%
  summarise(reactions = list(unique(pt)), .groups = "drop")

case_level <- drugs_by_case %>%
  left_join(reactions_by_case, by = "primaryid")

cat("Case-level dataset created.\n")


# Flag presence of selected drug + event


case_flagged <- case_level %>%
  mutate(
    has_drug  = map_lgl(drugs,     ~ any(str_detect(.x, DRUG_PATTERN))),
    has_event = map_lgl(reactions, ~ any(str_detect(.x, EVENT_PATTERN)))
  )

cat("Drug and event flagged.\n")



# Build the 2x2 table


# Raw summarized counts
two_by_two <- case_flagged %>%
  summarise(
    a = sum(has_drug & has_event, na.rm = TRUE),
    b = sum(has_drug & !has_event, na.rm = TRUE),
    c = sum(!has_drug & has_event, na.rm = TRUE),
    d = sum(!has_drug & !has_event, na.rm = TRUE)
  )

cat("Raw 1×4 tibble:\n")
print(two_by_two)

# 2×2 matrix 
two_by_two_matrix <- matrix(
  c(two_by_two$a, two_by_two$b, two_by_two$c, two_by_two$d),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Exposure = c("Drug Present", "Drug Absent"),
    Outcome  = c("Event Present", "Event Absent")
  )
)

cat("\nFormatted 2×2 Matrix:\n")
print(two_by_two_matrix)

# tibble version for readability or exporting for more visual appeal
two_by_two_pretty <- tibble::tibble(
  Exposure       = c("Drug Present", "Drug Absent"),
  Event_Present  = c(two_by_two$a, two_by_two$c),
  Event_Absent   = c(two_by_two$b, two_by_two$d)
)

cat("\nPretty 2×2 Table:\n")
print(two_by_two_pretty)

md_tbl <- suppressWarnings(
  knitr::kable(
    two_by_two_pretty,
    format = "markdown",
    caption = paste("2×2 Table:", DRUG_PATTERN, "vs", EVENT_PATTERN)
  )
)
writeLines(md_tbl, "results/2x2_table.md")



# ROR + PRR functions


ror_ci <- function(a, b, c, d) {
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  ror <- (a * d) / (b * c)
  se  <- sqrt(1/a + 1/b + 1/c + 1/d)
  tibble(
    measure  = "ROR",
    estimate = ror,
    lower    = exp(log(ror) - 1.96 * se),
    upper    = exp(log(ror) + 1.96 * se)
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
    measure  = "PRR",
    estimate = prr,
    lower    = exp(log(prr) - 1.96 * se),
    upper    = exp(log(prr) + 1.96 * se)
  )
}



# ROR + PRR
 
a <- as.numeric(two_by_two$a)
b <- as.numeric(two_by_two$b)
c <- as.numeric(two_by_two$c)
d <- as.numeric(two_by_two$d)

signal_results <- bind_rows(
  ror_ci(a, b, c, d),
  prr_ci(a, b, c, d)
)

cat("\nSignal Detection Results:\n")
print(signal_results)


# Outputs


if (!dir.exists("results")) dir.create("results")

# export pretty table as CSV
write.csv(
  two_by_two_pretty,
  file = "results/2x2_pretty_table.csv",
  row.names = FALSE
)

# export raw numeric matrix as CSV
write.csv(
  two_by_two_matrix,
  file = "results/2x2_matrix.csv",
  row.names = TRUE
)


# 2×2 cell counts


counts_long <- two_by_two %>%
  pivot_longer(cols = everything(), names_to = "cell", values_to = "count")

p <- ggplot(counts_long, aes(x = cell, y = count)) +
  geom_col(fill = "steelblue") +
  labs(
    title = paste("2x2 Cell Counts:", DRUG_PATTERN, "and", EVENT_PATTERN),
    x = "Cell (a, b, c, d)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/2x2_counts_plot.png", p, width = 6, height = 4, dpi = 300)

cat("\nPlot saved to 'figures/2x2_counts_plot.png'\n")
cat("\n==== ALL COMPLETE ====\n")
