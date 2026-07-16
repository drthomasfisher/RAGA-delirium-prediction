# ---- 1. Libraries & Setup ----
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(readxl)
  library(writexl)
  library(dplyr)
})

# Key 
api_key <- Sys.getenv("MS_TRANSLATOR_KEY")

`%||%` <- function(x, y) if (!is.null(x) && !is.na(x)) x else y

# ---- 2. Load Data & Dictionaries ----

# Main Data - Added .name_repair to ensure no NA or empty names
raga_data_chinese <- read_excel("data_raw/raga_data_chinese.xlsx", 
                                col_names = TRUE, 
                                .name_repair = "unique") # "unique" fixes "" and NA names

# Codebook (Primary Priority)
code_book <- read_excel("docs/code_book_custom.xlsx", col_names = TRUE)
names(code_book)[1:2] <- c("Chinese", "English")
code_book$Source <- "codebook"

# Headers
read_header_dict <- function(path) {
  if(!file.exists(path)) return(data.frame(Chinese=character(), English=character(), Source=character()))
  df <- read.csv(path, header = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  if (nrow(df) == 2) {
    headers_t <- as.data.frame(t(df), stringsAsFactors = FALSE)
    names(headers_t) <- c("Chinese", "English")
    headers_t$Source <- "headers"
    return(headers_t)
  }
  names(df)[1:2] <- c("Chinese", "English")
  df$Source <- "headers"
  df
}
header_dict <- read_header_dict("docs/headers_translated.csv") %>%
  filter(!is.na(Chinese) & Chinese != "" & !is.na(English) & English != "")

# Cache 
cache_file <- "data_interim/translation_cache.csv"
cached <- if (file.exists(cache_file)) {
  read.csv(cache_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
} else {
  data.frame(Chinese = character(), English = character(), Source = character(), stringsAsFactors = FALSE)
}

# ---- 3. Translator ----
ms_translate_batch <- function(texts, key, to = "en", from = "zh", region = "global", 
                               endpoint = "https://api.cognitive.microsofttranslator.com") {
  if (length(texts) == 0) return(character(0))
  
  out <- character(length(texts))
  batch_size <- 90 
  
  for (i in seq(1, length(texts), by = batch_size)) {
    idx <- i:min(i + batch_size - 1, length(texts))
    chunk <- texts[idx]
    body <- lapply(chunk, function(t) list(Text = t))
        
    res <- RETRY(
      verb = "POST",
      url = paste0(endpoint, "/translate"),
      add_headers(
        `Ocp-Apim-Subscription-Key` = key,
        `Ocp-Apim-Subscription-Region` = region,
        `Content-Type` = "application/json"
      ),
      query = list("api-version" = "3.0", "from" = from, "to" = to),
      body = toJSON(body, auto_unbox = TRUE),
      encode = "json",
      times = 3,      
      pause_base = 2  
    )
    
    if (status_code(res) != 200) {
      warning(paste("Batch starting at", i, "failed. Status:", status_code(res)))
      out[idx] <- NA_character_
      next
    }
    
    parsed <- fromJSON(content(res, "text", encoding = "UTF-8"))
    out[idx] <- vapply(parsed, function(item) item$translations[[1]]$text %||% NA_character_, character(1))
  }
  out
}

# ---- 4. Dictionary Building ----

# Identify Chinese content (Values + Column Headers)
all_values <- unlist(raga_data_chinese, use.names = FALSE)

unique_chinese <- unique(all_values[grepl("\\p{Han}", all_values, perl = TRUE)])
unique_chinese <- unique(c(unique_chinese, names(raga_data_chinese)))

# Current known translations
dict_prelim <- bind_rows(code_book, header_dict, cached) %>%
  filter(!duplicated(Chinese))

# Find new values
to_translate <- setdiff(unique_chinese, dict_prelim$Chinese)

if (length(to_translate) > 0) {
  cat("Translating", length(to_translate), "new values via Azure...\n")
  engs <- ms_translate_batch(to_translate, key = api_key)
  
  new_entries <- data.frame(
    Chinese = to_translate,
    English = engs,
    Source  = "microsoft",
    timestamp = Sys.time(), 
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(English) & English != "" & English != "-")
  
  cached <- bind_rows(cached, new_entries) %>% filter(!duplicated(Chinese))
write.csv(cached, cache_file, row.names = FALSE, fileEncoding = "UTF-8")
}

# Final  table
dict_table_final <- bind_rows(code_book, header_dict, cached) %>%
  filter(!duplicated(Chinese))
dict_lookup <- setNames(dict_table_final$English, dict_table_final$Chinese)

# ---- 5. Apply Translations ----

# 1: Fix column names of the raw data 
names(raga_data_chinese) <- make.names(names(raga_data_chinese), unique = TRUE)

# 2: Translate Data Values
raga_data_final <- raga_data_chinese %>%
  mutate(across(where(~ is.character(.) || is.factor(.)), function(col) {
    col <- as.character(col)
    recode_idx <- col %in% names(dict_lookup)
    col[recode_idx] <- dict_lookup[col[recode_idx]]
    col
  }))

# 3: Translate Column Headers
new_names <- names(raga_data_final)
match_idx <- match(new_names, names(dict_lookup))
names(raga_data_final) <- ifelse(!is.na(match_idx), dict_lookup[match_idx], new_names)

# ---- 6. Save & Cleanup ----
saveRDS(raga_data_final, "data_interim/raga_data_translated.rds")
write_xlsx(dict_table_final, "data_final/final_translation_dictionary.xlsx")

cat("Translation complete. RDS saved to data_interim.\n")
rm(list = ls()) 
gc()