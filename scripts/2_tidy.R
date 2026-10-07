# Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(writexl)
  library(here)
})

# Import  ----

rds_path <- here::here("data_interim", "raga_data_translated.rds")

# Load and flatten master list ----
if (file.exists(rds_path)) {
  # Load  rds if exists, else run translation script 
  data_translated <- readRDS(rds_path)
} else {
  source(here::here("scripts", "1_translation.R"))
  data_translated <- readRDS(rds_path)
}

data_translated <- data_translated |> janitor::clean_names()
data_translated <- data_translated %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%
  mutate(across(where(is.factor),    ~ na_if(as.character(.x), "NA")))

# Define functions ----

#Standardise date formats
excel_date <- function(x) {
  x_chr <- trimws(factor(x))
  x_chr[x_chr == ""] <- NA
  
  x_chr <- gsub("--+", "-", x_chr)
  
  # try Excel serials
  num <- suppressWarnings(as.numeric(x_chr))
  out <- as.Date(num, origin = "1899-12-30")
  
  # otherwise parse as text date
  bad <- is.na(out) & !is.na(x_chr)
  if (any(bad)) {
    out[bad] <- suppressWarnings(lubridate::parse_date_time(
      x_chr[bad],
      orders = c("Y-m-d", "d/m/Y", "m/d/Y", "Y/m/d")
    ))
  }
  
  as.Date(out)
}
block <- function(df, expr) {
  vars <- all.vars(expr)                  
  used <- intersect(vars, names(df))    
  out  <- eval(expr)                      
  attr(out, "used_cols") <- used          
  out
}
#Standardise time to mins
parse_duration_to_minutes <- function(x) {
  x_clean <- str_to_lower(x) |> str_trim()
  x_clean <- str_replace_all(x_clean, "mim", "min") 
  
  # Extract hours and minutes
  hours   <- str_extract(x_clean, "\\d+(?=h| hour)")
  minutes <- str_extract(x_clean, "\\d+(?=m|min| minutes)")
  
  # Convert to numeric
  hours   <- as.numeric(hours)
  minutes <- as.numeric(minutes)
  
  # Handle decimals 
  hours <- ifelse(
    is.na(hours) & str_detect(x_clean, "\\d+\\.\\d+h"),
    as.numeric(str_remove(x_clean, "h")),
    hours
  )
  
  # Total minutes
  total <- coalesce(hours * 60, 0) + coalesce(minutes, 0)
  
  # If plain number, treat as minutes
  total <- ifelse(
    is.na(total) & str_detect(x_clean, "^\\d+$"),
    as.numeric(x_clean),
    total
  )
  
  return(total)
}
normalise_conc <- function(x) {
  # convert to numeric
  val <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(val) & val <= 0.2, val * 100, val)
}
parse_conc_min <- function(x) {
  x_clean <- str_to_lower(x) |> str_trim()
  x_clean <- str_replace_all(x_clean, "[~–]", "-")
  x_clean <- str_remove_all(x_clean, "%")
  lo <- as.numeric(str_extract(x_clean, "^[0-9]+\\.?[0-9]*"))
  normalise_conc(lo)
}
parse_conc_max <- function(x) {
  x_clean <- str_to_lower(x) |> str_trim()
  x_clean <- str_replace_all(x_clean, "[~–]", "-")
  x_clean <- str_remove_all(x_clean, "%")
  hi <- as.numeric(str_extract(x_clean, "(?<=-)[0-9]+\\.?[0-9]*"))
  normalise_conc(hi)
}
parse_conc_mid <- function(x) {
  lo <- parse_conc_min(x)
  hi <- parse_conc_max(x)
  ifelse(is.na(hi), lo, (lo + hi) / 2)
}
#Standardise units
normalize_unit <- function(u) {
  u <- tolower(as.character(u))
  u <- trimws(u)
  u <- gsub("μ|µ", "u", u)
  dplyr::case_when(
    u %in% c("mg", "milligram", "milligrams") ~ "mg",
    u %in% c("g", "gram", "grams") ~ "g",
    u %in% c("ug", "mcg", "microgram", "micrograms") ~ "ug",
    TRUE ~ NA_character_
  )
}
to_mg <- function(dose, unit) {
  unit <- normalize_unit(unit)
  case_when(
    unit == "mg" ~ dose,
    unit == "g"  ~ dose * 1000,
    unit == "ug" ~ dose * 0.001,
    TRUE ~ NA_real_
  )
}
parse_dose_mg <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  x_clean <- gsub("µ|μ|Î¼", "u", x_clean)
  x_clean <- gsub("\\s+", "", x_clean)
  
  num <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", x_clean)))
  
  dplyr::case_when(
    is.na(x) | x_clean == ""          ~ NA_character_,   
    grepl("ug|mcg|u", x_clean)        ~ as.character(num / 1000),
    grepl("mg", x_clean)              ~ as.character(num),
    grepl("g", x_clean)               ~ as.character(num * 1000),
    TRUE                              ~ paste0("Other (", x_clean, ")")
  )
}
#Y/N to boolean 
normalise_yes_no <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  
  case_when(
    is.na(x_clean) | x_clean == "" ~ NA,                # missing stays NA
    x_clean %in% c("yes", "y", "true", "1", "1. yes", "1.y") ~ TRUE,
    x_clean %in% c("no", "n", "false", "0", "2. no", "2.n")  ~ FALSE,
    TRUE ~ NA
  )
}
normalise_fracture_type <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_clean %in% c("1", "1.ggj", "ggj", "femoral neck fracture") ~ "Femoral neck fracture",
    x_clean %in% c("2", "2.zzj", "zzj", "intertrochanteric") ~ "Intertrochanteric",
    x_clean %in% c("3", "3.zzx", "zzx", "subtrochanteric") ~ "Subtrochanteric",
    x_clean %in% c("4", "4.ggt", "ggt", "fracture of femoral head") ~ "Fracture of femoral head",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("Femoral neck fracture", "Intertrochanteric", "Subtrochanteric", "Fracture of femoral head"))
}
normalise_vas <- function(x) {
  x_clean <- tolower(trimws(factor(x)))
  dplyr::case_when(
    x_clean %in% c("1", "1.0-10", "0-10") ~ "1",
    x_clean %in% c("2.10-30", "10-30")    ~ "2",
    x_clean %in% c("3.40-60", "40-60")    ~ "3",
    x_clean %in% c("4.70-100", "70-100")  ~ "4",
    TRUE ~ NA_character_
  ) |> 
    factor(levels = c("1","2","3","4"), ordered = TRUE)
}
#Convert BMI categories to ordinal 
normalise_bmi_level <- function(x) {
  x_clean <- tolower(trimws(factor(x)))
  dplyr::case_when(
    x_clean %in% c("1.≤18") ~ "≤18",
    x_clean %in% c("2.>18")                    ~ ">18",
    TRUE ~ NA_character_
  ) |>
    factor(levels = c("≤18", ">18"), ordered = TRUE)
}
delirium_type <- function(x) {
  recode(x,
         "N - Other" = "Other",
         "H-High activity" = "Hyperactive",
         "L-Low activity" = "Hypoactive",
         "M-hybrid" = "Mixed") |> factor()
}
normalise_drug <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  
  case_when(
    is.na(x_clean) | x_clean == "" ~ NA_character_,
    
    # Volatiles
    str_detect(x_clean, "desflurane|deflurane") ~ "Desflurane",
    str_detect(x_clean, "sevoflurane") ~ "Sevoflurane",
    str_detect(x_clean, "heptaamine") ~ "Other (heptaamine)",
    
    # NSAIDs
    str_detect(x_clean, "parecoxib|paroxib|parexib") ~ "Parecoxib",
    str_detect(x_clean, "celecoxib|xilebao|celebao|semicoxib|celecoxib|celebao|celecobena|celai xibu") ~ "Celecoxib",
    str_detect(x_clean, "flurbiprofen|forbiprofen|flurprophenolate|flurprolone|flubinoshenol|flurbiprofen pointer|fen ester") ~ "Flurbiprofen",
    str_detect(x_clean, "ketorolac|ketoloic") ~ "Ketorolac",
    str_detect(x_clean, "ketoprofen|ketolotic|ketorolalate") ~ "Ketoprofen",
    str_detect(x_clean, "diclofenac") ~ "Diclofenac",
    str_detect(x_clean, "indomethacin|indomethymbab|indomethasimbuba") ~ "Indomethacin",
    str_detect(x_clean, "loxoprofen|loxolofena") ~ "Loxoprofen",
    str_detect(x_clean, "aspirin") ~ "Aspirin",
    str_detect(x_clean, "ibuprofen") ~ "Ibuprofen",
    str_detect(x_clean, "meloxican") ~ "Meloxicam",
    str_detect(x_clean, "flurbiprofen") ~ "Flurbiprofen",
    str_detect(x_clean, "phenol hempamine|phenol hemp methyl|tylenol") ~ "Paracetamol",
    
    # Opioids
    str_detect(x_clean, "remifentanil") ~ "Remifentanil",
    str_detect(x_clean, "sufentanil|sufentanyl") ~ "Sufentanil",
    str_detect(x_clean, "fentanyl") ~ "Fentanyl",
    str_detect(x_clean, "tramad|aminophenol tramadol|Tramado") ~ "Tramadol",
    str_detect(x_clean, "buprenorphine") ~ "Buprenorphine",
    str_detect(x_clean, "butorphanol|butor agreed|butophanol|butofino|butonofi|noyang") ~ "Butorphanol",
    str_detect(x_clean, "morphine") ~ "Morphine",
    str_detect(x_clean, "oxycodone") ~ "Oxycodone",
    str_detect(x_clean, "hydromorphone") ~ "Hydromorphone",
    str_detect(x_clean, "pentazocin") ~ "Pentazocin",
    
    # Psychotropics
    str_detect(x_clean, "diazepam") ~ "Diazepam",
    str_detect(x_clean, "midazolam|liyuexi") ~ "Midazolam",
    str_detect(x_clean, "estazolam") ~ "Estazolam",
    str_detect(x_clean, "alprazolam") ~ "Alprazolam",
    str_detect(x_clean, "olanzapine") ~ "Olanzapine",
    str_detect(x_clean, "paroxetine") ~ "Paroxetine",
    str_detect(x_clean, "duloxetine|loxetine") ~ "Duloxetine",
    str_detect(x_clean, "sertraline|shetralin") ~ "Sertraline",
    str_detect(x_clean, "fluoxetine|prozac") ~ "Fluoxetine",
    str_detect(x_clean, "citalopram") ~ "Citalopram",
    str_detect(x_clean, "clozapine") ~ "Clozapine",
    str_detect(x_clean, "paliperidone") ~ "Paliperidone",
    str_detect(x_clean, "carbamazepine") ~ "Carbamazepine",
    str_detect(x_clean, "donepezil") ~ "Donepezil",
    str_detect(x_clean, "flunarizine") ~ "Flunarizine",
    str_detect(x_clean, "citicoline") ~ "Citicoline",
    str_detect(x_clean, "theophylline") ~ "Theophylline",
    str_detect(x_clean, "olaracetam") ~ "Olaracetam",
    str_detect(x_clean, "chlorpromazine") ~ "Chlorpromazine",
    str_detect(x_clean, "haloperidol") ~ "Haloperidol",
    str_detect(x_clean, "clonazepam") ~ "Clonazepam",
    str_detect(x_clean, "lorazepam") ~ "Lorazepam",
    str_detect(x_clean, "oxazepam") ~ "Oxazepam",
    str_detect(x_clean, "venlafaxine") ~ "Venlafaxine",
    str_detect(x_clean, "gabapentin") ~ "Gabapentin",
    str_detect(x_clean, "oxcarthypine") ~ "Oxcarbazepine",
    
    
    # Hypnotics
    str_detect(x_clean, "propofol") ~ "Propofol",
    str_detect(x_clean, "ketamine") ~ "Ketamine",
    str_detect(x_clean, "dexmedetomidine|dextromedetomidine") ~ "Dexmedetomidine",
    
    # Anti-emetics
    str_detect(x_clean, "dexamethasone") ~ "Dexamethasone",
    str_detect(x_clean, "ondansetron|ondan siqiong|ondan sjones") ~ "Ondansetron",
    str_detect(x_clean, "tropane|tropanesetron") ~ "Tropane derivative",
    
    # Local anaesthetics
    str_detect(x_clean, "ropivacaine") ~ "Ropivacaine",
    str_detect(x_clean, "bupivacaine|boo ratio") ~ "Bupivacaine",
    str_detect(x_clean, "levobup|levuviviviaine|levobupivaine|levuviviviaine|lev-bup|dexbupivacaine|l-bup|l-boo") ~ "Levobupivacaine",
    str_detect(x_clean, "lidocaine|kedocaine|kedocain|lidocane") ~ "Lidocaine",
    
    # Catch noisy unknowns
    str_detect(x_clean, "noyan|nixon|right beauty|tolane stone|ondan siqiong") ~ paste0("Other (", x_clean, ")"),
    
    
    # Vasopressors
    str_detect(x_clean, "epinephrine") ~ "Adrenaline",
    str_detect(x_clean, "norepinephrine") ~ "Noradrenaline",
    str_detect(x_clean, "ephedrine") ~ "Ephedrine",
    str_detect(x_clean, "dopamine") ~ "Dopamine",
    str_detect(x_clean, "phenylephrine|neofulin") ~ "Phenylephrine",
    str_detect(x_clean, "m-hydroxylamine|metahydroxylamine|bitarrate|tartrate") ~ "Metahydroxylamine",
    str_detect(x_clean, "piperidine") ~ "Piperidine",
    str_detect(x_clean, "dezosin") ~ "Dezocine",
    
    
    # Other 
    str_detect(x_clean, "atropine") ~ "Atropine",
    str_detect(x_clean, "zinc") ~ "Zinc",
    str_detect(x_clean, "du lengding") ~ "Du Lengding",
    str_detect(x_clean, "dizoxin") ~ "Digoxin",
    str_detect(x_clean, "oryzanosin") ~ "Oryzanosin",
    str_detect(x_clean, "chlornitrate") ~ "Chlornitrate",
    str_detect(x_clean, "betastine") ~ "Betahistine",
    str_detect(x_clean, "promethazine") ~ "Promethazine",
    
    # Fluids
    str_detect(x_clean, "saline|sodium chloride") ~ "Saline",
    
    # Catch-all
    TRUE ~ paste0("Other (", x_clean, ")")
  )
}
drug_class <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  
  case_when(
    is.na(x_clean) | x_clean == "" ~ NA_character_,
    
    # Volatiles
    x_clean %in% c("desflurane", "sevoflurane") ~ "Volatile",
    
    # NSAIDs
    x_clean %in% c("Aspirin", "Diclofenac", "Celecoxib", "Parecoxib", "Flurbiprofen", "Ketorolac", "Indomethacin", "Loxoprofen", "Ibuprofen", "Paracetamol", "coxi") ~ "NSAID",
    
    # Opioids
    x_clean %in% c("Tramadol", "Fentanyl", "Butorphanol", "Buprenorphine", "Morphine", "Oxycodone", "Hydromorphone", "Pentazocin", "dezosin") ~ "Opioid",
    
    
    # Psychotropics
    x_clean %in% c("olanzapine", "paroxetine", "duloxetine", "sertraline", "fluoxetine", "citalopram",
                   "clozapine", "paliperidone", "carbamazepine", "donepezil", "flunarizine",
                   "citicoline", "theophylline", "olaracetam", "chlorpromazine", "haloperidol", "Carbamazepine", "venlafaxine", "gabapentin", "oxcarbazepine") ~ "Psychotropic",
    
    # Benzodiazepines 
    x_clean %in% c("clonazepam", "lorazepam", "oxazepam", "midazolam", "diazepam", "estazolam", "alprazolam") ~ "Psychotropic",
    
    # Local anaesthetics
    x_clean %in% c("lidocaine", "bupivacaine", "levobupivacaine", "ropivacaine") ~ "Anaesthetic",
    
    # Intravenous hypnotic agents
    x_clean %in% c("propofol", "ketamine") ~ "IV hypnotic",
    
    # Vasoactive
    x_clean %in% c("epinephrine", "norepinephrine", "dopamine", "phenylephrine",
                   "metahydroxylamine", "piperidine") ~ "Vasoactive",
    
    # Other knowns
    x_clean %in% c("atropine") ~ "Anticholinergic",
    x_clean %in% c("zinc") ~ "Supplement",
    x_clean %in% c("du lengding", "digoxin", "oryzanosin", "chlornitrate", "betahistine", "promethazine") ~ "Other",
    
    # Fluids
    x_clean %in% c("saline") ~ "Fluid",
    
    # Catch-all
    TRUE ~ "Other"
  )
}

extract_dfs <- function(x) {
  if (is.data.frame(x)) {
    list(x)
  } else if (is.list(x)) {
    purrr::map(x, extract_dfs) |> unlist(recursive = FALSE)
  } else {
    list()
  }
}

# Demographics ----
demographic_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id              = as.character(individual_id),
            dob             = as.Date(birth_date),
            sex             = factor(sex),
            height_cm       = as.numeric(height),
            weight_kg       = as.numeric(weight),
            bmi             = as.numeric(bmi),
            bmi_level       = normalise_bmi_level(bmi_level),
            centre          = factor(centres),
            age_years       = as.numeric(age),
            
            # Age group (recoded from age_2) ---
            age_group = factor(
              dplyr::recode(age_2,
                            "1:â‰¤80" = "≤80",
                            "1:≤80"  = "≤80",   
                            "2:>80"  = ">80",
                            .default = NA_character_
              ),
              levels = c("≤80", ">80")
            ),
            
            education_level = case_when(
              education_level == "illiterate"        ~ 1,
              education_level == "Elementary school" ~ 2,
              education_level == "≥Secondary school" ~ 3
            ),
            work_status     = factor(work),
            
            # Administrative identifiers
            hospital_number   = as.character(hospitalization_number),
            address           = as.character(address),
            id_card           = as.character(id_card),
            contact           = as.character(contact),
            patient_name      = as.character(name),
            remark            = as.character(remark),
            physician         = factor(physician),
            anaesthetist = factor(anesthesiologists)
  )
))
}
# Background health ----
background_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Fracture descriptors ---
            fracture_type    = normalise_fracture_type(fracture_type),
            fracture_type_2  = normalise_fracture_type(fracture_type_2),
            frac_type        = normalise_fracture_type(frac_type),
            
            # ASA ---
            asa_status = american_society_of_anesthesiologists_asa_physical_status |>
            as.character() |>
            toupper() |>        
            trimws() |>         
            dplyr::recode(
              "I"   = "1",
              "II"  = "2",
              "III" = "3",
              "IV"  = "4",
              "V"   = "5"
            ) |> 
            factor(levels = c("1","2","3","4","5"), ordered = FALSE),

            # Baseline comorbidities by system ---
            cardiovascular   = normalise_yes_no(cardiovascular_system),
            respiratory      = normalise_yes_no(respiratory_system),
            gastrointestinal = normalise_yes_no(gastrointestinal_system),
            urinary_tract    = normalise_yes_no(urinary_tract_system),
            cns              = normalise_yes_no(central_nervous_system),
            hematologic      = normalise_yes_no(hematologic_system),
            other_comorb     = normalise_yes_no(other),
            
            comorbidity_1 = as.character(name_1),
            comorbidity_2 = as.character(name_2),
            comorbidity_3 = as.character(name_3),
            comorbidity_4 = as.character(name_4),
            comorbidity_5 = as.character(name_5),
            comorbidity_6 = as.character(name_6),
            comorbidity_7 = as.character(name_7),
  )
))
}
# Consent ----
consent_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Initial consent ---
            consent_signed          = normalise_yes_no(whether_the_informed_consent_form_is_signed),
            consent_method          = how_informed_consent_is_obtained |> recode(
              "The patient himself" = "Patient") |> factor(),
            consent_date            = as.Date(date_of_informed_consent),
            second_consent_required = normalise_yes_no(is_a_second_informed_consent_required),
            
            # Withdrawal ---
            withdrawn                 = normalise_yes_no(whether_to_withdraw_from_the_study),
            withdrawal_date           = excel_date(date_of_withdrawal_from_the_study),
            withdrawal_decision_maker = factor(decision_maker_for_withdrawal_from_the_study),
            relationship_patient      = factor(determine_the_relationship_between_people_and_patients),
            agree_use_existing_data   = normalise_yes_no(if_you_withdraw_from_the_study_do_you_agree_to_use_the_existing_data),
            
            # Follow-up consent/exit flags ---
            exit_experiment_1_x = normalise_yes_no(do_you_want_to_exit_the_experiment_1_x),
            exit_experiment_2   = normalise_yes_no(do_you_want_to_exit_experiment_2),
            exit_experiment_3   = normalise_yes_no(do_you_want_to_exit_experiment_3),
            exit_experiment_1_y = normalise_yes_no(do_you_want_to_exit_the_experiment_1_y)
  )
))
}
# Repeated biological variables ----
repeated_bio_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Baseline (t0) ---
            temp_t0       = as.numeric(preoperative_body_temperature),
            sysbp_t0      = as.numeric(preoperative_systolic_blood_pressure),
            diabp_t0      = as.numeric(preoperative_diastolic_blood_pressure),
            hr_t0         = as.numeric(preoperative_heart_rate),
            spo2_t0 = case_when(
                as.numeric(preoperative_blood_oxygen_saturation) > 100 ~ NA_real_,   # impossible, set NA
                as.numeric(preoperative_blood_oxygen_saturation) < 1   ~ as.numeric(preoperative_blood_oxygen_saturation) * 100,  # proportion → percent
                TRUE ~ as.numeric(preoperative_blood_oxygen_saturation)
                        ),
            hb_t0 = case_when(
  as.numeric(preoperative_hemoglobin) < 22  ~ as.numeric(preoperative_hemoglobin) * 100,
  as.numeric(preoperative_hemoglobin) > 220 ~ as.numeric(preoperative_hemoglobin) / 10,
  TRUE                              ~ as.numeric(preoperative_hemoglobin)
),
            albumin_t0    = as.numeric(preoperative_albumin),
            creatinine_t0 = as.numeric(preoperative_creatinine),
            ph_t0         = as.numeric(preoperative_ph),
            po2_kpa_t0        = 0.133322*as.numeric(preoperative_po2),
            pco2_kpa_t0       = 0.133322*as.numeric(preoperative_pco2),
            hco3_t0       = as.numeric(preoperative_hco3),
            cl_t0         = as.numeric(preoperative_cl),
            na_t0         = if_else(as.numeric(before_surgery) < 80,
                                    as.numeric(before_surgery) + 100,
                                    as.numeric(before_surgery)),
            k_t0          = as.numeric(before_surgery_k),
            alt_t0        = as.numeric(preoperative_alanine_aminotransferase),
            ast_t0        = as.numeric(preoperative_aspartate_aminotransferase),
            urea_t0       = as.numeric(preoperative_urea_nitrogen),
            delirium_t0       = normalise_yes_no(before_surgery_whether_delirium_occurs),
            delirium_type_t0  = delirium_type(preoperative_delirium_classification),
            delirium_sev_t0   = as.numeric(preoperative_delirium_severity_score),
            delirium_tot_t0   = as.numeric(preoperative_delirium_total_score),
            pain_score_t0     = as.numeric(preoperative_pain_score),
            
            
            # Timepoint 1 (t1) ---
            temp_t1       = as.numeric(no_1_body_temperature),
            sysbp_t1      = as.numeric(no_1_systolic_blood_pressure),
            diabp_t1      = as.numeric(no_1_diastolic_blood_pressure),
            hr_t1         = as.numeric(no_1_heart_rate),
            spo2_t1       = as.numeric(no_1_blood_oxygen_saturation),
            hb_t1         = if_else(as.numeric(no_1_hemoglobin) < 22,
                                    as.numeric(no_1_hemoglobin)*100,
                                    as.numeric(no_1_hemoglobin)),
            albumin_t1    = as.numeric(no_1_albumin),
            creatinine_t1 = as.numeric(no_1_creatinine),
            ph_t1         = as.numeric(no1_ph),
            po2_kpa_t1        = 0.133322*as.numeric(no1_po2),
            pco2_kpa_t1       = 0.133322*as.numeric(no1_pco2),
            hco3_t1       = as.numeric(no1_hco3),
            na_t1         = if_else(as.numeric(no1_na) < 80,
                                    as.numeric(no1_na) + 100,
                                    as.numeric(no1_na)),
            k_t1          = as.numeric(no1_k),
            cl_t1         = as.numeric(no1_cl),
            alt_t1        = as.numeric(no_1_alanine_aminotransferase),
            ast_t1        = as.numeric(no_1_aspartate_aminotransferase),
            urea_t1       = as.numeric(no_1_urea_nitrogen),
            delirium_t1       = normalise_yes_no(no_1_whether_delirium_occurs),
            delirium_type_t1  = delirium_type(no_1_delirium_classification),
            delirium_sev_t1   = as.numeric(no_1_delirium_severity),
            delirium_tot_t1   = as.numeric(no_1_total_delirium_score),
            pain_score_t1     = as.numeric(no_1_pain_score),
            
            
            # Timepoint 2 (t2) ---
            temp_t2       = as.numeric(no_2_body_temperature),
            sysbp_t2      = as.numeric(no_2_systolic_blood_pressure),
            diabp_t2      = as.numeric(no_2_diastolic_blood_pressure),
            hr_t2         = as.numeric(no_2_heart_rate),
            spo2_t2       = as.numeric(no_2_blood_oxygen_saturation),
            hb_t2         = if_else(as.numeric(no_2_hemoglobin) < 22,
                                    as.numeric(no_2_hemoglobin)*100,
                                    as.numeric(no_2_hemoglobin)),
            albumin_t2    = as.numeric(no_2_albumin),
            creatinine_t2 = as.numeric(no_2_creatinine),
            ph_t2         = as.numeric(no2_ph),
            po2_kpa_t2        = 0.133322*as.numeric(no2_po2),
            pco2_kpa_t2       = 0.133322*as.numeric(no2_pco2),
            hco3_t2       = as.numeric(no2_hco3),
            na_t2         = if_else(as.numeric(no2_na) < 80,
                                    as.numeric(no2_na) + 100,
                                    as.numeric(no2_na)),
            k_t2          = as.numeric(no2_k),
            cl_t2         = as.numeric(no2_cl),
            alt_t2        = as.numeric(no_2_alanine_aminotransferase),
            ast_t2        = as.numeric(no_2_aspartate_aminotransferase),
            urea_t2       = as.numeric(no_2_urea_nitrogen),
            delirium_t2       = normalise_yes_no(no_2_whether_delirium_occurs),
            delirium_type_t2  = delirium_type(no_2_delirium_classification),
            delirium_sev_t2   = as.numeric(no_2_delirium_severity),
            delirium_tot_t2   = as.numeric(no_2_total_delirium_score),
            pain_score_t2     = as.numeric(no_2_pain_score),
            
            
            # Timepoint 3 (t3) ---
            temp_t3       = as.numeric(no_3_body_temperature),
            sysbp_t3      = as.numeric(no_3_systolic_blood_pressure),
            diabp_t3      = as.numeric(no_3_diastolic_blood_pressure),
            hr_t3         = as.numeric(no_3_heart_rate),
            spo2_t3       = as.numeric(no_3_blood_oxygen_saturation),
            hb_t3         = if_else(as.numeric(no_3_hemoglobin) < 22,
                                    as.numeric(no_3_hemoglobin)*100,
                                    as.numeric(no_3_hemoglobin)),
            albumin_t3    = as.numeric(no_3_albumin),
            creatinine_t3 = as.numeric(no_3_creatinine),
            ph_t3         = as.numeric(no3_ph),
            po2_kpa_t3        = 0.133322*as.numeric(no3_po2),
            pco2_kpa_t3       = 0.133322*as.numeric(no3_pco2),
            hco3_t3       = as.numeric(no3_hco3),
            na_t3         = if_else(as.numeric(no3_na) < 80,
                                    as.numeric(no3_na) + 100,
                                    as.numeric(no3_na)),
            k_t3          = as.numeric(no3_k),
            cl_t3         = as.numeric(no3_cl),
            alt_t3        = as.numeric(no_3_alanine_aminotransferase),
            ast_t3        = as.numeric(no_3_aspartate_aminotransferase),
            urea_t3       = as.numeric(no_3_urea_nitrogen),
            delirium_t3       = normalise_yes_no(no_3_whether_delirium_occurs),
            delirium_type_t3  = delirium_type(no_3_delirium_classification),
            delirium_sev_t3   = as.numeric(no_3_delirium_severity),
            delirium_tot_t3   = as.numeric(no_3_total_delirium_score),
            pain_score_t3     = as.numeric(no_3_pain_score),
            
            
            # Timepoint 4 (t4) ---
            temp_t4       = as.numeric(no_4_body_temperature),
            sysbp_t4      = as.numeric(no_4_systolic_blood_pressure),
            diabp_t4      = as.numeric(no_4_diastolic_blood_pressure),
            hr_t4         = as.numeric(no_4_heart_rate),
            spo2_t4       = as.numeric(no_4_blood_oxygen_saturation),
            hb_t4         = if_else(as.numeric(no_4_hemoglobin) < 22,
                                    as.numeric(no_4_hemoglobin)*100,
                                    as.numeric(no_4_hemoglobin)),
            albumin_t4    = as.numeric(no_4_albumin),
            creatinine_t4 = as.numeric(no_4_creatinine),
            ph_t4         = as.numeric(no4_ph),
            po2_kpa_t4        = 0.133322*as.numeric(no_4_po2),
            pco2_kpa_t4       = 0.133322*as.numeric(no4_pco2),
            hco3_t4       = as.numeric(no_4_hco3),
            na_t4         = if_else(as.numeric(no4_na) < 80,
                                    as.numeric(no4_na) + 100,
                                    as.numeric(no4_na)),
            k_t4          = as.numeric(no4_k),
            cl_t4         = as.numeric(no4_cl),
            alt_t4        = as.numeric(no_4_alanine_aminotransferase),
            ast_t4        = as.numeric(no_4_aspartate_aminotransferase),
            urea_t4       = as.numeric(no_4_urea_nitrogen),
            delirium_t4       = normalise_yes_no(no_4_whether_delirium_occurs),
            delirium_type_t4  = delirium_type(no_4_delirium_classification),
            delirium_sev_t4   = as.numeric(no_4_delirium_severity_score),
            delirium_tot_t4   = as.numeric(no_4_total_delirium_score),
            pain_score_t4     = as.numeric(no_4_pain_score),
            
            # Timepoint 5 (t5) ---
            temp_t5       = as.numeric(no_5_body_temperature),
            sysbp_t5      = as.numeric(no_5_systolic_blood_pressure),
            diabp_t5      = as.numeric(no_5_diastolic_blood_pressure),
            hr_t5         = as.numeric(no_5_heart_rate),
            spo2_t5       = as.numeric(no_5_blood_oxygen_saturation),
            hb_t5         = if_else(as.numeric(no_5_hemoglobin) < 22,
                                    as.numeric(no_5_hemoglobin)*100,
                                    as.numeric(no_5_hemoglobin)),
            albumin_t5    = as.numeric(no_5_albumin),
            creatinine_t5 = as.numeric(no_5_creatinine),
            ph_t5         = as.numeric(no5_ph),
            po2_kpa_t5        = 0.133322*as.numeric(no_5_po2),
            pco2_kpa_t5       = 0.133322*as.numeric(no5_pco2),
            hco3_t5       = as.numeric(no5_hco3),
            na_t5         = if_else(as.numeric(no5_na) < 80,
                                    as.numeric(no5_na) + 100,
                                    as.numeric(no5_na)),
            k_t5          = as.numeric(no_5_k),
            cl_t5         = as.numeric(no5_cl),
            alt_t5        = as.numeric(no_5_alanine_aminotransferase),
            ast_t5        = as.numeric(no_5_aspartate_aminotransferase),
            urea_t5       = as.numeric(no_5_urea_nitrogen),
            delirium_t5       = normalise_yes_no(no_5_whether_delirium_occurs),
            delirium_type_t5  = delirium_type(no_5_delirium_classification),
            delirium_sev_t5   = as.numeric(no_5_delirium_severity),
            delirium_tot_t5   = as.numeric(no_5_total_delirium_score),
            pain_score_t5     = as.numeric(no_5_pain_score),
            
            # Timepoint 6 (t6) ---
            temp_t6       = as.numeric(no_6_body_temperature),
            sysbp_t6      = as.numeric(no_6_systolic_blood_pressure),
            diabp_t6      = as.numeric(no_6_diastolic_blood_pressure),
            hr_t6         = as.numeric(no_6_heart_rate),
            spo2_t6       = as.numeric(no_6_blood_oxygen_saturation),
            hb_t6         = if_else(as.numeric(no_6_hemoglobin) < 22,
                                    as.numeric(no_6_hemoglobin)*100,
                                    as.numeric(no_6_hemoglobin)),
            albumin_t6    = as.numeric(no_6_albumin),
            creatinine_t6 = as.numeric(no_6_creatinine),
            ph_t6         = as.numeric(no6_ph),
            po2_kpa_t6        = 0.133322*as.numeric(no_6_po2),
            pco2_kpa_t6       = 0.133322*as.numeric(no6_pco2),
            hco3_t6       = as.numeric(no_6_hco3),
            na_t6         = if_else(as.numeric(no6_na) < 80,
                                    as.numeric(no6_na) + 100,
                                    as.numeric(no6_na)),
            k_t6          = as.numeric(no_6_k),
            cl_t6         = as.numeric(no6_cl),
            alt_t6        = as.numeric(no_6_alanine_aminotransferase),
            ast_t6        = as.numeric(no_6_aspartate_aminotransferase),
            urea_t6       = as.numeric(no_6_urea_nitrogen),
            delirium_t6       = normalise_yes_no(no_6_whether_delirium_occurs),
            delirium_type_t6  = delirium_type(no_6_delirium_classification),
            delirium_sev_t6   = as.numeric(no_6_delirium_severity_score),
            delirium_tot_t6   = as.numeric(no_6_total_delirium_score),
            pain_score_t6     = as.numeric(no_6_pain_score),
            
            # Timepoint 7 (t7) ---
            temp_t7       = as.numeric(no_7_body_temperature),
            sysbp_t7      = as.numeric(no_7_systolic_blood_pressure),
            diabp_t7      = as.numeric(no_7_diastolic_blood_pressure),
            hr_t7         = as.numeric(no_7_heart_rate),
            spo2_t7       = as.numeric(no_7_blood_oxygen_saturation),
            hb_t7         = if_else(as.numeric(no_7_hemoglobin) < 22,
                                    as.numeric(no_7_hemoglobin)*100,
                                    as.numeric(no_7_hemoglobin)),
            albumin_t7    = as.numeric(no_7_albumin),
            creatinine_t7 = as.numeric(no_7_creatinine),
            ph_t7         = as.numeric(no7_ph),
            po2_kpa_t7        = 0.133322*as.numeric(no7_po2),
            pco2_kpa_t7       = 0.133322*as.numeric(no7_pco2),
            hco3_t7       = as.numeric(no7_hco3),
            na_t7         = if_else(as.numeric(no7_na) < 80,
                                    as.numeric(no7_na) + 100,
                                    as.numeric(no7_na)),
            k_t7          = as.numeric(no_7_k),
            cl_t7         = as.numeric(no7_cl),
            alt_t7        = as.numeric(no_7_alanine_aminotransferase),
            ast_t7        = as.numeric(no_7_aspartate_aminotransferase),
            urea_t7       = as.numeric(no_7_urea_nitrogen),
            delirium_t7       = normalise_yes_no(no_7_whether_delirium_occurs),
            delirium_type_t7  = delirium_type(no_7_delirium_classification),
            delirium_sev_t7   = as.numeric(no_7_delirium_severity_score),
            delirium_tot_t7   = as.numeric(no_7_total_delirium_score),
            pain_score_t7     = as.numeric(no_7_pain_score)
            
  )
))
}
# Pre‑operative & Longitudinal Features----
preop_features_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Baseline pre‑op characteristics ---
            preop_char1 = factor(preoperative_characteristics_1),
            preop_char2 = factor(preoperative_characteristics_2),
            preop_char3 = factor(preoperative_characteristics_3),
            preop_char4 = factor(preoperative_characteristics_4),
            
            # Alternate baseline sets (duplicates / repeats) ---
            preop_char1_2 = factor(preoperative_characteristics_1_2),
            preop_char2_2 = factor(preoperative_characteristics_2_2),
            preop_char3_2 = factor(preoperative_characteristics_3_2),
            preop_char4_2 = factor(preoperative_characteristics_4_2),
            
            # Numbered peri‑op features ---
            no1_feat1 = factor(no_1_feature_1),
            no1_feat2 = factor(no_1_feature_2),
            no1_feat3 = factor(no_1_feature_3),
            no1_feat4 = factor(no_1_feature_4),
            
            no2_feat1 = factor(no_2_feature_1),
            no2_feat2 = factor(no_2_feature_2),
            no2_feat3 = factor(no_2_feature_3),
            no2_feat4 = factor(no_2_feature_4),
            
            no3_feat1 = factor(no_3_feature_1),
            no3_feat2 = factor(no_3_feature_2),
            no3_feat3 = factor(no_3_feature_3),
            no3_feat4 = factor(no_3_feature_4),
            
            no4_feat1 = factor(no_4_feature_1),
            no4_feat2 = factor(no_4_feature_2),
            no4_feat3 = factor(no_4_feature_3),
            no4_feat4 = factor(no_4_feature_4),
            
            no5_feat1 = factor(no_5_feature_1),
            no5_feat2 = factor(no_5_feature_2),
            no5_feat3 = factor(no_5_feature_3),
            no5_feat4 = factor(no_5_feature_4),
            
            no6_feat1 = factor(no_6_feature_1),
            no6_feat2 = factor(no_6_feature_2),
            no6_feat3 = factor(no_6_feature_3),
            no6_feat4 = factor(no_6_feature_4),
            
            no7_feat1 = factor(no_7_feature_1),
            no7_feat2 = factor(no_7_feature_2),
            no7_feat3 = factor(no_7_feature_3),
            no7_feat4 = factor(no_7_feature_4),
            
            # Daily features (day_1 … day_7) ---
            day1_feat1 = factor(day_1_feature_1),
            day1_feat2 = factor(day_1_feature_2),
            day1_feat3 = factor(day_1_feature_3),
            day1_feat4 = factor(day_1_feature_4),
            
            day2_feat1 = factor(day_2_feature_1),
            day2_feat2 = factor(day_2_feature_2),
            day2_feat3 = factor(day_2_feature_3),
            day2_feat4 = factor(day_2_feature_4),
            
            day3_feat1 = factor(day_3_feature_1),
            day3_feat2 = factor(day_3_feature_2),
            day3_feat3 = factor(day_3_feature_3),
            day3_feat4 = factor(day_3_feature_4),
            
            day4_feat1 = factor(day_4_feature_1),
            day4_feat2 = factor(day_4_feature_2),
            day4_feat3 = factor(day_4_feature_3),
            day4_feat4 = factor(day_4_feature_4),
            
            day5_feat1 = factor(day_5_feature_1),
            day5_feat2 = factor(day_5_feature_2),
            day5_feat3 = factor(day_5_feature_3),
            day5_feat4 = factor(day_5_feature_4),
            
            day6_feat1 = factor(day_6_feature_1),
            day6_feat2 = factor(day_6_feature_2),
            day6_feat3 = factor(day_6_feature_3),
            day6_feat4 = factor(day_6_feature_4),
            
            day7_feat1 = factor(day_7_feature_1),
            day7_feat2 = factor(day_7_feature_2),
            day7_feat3 = factor(day_7_feature_3),
            day7_feat4 = factor(day_7_feature_4)
  )
))
}
# Perioperative variables ----
periop_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Surgery details ---
            surgery_type      = factor(surgery_type),
            surgery_type_alt  = dplyr::recode(surgery_type_2,
                                              "1.closed" = "Closed",
                                              "2.open"   = "Open",
                                              "2.open surgery_type_alt" = "Open") |>
              factor(levels = c("Closed","Open")),
            open_surgery_name = as.character(open_surgery_name),
            closure_surgery   = as.character(name_of_closure_surgery),
            surgeon_level     = factor(surgeon_level),
            
            # Anaesthesia details ---
            anaesth_level      = factor(anesthesiologist_level),
            anaesthesia_method = factor(anesthesia_method),
            anaesthesia_method_alt = factor(anesthesia_method_x),
            anaesthesia_start  = as.POSIXct(anesthesia_start_time),
            anaesthesia_end    = as.POSIXct(end_of_anesthesia),
            duration_hrs = as.numeric(mztime_int),
            
            # Blood loss / transfusion ---
            estimated_blood_loss   = as.numeric(estimated_blood_loss),
            blood_transfusion      = normalise_yes_no(blood_transfusion_or_not),
            total_transfusion      = as.numeric(total_blood_transfusion),
            autologous_transfusion = as.numeric(autologous_blood_transfusion),
            
            # Method detail (folded in) ---
            general_anaesthesia_1 = factor(general_anesthesia_1),
            general_anaesthesia_adjunct = case_when(
              general_anesthesia_2 == "GA" ~ "None",
              TRUE ~ str_remove(general_anesthesia_2, "^GA\\+")
            ) |> factor(),
            regional_method_1    = factor(local_anesthesia_method_1),
            regional_method_2    = factor(local_anesthesia_method_2),
            regional_method_3    = factor(local_anesthesia_method_3),
            spinal_level         = factor(puncture_gap),
            other_method         = as.character(other_please_specify),
            
            #Other
            anaesthesia_change = normalise_yes_no(whether_to_change_the_anesthesia_method),
            anaesthesia_change_reason = factor(reason_for_changing_anesthesia),
            anaesthesia_change_comment = as.character(brief_description)
  )
))
}
# Perioperative treatment ----
data_translated_fixed_volatiles <- {data_translated %>%
    rowwise() %>%
    mutate(
      inhalation_agent_1 = {
        vals <- na.omit(c(inhalation_anesthetics1, inhalation_anesthesia1))
        vals <- vals[vals != ""]
        if (length(unique(vals)) == 1) unique(vals) else NA_character_
      },
      inhalation_agent_2 = {
        vals <- na.omit(c(inhalation_anesthetics2, inhalation_anesthesia2))
        vals <- vals[vals != ""]
        if (length(unique(vals)) == 1) unique(vals) else NA_character_
      }
    ) %>%
    ungroup() %>%
    mutate(
      # Blank agent2 if identical to agent1
      inhalation_agent_2 = ifelse(
        !is.na(inhalation_agent_1) & inhalation_agent_1 == inhalation_agent_2,
        NA_character_,
        inhalation_agent_2
      ),
      # Normalise spelling/capitalisation of Desflurane
      inhalation_agent_1 = case_when(
        tolower(inhalation_agent_1) %in% c("deflurane", "desflurane") ~ "Desflurane",
        TRUE ~ inhalation_agent_1
      ),
      inhalation_agent_2 = case_when(
        tolower(inhalation_agent_2) %in% c("deflurane", "desflurane") ~ "Desflurane",
        TRUE ~ inhalation_agent_2
      )
    )
}
periop_treatment_vars <- {data_translated_fixed_volatiles  |> 
    transmute(
      id = as.character(individual_id),
      
      # Planned anaesthetic drugs 1–13 (tidy)
      anaesthetic_drug_tidy_1  = normalise_drug(anesthetic_drugs_1),
      anaesthetic_drug_tidy_2  = normalise_drug(anesthetic_drugs_2),
      anaesthetic_drug_tidy_3  = normalise_drug(anesthetic_drugs_3),
      anaesthetic_drug_tidy_4  = normalise_drug(anesthetic_drugs_4),
      anaesthetic_drug_tidy_5  = normalise_drug(anesthetic_drugs_5),
      anaesthetic_drug_tidy_6  = normalise_drug(anesthetic_drugs_6),
      anaesthetic_drug_tidy_7  = normalise_drug(anesthetic_drugs7),
      anaesthetic_drug_tidy_8  = normalise_drug(anesthetic_drugs8),
      anaesthetic_drug_tidy_9  = normalise_drug(anesthetic_drugs_9),
      anaesthetic_drug_tidy_10 = normalise_drug(anesthetic_drugs_10),
      anaesthetic_drug_tidy_11 = normalise_drug(anesthetic_drugs11),
      anaesthetic_drug_tidy_12 = normalise_drug(anesthetic_drugs_12),
      anaesthetic_drug_tidy_13 = normalise_drug(anesthetic_drugs_13),
      
      # Planned anaesthetic drugs 1–13 (raw)
      anaesthetic_drug_raw_1  = as.character(anesthetic_drugs_1),
      anaesthetic_drug_raw_2  = as.character(anesthetic_drugs_2),
      anaesthetic_drug_raw_3  = as.character(anesthetic_drugs_3),
      anaesthetic_drug_raw_4  = as.character(anesthetic_drugs_4),
      anaesthetic_drug_raw_5  = as.character(anesthetic_drugs_5),
      anaesthetic_drug_raw_6  = as.character(anesthetic_drugs_6),
      anaesthetic_drug_raw_7  = as.character(anesthetic_drugs7),
      anaesthetic_drug_raw_8  = as.character(anesthetic_drugs8),
      anaesthetic_drug_raw_9  = as.character(anesthetic_drugs_9),
      anaesthetic_drug_raw_10 = as.character(anesthetic_drugs_10),
      anaesthetic_drug_raw_11 = as.character(anesthetic_drugs11),
      anaesthetic_drug_raw_12 = as.character(anesthetic_drugs_12),
      anaesthetic_drug_raw_13 = as.character(anesthetic_drugs_13),
      
      #Doseage in mg
      
      anaesthetic_dose_mg_1  = to_mg(as.numeric(anesthetic_dose_1),  anesthetic_dosage_unit_1),
      anaesthetic_dose_mg_2  = to_mg(as.numeric(anesthetic_dose_2),  anesthetic_dosage_unit_2),
      anaesthetic_dose_mg_3  = to_mg(as.numeric(anesthetic_dose_3),  anesthetic_dosage_unit_3),
      anaesthetic_dose_mg_4  = to_mg(as.numeric(anesthetic_dose_4),  anesthetic_dosage_unit_4),
      anaesthetic_dose_mg_5  = to_mg(as.numeric(anesthetic_dose_5),  anesthesia_dosage_unit_5),
      anaesthetic_dose_mg_6  = to_mg(as.numeric(anesthetic_dose_6),  anesthesia_dosage_unit_6),
      anaesthetic_dose_mg_7  = to_mg(as.numeric(anesthetic_dose_7),  anesthesia_dosage_unit_7),
      anaesthetic_dose_mg_8  = to_mg(as.numeric(anesthetic_dose_8),  anesthesia_dosage_unit_8),
      anaesthetic_dose_mg_9  = to_mg(as.numeric(anesthetic_dose_9),  anesthesia_dosage_unit_9),
      anaesthetic_dose_mg_10 = to_mg(as.numeric(anesthetic_dose_10), anesthesia_dosage_unit_10),
      anaesthetic_dose_mg_11 = to_mg(as.numeric(anesthetic_dose_11), anesthesia_dosage_unit_11),
      anaesthetic_dose_mg_12 = to_mg(as.numeric(anesthetic_dose_12), anesthesia_dosage_unit_12),
      anaesthetic_dose_mg_13 = to_mg(as.numeric(anesthetic_dose_13), anesthesia_dosage_unit_13),
      
      # Administration methods
      method_1  = as.character(anesthesia_administration_method_1),
      method_2  = as.character(anesthesia_administration_method_2),
      method_3  = as.character(anesthesia_administration_method_3),
      method_4  = as.character(anesthesia_administration_method_4),
      method_5  = as.character(anesthesia_administration_method_5),
      method_6  = as.character(anesthesia_administration_method6),
      method_7  = as.character(anesthesia_administration_method7),
      method_8  = as.character(anesthesia_administration_method8),
      method_9  = as.character(anesthesia_administration_method9),
      method_10 = as.character(anesthesia_administration_method_10),
      method_11 = as.character(anesthesia_administration_method11),
      method_12 = as.character(anesthesia_administration_method12),
      method_13 = as.character(anesthesia_administration_method_13),
      
      # Notes
      note_1  = as.character(anesthesia_note_1),
      note_2  = as.character(anesthesia_note_2),
      note_3  = as.character(anesthesia_note_3),
      note_4  = as.character(anesthesia_note_4),
      note_5  = as.character(anesthesia_note_5),
      note_6  = as.character(anesthesia_note_6),
      note_7  = as.character(anesthesia_note_7),
      note_8  = as.character(anesthesia_note_8),
      note_11 = as.character(anesthesia_note_11),
      note_12 = as.character(anesthesia_note_12),
      
      # Primary inhalation agents (x-series) ---
      inhalation_1 = normalise_drug(inhalation_anesthetics1),
      inhalation_2 = normalise_drug(inhalation_anesthetics2),

      # Secondary inhalation agents (y-series, if present) ---
      inhalation_y1 = factor(inhalation_agent_1),
      inhalation_y2 = factor(inhalation_agent_2),
      
      # Volatile vs TIVA logic 

      has_volatile = (
      (!is.na(inhalation_1)  & inhalation_1 != "None"  & inhalation_1 != "") |
      (!is.na(inhalation_2)  & inhalation_2 != "None"  & inhalation_2 != "") |
      (!is.na(inhalation_y1) & inhalation_y1 != "None" & inhalation_y1 != "") |
      (!is.na(inhalation_y2) & inhalation_y2 != "None" & inhalation_y2 != "")
    ),
    # Anaesthesia type

      anaesthesia_phenotype = case_when(
        str_detect(anesthesia_method_x, "(?i)general anesthesia") & has_volatile  ~ "Volatile_GA",
        str_detect(anesthesia_method_x, "(?i)general anesthesia") & !has_volatile ~ "TIVA_GA",
        str_detect(anesthesia_method_x, "(?i)local anesthesia")                   ~ "Regional",
        TRUE ~ NA_character_
      ),

      # 2. Binary coding for LASSO
      anaes_volatile = ifelse(anaesthesia_phenotype == "Volatile_GA", 1, 0),
      anaes_tiva     = ifelse(anaesthesia_phenotype == "TIVA_GA", 1, 0),
      anaes_regional = ifelse(anaesthesia_phenotype == "Regional", 1, 0),
            
      
      conc1_min = parse_conc_min(maintain_concentration_range_1_x),
      conc1_max = parse_conc_max(maintain_concentration_range_1_x),
      conc1_mid = parse_conc_mid(maintain_concentration_range_1_x),
      conc1_duration = parse_duration_to_minutes(maintain_concentration_1),
      
      conc2_min = parse_conc_min(maintain_concentration_range_2_x),
      conc2_max = parse_conc_max(maintain_concentration_range_2_x),
      conc2_mid = parse_conc_mid(maintain_concentration_range_2_x),
      conc2_duration = parse_duration_to_minutes(maintain_concentration_2),
      
      conc1_y_min = parse_conc_min(maintain_concentration_range_1_y),
      conc1_y_max = parse_conc_max(maintain_concentration_range_1_y),
      conc1_y_mid = parse_conc_mid(maintain_concentration_range_1_y),
      
      conc2_y_min = parse_conc_min(maintain_concentration_range_2_y),
      conc2_y_max = parse_conc_max(maintain_concentration_range_2_y),
      conc2_y_mid = parse_conc_mid(maintain_concentration_range_2_y),
      
      
      # Maintenance times (numeric minutes)
      maintenance_time_1 = parse_duration_to_minutes(maintenance_time_1),
      maintenance_time_2 = parse_duration_to_minutes(maintenance_time_2),
      
      #Unplanned drugs
      intraop_hypotension   = normalise_yes_no(if_intraoperative_hypotension_happened),
      hypotension_medicated = normalise_yes_no(if_medicated_due_to_hypotension),
      
      # Drug treatments 
      drug_treatment_1       = normalise_drug(drug_treatment_1),
      drug_treatment_2       = normalise_drug(drug_treatment_2),
      drug_treatment_3       = normalise_drug(drug_treatment_3),
      drug_treatment_4       = normalise_drug(drug_treatment_4),
      drug_treatment_5       = normalise_drug(drug_treatment_5),
      drug_treatment_6       = normalise_drug(drug_treatment_6),
      drug_treatment_7       = normalise_drug(drug_treatment_7),
      
      # Standardised doses in mg
      dose_mg_1 = to_mg(as.numeric(treatment_dose_1), treatment_dose_unit_1),
      dose_mg_2 = to_mg(as.numeric(treatment_dose_2), treatment_dose_unit_2),
      dose_mg_3 = to_mg(as.numeric(treatment_dose_3), treatment_dose_unit_3),
      dose_mg_4 = to_mg(as.numeric(treatment_dose_4), treatment_dose_unit_4),
      dose_mg_5 = to_mg(as.numeric(treatment_dose_5), treatment_dose_unit_5),
      
      # Administration methods
      method1               = factor(treatment_and_administration_method_1),
      method2               = factor(treatment_and_administration_method_2),
      method3               = factor(treatment_and_administration_method_3),
      method4               = factor(treatment_and_administration_method_4),
      method5               = factor(treatment_and_administration_method_5),
      method6               = factor(treatment_and_administration_method_6),
      method7               = factor(treatment_and_administration_method_7),
      
      drug_treatment_1_note = factor(processing_notes_1),
      drug_treatment_2_note = factor(processing_note_2),
      drug_treatment_3_note = factor(processing_note_3),
      
      #Blood transfusion
      blood_transfusion_units = factor(total_blood_transfusion_unit),
      blood_transfusion_autologous_units  = factor(autologous_blood_transfusion_unit)
    )
  
}
# Misc perioperative drugs ----

periop_treatment_vars_misc <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Drug names tidy ---
            drug_name_tidy_1 = normalise_drug(drug_name_1_x),
            drug_name_tidy_2 = normalise_drug(drug_name_2_x),
            drug_name_tidy_3 = normalise_drug(drug_name_3_x),
            drug_name_tidy_4 = normalise_drug(drug_name_4_x),
            drug_name_tidy_5 = normalise_drug(drug_name_5_x),
            drug_name_tidy_6 = normalise_drug(drug_name_6_x),
            drug_name_tidy_7 = normalise_drug(drug_name_7_x),
            drug_name_tidy_8 = normalise_drug(drug_name_8_x),
            drug_name_tidy_9 = normalise_drug(drug_name_9_x),
            
            # Drug names raw ---
            drug_name_raw_1 = factor(drug_name_1_x),
            drug_name_raw_2 = factor(drug_name_2_x),
            drug_name_raw_3 = factor(drug_name_3_x),
            drug_name_raw_4 = factor(drug_name_4_x),
            drug_name_raw_5 = factor(drug_name_5_x),
            drug_name_raw_6 = factor(drug_name_6_x),
            drug_name_raw_7 = factor(drug_name_7_x),
            drug_name_raw_8 = factor(drug_name_8_x),
            drug_name_raw_9 = factor(drug_name_9_x),
            
            # Dosages / totals ---
            total_dosage_1 = as.numeric(total_medication_amount_1),
            total_dosage_2 = as.numeric(total_medication_amount_2),
            total_dosage_3    = as.numeric(total_dosage_3),
            total_dosage_4    = as.numeric(total_dosage_4),
            total_dosage_5    = as.numeric(total_dosage_5),
            total_dosage_6     = as.numeric(drug_dosage_6),
            total_dosage_7     = as.numeric(drug_dosage7),
            dosage_1          = factor(dosage_1),
            dosage_2          = factor(dosage_2),
            dosage_3          = factor(dosage_3),
            dosage_4          = factor(dosage_4),
            dosage_5          = factor(dosage_5),
            dosage_6          = factor(dosage_method_6),
            dosage_7          = factor(dosage_7),
            dosage_8          = factor(dosage_8),
            dosage_9          = factor(dosage_9),
            
            # Drug names  ---
            drug_name_y1 = normalise_drug(drug_name_1_y),
            drug_name_y2 = normalise_drug(drug_name_2_y),
            drug_name_y3 = normalise_drug(drug_name_3_y),
            drug_name_y4 = normalise_drug(drug_name_4_y),
            drug_name_y5 = normalise_drug(drug_name_5_y),
            drug_name_y6 = normalise_drug(drug_name_6_y),
            drug_name_y7 = normalise_drug(drug_name_7_y),
            
            # Usage / administration ---
            usage_mg_1 = parse_dose_mg(usage_1),
            usage_mg_2 = parse_dose_mg(usage_2),
            usage_mg_3 = parse_dose_mg(usage_3),
            usage_mg_4 = parse_dose_mg(usage_4),
            usage_mg_5 = parse_dose_mg(usage_5),
            usage_mg_6 = parse_dose_mg(usage_6),
            usage_mg_7 = parse_dose_mg(usage_7),
            
            # Pathways / routes ---
            path_1 = factor(path_1),
            path_2 = factor(path_2),
            path_3 = factor(path_3),
            path_4 = factor(path_4),
            path_5 = factor(path_5),
            path_6 = factor(pathway_6),
            path_7 = factor(path_7),
            
            # Frequency per day ---
            times_per_day1 = factor(x1_times_per_day),
            times_per_day2 = factor(x2_times_per_day),
            times_per_day3 = factor(x3_times_per_day),
            times_per_day4 = factor(x4_times_per_day),
            times_per_day5 = factor(x5_times_a_day),
            times_per_day6 = factor(x6_times_a_day),
            times_per_day7 = factor(x7_times_a_day),
            
            # Start / stop times ---
            start_time_1 = excel_date(start_time_1_2),
            start_time_2 = excel_date(start_time_2_2),
            start_time_3 = excel_date(start_time_3),
            start_time_4 = excel_date(start_time_4),
            start_time_5 = excel_date(start_time_5),
            start_time_6 = excel_date(start_time_6),
            start_time_7 = excel_date(start_time_7),
            
            stop_time_1  = excel_date(stop_time_1),
            stop_time_2  = excel_date(stop_time_2),
            stop_time_3  = excel_date(stop_time_3),
            stop_time_4  = excel_date(stop_time_4),
            stop_time_5  = excel_date(stop_time_5),
            stop_time_6  = excel_date(stop_time_6),
            stop_time_7  = excel_date(stop_time_7),
            
            # Cycle ---
            cycle = factor(cycle_y)
  )
))}

# Postoperative treatment (including analgesics, sedatives) ----
analgesic_sedative_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            
            # Analgesics 
            analgesic_1   = normalise_drug(analgesics_1),
            analgesic_2   = normalise_drug(analgesics_2),
            analgesic_3   = normalise_drug(analgesics_3),
            analgesic_4  = normalise_drug(analgesics_1_2),
            analgesic_5  = normalise_drug(painkillers_2),
            analgesic_6  = normalise_drug(painkillers_3),
            analgesic_7  = normalise_drug(painkillers_4),
            
            analgesic_dose_1_mg = parse_dose_mg(analgesic_dosage_1),
            analgesic_dose_2_mg = parse_dose_mg(analgesic_dosage_2),
            analgesic_dose_3_mg = parse_dose_mg(analgesic_dosage_3),
            analgesic_dose_4_mg = parse_dose_mg(analgesic_dosage_4),
            
            # PCA / pump settings
            total_analgesia_pump_mg = as.numeric(total_amount_of_analgesia_pump),
            basal_infusion_rate_mg  = as.numeric(basal_infusion_rate),
            pca                  = factor(pca),
            pca_alt              = factor(pca_1),
            
            # Sedatives
            sedative_1 = normalise_drug(sedative_1),
            sedative_2 = normalise_drug(sedative_2),
            sedative_3 = normalise_drug(sedative_3),
            
            sedative_dose_1_unknown_unit = factor(sedation_dose_1),
            sedative_dose_2_unknown_unit = factor(sedative_dose_2),
            sedative_dose_3_unknown_unit = factor(sedative_dose_3),
            
            # Post‑operative analgesia ---
            postop_analgesia   = postoperative_analgesia |> recode(
              "vein" = "Intravenous") |> factor(),
            analgesic_dose_1   = factor(analgesic_dose_1),
            analgesic_dose_2   = factor(analgesic_dose_2),
            analgesic_dose_3   = factor(analgesic_dose_3),
            
            # Analgesic / sedative classes ---
            opioid_analgesics        = normalise_yes_no(opioid_analgesics),
            nonsteroidal_analgesics  = normalise_yes_no(nonsteroidal_analgesics),
            benzodiazepine_sedatives = normalise_yes_no(benzodiazepine_sedatives),
            fluphenazine             = normalise_yes_no(fluphenazine),
            other_psychotropic_drugs = normalise_yes_no(other_psychotropic_drugs),
            
            # Sedation & pain scores ---
            sedation   = normalise_yes_no(sedation),
            vas_category_pre    = normalise_vas(vas_pre),
            vas_category_post   = normalise_vas(vas_post),
            vas_max_7d = as.numeric(vas_max_7d)
            
  )
))
}
# Adverse events ----
adverse_event_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Overall flag
            ae_any_flag = normalise_yes_no(are_there_any_adverse_events),
            
            # Event 1 ---
            ae_name_1     = factor(adverse_event_name_1_x),
            ae_start_1    = excel_date(start_time_x),
            ae_severity_1 = adverse_event_severity_1_x |> recode(
                             "light" = "Mild", "middle" = "Moderate", "heavy" = "Severe") |> factor(),
            ae_rel_anes_1 = factor(relationship_with_anesthesia_1_x),
            ae_outcome_1  = factor(outcome_of_adverse_events_1_x),
            ae_serious_1  = factor(is_it_a_serious_adverse_event, levels = c("No","Yes")),
            ae_mitig_date_1 = excel_date(mitigation_date_1_x),
            ae_remis_date_1 = excel_date(remission_date_1_y),
            
            # Event 2 ---
            ae_name_2     = factor(adverse_event_name_2),
            ae_start_2    = excel_date(start_time_1), 
            ae_severity_2 = adverse_event_severity_2 |> recode(
                             "light" = "Mild", "middle" = "Moderate", "heavy" = "Severe") |> factor(),
            ae_rel_anes_2 = factor(relationship_with_anesthesia2),
            ae_outcome_2  = factor(outcome_of_adverse_events2),
            ae_serious_2  = factor(is_it_a_serious_adverse_event_2, levels = c("No","Yes")),
            ae_mitig_date_2 = excel_date(remission_date_2),
            
            # Event 3 ---
            ae_name_3     = factor(adverse_event_name_3),
            ae_start_3    = excel_date(start_time_2),
            ae_severity_3 = adverse_event_severity_3 |> recode(
                             "light" = "Mild", "middle" = "Moderate", "heavy" = "Severe") |> factor(),
            ae_rel_anes_3 = factor(relationship_with_anesthesia_3),
            ae_outcome_3  = factor(outcome_of_adverse_events3),
            ae_serious_3  = factor(is_it_a_serious_adverse_event_3, levels = c("No","Yes")),
            ae_mitig_date_3 = excel_date(remission_date_3),
            
            # Event 4 ---
            ae_name_4     = factor(adverse_event_name_4),
            
            # “_y” variant (Wave 2) ---
            ae_any_1y     = normalise_yes_no(are_there_any_adverse_events_2), 
            ae_name_1y    = factor(adverse_event_name_1_y),
            ae_date_1y    = excel_date(start_time_y),
            ae_severity_1y = severity_of_adverse_events_1_y |> recode(
                             "light" = "Mild", "middle" = "Moderate", "heavy" = "Severe") |> factor(),
            ae_rel_anes_1y = factor(relationship_with_anesthesia_1_y),
            ae_outcome_1y  = factor(outcome_of_adverse_events1_y),
            
            # Extra names 
            ae_extra_name_1 = factor(a_ename),
            ae_extra_name_2 = factor(a_ename2),
            ae_extra_name_3 = factor(a_ename3)
  )
))}

# Cognitive outcomes ----
cognitive_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # 1. Pre‑operative Delirium ---
            cog_pre_delirium_exists  = normalise_yes_no(before_surgery_whether_delirium_occurs),
            cog_pre_delirium_type    = delirium_type(preoperative_delirium_classification),
            
            # 2. PRE-EXISTING DEMENTIA (baseline fields only) ---
            # follow-up dementia fields left out (they flagged 13 new post-operative diagnoses)
            cog_pre_dementia_flag = case_when(
              normalise_yes_no(is_there_cognitive_impairment_alzheimers_disease) == "TRUE" ~ "TRUE",
              normalise_yes_no(dementia_preop) == "TRUE" ~ "TRUE",
              TRUE ~ "FALSE"
            ),
            
            # 3. Post‑operative Delirium (Outcome) ---
            cog_delirium_7d_binary   = normalise_yes_no(delirium_7d),
            
            # 4. MMSE Scores ---
            cog_mmse_pre  = readr::parse_number(as.character(mmse_pre)),
            cog_mmse_6m   = readr::parse_number(as.character(mmse_6m)),
            cog_mmse_12m  = readr::parse_number(as.character(mmse_12m))
  )
))}
# Cognitive outcomes: MMSE ----
mmse_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Global scores ---
            mmse            = as.numeric(mmse),
            new_mmse_pre    = as.numeric(new_mmse_pre),
            new_mmse_fo     = as.numeric(new_mmse_fo),
            
            # Pre-op subdomains ---
            pre_mmse_temporal   = as.numeric(pre_mmse_temporal),
            pre_mmse_sp         = as.numeric(pre_mmse_sp),
            pre_mmse_memory     = as.numeric(pre_mmse_memory),
            pre_mmse_attention  = as.numeric(pre_mmse_attention),
            pre_mmse_only       = as.numeric(pre_mmse_only),
            pre_mmse_design     = as.numeric(pre_mmse_design),
            
            # Follow-up subdomains ---
            fo_mmse_temporal    = as.numeric(fo_mmse_temporal),
            fo_mmse_sp          = as.numeric(fo_mmse_sp),
            fo_mmse_memory      = as.numeric(fo_mmse_memory),
            fo_mmse_attention   = as.numeric(fo_mmse_attention),
            fo_mmse_lang        = as.numeric(fo_mmse_lang),
            fo_mmse_design      = as.numeric(fo_mmse_design),
            
            # Differences ---
            diff_mmse_temporal  = as.numeric(diff_mmse_temporal),
            diff_mmse_sp        = as.numeric(diff_mmse_sp),
            diff_mmse_memory    = as.numeric(diff_mmse_memory),
            diff_mmse_attention = as.numeric(diff_mmse_attention),
            diff_mmse_lang      = as.numeric(diff_mmse_lang),
            diff_mmse_design    = as.numeric(diff_mmse_design),
            mms_ediff           = as.numeric(mms_ediff),
            
            # Derived outcomes ---
            mci                 = factor(mci),
            mmse_sd             = as.numeric(mmse_sd),
            
            # Standard deviations ---
            diff_mmse_temporal_sd  = as.numeric(diff_mmse_temporal_sd),
            diff_mmse_sp_sd        = as.numeric(diff_mmse_sp_sd),
            diff_mmse_memory_sd    = as.numeric(diff_mmse_memory_sd),
            diff_mmse_attention_sd = as.numeric(diff_mmse_attention_sd),
            diff_mmse_lang_sd      = as.numeric(diff_mmse_lang_sd),
            diff_mmse_design_sd    = as.numeric(diff_mmse_design_sd),
            
            # NCD flags ---
            diff_mmse_temporal_ncd  = factor(diff_mmse_temporal_ncd),
            diff_mmse_sp_ncd        = factor(diff_mmse_sp_ncd),
            diff_mmse_memory_ncd    = factor(diff_mmse_memory_ncd),
            diff_mmse_attention_ncd = factor(diff_mmse_attention_ncd),
            diff_mmse_lang_ncd      = factor(diff_mmse_lang_ncd),
            diff_mmse_design_ncd    = factor(diff_mmse_design_ncd)
  )
))}
# Psychological outcomes / scales ----
psych_outcomes_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Single-item scales ---
            ksrq4_item   = factor(ksrq4),
            yzcd4_item   = factor(yzcd4),
            mzgx4_item   = factor(mzgx4),
            yzblsj4_item = factor(yzblsj4),
            jj4_item     = factor(jj4),
            hjrq4_item   = factor(hjrq4),
            tcsy4_item   = factor(tcsy4),
            
            # Other psychometric items ---
            fzl_zty_item   = normalise_yes_no(fzl_zty),
            apl_zty_item   = normalise_yes_no(apl_zty),
            bedzl_zjy_item = normalise_yes_no(bedzl_zjy),
            else_jsl_item  = normalise_yes_no(else_jsl),
            ca_mcat_item   = factor(ca_mcat),
            icdat_item     = factor(icdat),
            
            # Structured sub-scales (Form X) ---
            x1_1_x = factor(x1_1_x),
            x1_2_x = factor(x1_2_x),
            x1_3_x = factor(x1_3_x),
            x1_4_x = factor(x1_4_x),
            x1_5_x = factor(x1_5_x),
            x2_1_x = factor(x2_1_x),
            x2_2_x = factor(x2_2_x),
            x2_3_x = factor(x2_3_x),
            x2_4_x = factor(x2_4_x),
            x2_5_x = factor(x2_5_x),
            x3_x   = factor(x3_x),
            x4_x   = factor(x4_x),
            x5_x   = factor(x5_x),
            x6_1_x = factor(x6_1_x),
            x6_2_x = factor(x6_2_x),
            x7_x   = factor(x7_x),
           
  )
))}
# Quality of life ----
qol_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # SF-36 Completion Flags ---
            sf36_x_completed = normalise_yes_no(sf_36_x),
            sf36_y_completed = normalise_yes_no(sf_36_y),
            
            # SF-36 Subscales: Wave X (first follow-up) ---
            physical_functioning_x       = as.numeric(pf_x),
            role_physical_x              = as.numeric(rp_x),
            bodily_pain_x                = as.numeric(bp_x),
            general_health_x             = as.numeric(gh_x),
            vitality_x                   = as.numeric(vt_x),
            social_functioning_x         = as.numeric(sf_x),
            role_emotional_x             = as.numeric(re_x),
            mental_health_x              = as.numeric(mh_x),
            health_transition_x          = as.numeric(ht_x),
            
            # SF-36 Subscales: Wave Y (second follow-up) ---
            physical_functioning_y       = as.numeric(pf_y),
            role_physical_y              = as.numeric(rp_y),
            bodily_pain_y                = as.numeric(bp_y),
            general_health_y             = as.numeric(gh_y),
            vitality_y                   = as.numeric(vt_y),
            social_functioning_y         = as.numeric(sf_y),
            mental_health_y              = as.numeric(mh_y),
            health_transition_y          = as.numeric(ht_y),
            king_score                   = as.numeric(king),   
            
            # SF-36 Derived Scores (Standardised) ---
            physical_functioning_sd      = as.numeric(pf_sd),
            role_physical_sd             = as.numeric(rp_sd),
            bodily_pain_sd               = as.numeric(bp_sd),
            general_health_sd            = as.numeric(gh_sd),
            vitality_sd                  = as.numeric(vt_sd),
            social_functioning_sd        = as.numeric(sf_sd),
            role_emotional_sd            = as.numeric(re_sd),
            mental_health_sd             = as.numeric(mh_sd),
            
            # SF-36 Derived Scores (Component Metrics) ---
            physical_functioning_cm      = as.numeric(pf_cm),
            role_physical_cm             = as.numeric(rp_cm),
            bodily_pain_cm               = as.numeric(bp_cm),
            general_health_cm            = as.numeric(gh_cm),
            vitality_cm                  = as.numeric(vt_cm),
            social_functioning_cm        = as.numeric(sf_cm),
            role_emotional_cm            = as.numeric(re_cm),
            mental_health_cm             = as.numeric(mh_cm)
  )
))
}
# Survival ----
survival_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Early survival ---
            survival_30d        = normalise_yes_no(survival_30_days_after_surgery),
            
            # Survival status at follow-up ---
            survival_status_wave_x   = normalise_yes_no(survival_status_x),
            survival_status_wave_y   = normalise_yes_no(survival_status_y),
            
            # Time of death / survival time ---
            date_of_death            = excel_date(time_of_death),
            death_time_recorded_flag = normalise_yes_no(whether_there_is_a_time_of_death),
            survival_time_days       = as.numeric(survtime),
            
            # Mortality intervals ---
            mortality_0to1m          = normalise_yes_no(dead0to1m),
            mortality_1mto6m         = normalise_yes_no(dead1mto6m),
            mortality_6mto12m        = normalise_yes_no(dead6mto12m),
            mortality_0to12m         = normalise_yes_no(dead0mto12m),
            mortality_12mto3y        = normalise_yes_no(dead12mto3y),
            mortality_0to3y          = normalise_yes_no(dead0to3y),
            
            # Alternative mortality encodings ---
            mortality_0to1_flag      = normalise_yes_no(mort0to1),
            mortality_0to6_flag      = normalise_yes_no(mort0to6),
            mortality_0to12_flag     = normalise_yes_no(mort0to12),
            mortality_3yr_flag       = normalise_yes_no(mort3yr),
            
            # Study attrition metadata ---
            divergence_time          = as.numeric(divergence_time),
            deviating_content        = factor(deviating_content)
  )
))
}
# Costs, discharge and length of stay ----
costs_resource_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Direct costs ---
            total_medical_expenses   = as.numeric(total_medical_expenses),
            anaesthesia_cost          = as.numeric(cost_of_anesthesia),
            
            # Discharge & length of stay ---
            discharge_status         = factor(discharge),
            discharge_date           = excel_date(date_of_discharge),
            time_out_of_pacu         = as.numeric(time_out_of_pacu),
            hospital_days            = as.numeric(zhuyuan_days),
            hospitalised_7d_flag     = normalise_yes_no(hos7d),
            
            # Medication counts ---
            total_medication_count   = as.numeric(med_counts),
            number_of_drug_types     = as.numeric(number_of_drug_types),
            med_counts_zty           = as.numeric(zty_medcounts),
            
            # Resource classifications ---
            blood_loss_category      = factor(blood_loss_classification),
            transfusion_category     = factor(blood_transfusion_classification),
            sxzl_code                = factor(sxzl),
            shuxue_shixue_code       = factor(shuxue_shixue),
            sxzl_category            = factor(sxzl_cat),
            shixue_category          = factor(shixue_cat)
  )
))
}
# Follow‑up & Contact ----
followup_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Follow‑up dates ---
            follow_up_date_x = excel_date(follow_up_date_x),
            follow_up_date_y = excel_date(follow_up_date_y),
            
            # Contact status ---
            can_contact      = normalise_yes_no(can_the_patient_be_contacted),
            contacted        = normalise_yes_no(whether_the_patient_was_contacted),
            
            # Reporting / period ---
            period_x         = factor(period_x),
            report_date_1    = excel_date(report_date_1),
            report_date_2    = excel_date(report_date_2)
  )
))
}
# Trial / Randomisation ----
trial_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            # Randomisation process ---
            randomisation_datetime = excel_date(date_time_for_randomisation),
            treatment_allocation       = factor(random_results),
            protocol_violation     = normalise_yes_no(whether_it_violates),
            
            # Misc trial flags ---
            and_flag               = factor(and)
  )
))
}
# Redundant variables ----
redundant_vars <- {block(data_translated, quote(
  transmute(data_translated,
            id = as.character(individual_id),
            
            sex_2 = factor(sex_2),
            no    = factor(no),
            ytx   = factor(ytx),
            all_symptom = factor(all_symptom),
            anaesthesia_start_2 = as.POSIXct(anesthesia_start), #same as another column
            anaesthesia_end_2    = as.POSIXct(anesthesia_end),
            mortality_30d_flag       = normalise_yes_no(mort_30d), # same as survival
            random                 = factor(random) #same as treatment allocation
  )
))}

# Master list ----
master_list <- list(
  demographics = demographic_vars,
  background   = background_vars,
  bio_vars     = repeated_bio_vars,
  preop_char   = preop_features_vars,
  periop_vars  = periop_vars,
  periop_treatment = list(
    core = periop_treatment_vars,
    misc = periop_treatment_vars_misc),
  postop_vars  = analgesic_sedative_vars,
  adverse_vars = adverse_event_vars,
  outcomes = list(
    cognitive = list(
      core = cognitive_vars,
      mmse = mmse_vars,
      other = psych_outcomes_vars),
    survival = survival_vars,
    costs = costs_resource_vars),
  admin = list(
    consent      = consent_vars,
    followup     = followup_vars,
    trial_process = trial_vars
  ),
  misc = redundant_vars)

# Identify withdrawn patients
withdrawn_ids <- master_list$admin$consent |>
  filter(withdrawn == TRUE) |>
  pull(id)

#  extract all data frames
extract_dfs <- function(x) {
  if (is.data.frame(x)) list(x)
  else if (is.list(x)) map(x, extract_dfs) |> unlist(recursive = FALSE)
  else list()
}

# Flatten only valid IDs
valid_ids <- master_list$demographics$id
all_dfs <- extract_dfs(master_list) |> map(~ filter(.x, id %in% valid_ids))
master_flat <- reduce(all_dfs, full_join, by = "id")

# Variable list 

library(dplyr)
library(stringr)

clean_with_counts <- function(vars_string) {
  if (is.na(vars_string) || vars_string == "-" || vars_string == "") return("-")
  
  # Split the long string into a vector of names
  vec <- str_trim(unlist(str_split(vars_string, "\\|")))
  vec <- vec[vec != ""]
  
  # 1. Handle Time Series (always t0-t7)
  is_t <- str_detect(vec, "_t[0-7]$")
  t_bases <- unique(str_remove(vec[is_t], "_t[0-7]$"))
  t_part <- if(length(t_bases) > 0) paste0(t_bases, " [t0-7]") else NULL
  
  # 2. Handle Numbered Sequences (finding max N)
  others <- vec[!is_t]
  
  if(length(others) > 0) {
    df_counts <- data.frame(
      original = others,
      base = str_remove(others, "_?[0-9]{1,2}$"),
      num = as.numeric(str_extract(others, "[0-9]{1,2}$"))
    )
    
    # Only collapse if a base appears more than once
    seq_summary <- df_counts %>%
      filter(!is.na(num)) %>%
      group_by(base) %>%
      summarise(max_n = max(num), count = n(), .groups = 'drop') %>%
      filter(count > 1)
    
    final_vec <- c(t_part)
    processed_bases <- c()
    
    for (v in others) {
      b <- str_remove(v, "_?[0-9]{1,2}$")
      if (b %in% seq_summary$base) {
        if (!(b %in% processed_bases)) {
          m_n <- seq_summary$max_n[seq_summary$base == b]
          style <- if(str_detect(v, "_[0-9]")) " [1-" else "[1-"
          final_vec <- c(final_vec, paste0(b, style, m_n, "]"))
          processed_bases <- c(processed_bases, b)
        }
      } else {
        final_vec <- c(final_vec, v)
      }
    }
  } else {
    final_vec <- t_part
  }
  
  # 3. Clean up the list
  final_vec <- unique(final_vec)
  if ("id" %in% final_vec) {
    final_vec <- c("id", final_vec[final_vec != "id"])
  }
  
  return(paste(final_vec, collapse = " | "))
}

# Build flat_metadata from master_list ----
extract_metadata <- function(lst, section = "", sub_section = "") {
  results <- list()
  
  for (nm in names(lst)) {
    item <- lst[[nm]]
    if (is.data.frame(item)) attributes(item) <- attributes(item)[c("names", "row.names", "class")]
    
    current_section    <- if (section == "") nm else section
    current_subsection <- if (section == "") "" else if (sub_section == "") nm else paste(sub_section, nm, sep = " > ")
    
    if (is.data.frame(item)) {
      results <- append(results, list(data.frame(
        Section     = current_section,
        Sub_Section = current_subsection,
        Vars        = paste(names(item), collapse = " | "),
        Count       = ncol(item),
        stringsAsFactors = FALSE
      )))
    } else if (is.list(item)) {
      results <- append(results, extract_metadata(item, current_section, current_subsection))
    }
  }
    do.call(rbind, results)
}

flat_metadata <- extract_metadata(master_list)

flat_metadata <- flat_metadata %>%
  mutate(Vars_Display = sapply(Vars, clean_with_counts))

# Append tidy-defined name mappings back to the translation cache
cache_file <- here::here("data_interim", "translation_cache.csv")

# Read existing cache 
existing_cache <- if (file.exists(cache_file)) {
  read.csv(cache_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
} else {
  data.frame(Chinese = character(), English = character(), 
             Source = character(), stringsAsFactors = FALSE)
}

# Build a mapping from the translated (pre-tidy) names -> tidy final names
translated_names <- names(readRDS(here::here("data_interim", "raga_data_translated.rds")))
tidy_names       <- names(master_flat)   # or janitor::clean_names equivalent

# Only include columns that survived into master_flat
shared <- intersect(
  janitor::make_clean_names(translated_names),  
  tidy_names
)

tidy_map <- data.frame(
  Chinese = translated_names[janitor::make_clean_names(translated_names) %in% shared],
  English = tidy_names[tidy_names %in% shared],
  Source  = "tidy_rename",
  stringsAsFactors = FALSE
)

# Merge into cache, deduplicating on Chinese key
updated_cache <- dplyr::bind_rows(existing_cache, tidy_map) |>
  dplyr::filter(!duplicated(Chinese))

write.csv(updated_cache, cache_file, row.names = FALSE, fileEncoding = "UTF-8")
cat("Cache updated with", nrow(tidy_map), "tidy-defined name mappings.\n")

# Render the table
library(kableExtra)
flat_metadata %>%
  select(Section, Sub_Section, Vars_Display, Count) %>%
  kbl(col.names = c("Master Section", "Sub-Level", "Variables (Grouped)", "Total Count")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12) %>%
  column_spec(1, bold = T, width = "10em") %>%
  column_spec(3, width = "45em") %>% 
  collapse_rows(columns = 1, valign = "top") %>%
  row_spec(0, background = "#2c3e50", color = "white", bold = T)

# Export
write_xlsx(master_flat, here::here("data_interim", "master_table.xlsx"))
saveRDS(master_list, here::here("data_interim", "master_list.rds"))
saveRDS(master_flat, here::here("data_interim", "master_flat_raw.rds"))  
saveRDS(flat_metadata, here::here("data_interim", "metadata_dictionary.rds"))
rm(
  data_translated,      
  rds_path,             
  demographic_vars,     
  background_vars,      
  consent_vars,         
  repeated_bio_vars,
  periop_treatment_vars,
  periop_treatment_vars_misc,
  periop_vars,
  adverse_event_vars
)

gc()
