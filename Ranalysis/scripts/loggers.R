library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(tidyr)
library(patchwork)
library(ggbreak)
library(emmeans)
library(ggpattern)
library(anytime)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

setwd("/home/gospozha/haifa/O.patagonica/loggers/data/")

# prepare each dataset from different depths

m10=read.csv("./Licor general code/10m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m10$Depth="10"
m10$Depth=as.factor(m10$Depth)
m10$Local_timestamp=as.POSIXct(m10$Unix_Timestamp)
str(m10)

m25=read.csv("./Licor general code/25m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m25$Depth="25"
m25$Depth=as.factor(m25$Depth)
m25$Local_timestamp=as.POSIXct(m25$Unix_Timestamp)
str(m25)

m45=read.csv("./Licor general code/45m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m45$Depth="45"
m45$Depth=as.factor(m45$Depth)
m45$Local_timestamp=as.POSIXct(m45$Unix_Timestamp)
str(m45)

# merge all depths

m_all=rbind(m10, m25, m45)
str(m_all)

# cat datasets so they start and end at the same date
df <- m_all %>%
  filter(Local_timestamp >= "2024-06-20 11:19:00",
         Local_timestamp <="2024-09-29 13:39:00") %>%
  select(Local_timestamp, Depth, Temperature, PAR)%>%
  mutate(
    datetime = ymd_hms(Local_timestamp),
    date = as.Date(datetime)
  )


# Temperature - as one line

# Temperature: daily min and max per depth
daily_temp <- df %>%
  group_by(date, Depth) %>%
  summarise(
    min_temp = min(Temperature, na.rm = TRUE),
    max_temp = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # reshape to long format for plotting as one zig-zag line
  tidyr::pivot_longer(cols = c(min_temp, max_temp),
                      names_to = "stat",
                      values_to = "temp") %>%
  arrange(Depth, date, stat)

# PAR: daily max per depth, only until July 20
daily_par <- df %>%
  filter(date <= as.Date("2024-07-20")) %>%
  group_by(date, Depth) %>%
  summarise(
    max_PAR = max(PAR, na.rm = TRUE),
    .groups = "drop"
  )

# Temp plot
temp_plot <- ggplot(daily_temp, aes(x = date, y = temp, color = Depth, group = Depth)) +
  geom_line(size = 0.4) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  labs(
    title = "Daily temperature",
    x = "Date (DD/MM)",
    y = "Temperature (°C)",
    color = "Depth (m)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# PAR plot
par_plot <- ggplot(daily_par, aes(x = date, y = max_PAR, color = Depth)) +
  geom_line(size = 0.4) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  labs(
    title = "Daily maximum PAR",
    x = "Date (DD/MM)",
    y = "PAR (µmol/(s·m²))",
    color = "Depth (m)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined <- (par_plot / temp_plot) +
  plot_annotation(tag_levels = "A")

combined

ggsave("combined_loggers.jpg", combined, width = 8, height = 7)


# Temperature: as ribbon
daily_temp <- df %>%
  group_by(date, Depth) %>%
  summarise(
    min_temp = min(Temperature, na.rm = TRUE),
    max_temp = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  )

# Temperature plot (min and max as ribbons )
temp_plot <- ggplot(daily_temp, aes(x = date)) +
  geom_ribbon(aes(ymax = max_temp, ymin = min_temp,fill = Depth), alpha = 0.7) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  labs(
    title = "Daily temperature",
    x = "Date (MM/YY)",
    y = "Temperature (°C)",
    color = "Depth (m)",
    linetype = "Statistic"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined <- (par_plot / temp_plot) +
  plot_annotation(tag_levels = "A")

combined

ggsave("combined_loggers_ribbon.jpg", combined, width = 8, height = 7)
