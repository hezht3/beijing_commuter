require(tidyverse)
require(chron)

setwd("F:/Box Sync/micro-environment/beijing_commuter/beijing_commuter (zh133@duke.edu)/data_verify")
data <- read_csv("beijing.csv")

# verify weekdays
data <- data %>%
            separate(start_time, c("start_date", "start_time"), sep = " ") %>%
                separate(end_time, c("end_date", "end_time"), sep = " ") %>%
                    separate(start_date, c("month", "day", "year"), sep = "/") %>%
                        unite(col = date_ver, year, month, day, sep = "/") %>%
                            mutate(date_ver = as.Date(date_ver)) %>%
                                mutate(weekdays_ver = ifelse(is.weekend(date_ver) == FALSE, "weekday", "weekend")) # no differences detected

# verify pattern
# peak is defined as 7:00 - 9:00 (420 minutes - 540 minutes) as well as 17:00 - 20:00 (1020 minutes - 1200 minutes)
data <- data %>%
            separate(start_time, c("start_hour", "start_minute"), sep = ":") %>%
                separate(end_time, c("end_hour", "end_minute"), sep = ":") %>%
                    mutate(start_hour = as.numeric(start_hour),
                           start_minute = as.numeric(start_minute),
                           end_hour = as.numeric(end_hour),
                           end_minute = as.numeric(end_minute)) %>%
                        mutate(start_min = start_hour * 60 + start_minute) %>%
                            mutate(end_min = end_hour * 60 + end_minute) %>%
                                mutate(mid_min = (end_min + start_min) / 2) %>%
                                    mutate(pattern_ver = case_when(mid_min < 420 ~ "off.peak",
                                                                   mid_min >= 420 & mid_min <= 540 ~ "peak",
                                                                   mid_min > 540 & mid_min < 1020 ~ "off.peak",
                                                                   mid_min >= 1020 & mid_min <= 1200 ~ "peak",
                                                                   mid_min > 1200 ~ "off.peak")) %>%
                                        mutate(pattern_ver = ifelse(weekdays == "weekend", "off.peak", pattern_ver)) # no differences detected

# import and clean meteorological factors dataset
data.met1 <- read_csv("met_data(1).csv")
data.met2 <- read_csv("wanliu_2015.csv")

data.met1 <- data.met1 %>%
                 separate(hour, c("date", "time"), sep = " ") %>%
                     separate(date, c("month", "day", "year"), sep = "/") %>%
                         unite(col = date, year, month, day, sep = "/") %>%
                             mutate(date = as.Date(date)) %>%
                                 mutate(time = str_replace(time, ":00", "")) %>%
                                     mutate(time = as.numeric(time)) %>%
                                         mutate(drop = NULL) %>%
                                             filter(date >= min(data$date_ver) & date <= max(data$date_ver))

data.met2 <- data.met2 %>%
                 separate(date, c("month", "day", "year"), sep = "/") %>%
                     unite(col = date, year, month, day, sep = "/") %>%
                         mutate(date = as.Date(date)) %>%
                             mutate(time = str_replace(hour, ":00:00", "")) %>%
                                 mutate(time = as.numeric(time)) %>%
                                     mutate(hour = NULL) %>%
                                         mutate(PM25_wanliu = pm2.5) %>%
                                             filter(date >= min(data$date_ver) & date <= max(data$date_ver))

data.met <- full_join(data.met1, data.met2, by = c("date", "time"))
data.met <- data.met %>% select(date, time, PM25_wanliu, TEMP, RH, BP, WD, WS, no2)

data.met <- data.met %>%
                mutate(time_start = time * 60) %>%
                    mutate(time_end = case_when(time == 23 ~ 24 * 60,
                                               TRUE ~ (time + 1) * 60))

# merge data with meteorological dataset
data <- data %>% mutate(period = end_min - start_min) # summary(data$period), min = 21, max = 49; ceiling(start_min/60)*60 = floor(end_min/60)*60

data <- data %>%
            mutate(time_point1 = floor(start_min / 60) * 60,
                   time_point2 = floor(end_min / 60) * 60,
                   time_point3 = ceiling(end_min / 60) * 60) %>%
                mutate(time_start = ifelse(mid_min < time_point2, time_point1, time_point2),
                       time_end = ifelse(mid_min < time_point2, time_point2, time_point3))

data.met <- data.met %>%
                mutate(date_ver = date,
                       PM25_wanliu_ver = PM25_wanliu,
                       TEMP_ver = TEMP,
                       RH_ver = RH,
                       BP_ver = BP,
                       WD_ver = WD,
                       WS_ver = WS,
                       NO2_ver = no2)

data <- data %>%
            left_join(data.met[, c("date_ver", "time_start", "time_end", "PM25_wanliu_ver", "TEMP_ver", "RH_ver", "BP_ver", "WD_ver", "WS_ver", "NO2_ver")],
                      by = c("date_ver", "time_start", "time_end"))

data_na <- data %>%
               filter(is.na(PM25_wanliu_ver) == TRUE | is.na(TEMP_ver) == TRUE | is.na(RH_ver) == TRUE | is.na(BP_ver) == TRUE |
                      is.na(WD_ver) == TRUE | is.na(WS_ver) == TRUE | is.na(NO2_ver) == TRUE) %>%
                   right_join(data.met, by = c("date_ver", "time_start", "time_end")) %>%
                       arrange(date, time) %>%
                           mutate(time_point0 = time_point1 - 60,
                                  time_point4 = time_point3 + 60) %>%
                               mutate(time_start = case_when(is.na(start_min) == FALSE & mid_min < time_point2 & (time_point3 - end_min) < (start_min - time_point0) ~ time_point2,
                                                             is.na(start_min) == FALSE & mid_min < time_point2 & (time_point3 - end_min) >= (start_min - time_point0) ~ time_point0,
                                                             is.na(start_min) == FALSE & mid_min >= time_point2 & (time_point4 - end_min) < (start_min - time_point1) ~ time_point3,
                                                             is.na(start_min) == FALSE & mid_min >= time_point2 & (time_point4 - end_min) >= (start_min - time_point1) ~ time_point1,
                                                             TRUE ~ time_start),
                                      time_end = case_when(is.na(start_min) == FALSE & mid_min < time_point2 & (time_point3 - end_min) < (start_min - time_point0) ~ time_point2 + 60,
                                                           is.na(start_min) == FALSE & mid_min < time_point2 & (time_point3 - end_min) >= (start_min - time_point0) ~ time_point0 + 60,
                                                           is.na(start_min) == FALSE & mid_min >= time_point2 & (time_point4 - end_min) < (start_min - time_point1) ~ time_point3 + 60,
                                                           is.na(start_min) == FALSE & mid_min >= time_point2 & (time_point4 - end_min) >= (start_min - time_point1) ~ time_point1 + 60,
                                                           TRUE ~ time_end)) %>%
                                          select(colnames(data)) %>%
                                              right_join(data.met, by = c("date_ver", "time_start", "time_end")) %>%
                                                  filter(is.na(start_min) == FALSE) %>%
                                                      select(colnames(data), time_start, time_end, PM25_wanliu_ver, TEMP_ver, RH_ver, BP_ver, WD_ver, WS_ver, NO2_ver)

data <- data %>%
            filter(is.na(PM25_wanliu_ver) == TRUE | is.na(TEMP_ver) == TRUE | is.na(RH_ver) == TRUE | is.na(BP_ver) == TRUE |
                       is.na(WD_ver) == TRUE | is.na(WS_ver) == TRUE | is.na(NO2_ver) == TRUE) %>%
                bind_rows(data_na)
# differences
# nrow = 16, PM25_wanliu_ver = NA --> time_start = 1080, time_end = 1140, PM25_wanliu_Ver = 19.0, correct
# nrow = 25, NO2_ver = NA --> start_min - 1080 < 1140 - end_min --> time_start = 1020, time_end = 1080, NO2 = 39.6, correct
# nrow = 61, NO2_ver = NA --> 540 - start_min < 600 - end_min --> time_start =480, time_end = 540, NO2 = 55.0, correct

# final conclusion: all cells are correct

saveRDS(data.met[, c("date", "time", "PM25_wanliu", "TEMP", "RH", "BP", "WD", "WS", "no2")],
        file = "F:/Box Sync/micro-environment/beijing_commuter/beijing_commuter (zh133@duke.edu)/data_met.rds")