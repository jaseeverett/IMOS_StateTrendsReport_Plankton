# Code to analyse and plot the Integrated Marine Observing System (IMOS) Phytoplankton and Zooplankton data
# for the State and Trends Report.
#
# Jason Everett (UQ/CSIRO/UNSW) and Anthony Richardson (UQ/CSIRO)
# Last Updated: November 2019

library(tidyverse)
library(ggpubr)
library(lubridate)

re_down <- 0 # Should we download the files again

## Plot the measured Phytoplankton Biomass
p_file = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=anmn_nrs_bgc_plankton_phytoplankton_data&SERVICE=WFS&outputFormat=csv&REQUEST=GetFeature&VERSION=1.0.0";
p_out_file = paste0('Data',.Platform$file.sep,'IMOS_National_Reference_Station_(NRS)_-_Phytoplankton_Abundance_and_Biovolume.csv')

if (re_down==1){
  download.file(p_file, p_out_file, method = "auto", quiet=FALSE)
}

z_file = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=anmn_nrs_bgc_plankton_biomass_data&SERVICE=WFS&outputFormat=csv&REQUEST=GetFeature&VERSION=1.0.0";
z_out_file = paste0('Data',.Platform$file.sep,'IMOS National Reference Station (NRS) - Zooplankton Biomass.csv')

if (re_down==1){
  download.file(z_file,z_out_file, method = "auto", quiet=FALSE);
}

## Phytoplankton
pdat <- read_csv(p_out_file)
pdat2 <- pdat %>% group_by(NRS_TRIP_CODE) %>%
  summarise(Station = STATION_NAME[1],
            LocalTime = as.Date(LOCAL_TRIP_START_TIME[1]),
            IMOS_SiteCode = IMOS_SITE_CODE[1],
            NRS_TripCode = NRS_TRIP_CODE[1],
            BioVol_m3 = log10(sum(BIOVOLUME_UM3_PER_L, na.rm = TRUE)/1e3),
            Abundance_m3 = log10(sum(CELL_PER_LITRE, na.rm = TRUE)/1e3),
            Month = as.factor(month(LOCAL_TRIP_START_TIME[1])))
rm(pdat)

## Zooplankton
zdat <- read_csv(z_out_file)
zdat2 <- zdat %>% group_by(NRS_TRIP_CODE) %>%
  summarise(Station = STATION_NAME[1],
            LocalTime = as.Date(LOCAL_TRIP_START_TIME[1]),
            IMOS_SiteCode = IMOS_SITE_CODE[1],
            NRS_TripCode = NRS_TRIP_CODE[1],
            ZooBiomass_m3 = log10(MG_PER_M3),
            Month = as.factor(month(LOCAL_TRIP_START_TIME[1])))
rm(zdat)

# There is currently a missing pt in the Zoop Biomass Data. Claire is checking. For the moment, lets remove inf values
zdat2 <- drop_na(zdat2)
# zdat2 <- zdat2[is.finite(zdat2$ZooBiomass_m3),]

## Plotting        
sites <- c("NRSDAR", "NRSYON", "NRSNSI", "NRSPHB",  "NRSROT", "NRSKAI", "NRSMAI")
site_names <- c("a) Darwin Harbour","b) Darwin Harbour",
                "c) Yongala","d) Yongala",
                "e) North Stradbroke Island", "f) North Stradbroke Island",
                "g) Port Hacking","h) Port Hacking",
                "i) Rottnest Island","j) Rottnest Island", 
                "k) Kangaroo Island", "l) Kangaroo Island", 
                "m) Maria Island", "n) Maria Island")

myplots <- list() # Initialise plot list
p_Bio_slope <- NULL
z_Bio_slope <- NULL
p_Bio_decade_change <- NULL
z_Bio_decade_change <- NULL
p_BioPerc_decade_change <- NULL
z_BioPerc_decade_change <- NULL

counter <- 1

for (i in 1:length(sites)) {
  pm_dat <- pdat2 %>% filter(IMOS_SiteCode==sites[i])
  pmdl <- lm(BioVol_m3 ~ LocalTime + Month, data = pm_dat)
  p_pval <- anova(pmdl)$'Pr(>F)'[1]
  p_Bio_slope[i] <- pmdl$coefficients[[2]] # Save the model slope
  
  zm_dat <- zdat2 %>% filter(IMOS_SiteCode==sites[i])
  zmdl <- lm(ZooBiomass_m3 ~ LocalTime + Month, data = zm_dat)
  z_pval <- anova(zmdl)$'Pr(>F)'[1]
  z_Bio_slope[i] <- zmdl$coefficients[[2]] # Save the model slope
  
  # I have checked and the slopes are identical if I use date (above) or days (below).
  # Obviously R correctly deals with the data in the lm.
  ref_date = dmy('01-1-2009')
  # pm_dat$day_number <- as.numeric(difftime(pm_dat$LocalTime, ref_date)) 
  # pmdl2 <- lm(BioVol_m3 ~ day_number + Month, data = pm_dat)
  
  # So this means that the slope is the change per day.
  # Now I calculate the days over a decade and calculate the change
  decade_length <- as.numeric(difftime(dmy("01-01-2019"), ref_date))
  
  ## Phytoplankton decade trends
  pdata_length <- as.numeric(difftime(pm_dat$LocalTime[length(pm_dat$LocalTime)], pm_dat$LocalTime[1]))
  ptemp <- predict(pmdl) # Predict on existing dates
  ppred <- ptemp[c(1,length(ptemp))] # First and last predicted answer to get overall change
  p_Bio_decade_change[i] <- ((diff(10^ppred))/pdata_length)*decade_length # Correct to decade
  p_BioPerc_decade_change[i] <- p_Bio_decade_change[i]/(10^ppred[1])
  
  ## Zooplankton decade trends
  zdata_length <- as.numeric(difftime(zm_dat$LocalTime[length(zm_dat$LocalTime)], zm_dat$LocalTime[1]))
  ztemp <- predict(zmdl) # Predict on existing dates
  zpred <- ztemp[c(1,length(ztemp))] # First and last predicted answer to get overall change
  z_Bio_decade_change[i] <- ((diff(10^zpred))/zdata_length)*decade_length # Correct to decade
  z_BioPerc_decade_change[i] <- z_Bio_decade_change[i]/(10^zpred[1])
  
  if (p_Bio_slope[i]>0 & p_pval <= 0.05){
    pclr <- "red"
  } else if (p_Bio_slope[i]<0 & p_pval <= 0.05){
    pclr <- "blue"
  } else {
    pclr <- "black"
  }
  
  if (z_Bio_slope[i]>0 & z_pval <= 0.05){
    zclr <- "red"
  } else if (z_Bio_slope[i]<0 & z_pval <= 0.05){
    zclr <- "blue"
  } else {
    zclr <- "black"
  }
  
  gg <-  ggplot() +
    geom_point(data = pm_dat, aes(x=LocalTime, y=BioVol_m3)) +
    stat_smooth(data = pmdl, method = "lm", col = pclr, fill = pclr, aes_string(x = names(pmdl$model)[2], y = names(pmdl$model)[1])) +
    theme_bw() + 
    scale_y_continuous(limits=c(4, 8), expand = c(0, 0)) +
    annotate("text", x = as.Date("2009-2-1"), y = 7.6, label = site_names[counter], hjust = 0)
  
  if (i == 7) {
    gg <- gg + 
      scale_x_date(date_labels = "%Y", date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # axis.text.x = element_text(angle = 45, hjust = 1),
  } else {
    gg <- gg +
      scale_x_date(labels = NULL, date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  myplots[[counter]] <- gg
  counter <- counter + 1
  rm(gg)
  
  gg <-  ggplot() +
    geom_point(data = zm_dat, aes(x=LocalTime, y=ZooBiomass_m3)) +
    stat_smooth(data = zmdl, method = "lm", col = zclr, fill = zclr, aes_string(x = names(zmdl$model)[2], y = names(zmdl$model)[1])) +
    theme_bw() + 
    annotate("text", x = as.Date("2009-2-1"), y = 2.6, label = site_names[counter], hjust = 0) +
    scale_y_continuous(position="right",limits=c(-1, 3), expand = c(0, 0)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  if (i == 7) {
    gg <- gg + 
      scale_x_date(labels = NULL, date_labels = "%Y", date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  } else {
    gg <- gg +
      scale_x_date(labels = NULL, date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  }
  
  myplots[[counter]] <- gg
  counter <- counter + 1
  rm(gg, pmdl, zmdl)
}

graphics.off()
figure <- ggarrange(plotlist=myplots, ncol = 2, nrow = 7,widths = 1, heights = 1)
annotate_figure(figure, 
                left = text_grob(expression(paste('Phytoplankton Biovolume log'[10],'(',mu,'m m'^{-3},')')), 
                                         color = "black", rot = 90, size = 12),
                right = text_grob(expression(paste('Zooplankton Biomass log'[10],'(mg m'^{-3},')')), 
                                    color = "black", rot = 270, size = 12))
ggsave(paste0('Figures',.Platform$file.sep,'STAR_Plankton_Biomass.png'), dpi=300)

rm(zdat2)

##################################################
##################################################
##################################################

## Phyto abundance has already been done above....

## Just do Zoop
z_file = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=anmn_nrs_bgc_plankton_zooplankton_data&SERVICE=WFS&outputFormat=csv&REQUEST=GetFeature&VERSION=1.0.0";
z_out_file = paste0('Data',.Platform$file.sep,'IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance.csv')

if (re_down==1){
  download.file(z_file,z_out_file, method = "auto", quiet=FALSE);
}

## Zooplankton
zdat <- read_csv(z_out_file)

# Remove Noctiluca, eggs and plastics
zdat <- zdat %>% 
  filter(TAXON_ECO_GROUP != "Noctilucaceae" & 
           TAXON_NAME != "Fish egg" & 
           TAXON_ECO_GROUP != "Egg" & 
           TAXON_GROUP != "NON BIOLOGICAL")

zdat2 <- zdat %>% group_by(NRS_TRIP_CODE) %>%
  summarise(Station = STATION_NAME[1],
            LocalTime = as.Date(LOCAL_TRIP_START_TIME[1]),
            IMOS_SiteCode = IMOS_SITE_CODE[1],
            NRS_TripCode = NRS_TRIP_CODE[1],
            Abundance_m3 = log10(sum(TAXON_PER_M3, na.rm = TRUE)),
            Month = as.factor(month(LOCAL_TRIP_START_TIME[1])))
rm(zdat)

myplots <- list() # Initialise plot list
p_Abun_slope <- NULL
z_Abun_slope <- NULL
p_Abun_decade_change <- NULL
z_Abun_decade_change <- NULL

p_AbunPerc_decade_change <- NULL
z_AbunPerc_decade_change <- NULL

counter <- 1

for (i in 1:length(sites)) {
  pm_dat <- pdat2 %>% filter(IMOS_SiteCode==sites[i])
  pmdl <- lm(Abundance_m3 ~ LocalTime + Month, data = pm_dat)
  p_pval <- anova(pmdl)$'Pr(>F)'[1]
  p_Abun_slope[i] <- pmdl$coefficients[[2]] # Save the model slope
  
  zm_dat <- zdat2 %>% filter(IMOS_SiteCode==sites[i])
  zmdl <- lm(Abundance_m3 ~ LocalTime + Month, data = zm_dat)
  z_pval <- anova(zmdl)$'Pr(>F)'[1]
  z_Abun_slope[i] <- zmdl$coefficients[[2]]
  
  # Now I calculate the days over a decade and calculate the change
  decade_length <- as.numeric(difftime(dmy("01-01-2019"), ref_date))
  
  pdata_length <- as.numeric(difftime(pm_dat$LocalTime[length(pm_dat$LocalTime)], pm_dat$LocalTime[1]))
  ptemp <- predict(pmdl)
  ppred <- ptemp[c(1,length(ptemp))]
  p_Abun_decade_change[i] <- ((diff(10^ppred))/pdata_length)*decade_length
  p_AbunPerc_decade_change[i] <- p_Abun_decade_change[i]/(10^ppred[1])
  
  zdata_length <- as.numeric(difftime(zm_dat$LocalTime[length(zm_dat$LocalTime)], zm_dat$LocalTime[1]))
  ztemp <- predict(zmdl)
  zpred <- ztemp[c(1,length(ztemp))]
  z_Abun_decade_change[i] <- ((diff(10^zpred))/zdata_length)*decade_length
  z_AbunPerc_decade_change[i] <- z_Abun_decade_change[i]/(10^zpred[1])
  
  if (p_Abun_slope[i]>0 & p_pval <= 0.05){
    pclr <- "red"
  } else if (p_Abun_slope[i]<0 & p_pval <= 0.05){
    pclr <- "blue"
  } else {
    pclr <- "black"
  }
  
  if (z_Abun_slope[i]>0 & z_pval <= 0.05){
    zclr <- "red"
  } else if (z_Abun_slope[i]<0 & z_pval <= 0.05){
    zclr <- "blue"
  } else {
    zclr <- "black"
  }
  
  gg <-  ggplot() +
    geom_point(data = pm_dat, aes(x=LocalTime, y=Abundance_m3)) +
    stat_smooth(data = pmdl, method = "lm", col = pclr, fill = pclr, aes_string(x = names(pmdl$model)[2], y = names(pmdl$model)[1])) +
    theme_bw() + 
    scale_y_continuous(limits=c(0, 4), expand = c(0, 0)) +
    annotate("text", x = as.Date("2009-2-1"), y = 3.6, label = site_names[counter], hjust = 0) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  if (i == 7) {
    gg <- gg + 
      scale_x_date(labels = NULL, date_labels = "%Y", date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  } else {
    gg <- gg +
      scale_x_date(labels = NULL, date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  }
  myplots[[counter]] <- gg
  counter <- counter + 1
  rm(gg)
  
    gg <-  ggplot() +
    geom_point(data = zm_dat, aes(x=LocalTime, y=Abundance_m3)) +
    stat_smooth(data = zmdl, method = "lm", col = zclr, fill = zclr, aes_string(x = names(zmdl$model)[2], y = names(zmdl$model)[1])) +
    theme_bw() + 
    annotate("text", x = as.Date("2009-2-1"), y = 4.6, label = site_names[counter], hjust = 0) +
    scale_y_continuous(position="right",limits=c(1, 5), expand = c(0, 0)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  if (i == 7) {
    gg <- gg + 
      scale_x_date(labels = NULL, date_labels = "%Y", date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  } else {
    gg <- gg +
      scale_x_date(labels = NULL, date_breaks = "years", limits = c(as.Date("2009-1-1"), as.Date("2018-12-31")), expand = c(0, 0))
  }
  
  myplots[[counter]] <- gg
  counter <- counter + 1
  rm(gg, pmdl, zmdl)
}

graphics.off()
figure <- ggarrange(plotlist=myplots, ncol = 2, nrow = 7)
annotate_figure(figure, 
                left = text_grob(expression(paste('Phytoplankton Abundance log'[10],'(ind. m'^{-3},')')), 
                                 color = "black", rot = 90, size = 12),
                right = text_grob(expression(paste('Zooplankton Abundance log'[10],'(ind. m'^{-3},')')), 
                                  color = "black", rot = 270, size = 12))
ggsave(paste0('Figures',.Platform$file.sep,'STAR_Plankton_Abundance.png'), dpi=300)


