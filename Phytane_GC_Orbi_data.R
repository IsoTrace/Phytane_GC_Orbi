#MEASURING INTRAMOLEUCLAR ISOTOPIC SIGNATURES OF PHYTANE USING ULTRA HIGH RESOLUTION ORBITRAP MASS SPECTROMETRY#
#Direct injections
#------------------
# created by Merve Oeztoprak
# m.oeztoprak@gmail.com
#------------------
#Load required packages-------
packages = c("plyr","dplyr","ggplot2","ggpubr","slider","psych","tidyverse","scales","cowplot","tidyquant","data.table","zoo","patchwork",
             "gtable","ggExtra","ggstance","gridExtra","grid","lubridate","datasets","gapminder","ggsignif","reshape2", "viridis","Thermimage")
# Now load or install&load all
package.check <- lapply(packages,
                        FUN = function(x) {
                          if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)
                          }
                        }
)
#
# Data import------------
# store all .isox files of mixed sequence same day acquisitions in a folder inside your working directory named mixed
filenames <- list.files("isox_files", pattern="*.isox", full.names=TRUE)
list2env(lapply(setNames(filenames, make.names(gsub("*.isox$", "", filenames))),
                read.table, header=T),  envir = .GlobalEnv)# will give "mixed.filename" in global R environment

#for exact masses use isotopologs.tsv same file as IsoX stored inside your working directory
masses= read.table("isotopologs.tsv", header=F) # this isotopologs.tsv is also used by isoX to extract the target masses from the raw thermos output files.
colnames(masses)=c("compound","isotopolog","m/z", "tolerance [mmu]", "z")
# mass tolerance for mass assignment in daltons (currently same as in isox exporting file but may be useful to test things in the future)
tolerance= 0.001
#supply list of isotopologs of interest as defined in isox file
ilogs=c("13C","0U","2H")
#define ions of interest
ions=c('C6H13')
#define fragments of interest
frags=c('C6H13')
#
# ----------------
#Load functions------
#imports data, keeping only compounds and isotopologs of interest, expands columns ranking satellite peaks
I_df <- function(x,name,sampleID,AQdate,masses,ions,ilogs){
  #reformatting
  data= x
  colnames(masses)=c("compound","isotopolog","m/z", "tolerance [mmu]", "z")
  data= merge(data,masses[,c(1:3)], by=c("compound","isotopolog")) #append data.frame with exact masses for fragments\
  #add additional data quality columns
  data= data %>%
    mutate("massdefect[da]" = abs(data$mzMeasured-data$`m/z`)) %>%#append column mass defect = absolute amount deviation from exact mass in daltons
    filter(data$compound %in% ions & #only include scans with fragments/ions of interest i.e. C6H13 fragment
             data$isotopolog %in% ilogs #only include scans with isotopologs of interest i.e. C6H13 fragment with 12C (unsubstituted) and 1data13C
    ) %>%
    ungroup() %>%
    mutate(filename=name)%>%
    group_by(filename, scan.no, time.min, compound, isotopolog) %>% # rank duplicate isotopolog (assigned b/c mass match) in same scan based on signal intensity
    mutate(num_dups = n(), 
           rank= rank(-intensity),
           is_duplicated = num_dups > 1, # logical column identifying if there are duplicate isotopologs in single scan
           "ticxit"= tic * it.ms,
           element = str_extract(isotopolog, "(?<=\\d)\\p{L}+"),
           ID= sampleID,
           date= AQdate)%>%
    ungroup()%>%
    group_by(filename,compound, isotopolog) %>%
    ungroup()%>%
    mutate(n.H=case_when(compound=="C6H9"~ 9,
                         compound=="C6H11"~ 11,
                         compound=="C6H13"~ 13),
           n.C=case_when(compound=="C6H9"|compound=="C6H11"|compound=="C6H13"~ 6))
}
#defines background dataframe calculating mean value for ions of interest
b_df <- function (I_data,start.min){
  #determine intensity of backgorund
  data.background=I_data%>%
    filter(time.min < start.min - 1)%>%
    group_by(filename,compound,element,isotopolog)%>%
    mutate(major_mz = if (any(is_duplicated == "FALSE")){mean(mzMeasured[is_duplicated == FALSE])} 
           else {mean(mzMeasured[rank == 1])})%>% # calculate mean mass measured for non duplicate scans for each compound
    ungroup()%>%
    group_by(filename,compound,element,isotopolog,n.H,n.C, scan.no)%>%
    mutate(massdef.major_mz = abs(major_mz-mzMeasured))%>%
    mutate(main_mz=ifelse(massdef.major_mz == min(massdef.major_mz),TRUE,FALSE))%>%
    ungroup()%>%
    filter(main_mz == TRUE)%>%
    group_by(filename,compound,element,isotopolog)%>%
    summarise(b.int=median(intensity, na.rm = T))
}
#defines sample dataframe after valve turning excluding scans below background level, 
#gives percent of scans of outliers, satellite peaks and below background
s_df <- function(I_data, b_data, outlier_factor) {
  data.sample <- I_data %>%
    filter(time.min > start.min) %>%
    group_by(filename, compound, element, isotopolog) %>%
    mutate(major_mz = if (any(is_duplicated == "FALSE")) {
      mean(mzMeasured[is_duplicated == FALSE])
    } else {
      mean(mzMeasured[rank == 1])
    }) %>% # calculate mean mass measured for non duplicate scans for each compound
    ungroup() %>%
    group_by(filename, compound, element, isotopolog, n.H, n.C, scan.no) %>%
    mutate(massdef.major_mz = abs(major_mz - mzMeasured)) %>%
    mutate(main_mz = ifelse(massdef.major_mz == min(massdef.major_mz), TRUE, FALSE)) %>%
    ungroup() %>%
    filter(main_mz == TRUE) # remove scans with duplicate/satellite peak
  
  s_data.out <- left_join(data.sample, b_data, by = c("filename", "compound", "element", "isotopolog")) %>%
    group_by(filename, compound, element, isotopolog) %>%
    mutate(
      below.background = ifelse(intensity < b.int, TRUE, FALSE),
      prcnt.below.background.scans = 100 - (length(which(below.background == FALSE))) / length(scan.no) * 100
    ) %>%
    ungroup() %>%
    group_by(filename, compound) %>%
    mutate(
      # Calculate the 10th percentile of TIC across all scans for the compound
      tenth_percentile_tic = quantile(tic, probs = 0.1, na.rm = TRUE),
      # Identify scans below 10% of TIC distribution
      below.ten.percent = ifelse(tic < tenth_percentile_tic, TRUE, FALSE),
      prcnt.below.10.max = 100 * sum(below.ten.percent, na.rm = TRUE) / length(unique(scan.no))
    ) %>%
    # identify outliers
    group_by(filename, compound) %>%
    mutate(
      medianTICxIT = median(ticxit),
      sdTICxIT = sd(ticxit),
      upsd = medianTICxIT + outlier_factor * sdTICxIT,
      lowsd = medianTICxIT - outlier_factor * sdTICxIT,
      outlier = ifelse(ticxit < upsd & ticxit > lowsd, FALSE, TRUE)
    ) %>%
    ungroup() %>%
    group_by(filename, compound) %>% 
    mutate(prcnt.outlier.scans = 100 - (length(which(outlier == FALSE))) / length(scan.no) * 100) %>%
    ungroup() %>%
    #percent of scnas missing 0U,13C or 2H
    group_by(filename, compound, scan.no) %>%
    mutate(
      missing.0U   = !any(isotopolog == "0U"),
      missing.13C  = !any(isotopolog == "13C"),
      missing.2H   = !any(isotopolog == "2H")
    ) %>%
    ungroup() %>%
    group_by(filename, compound) %>%
    mutate(n.total.scans = length(unique(scan.no))) %>% # length unique scan.no since multiple ilog per scan separated into multiple rows
    ungroup()%>%
  return(s_data.out)
}
#generates final dataframe with calculated Nio values, outliers removed
nio_df <- function(s_data){
  data_out=s_data %>%
    filter(!below.ten.percent,!below.background,!outlier,!missing.0U,!missing.13C)%>% #filtering
    ungroup()%>%
    group_by(filename, compound, scan.no, time.min)%>%
    filter(all(c("13C","0U") %in% isotopolog))%>% # filter scans where "13C","0U"
    ungroup()%>%
    #calculate Nio
    group_by(filename, compound, scan.no, time.min, isotopolog) %>%
    mutate(Nio = (intensity/peakNoise)*(4.4 /1)*((120000/resolution)^0.5)*(microscans^0.5))%>% #calculated Nio values (Eiler 2017) CN=4.4 ref.resolution=120000, z=1)
    ungroup()%>%
    group_by(filename,compound,element,isotopolog) %>%
    mutate(n.total.scans =length(unique(scan.no)))%>% #length unique scan.no since multiple ilog per scan separated into multiple rows
    ungroup() %>%
    select(filename, ID, type, date, compound, isotopolog, element, n.H, n.C, scan.no, time.min, n.total.scans, 
           tic, ticxit,intensity,b.int,ions.incremental, peakResolution, peakNoise, 
           Nio, missing.2H)%>%
    arrange(filename,isotopolog,n.C,n.H)%>%
  #shot noise calculations
  group_by(filename, compound, scan.no)%>%
    mutate(
      Nio_13C = Nio[isotopolog == "13C"],
      Nio_0U = Nio[isotopolog == "0U"],
      Nio.C.ratio = Nio_13C / Nio_0U
    )%>%
    ungroup()
}
#calculate fractional abundances
frac.output.weighted.avg.ratio = function(all_nio) {
  data_prep <- all_nio %>%
    complete(nesting(filename, compound, scan.no), isotopolog, 
             fill = list(Nio = 0, intensity = 0), explicit = FALSE) %>%
    group_by(filename, compound, isotopolog) %>%
    split(.$isotopolog) %>%
    lapply(function(i) i[order(i$filename, i$scan.no), ]) %>%
    as.data.frame() %>%
    rename(
      filename = X0U.filename, ID = X0U.ID, date = X0U.date,
      compound = X0U.compound, scan.no = X0U.scan.no, time.min = X0U.time.min,
      n.H = X0U.n.H, n.C = X0U.n.C, tic = X0U.tic, missing.2H = X0U.missing.2H,
      intensity.0U = X0U.intensity, intensity.13C = X13C.intensity, intensity.2H = X2H.intensity,
      Nio.0U = X0U.Nio, Nio.13C = X13C.Nio, Nio.2H = X2H.Nio
    ) %>%
    select(
      filename, ID, date, compound, scan.no, time.min, n.H, n.C, 
      tic, missing.2H, intensity.0U, intensity.13C, intensity.2H,
      Nio.0U, Nio.13C, Nio.2H
    ) %>%
    group_by(filename, ID, date, compound, n.H, n.C) %>%
    summarise(
      # Calculate basic statistics
      n_scans = n(),
      n_zero_2H = sum(intensity.2H == 0),
      prcnt.missing.2H = n_zero_2H/n_scans*100,
      
      # Calculate weighted Nio values
      Nio.w.0U = sum(Nio.0U * intensity.0U, na.rm = TRUE) / sum(intensity.0U, na.rm = TRUE),
      Nio.w.13C = sum(Nio.13C * intensity.13C, na.rm = TRUE) / sum(intensity.13C, na.rm = TRUE),
      
      # Variance
      var.Nio.w.0U = sum(intensity.0U * (Nio.0U - Nio.w.0U)^2, na.rm = TRUE) /
                              sum(intensity.0U, na.rm = TRUE),
      var.Nio.w.13C = sum(intensity.13C * (Nio.13C - Nio.w.13C)^2, na.rm = TRUE) /
                               sum(intensity.13C, na.rm = TRUE),
      
      # Effective sample size
      N_eff_0U = (sum(intensity.0U)^2) / sum(intensity.0U^2),
      N_eff_13C = (sum(intensity.13C)^2) / sum(intensity.13C^2),
      
      # Standard error
      se.Nio.w.0U = sqrt(var.Nio.w.0U / N_eff_0U),
      se.Nio.w.13C = sqrt(var.Nio.w.13C / N_eff_13C),
      
      # Relative standard error 
      rse.Nio.w.0U = se.Nio.w.0U / Nio.w.0U * 1000,
      rse.Nio.w.13C = se.Nio.w.13C / Nio.w.13C * 1000,
      
      # Total and per carbon fractional abundance
      ro = Nio.w.13C / (Nio.w.0U + Nio.w.13C),
      R13C = ro / n.C,
      
      # Propagated error for fractional abundance
      ro_variance = (Nio.w.0U^2 * se.Nio.w.13C^2 + Nio.w.13C^2 * se.Nio.w.0U^2) / 
        (Nio.w.0U + Nio.w.13C)^4,
      ro_se = sqrt(ro_variance), 
      R13C_se = ro_se / n.C,
      R13C_rse = R13C_se / R13C * 1000,
    ) %>%
    ungroup()%>%
    distinct()
  
  return(data_prep)
}
#RSE vs effective number of counts. Shot-noise limit and measurement
calc.shotnoise <- function(d.neff){
  df.out <- data.frame(scans =1:length(d.neff$scan.no),
                       neff = 1:length(d.neff$scan.no),
                       rse= 1:length(d.neff$scan.no),
                       shotnoise = 1:length(d.neff$scan.no),
                       shotnoise.neff = 1:length(d.neff$scan.no))
  #Geometric mean, geometric SD and SE
  gmean <- function(x) exp(mean(log(x))) #define geometric mean
  gsd <- function(x) exp(mean(log(x)) + sd(log(x))) - exp(mean(log(x)))
  gse <- function(x) (exp(mean(log(x)) + sd(log(x))) - exp(mean(log(x)))) /sqrt(length(x))
  
  for (l in 1:length(d.neff$scan.no)) {
    df <- d.neff[1:l,]
    neff <- (sum(df$Nio_0U) * sum(df$Nio_13C))/ (sum(df$Nio_0U) + sum(df$Nio_13C))
    df.out$neff[l] <- neff
    shotnoise.neff <- neff^-0.5
    df.out$shotnoise.neff[l] <- shotnoise.neff
    shotnoise <- ((1/sum(df$Nio_13C) + (1/sum(df$Nio_0U))))^0.5
    df.out$shotnoise[l] <- shotnoise
    rse <- gse(d.neff[1:l,]$Nio.C.ratio) / gmean(d.neff[1:l,]$Nio.C.ratio)
    df.out$rse[l] <- rse
  }
  df.out <- df.out %>% na.omit() %>% select(-shotnoise, -shotnoise.neff)
  df.out <- df.out %>% gather(key, value, -scans, -neff)
  return(df.out)
}

#Data import----
I.mixed.std.1 = I_df(mixed.Std_Phytane_1,"mixed.Std_Phytane_1","syn","24/07/2019",masses,ions,ilogs)
I.mixed.std.2 = I_df(mixed.Std_Phytane_2,"mixed.Std_Phytane_2","syn","25/07/2019",masses,ions,ilogs)
I.mixed.std.3 = I_df(mixed.Std_Phytane_3,"mixed.Std_Phytane_3","syn","26/07/2019",masses,ions,ilogs)

I.mixed.halo.1 = I_df(mixed.Halo_Phytane_1,"mixed.Halo_Phytane_1","MVA","24/07/2019",masses,ions,ilogs)
I.mixed.halo.2 = I_df(mixed.Halo_Phytane_2,"mixed.Halo_Phytane_2","MVA","25/07/2019",masses,ions,ilogs)
I.mixed.halo.3 = I_df(mixed.Halo_Phytane_3,"mixed.Halo_Phytane_3","MVA","26/07/2019",masses,ions,ilogs)

I.mixed.cyano.1 = I_df(mixed.Cyano_Phytane_1,"mixed.Cyano_Phytane_1","MEP","24/07/2019",masses,ions,ilogs)
I.mixed.cyano.2 = I_df(mixed.Cyano_Phytane_2,"mixed.Cyano_Phytane_2","MEP","25/07/2019",masses,ions,ilogs)
I.mixed.cyano.3 = I_df(mixed.Cyano_Phytane_3,"mixed.Cyano_Phytane_3","MEP","26/07/2019",masses,ions,ilogs)

I.all= rbind(I.mixed.std.1,I.mixed.std.2,I.mixed.std.3,
             I.mixed.halo.1,I.mixed.halo.2,I.mixed.halo.3,
             I.mixed.cyano.1,I.mixed.cyano.2,I.mixed.cyano.3)
#
#Background data frame definition----
start.min=9.6
b.mixed.std.1 = b_df(I.mixed.std.1,start.min)
b.mixed.std.2 = b_df(I.mixed.std.2,start.min)
b.mixed.std.3 = b_df(I.mixed.std.3,start.min)

b.mixed.halo.1 = b_df(I.mixed.halo.1,start.min)
b.mixed.halo.2 = b_df(I.mixed.halo.2,start.min)
b.mixed.halo.3 = b_df(I.mixed.halo.3,start.min)

b.mixed.cyano.1 = b_df(I.mixed.cyano.1,start.min)
b.mixed.cyano.2 = b_df(I.mixed.cyano.2,start.min)
b.mixed.cyano.3 = b_df(I.mixed.cyano.3,start.min)

b.all= rbind(b.mixed.std.1,b.mixed.std.2,b.mixed.std.3,
             b.mixed.halo.1,b.mixed.halo.2,b.mixed.halo.3,
             b.mixed.cyano.1,b.mixed.cyano.2,b.mixed.cyano.3)
#
#Sample dataframe----
outlier_factor=2
s.mixed.std.1 = s_df(I.mixed.std.1,b.mixed.std.1,outlier_factor)
s.mixed.std.2 = s_df(I.mixed.std.2,b.mixed.std.2,outlier_factor)
s.mixed.std.3 = s_df(I.mixed.std.3,b.mixed.std.3,outlier_factor)

s.mixed.halo.1 = s_df(I.mixed.halo.1,b.mixed.halo.1,outlier_factor)
s.mixed.halo.2 = s_df(I.mixed.halo.2,b.mixed.halo.2,outlier_factor)
s.mixed.halo.3 = s_df(I.mixed.halo.3,b.mixed.halo.3,outlier_factor)

s.mixed.cyano.1 = s_df(I.mixed.cyano.1,b.mixed.cyano.1,outlier_factor)
s.mixed.cyano.2 = s_df(I.mixed.cyano.2,b.mixed.cyano.2,outlier_factor)
s.mixed.cyano.3 = s_df(I.mixed.cyano.3,b.mixed.cyano.3,outlier_factor)
#
#Meta data ----
s.all=do.call(rbind, mget(ls(pattern="s.mixed.")))%>%
  mutate(type = "mixed")

# culled data
cull.data=s.all%>%
  select("filename","date","compound","isotopolog",
         "prcnt.below.10.max", "prcnt.below.background.scans", "prcnt.outlier.scans")%>%
  unique()

#calculate percent of missing isotopologues
missing.i.log=s.all%>%
  group_by(filename,compound)%>%
  summarise(
    prcnt.missing.U= (length(which(missing.0U == TRUE))) / length(unique(scan.no)) * 100,
    prcnt.missing.13C= (length(which(missing.13C == TRUE))) / length(unique(scan.no)) * 100,
    prcnt.missing.2H= (length(which(missing.2H == TRUE))) / length(unique(scan.no)) * 100,
  )
missing.i.log.filtered=s.all%>%
  filter(!below.ten.percent,!below.background,!outlier)%>%
  group_by(filename,compound)%>%
  summarise(
    prcnt.missing.U= (length(which(missing.0U == TRUE))) / length(unique(scan.no)) * 100,
    prcnt.missing.13C= (length(which(missing.13C == TRUE))) / length(unique(scan.no)) * 100,
    prcnt.missing.2H= (length(which(missing.2H == TRUE))) / length(unique(scan.no)) * 100,
  )

#data output quality check
q.all = s.all%>%
  select(filename, date, ID, compound, element, isotopolog, intensity, ions.incremental, tic, it.ms,
         prcnt.below.10.max, prcnt.below.background.scans, prcnt.outlier.scans, scan.no)%>%
  group_by(filename,date,ID,compound)%>%
  summarise(mean.intensity=mean(intensity),
            mean.ions.incremental=mean(ions.incremental),
            max.tic=max(tic),
            mean.it=mean(it.ms),
            n.inc.scans=length(unique(scan.no)),#n.total.scans minus scans satellite, below background and outliers
            prcnt.below.10.max=mean(prcnt.below.10.max),
            prcnt.below.background=mean(prcnt.below.background.scans),
            prcnt.outlier=mean(prcnt.outlier.scans))
write.csv(path = "output.tables",q.all,"quality.all.csv",row.names = F)
#
#compute Nio data -----
all.nio=nio_df(s.all)

# Compute rel.tic.int per
rel.tic <- all.nio %>%
  group_by(filename, ID, date, n.H, n.C) %>%
  filter(all(isotopolog %in% c("13C","0U","2H"))) %>%  # group must contain all isotopologues
  ungroup()%>%
  group_by(filename, ID, date, n.H, n.C, isotopolog) %>%
  summarise(rel.tic.int = sum(intensity)/sum(tic)*100, .groups = "drop")%>%
  ungroup()%>%
  group_by(ID, n.H, n.C, isotopolog) %>%
  mutate(m.rel.tic.int=mean(rel.tic.int),
         sd.rel.tic.int=sd(rel.tic.int))%>%
  ungroup()%>%
  group_by(filename, ID, date, n.H, n.C)%>%
  mutate(rel.C6H13.int=sum(rel.tic.int))%>%
  ungroup()%>%
  group_by(ID, n.H, n.C, isotopolog) %>%
  mutate(m.rel.C6H13.int=mean(rel.C6H13.int),
         sd.rel.C6H13.int=sd(rel.C6H13.int))%>%
  ungroup()

all.nio.m <- all.nio %>%
  select(-c("type","isotopolog","n.total.scans","element","intensity","b.int","ions.incremental","peakResolution","peakNoise","Nio"))%>%
  distinct()%>%
  group_by(filename, ID, date, n.H, n.C) %>%
  mutate(
    neff = (sum(Nio_13C) * sum(Nio_0U)) / (sum(Nio_13C) + sum(Nio_0U))
    )%>%
  ungroup()%>%
  group_by(filename, ID, date, n.H, n.C) %>%
  summarise(
    n.scans=length(scan.no),
    max.TIC=max(tic),
    m.Nio_0U = mean(Nio_0U),
    sd.Nio_0U = sd(Nio_0U),
    se.Nio_0U = sd(Nio_0U)/sqrt(n()),
    rse.Nio_0U = se.Nio_0U/m.Nio_0U*1000,
    m.Nio_13C = mean(Nio_13C),
    sd.Nio_13C = sd(Nio_13C),
    se.Nio_13C = sd(Nio_13C)/sqrt(n()),
    rse.Nio_13C = se.Nio_13C/m.Nio_13C*1000,
    m.ratio = mean(Nio.C.ratio),
    sd.ratio = sd(Nio.C.ratio),
    se.ratio = sd(Nio.C.ratio)/sqrt(n()),
    rse.ratio = se.ratio/m.ratio*1000,
    m.neff=max(neff),
    SNL=max(neff)^-0.5*1000,
    rse.snl=rse.ratio/SNL,
    .groups = "drop"
  ) %>%
  ungroup()

# Table 1 output ----
out.Nio.ratio= all.nio.m %>%
  select("filename", "ID", "date","n.scans","max.TIC",
         "m.Nio_0U","sd.Nio_0U","se.Nio_0U","rse.Nio_0U",
         "m.Nio_13C","sd.Nio_13C","se.Nio_13C","rse.Nio_13C",
         "m.ratio","sd.ratio","se.ratio","rse.ratio",
         "SNL","rse.snl")
write.csv(path = "output.tables", out.Nio.ratio, "Table1.csv", row.names = F)
#quality check output data-----
quality_figure <- function(x, y, file, pname) {
  my_colors1 <- c("orange", "darkorange", "darkorange3")
  
  # Calculate the 10th percentile used for filtering
  tic_data <- x %>% filter(filename == file)
  tenth_percentile_cutoff <- unique(tic_data$tenth_percentile_tic)[1]  # Get the 10th percentile value
  
  # Figure 1: Signal intensities 
  q_fig.1 <- ggplot(x %>% filter(filename == file)) +
    # All TIC points (blue)
    geom_point(
      data = subset(x, filename == file),
      aes(x = time.min, y = tic, colour = "TIC"),
      size = 3, alpha = 0.6, shape = 1
    ) +
    # below.background points (darkred) - include in legend
    geom_point(
      data = subset(x, filename == file & below.background == TRUE),
      aes(x = time.min, y = tic, colour = "below.background"),
      size = 3, alpha = 0.6
    ) +
    # below.ten.percent points (red) - include in legend
    geom_point(
      data = subset(x, filename == file & below.ten.percent == TRUE),
      aes(x = time.min, y = tic, colour = "below.ten.percent"),
      size = 3, alpha = 0.6
    ) +
    # Add the 10th percentile cutoff line 
    geom_hline(yintercept = tenth_percentile_cutoff, 
               color = "red", linetype = "dashed", size = 1) +
    scale_colour_manual(
      name = "Signal type",
      values = c(
        "TIC" = "blue",
        "below.background" = "darkred",
        "below.ten.percent" = "red"
      ),
      labels = c(
        "TIC" = "TIC",
        "below.ten.percent" = "< 10th percentile signal",
        "below.background" = "< background signal"
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      text = element_text(size = 22),
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", colour = "black", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA)  # White background
    ) +
    ylab("Signal (S)") +
    xlab("Time (min)")
  
  # Figure 2: TIC x IT
  q_fig.2 <- ggplot(x %>%
                      filter(filename %in% c(file)) %>%
                      filter(below.background == FALSE & below.ten.percent == FALSE),
                    aes(x = time.min, y = tic * it.ms)) +
    geom_point(size = 3, alpha = 0.6, shape = 1, colour = 'blue') +
    geom_point(data = x %>%
                 filter(filename %in% c(file)) %>%
                 filter(below.background == FALSE & below.ten.percent == FALSE & outlier == TRUE),
               size = 3, alpha = 0.6, colour = 'red') +
    geom_hline(aes(yintercept = medianTICxIT), color = "black", size = 1) +
    geom_hline(aes(yintercept = upsd), color = "red", linetype = "dashed") +
    geom_hline(aes(yintercept = lowsd), color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 22),
          plot.background = element_rect(fill = "white", color = NA)) +  # White background
    ylab("TICxIT") +
    xlab("Time (min)")
  
  # Calculate statistics for the isotope ratio plot
  y_data <- y %>% filter(filename %in% c(file)) %>% pull(Nio.C.ratio)
  median_val <- median(y_data, na.rm = TRUE)
  sd_val <- sd(y_data, na.rm = TRUE)
  
  # Create main isotope ratio plot
  main_plot <- ggplot(y %>% filter(filename %in% c(file)), 
                      aes(x = time.min, y = Nio.C.ratio)) +
    geom_point(size = 3, alpha = 0.6, shape = 1, colour = "blue") +
    geom_ma(ma_fun = SMA, n = 228, color = "red", size = 3) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 22),
          plot.background = element_rect(fill = "white", color = NA)) +  # White background
    labs(y = bquote({""^13~C~italic((N[io]))}/{""^12~C~italic((N[io]))}),
         x = bquote("Time (min)"))
  
  # Create marginal histogram with proper alignment
  marginal_hist <- ggplot(y %>% filter(filename %in% c(file)), aes(x = Nio.C.ratio)) +
    geom_histogram(aes(y = after_stat(density)), 
                   fill = "blue", 
                   color = "black",
                   size = 0.3,
                   alpha = 0.6, 
                   bins = 20) +
    geom_vline(xintercept = median_val, color = "black", linetype = "solid", size = 1) +
    geom_vline(xintercept = c(median_val - sd_val, median_val + sd_val), 
               color = "red", linetype = "dashed", size = 0.8) +
    geom_vline(xintercept = c(median_val - 2*sd_val, median_val + 2*sd_val), 
               color = "red", linetype = "dashed", size = 0.8) +
    geom_vline(xintercept = c(median_val - 3*sd_val, median_val + 3*sd_val), 
               color = "red", linetype = "dashed", size = 0.8) +
    coord_flip() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", color = NA))  # White background
  
  # Combine main plot and marginal histogram using grid.arrange
  library(gridExtra)
  q_fig.4_combined <- arrangeGrob(main_plot, marginal_hist, 
                                  ncol = 2, widths = c(4, 1))
  
  # Shot noise plot
  data.SN <- y %>% filter(filename %in% c(file)) %>%
    do(calc.shotnoise(.))
  q_fig.4_c <- ggplot(data = data.SN, aes(x = neff, y = value)) + 
    geom_point(size = 3, alpha = 0.5, colour = 'blue') +
    geom_abline(intercept = 0, slope = -0.5, col = "black") +
    scale_color_manual(values = my_colors1) +
    labs(y = "RSE", x = "effective number of ions") +
    annotate(geom = 'text', size = 8, label = 'shot noise \nlimit', 
             color = "black", x = 0.005 * median(data.SN$neff), y = 2 * median(data.SN$value)) + 
    theme(rect = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "white", color = NA),  # White background
          panel.border = element_blank(),
          legend.position = "none",
          text = element_text(size = 22)) +
    scale_x_continuous(trans = log10_trans(), 
                       breaks = c(10000, 100000, 1000000, 10000000, 100000000), 
                       labels = c("1E4", "1E5", "1E6", "1E7", "1E8")) +
    scale_y_continuous(trans = log10_trans(), 
                       breaks = c(0.05, 0.01, 0.002, 0.001, 0.0005, 0.0002, 0.0001), 
                       labels = c("50 permil", "10 permil", "2 permil", "1 permil", 
                                  "0.5 permil", "0.2 permil", "0.1 permil"))
  
  # Final figure assembly using grid.arrange with proper labels
  # Add labels to individual plots
  q_fig.1_labeled <- q_fig.1 + labs(tag = "A") + theme(plot.tag = element_text(size = 24, face = "bold"))
  q_fig.2_labeled <- q_fig.2 + labs(tag = "B") + theme(plot.tag = element_text(size = 24, face = "bold"))
  
  # Create the layout
  top_row <- plot_grid(q_fig.1_labeled, q_fig.2_labeled, ncol = 2, align = 'h', rel_widths = c(1, 1))
  
  # For the bottom row
  bottom_left <- ggplot() + 
    annotation_custom(
      grob = arrangeGrob(q_fig.4_combined),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    labs(tag = "C") +
    theme(plot.tag = element_text(size = 24, face = "bold"),
          plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", color = NA))  # White background
  
  bottom_right <- q_fig.4_c + labs(tag = "D") + theme(plot.tag = element_text(size = 24, face = "bold"))
  bottom_row <- plot_grid(bottom_left, bottom_right, ncol = 2, align = 'h', rel_widths = c(2, 1))
  
  # Combine everything
  q_figure <- plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1, 1))
  
  # Save the plot with white background
  ggsave(path = "qual.figs", filename = paste0(pname, ".pdf"), 
         plot = q_figure, width = 20, height = 15, bg = "white")  # Set background to white
}

quality_figure(s.all, all.nio, "mixed.Std_Phytane_1", "Std_Phytane_1.test")
quality_figure(s.all,all.nio,"mixed.Std_Phytane_1" ,"Std_Phytane_1")
quality_figure(s.all,all.nio,"mixed.Std_Phytane_2" ,"Std_Phytane_2")
quality_figure(s.all,all.nio,"mixed.Std_Phytane_3" ,"Std_Phytane_3")
quality_figure(s.all,all.nio,"mixed.Cyano_Phytane_1" ,"Cyano_Phytane_1")
quality_figure(s.all,all.nio,"mixed.Cyano_Phytane_2" ,"Cyano_Phytane_2")
quality_figure(s.all,all.nio,"mixed.Cyano_Phytane_3" ,"Cyano_Phytane_3")
quality_figure(s.all,all.nio,"mixed.Halo_Phytane_1" ,"Halo_Phytane_1")
quality_figure(s.all,all.nio,"mixed.Halo_Phytane_2" ,"Halo_Phytane_2")
quality_figure(s.all,all.nio,"mixed.Halo_Phytane_3" ,"Halo_Phytane_3")
#
#summary-------
frac.output.weighted.all=frac.output.weighted.avg.ratio(all.nio)
output.nio=all.nio%>%
  select(-c("type","isotopolog","n.total.scans","element","intensity","b.int","ions.incremental","peakResolution","peakNoise","Nio"))%>%
  distinct()%>%
  left_join(frac.output.weighted.all, by = c("filename","ID","date", "compound", "n.H", "n.C"))
#reproducibility tests-----
ggplot(output.nio,aes(x=ID,y=R13C,colour=filename))+
  geom_point(size=3, position = position_dodge(width = 0.90))+
  geom_errorbar(aes(ymin=R13C-R13C_se, ymax=R13C+R13C_se), width=.3,
                position=position_dodge(0.90))+
  theme_minimal()+
  theme(legend.position="none",text = element_text(size = 15))+
  ggtitle("C6H13")+
  ylab("weighted 13R")

#RSE vs effective number of counts. Shot-noise limit and measurement-----
#make the data frame:
output.m.syn.C6H13=output.nio%>%
  filter(ID %in% c("syn"))
output.m.MEP.C6H13=output.nio%>%
  filter(ID %in% c("MEP"))
output.m.MVA.C6H13=output.nio%>%
  filter(ID %in% c("MVA"))

data.SN.syn.C6H13 <-output.m.syn.C6H13%>%
  group_by(date)%>%
  do(calc.shotnoise(.))
data.SN.MEP.C6H13 <- output.m.MEP.C6H13 %>%
  group_by(date)%>%
  do(calc.shotnoise(.))
data.SN.MVA.C6H13 <- output.m.MVA.C6H13 %>%
  group_by(date)%>%
  do(calc.shotnoise(.))

my_colors1 <- c("orange","darkorange","darkorange3")
my_colors3 <- c("mediumpurple1","mediumpurple3","mediumpurple4")
my_colors2 <- c("#74AC64","#4C7434","#053E21")
#syn#####
sp.syn.C6H13 <- ggplot(output.m.syn.C6H13, aes(x = time.min , y = Nio.C.ratio,color = as.factor(date))) +
  geom_point(size = 3, alpha = 0.6, shape=1)+
  scale_color_manual(values=my_colors1)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  geom_ma(ma_fun = SMA, n = 100, color = "red")+ #100scans moving average
  labs( x ="Time (min)", y = expression("13R"))+ 
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.margin = margin())+
  theme(text = element_text(size = 20))     
# Define marginal histogram
hist_syn.C6H13 <- ggMarginal(sp.syn.C6H13 , groupColour=TRUE,groupFill=TRUE, margins="y")
#align with boxplot
#boxplot with error bar on distribution
b_syn.C6H13=ggplot(output.m.syn.C6H13, aes(x = date, y = Nio.C.ratio,color = as.factor(date))) + 
  geom_boxplot(outlier.shape=NA, size=1, show.legend = FALSE)+
  scale_color_manual(values=my_colors1)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  theme(plot.margin = margin(),panel.grid.major=element_blank(),panel.background=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  scale_x_discrete(limits = rev(levels(output.m.syn.C6H13$date)))
a1.syn.C6H13 <- plot_grid(hist_syn.C6H13,b_syn.C6H13,ncol=2,rel_widths = c(1,0.1), rel_heights = c(1,1), align = "h",axis="bt")
#make the plot
SN.syn.C6H13 <- ggplot(data = data.SN.syn.C6H13, aes(x= neff, y=  value, color=date)) + 
  geom_point(size=3, alpha=0.5) +
  geom_abline(intercept = 0, slope = -0.5, col="blue") +
  scale_color_manual(values=my_colors1)+
  labs(y="RSE", x="effective number of ions")  +
  annotate(geom = 'text', size= 6, label = 'shot noise \nlimit', color="blue", x=0.001*median(data.SN.syn.C6H13$neff),y=2*median(data.SN.syn.C6H13$value))+ 
  theme(rect = element_rect(fill = "transparent"),plot.background = element_rect(color = "white"),
        legend.position="none",text = element_text(size = 20))+
  scale_x_continuous(trans = log10_trans(), breaks = c(10000, 100000, 1000000,10000000, 100000000), 
                     labels=c("1E4","1E5", "1E6", "1E7", "1E8")) +
  scale_y_continuous(trans = log10_trans(), breaks = c( 0.05, 0.01,  0.002, 0.001, 0.0005,0.0002, 0.0001), 
                     labels=c("50 permil", "10 permil", "2 permil" , "1 permil","0.5 permil", "0.2 permil", "0.1 permil"))
Std_plot.C6H13=plot_grid(a1.syn.C6H13
                   , SN.syn.C6H13 
                   , align = "v"
                   , ncol = 2
                   , rel_widths = c(1.2,0.5)
)
#MEP#####
sp.MEP.C6H13 <- ggplot(output.m.MEP.C6H13, aes(x = time.min , y = Nio.C.ratio,color = as.factor(date))) +
  geom_point(size = 3, alpha = 0.6, shape=1)+
  scale_color_manual(values=my_colors2)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  geom_ma(ma_fun = SMA, n = 100, color = "red")+ #100scans moving average
  labs( x ="Time (min)", y = expression("13R"))+ 
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.margin = margin())+
  theme(text = element_text(size = 20))     
# Define marginal histogram
hist_MEP.C6H13 <- ggMarginal(sp.MEP.C6H13 , groupColour=TRUE,groupFill=TRUE, margins="y")
#align with boxplot
#boxplot with error bar on distribution
b_MEP.C6H13=ggplot(output.m.MEP.C6H13, aes(x = date, y = Nio.C.ratio,color = as.factor(date))) + 
  geom_boxplot(outlier.shape=NA, size=1, show.legend = FALSE)+
  scale_color_manual(values=my_colors2)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  theme(plot.margin = margin(),panel.grid.major=element_blank(),panel.background=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  scale_x_discrete(limits = rev(levels(output.m.MEP.C6H13$date)))
a1.MEP.C6H13 <- plot_grid(hist_MEP.C6H13,b_MEP.C6H13,ncol=2,rel_widths = c(1,0.1), rel_heights = c(1,1), align = "h",axis="bt")
#make the plot
SN.MEP.C6H13 <- ggplot(data = data.SN.MEP.C6H13, aes(x= neff, y=  value, color=date)) + 
  geom_point(size=3, alpha=0.5) +
  geom_abline(intercept = 0, slope = -0.5, col="blue") +
  scale_color_manual(values=my_colors2)+
  labs(y="RSE", x="effective number of ions")  +
  annotate(geom = 'text', size= 6, label = 'shot noise \nlimit', color="blue", x=0.001*median(data.SN.MEP.C6H13$neff),y=2*median(data.SN.MEP.C6H13$value))+ 
  theme(rect = element_rect(fill = "transparent"),plot.background = element_rect(color = "white"),
        legend.position="none",text = element_text(size = 20))+
  scale_x_continuous(trans = log10_trans(), breaks = c(10000, 100000, 1000000,10000000, 100000000), 
                     labels=c("1E4","1E5", "1E6", "1E7", "1E8")) +
  scale_y_continuous(trans = log10_trans(), breaks = c( 0.05, 0.01,  0.002, 0.001, 0.0005,0.0002, 0.0001), 
                     labels=c("50 permil", "10 permil", "2 permil" , "1 permil","0.5 permil", "0.2 permil", "0.1 permil"))
MEP_plot.C6H13=plot_grid(a1.MEP.C6H13
                   , SN.MEP.C6H13 
                   , align = "v"
                   , ncol = 2
                   , rel_widths = c(1.2,0.5)
)
#MVA#####
sp.MVA.C6H13 <- ggplot(output.m.MVA.C6H13, aes(x = time.min , y = Nio.C.ratio,color = as.factor(date))) +
  geom_point(size = 3, alpha = 0.6, shape=1)+
  scale_color_manual(values=my_colors3)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  geom_ma(ma_fun = SMA, n = 100, color = "red")+ #100scans moving average
  labs( x ="Time (min)", y = expression("13R"))+ 
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.margin = margin())+
  theme(text = element_text(size = 20))     
# Define marginal histogram
hist_MVA.C6H13 <- ggMarginal(sp.MVA.C6H13 , groupColour=TRUE,groupFill=TRUE, margins="y")
#align with boxplot
#boxplot with error bar on distribution
b_MVA.C6H13=ggplot(output.m.MVA.C6H13, aes(x = date, y = Nio.C.ratio,color = as.factor(date))) + 
  geom_boxplot(outlier.shape=NA, size=1, show.legend = FALSE)+
  scale_color_manual(values=my_colors3)+
  scale_y_continuous(limits = c(0.045,0.09), breaks = seq(0.045,0.085,by=0.01))+
  theme(plot.margin = margin(),panel.grid.major=element_blank(),panel.background=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  scale_x_discrete(limits = rev(levels(output.m.MVA.C6H13$date)))
a1.MVA.C6H13 <- plot_grid(hist_MVA.C6H13,b_MVA.C6H13,ncol=2,rel_widths = c(1,0.1), rel_heights = c(1,1), align = "h",axis="bt")
#make the plot
SN.MVA.C6H13 <- ggplot(data = data.SN.MVA.C6H13, aes(x= neff, y=  value, color=date)) + 
  geom_point(size=3, alpha=0.5) +
  geom_abline(intercept = 0, slope = -0.5, col="blue") +
  scale_color_manual(values=my_colors3)+
  labs(y="RSE", x="effective number of ions")  +
  annotate(geom = 'text', size= 6, label = 'shot noise \nlimit', color="blue", x=0.001*median(data.SN.MVA.C6H13$neff),y=2*median(data.SN.MVA.C6H13$value))+ 
  theme(rect = element_rect(fill = "transparent"),plot.background = element_rect(color = "white"),
        legend.position="none",text = element_text(size = 20))+
  scale_x_continuous(trans = log10_trans(), breaks = c(10000, 100000, 1000000,10000000, 100000000), 
                     labels=c("1E4","1E5", "1E6", "1E7", "1E8")) +
  scale_y_continuous(trans = log10_trans(), breaks = c( 0.05, 0.01,  0.002, 0.001, 0.0005,0.0002, 0.0001), 
                     labels=c("50 permil", "10 permil", "2 permil" , "1 permil","0.5 permil", "0.2 permil", "0.1 permil"))
MVA_plot.C6H13=plot_grid(a1.MVA.C6H13
                   , SN.MVA.C6H13 
                   , align = "v"
                   , ncol = 2
                   , rel_widths = c(1.2,0.5)
)
#
#all####
title <- textGrob('C6H13', gp = gpar(fontsize = 22, fontface = 'bold'))
STD <- textGrob('STD', gp = gpar(fontsize = 20, fontface = 'bold'))
MEP <- textGrob('MEP', gp = gpar(fontsize = 20, fontface = 'bold'))
MVA <- textGrob('MVA', gp = gpar(fontsize = 20, fontface = 'bold'))
C6H13=plot_grid(title,
          STD,
          Std_plot.C6H13,
          MEP,
          MEP_plot.C6H13,
          MVA,
          MVA_plot.C6H13,
          nrow=7,
          rel_heights = c(0.3,0.1,1.5,0.1,1.5,0.1,1.5))
C6H13
#
#----save----
C6H13
ggsave(filename="C6H13.png",plot=C6H13, width = 17.5, height=10)
# Sample-STD---- 
view(frac.output.weighted.all)
sample=frac.output.weighted.all%>%
  filter(ID %in% c("MEP","MVA"))%>%
  select(date,ID,R13C,R13C_se,R13C_rse)%>%
  mutate(IRMS=case_when(ID == "MEP"~-43.679,
                        ID == "MVA"~-24.675),
         IRMS_se=case_when(ID == "MEP"~0.034908,
                           ID == "MVA"~0.046263),
         IRMS_rse=abs(IRMS_se/IRMS)*1000)
STD=frac.output.weighted.all%>%
  filter(ID == "syn")%>%
  select(date,ID,R13C,R13C_se,R13C_rse)%>%
  mutate(IRMS=case_when(ID == "syn"~-30.017),
         IRMS_se=case_when(ID == "syn"~0.015690),
         IRMS_rse=abs(IRMS_se/IRMS)*1000)%>%
  pivot_wider(names_from = c(ID), values_from = c(R13C,R13C_se,R13C_rse,IRMS,IRMS_se,IRMS_rse))

STDsample = left_join(sample, STD, by = "date") %>%
  mutate(
    IRMS_13R = (((IRMS/1000)+1)*0.01123372),
    IRMS_13R_syn = (((IRMS_syn/1000)+1)*0.01123372),
    offset = ((R13C/6)-IRMS_13R)*1000,
    
    sample_STD.IRMS = IRMS - IRMS_syn,
    sample_STD.IRMS.err.prop = sqrt((IRMS^2 * IRMS_se^2 + IRMS_syn^2 * IRMS_se_syn^2)/
                                      (IRMS+IRMS_syn)^4),
    sample_STD.IRMS.RSE = abs(sample_STD.IRMS.err.prop / sample_STD.IRMS) * 1000,  # in permil
    
    # Delta value in per mil
    sample_STD.ratio = ((R13C / R13C_syn)-1)*1000,# in permil
    # Error propagation for delta 
    sample_STD.ratio.err.prop =  sqrt((R13C^2 * R13C_se^2 + R13C_syn^2 * R13C_se_syn^2)/
                                        (R13C+R13C_syn)^4)*1000,# in permil
    sample_STD.ratio.RSE = abs(sample_STD.ratio.err.prop / sample_STD.ratio) ,
    delt.sample.STD.IRMS = sample_STD.ratio - sample_STD.IRMS
  )

# table 2 output ----
out.STDsample= STDsample%>%
  select("date","ID",
         "IRMS","IRMS_se","IRMS_syn","IRMS_se_syn",
         "R13C","R13C_se","R13C_rse",
         "R13C_syn","R13C_se_syn","R13C_rse_syn",
         "sample_STD.ratio","sample_STD.ratio.err.prop","sample_STD.ratio.RSE")
write.csv(path = "output.tables",out.STDsample, "Table2.csv", row.names = F)

#syn vs. sample plot-----
ggplot(STDsample, aes(y=sample_STD.ratio,x=sample_STD.IRMS, colour=ID, shape=date))+
  geom_point(size=3)+
  scale_x_reverse()+
  geom_errorbar(aes(ymin = sample_STD.ratio - sample_STD.ratio.err.prop,
                    ymax = sample_STD.ratio + sample_STD.ratio.err.prop),
                width=0.3)+
  geom_errorbarh(aes(xmin = sample_STD.IRMS - sample_STD.IRMS.err.prop,
                     xmax = sample_STD.IRMS + sample_STD.IRMS.err.prop),
                 height=0.2)+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  labs(x = expression(delta^{13}*C[paste("VPDB"," ", "sample-STD")]),
       y = expression({}^{13}*R[sample-STD]))+
  theme_bw()

ggplot(STDsample, aes(y=delt.sample.STD.IRMS,x=sample_STD.IRMS, colour=ID, shape=date))+
  geom_point(size=3)+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  labs(x = expression("sample_STD.IRMS"),
       y = expression("delt.sample.STD.IRMS"))+
  theme_bw()

dat=frac.output.weighted.all%>%
  select(date,ID,R13C,R13C_se,R13C_rse)%>%
  mutate(IRMS=case_when(ID == "MEP"~-43.679,
                        ID == "MVA"~-24.675,
                        ID == "syn"~-30.017),
         IRMS_se=case_when(ID == "MEP"~0.034908,
                           ID == "MVA"~0.046263,
                           ID == "syn"~0.015690))
ggplot(dat, aes(y=R13C*1000,x=IRMS, colour=ID, shape=date))+
  geom_point(size=3)+
  scale_x_reverse()+
  geom_errorbar(aes(ymin = (R13C*1000) - (R13C_se*1000),
                    ymax = (R13C*1000) + (R13C_se*1000)),
                width=0.3)+
  geom_errorbarh(aes(xmin = IRMS - IRMS_se,
                     xmax = IRMS + IRMS_se),
                 height=0.2)+
  scale_y_continuous() +
  labs(x = expression(delta^{13}*C[paste("VPDB")]),
       y = expression(bar({}^{13}*R[C6H13])))+
  theme_bw()
#
fig4=ggplot(STDsample, aes(y=sample_STD.ratio,x=sample_STD.IRMS, colour=ID))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin = sample_STD.ratio - sample_STD.ratio.err.prop,
                    ymax = sample_STD.ratio + sample_STD.ratio.err.prop),
                width=1,
                colour="black")+
  geom_errorbarh(aes(xmin = sample_STD.IRMS - sample_STD.IRMS.err.prop,
                     xmax = sample_STD.IRMS + sample_STD.IRMS.err.prop),
                 height=0.1,
                 colour="black")+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  labs(x =expression("mol.avrg"[sample-STD]))+
  labs(y =expression("C6H13"[sample-STD]))+
  scale_x_continuous(limits=c(-17,7)) +
  scale_y_continuous(limits=c(-17,7)) +
  theme_bw()+
  theme(text = element_text(size=15))
fig4

ggsave(filename="fig4.pdf",plot=fig4, width = 5, height=5)
