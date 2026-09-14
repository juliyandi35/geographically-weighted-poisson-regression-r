# Library
library(GWmodel)      ### GW models
library(sp)           ## Data management
library(spdep)        ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(raster)       ## spatial data
library(grid)         # plot
library(gridExtra)    # Multiple plot
library(ggplot2)      # Multiple plot
library(gtable)
library(readxl)       # Read excel file
library(energy)       # Poisson Distribution Test
library(ape)          # Moran I test
library(dplyr)        # Dataframe mutation
library(MASS)         # Fit Distribution

# Import dataset
Dataset <- read_excel("Dataset.xlsx")

# Pre-processing
colSums(is.na(Dataset))

# Analisis deskriptif
# Keragaman antar individu
summary(Dataset)

# Scale co-variates
Dataset[, 3:5] = scale(Dataset[, 3:5])

# Poisson Distribution test
# By Histogram
hist(Dataset$Covid19, col = 'lightblue', main = 'Histogram', 
     xlab = 'Number of Events', ylab = 'Frequency')

# By Test
distFit<-fitdistr(Dataset$Covid19,"Poisson") # Mendapatkan lambda yang sesuai
ks.test(Dataset$Covid19, "ppois", lambda = distFit$estimate) # Kolmogorov-Smirnov for Poisson Distribution

# Multicolinierity test
library(car)
vif(lm(Covid19~IE + BCR + SPU,data=Dataset))

# Plot peta
Ind_map = read_sf('BATAS WILAYAH KELURAHAN-DESA 10K/Batas_Wilayah_KelurahanDesa_10K_AR.shp')
map <- subset(Ind_map, WADMPR %in% "Daerah Istimewa Yogyakarta")
map <- subset(map, NAMOBJ %in% Dataset$Kelurahan)
map <- map[,-c(1,3:26)]
colnames(map) <- c("Kelurahan","geometry")
map <- map[-46,]
map

# Membuat plot peta
ggplot() +
  geom_sf(data = map) +
  labs(title = "Peta Kelurahan di DI Yogyakarta")

Dataset <- merge(map,Dataset,by="Kelurahan")
Dataset <- st_zm(Dataset)
Dataset

# moran test untuk melihat heterogenitas dan autokorelasi spasial
# Buat matriks ketetanggaan (neighbor matrix)
Longitude <- coordinates(as(Dataset,"Spatial"))[,1]
Latitude <- coordinates(as(Dataset,"Spatial"))[,2]
dists <- as.matrix(dist(cbind(Longitude, Latitude)))

dists.inv <- 1/dists
diag(dists.inv) <- 0

# Hitung Moran's I
Moran.I(Dataset$Covid19,dists.inv)

# Bandwidth selection
Dataset <- as(Dataset,"Spatial")
DM<-gw.dist(dp.locat=coordinates(Dataset))
bw.gwr <- bw.ggwr(Covid19 ~ IE+BCR+SPU,  
                  data = Dataset,
                  family = "poisson",
                  approach = "AICc",
                  kernel = "bisquare", 
                  adaptive = TRUE,
                  dMat = DM )
bw.gwr

# Fit the model
bgwr.res <- ggwr.basic(Covid19 ~ IE+BCR+SPU,  
                       data = Dataset,
                       family = "poisson",
                       bw = bw.gwr, 
                       kernel = "bisquare", 
                       adaptive = TRUE,
                       dMat = DM)
bgwr.res

# Extract GWPR results
### Create spatial data frame

Dataset@data$y<-bgwr.res$SDF$y
Dataset@data$yhat<-bgwr.res$SDF$yhat
Dataset@data$residual<-bgwr.res$SDF$residual
rsd=sd(Dataset@data$residual)
Dataset@data$stdRes<-(Dataset@data$residual)/sd(Dataset@data$residual)
Dataset@data$LLN=Dataset@data$yhat-1.645*rsd # Lower Limit of Normal Range
Dataset@data$ULN=Dataset@data$yhat+1.645*rsd # Upper Limit of Normal Range

# Intercept
Dataset@data$Intercept<-bgwr.res$SDF$Intercept
Dataset@data$est_IE<-bgwr.res$SDF$IE
Dataset@data$est_BCR<-bgwr.res$SDF$BCR
Dataset@data$est_SPU<-bgwr.res$SDF$SPU

# T-values
Dataset@data$t_Intercept<-bgwr.res$SDF$Intercept_TV
Dataset@data$t_IE<-bgwr.res$SDF$IE_TV
Dataset@data$t_BCR<-bgwr.res$SDF$BCR_TV
Dataset@data$t_SPU<-bgwr.res$SDF$SPU_TV

# Calculate psudo-t values
Dataset@data$p_IE<-2*pt(-abs(bgwr.res$SDF$IE_TV),df=44)
Dataset@data$p_BCR<-2*pt(-abs(bgwr.res$SDF$BCR_TV),df=44)
Dataset@data$p_SPU<-2*pt(-abs(bgwr.res$SDF$SPU_TV),df=44)

Dataset$sig_IE <-ifelse(Dataset@data$est_IE > 0 &
                           Dataset@data$p_IE <= 0.05 , "Significant", "Not Significant")
Dataset$sig_BCR <-ifelse(Dataset@data$est_BCR > 0 &
                          Dataset@data$p_BCR <= 0.05 , "Significant", "Not Significant")
Dataset$sig_SPU <-ifelse(Dataset@data$est_SPU > 0 &
                           Dataset@data$p_SPU <= 0.05 , "Significant", "Not Significant")

# Plot GWRP Statistics
polys<- list("sp.lines", as(as(st_zm(map),"Spatial"), "SpatialLines"), col="grey", lwd=.8,lty=1)
col.palette<-colorRampPalette(c("skyblue", "green","yellow"),space="rgb",interpolate = "linear")

# Plot Local Estimates
col.palette<-colorRampPalette(c("lightcyan","cyan","cyan1", "cyan2","cyan3","cyan4", "darkblue"),space="rgb",interpolate = "linear") 
est_IE<-spplot(Dataset,"est_IE", main = "Index Entropy", 
                 sp.layout=list(polys),
                 col="transparent",
                 col.regions=col.palette(100))

est_BCR<-spplot(Dataset,"est_BCR", main = "Building Coverage Ratio (%)", 
                sp.layout=list(polys),
                col="transparent",
                col.regions=col.palette(100))

est_SPU<-spplot(Dataset,"est_SPU", main = "Klasifikasi Sarana Pelayanan Umum", 
                 sp.layout=list(polys),
                 col="transparent",
                 col.regions=col.palette(100))
grid.arrange(est_IE, est_BCR,est_SPU,ncol= 3, heights = c(30,6), top = textGrob("Local Estimates",gp=gpar(fontsize=25)))

# Plot Local t-values
col.palette.t<-colorRampPalette(c("skyblue", "green","yellow"),space="rgb",interpolate = "linear") 

t_IE<-spplot(Dataset,"t_IE", main = "Index Entropy", 
               sp.layout=list(polys),
               col="transparent",
               col.regions=rev(col.palette.t(100)))

t_BCR<-spplot(Dataset,"t_BCR", main = "Building Coverage Ratio (%)", 
              sp.layout=list(polys),
              col="transparent",
              col.regions=rev(col.palette.t(100)))

t_SPU<-spplot(Dataset,"t_SPU", main = "Klasifikasi Sarana Pelayanan Umum", 
               sp.layout=list(polys),
               col="transparent",
               col.regions=rev(col.palette.t(100)))

grid.arrange(t_IE, t_BCR,t_SPU,ncol=3, heights = c(30,6), top = textGrob("Local t-values",gp=gpar(fontsize=25)))

# Significancy Plot
Dataset <- st_as_sf(Dataset)
ggplot(data=Dataset) +
  geom_sf(mapping=aes(fill =sig_IE)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="Signifikansi")+
  geom_text(
    aes(label = Kelurahan, x = coordinates(as(Dataset,"Spatial"))[,1], y = coordinates(as(Dataset,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi Index Entropy")+xlab("Longitude")+ylab("Latitude")


ggplot(data=Dataset) +
  geom_sf(mapping=aes(fill =sig_BCR)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="Signifikansi")+
  geom_text(
    aes(label = Kelurahan, x = coordinates(as(Dataset,"Spatial"))[,1], y = coordinates(as(Dataset,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi Building Coverage Ratio")+xlab("Longitude")+ylab("Latitude")

ggplot(data=Dataset) +
  geom_sf(mapping=aes(fill =sig_SPU)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="Signifikansi")+
  geom_text(
    aes(label = Kelurahan, x = coordinates(as(Dataset,"Spatial"))[,1], y = coordinates(as(Dataset,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi Klasifikasi Sarana Pelayanan Umum")+xlab("Longitude")+ylab("Latitude")

# Buat kolom baru untuk kombinasi signfikansi
Dataset <- Dataset %>%
  mutate(Variabel_Signifikan = case_when(
    # Utuh
    sig_IE == "Significant" & sig_BCR == "Significant" &
      sig_SPU == "Significant" ~ "IE,BCR,SPU",
    
    # Eliminasi 1
    sig_IE == "Significant" & sig_BCR == "Significant" ~ "IE,BCR",
    sig_IE == "Significant" & sig_SPU == "Significant" ~ "IE,SPU",
    sig_BCR == "Significant" & sig_SPU == "Significant" ~ "BCR,SPU",
    
    # Eliminasi 2
    sig_IE == "Significant" ~ "IE",
    sig_BCR == "Significant" ~ "BCR",
    sig_SPU == "Significant" ~ "SPU",
    
    TRUE ~ "Tidak Signifikan"
  ))

# Buat skema warna kustom
warna_custom <- c(
  "IE,BCR,SPU" = "lavender",
  "IE,BCR" = "tomato",
  "IE,SPU" = "skyblue",
  "BCR,SPU" = "aquamarine",
  "IE" = "cyan",
  "BCR" = "magenta",
  "SPU" = "darkgray",
  "Tidak Signifikan" = "white"
)

ggplot(data = Dataset) +
  geom_sf(mapping=aes(geometry = geometry,fill = Variabel_Signifikan)) +
  scale_fill_manual(values = warna_custom)+
  geom_text(
    aes(label = Kelurahan, x = coordinates(as(Dataset,"Spatial"))[,1], y = coordinates(as(Dataset,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+
  labs(fill="Variabel Signifikan Tiap Kelurahan")+xlab("Longitude")+ylab("Latitude")
