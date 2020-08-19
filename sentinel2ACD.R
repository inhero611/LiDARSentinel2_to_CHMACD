setwd("...")

#install.packages(c("raster","sp",'mapview',"rgdal","rgl","glcm"))
#install.packages("randomForest")

# Read data
library(sp)           
library(raster)       
library(mapview)
library(rgdal)
library(rgl)          
library(glcm)         
library(randomForest) 


# 2016 sentinel2 acd 
# preparing sentinel2 data
# b2(B) b3(G) b4(R) b5(RE1) b6(RE2) b7(NIR) b8A(NIRn2) b11(SWIR1) b12(SWIR2)
# original data
#sentinelb2 <- raster("sentinel2016_B2_50.tif")
#sentinelb3 <- raster("sentinel2016_B3_50.tif")
#sentinelb4 <- raster("sentinel2016_B4_50.tif")
#sentinelb5 <- raster("sentinel2016_B5_50.tif")
#sentinelb6 <- raster("sentinel2016_B6_50.tif")
#sentinelb7 <- raster("sentinel2016_B7_50.tif")
#sentinelb8A <- raster("sentinel2016_B8A_50.tif")
#sentinelb11 <- raster("sentinel2016_B11_50.tif")
#sentinelb12 <- raster("sentinel2016_B12_50.tif")

#no cloud data
sentinelb2 <- raster("sentinel2_2016_B2_50.tif")
sentinelb3 <- raster("sentinel2_2016_B3_50.tif")
sentinelb4 <- raster("sentinel2_2016_B4_50.tif")
sentinelb5 <- raster("sentinel2_2016_B5_50.tif")
sentinelb6 <- raster("sentinel2_2016_B6_50.tif")
sentinelb7 <- raster("sentinel2_2016_B7_50.tif")
sentinelb8A <- raster("sentinel2_2016_B8A_50.tif")
sentinelb11 <- raster("sentinel2_2016_B11_50.tif")
sentinelb12 <- raster("sentinel2_2016_B12_50.tif")

hist(sentinelb2, maxpixels=100000000, breaks=1000)
sentinelb2[sentinelb2 == 0] = NA
hist(sentinelb2, maxpixels=100000000, breaks=1000)
sentinelb3[sentinelb3 == 0] = NA
sentinelb4[sentinelb4 == 0] = NA
sentinelb5[sentinelb5 == 0] = NA
sentinelb6[sentinelb6 == 0] = NA
sentinelb7[sentinelb7 == 0] = NA
sentinelb8A[sentinelb8A == 0] = NA 
sentinelb11[sentinelb11 == 0] = NA
sentinelb12[sentinelb12 == 0] = NA
crs(sentinelb2)
res(sentinelb2)

# Spectral Indices
ARVI = (sentinelb7 - sentinelb4 - 0.106 * (sentinelb4 - sentinelb2)) / (sentinelb7 + sentinelb4 - 0.106 * (sentinelb4 - sentinelb2))
CVI = (sentinelb7 * sentinelb4)/(sentinelb3 ^2)
CI_G = sentinelb7/sentinelb3 -1
CI_RE1 = sentinelb7/sentinelb5 -1
DVI = sentinelb7 - sentinelb4
EVI_RE1 = 2.5 * (sentinelb5 - sentinelb4) / ((sentinelb5 + 6 * sentinelb4 - 7.5 * sentinelb2) + 1)
EVI_RE2 = 2.5 * (sentinelb6 - sentinelb4) / ((sentinelb6 + 6 * sentinelb4 - 7.5 * sentinelb2) + 1)
EVI = 2.5 * (sentinelb7 - sentinelb4) / ((sentinelb7 + 6 * sentinelb4 - 7.5 * sentinelb2) + 1)
EVI_NIR2 = 2.5 * (sentinelb8A - sentinelb4) / ((sentinelb8A + 6 * sentinelb4 - 7.5 * sentinelb2) + 1)
ENDVI = ((sentinelb7 + sentinelb3) - (2 * sentinelb2))/((sentinelb7 + sentinelb3) + (2 * sentinelb2))
IRECI = (sentinelb7 - sentinelb4) / (sentinelb5 / sentinelb6)
GNDVI = (sentinelb7 - sentinelb3) / (sentinelb7 + sentinelb3)
GDVI = sentinelb7 - sentinelb3
GARI = (sentinelb7 - (sentinelb3 - 1.7 * (sentinelb2 - sentinelb4))) / (sentinelb7 + (sentinelb3 - 1.7 * (sentinelb2 - sentinelb4)))
MSI =  sentinelb11 / sentinelb7
MSR_RE1 = (sentinelb5/sentinelb4 - 1)/(sqrt(sentinelb5/sentinelb4 + 1))
MSR_RE2 = (sentinelb6/sentinelb4 - 1)/(sqrt(sentinelb6/sentinelb4 + 1))
MSR = (sentinelb7/sentinelb4 - 1)/(sqrt(sentinelb7/sentinelb4 + 1))
MSR_NIR2 = (sentinelb8A/sentinelb4 - 1)/(sqrt(sentinelb8A/sentinelb4 + 1))
MCARI = ((sentinelb5 - sentinelb4) - 0.2 * (sentinelb5 - sentinelb3)) * (sentinelb5 / sentinelb4)
MARI = ((1 / sentinelb3) - (1 / sentinelb5)) * sentinelb7
NG = sentinelb3/(sentinelb7 + sentinelb4 + sentinelb3)
NR = sentinelb4/(sentinelb7 + sentinelb4 + sentinelb3)
NLI_RE1 = (sentinelb5 ^2 - sentinelb4)/ (sentinelb5 ^2 + sentinelb4)
NLI_RE2 = (sentinelb6 ^2 - sentinelb4)/ (sentinelb6 ^2 + sentinelb4)
NLI = (sentinelb7 ^2 - sentinelb4)/ (sentinelb7 ^2 + sentinelb4)
NLI_NIR2 = (sentinelb8A ^2 - sentinelb4)/ (sentinelb8A ^2 + sentinelb4)
NNIR = sentinelb7 / (sentinelb7 + sentinelb4 + sentinelb3)
NDVI = (sentinelb7 - sentinelb4) / (sentinelb7 + sentinelb4)
NDVI705 = (sentinelb5 - sentinelb4) / (sentinelb5 + sentinelb4)
NDII = (sentinelb7 - sentinelb11) / (sentinelb7 + sentinelb11)
NBR = (sentinelb7 - sentinelb12) / (sentinelb7 + sentinelb12)
NDWI = (sentinelb8A - sentinelb12) / (sentinelb8A + sentinelb12)
PSSR = sentinelb8A / sentinelb4
PSRI = (sentinelb4 - sentinelb3) / sentinelb8A
RDVI = (sentinelb7 - sentinelb4) / (sqrt(sentinelb7 + sentinelb4))
SAVI = (sentinelb7 - sentinelb4) / (sentinelb7 + sentinelb4 + 0.428) * (1 + 0.428);
TSAVI = (1.2 * (sentinelb7 - 1.2 * sentinelb4 - 0.04)) / (1.2 * sentinelb7 + sentinelb4 - (1.2 * 0.04) + 0.8 * (1 + 1.2 ^ 2) )
WDRVI_RE = (0.1 * sentinelb7 - sentinelb5) / (0.1 * sentinelb7 + sentinelb5) + (1 - 0.01) / (1 + 0.01)
WDRVI = (0.1 * sentinelb7 - sentinelb4) / (0.1 * sentinelb7 + sentinelb4) + (1 - 0.01) / (1 + 0.01)
VARIg = (sentinelb3 - sentinelb4) / (sentinelb3 + sentinelb4 - sentinelb2)

# Texture variable
# mean
mean_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")
mean_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics ="mean", min_x=NULL, max_x=NULL, na_opt="any")

# variance
var_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")
var_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "variance", min_x=NULL, max_x=NULL, na_opt="any")

# homogeneity
homo_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")
homo_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "homogeneity", min_x=NULL, max_x=NULL, na_opt="any")

# contrast
cont_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")
cont_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "contrast", min_x=NULL, max_x=NULL, na_opt="any")

# dissimilarity
diss_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")
diss_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "dissimilarity", min_x=NULL, max_x=NULL, na_opt="any")

# entropy
et_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")
et_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "entropy", min_x=NULL, max_x=NULL, na_opt="any")

# second_moment
sm_glcmb2 <- glcm(sentinelb2, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb3 <- glcm(sentinelb3, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb4 <- glcm(sentinelb4, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb5 <- glcm(sentinelb5, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb6 <- glcm(sentinelb6, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb7 <- glcm(sentinelb7, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb8A <- glcm(sentinelb8A, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb11 <- glcm(sentinelb11, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")
sm_glcmb12 <- glcm(sentinelb12, n_grey = 64, window = c(3, 3), shift = c(1, 1), statistics = "second_moment", min_x=NULL, max_x=NULL, na_opt="any")

# make training data 
# extract corresponding sentinel2 data
acd <- raster("acd_2016_50b.tif")
acd[acd == 0] = NA
acdp <- rasterToPoints(acd)
acdp_sp <- SpatialPoints(acdp, crs(acd))
acdp <- data.frame(acdp)[3]
head(acdp)


# if I can write them into a loop that will be better....
acdpb2 <- extract(sentinelb2,acdp_sp)
acdpb3 <- extract(sentinelb3,acdp_sp)
acdpb4 <- extract(sentinelb4,acdp_sp)
acdpb5 <- extract(sentinelb5,acdp_sp)
acdpb6 <- extract(sentinelb6,acdp_sp)
acdpb7 <- extract(sentinelb7,acdp_sp)
acdpb8A <- extract(sentinelb8A,acdp_sp)
acdpb11 <- extract(sentinelb11,acdp_sp)
acdpb12 <- extract(sentinelb12,acdp_sp)
acdpARVI <- extract(ARVI,acdp_sp)
acdpCVI <- extract(CVI,acdp_sp)
acdpCI_G <- extract(CI_G,acdp_sp)
acdpCI_RE1 <- extract(CI_RE1,acdp_sp)
acdpDVI <- extract(DVI,acdp_sp)
acdpEVI_RE1 <- extract(EVI_RE1,acdp_sp)
acdpEVI_RE2 <- extract(EVI_RE2,acdp_sp)
acdpEVI  <- extract(EVI,acdp_sp)
acdpEVI_NIR2 <- extract(EVI_NIR2,acdp_sp)
acdpENDVI <- extract(ENDVI,acdp_sp)
acdpIRECI <- extract(IRECI,acdp_sp)
acdpGNDVI <- extract(GNDVI,acdp_sp)
acdpGDVI <- extract(GDVI,acdp_sp)
acdpGARI <- extract(GARI,acdp_sp)
acdpMSI <- extract(MSI,acdp_sp)
acdpMSR_RE1 <- extract(MSR_RE1,acdp_sp)
acdpMSR_RE2 <- extract(MSR_RE2,acdp_sp)
acdpMSR <- extract(MSR,acdp_sp)
acdpMSR_NIR2 <- extract(MSR_NIR2,acdp_sp)
acdpMCARI <- extract(MCARI,acdp_sp)
acdpMARI <- extract(MARI,acdp_sp)
acdpNG <- extract(NG,acdp_sp)
acdpNR <- extract(NR,acdp_sp)
acdpNLI_RE1 <- extract(NLI_RE1,acdp_sp)
acdpNLI_RE2 <- extract(NLI_RE2,acdp_sp)
acdpNLI <- extract(NLI,acdp_sp)
acdpNLI_NIR2 <- extract(NLI_NIR2,acdp_sp)
acdpNNIR <- extract(NNIR,acdp_sp)
acdpNDVI <- extract(NDVI,acdp_sp)
acdpNDVI705 <- extract(NDVI705,acdp_sp)
acdpNDII <- extract(NDII,acdp_sp)
acdpNBR <- extract(NBR,acdp_sp)
acdpNDWI <- extract(NDWI,acdp_sp)
acdpPSSR <- extract(PSSR,acdp_sp)
acdpPSRI <- extract(PSRI,acdp_sp)
acdpRDVI <- extract(RDVI,acdp_sp)
acdpSAVI <- extract(SAVI,acdp_sp)
acdpTSAVI <- extract(TSAVI,acdp_sp)
acdpWDRVI_RE <- extract(WDRVI_RE,acdp_sp)
acdpWDRVI <- extract(WDRVI,acdp_sp)
acdpVARIg <- extract(VARIg,acdp_sp)

acdpmean_glcmb2 <- extract(mean_glcmb2,acdp_sp)
acdpmean_glcmb3 <- extract(mean_glcmb3,acdp_sp)
acdpmean_glcmb4 <- extract(mean_glcmb4,acdp_sp)
acdpmean_glcmb5 <- extract(mean_glcmb5,acdp_sp)
acdpmean_glcmb6 <- extract(mean_glcmb6,acdp_sp)
acdpmean_glcmb7 <- extract(mean_glcmb7,acdp_sp)
acdpmean_glcmb8A <- extract(mean_glcmb8A,acdp_sp)
acdpmean_glcmb11 <- extract(mean_glcmb11,acdp_sp)
acdpmean_glcmb12 <- extract(mean_glcmb12,acdp_sp)

acdpvar_glcmb2 <- extract(var_glcmb2,acdp_sp)
acdpvar_glcmb3 <- extract(var_glcmb3,acdp_sp)
acdpvar_glcmb4 <- extract(var_glcmb4,acdp_sp)
acdpvar_glcmb5 <- extract(var_glcmb5,acdp_sp)
acdpvar_glcmb6 <- extract(var_glcmb6,acdp_sp)
acdpvar_glcmb7 <- extract(var_glcmb7,acdp_sp)
acdpvar_glcmb8A <- extract(var_glcmb8A,acdp_sp)
acdpvar_glcmb11 <- extract(var_glcmb11,acdp_sp)
acdpvar_glcmb12 <- extract(var_glcmb12,acdp_sp)

acdphomo_glcmb2 <- extract(homo_glcmb2,acdp_sp)
acdphomo_glcmb3 <- extract(homo_glcmb3,acdp_sp)
acdphomo_glcmb4 <- extract(homo_glcmb4,acdp_sp)
acdphomo_glcmb5 <- extract(homo_glcmb5,acdp_sp)
acdphomo_glcmb6 <- extract(homo_glcmb6,acdp_sp)
acdphomo_glcmb7 <- extract(homo_glcmb7,acdp_sp)
acdphomo_glcmb8A <- extract(homo_glcmb8A,acdp_sp)
acdphomo_glcmb11 <- extract(homo_glcmb11,acdp_sp)
acdphomo_glcmb12 <- extract(homo_glcmb12,acdp_sp)

acdpcont_glcmb2 <- extract(cont_glcmb2,acdp_sp)
acdpcont_glcmb3 <- extract(cont_glcmb3,acdp_sp)
acdpcont_glcmb4 <- extract(cont_glcmb4,acdp_sp)
acdpcont_glcmb5 <- extract(cont_glcmb5,acdp_sp)
acdpcont_glcmb6 <- extract(cont_glcmb6,acdp_sp)
acdpcont_glcmb7 <- extract(cont_glcmb7,acdp_sp)
acdpcont_glcmb8A <- extract(cont_glcmb8A,acdp_sp)
acdpcont_glcmb11 <- extract(cont_glcmb11,acdp_sp)
acdpcont_glcmb12 <- extract(cont_glcmb12,acdp_sp)

acdpdiss_glcmb2 <- extract(diss_glcmb2,acdp_sp)
acdpdiss_glcmb3 <- extract(diss_glcmb3,acdp_sp)
acdpdiss_glcmb4 <- extract(diss_glcmb4,acdp_sp)
acdpdiss_glcmb5 <- extract(diss_glcmb5,acdp_sp)
acdpdiss_glcmb6 <- extract(diss_glcmb6,acdp_sp)
acdpdiss_glcmb7 <- extract(diss_glcmb7,acdp_sp)
acdpdiss_glcmb8A <- extract(diss_glcmb8A,acdp_sp)
acdpdiss_glcmb11 <- extract(diss_glcmb11,acdp_sp)
acdpdiss_glcmb12 <- extract(diss_glcmb12,acdp_sp)

acdpet_glcmb2 <- extract(et_glcmb2,acdp_sp)
acdpet_glcmb3 <- extract(et_glcmb3,acdp_sp)
acdpet_glcmb4 <- extract(et_glcmb4,acdp_sp)
acdpet_glcmb5 <- extract(et_glcmb5,acdp_sp)
acdpet_glcmb6 <- extract(et_glcmb6,acdp_sp)
acdpet_glcmb7 <- extract(et_glcmb7,acdp_sp)
acdpet_glcmb8A <- extract(et_glcmb8A,acdp_sp)
acdpet_glcmb11 <- extract(et_glcmb11,acdp_sp)
acdpet_glcmb12 <- extract(et_glcmb12,acdp_sp)

acdpsm_glcmb2 <- extract(sm_glcmb2,acdp_sp)
acdpsm_glcmb3 <- extract(sm_glcmb3,acdp_sp)
acdpsm_glcmb4 <- extract(sm_glcmb4,acdp_sp)
acdpsm_glcmb5 <- extract(sm_glcmb5,acdp_sp)
acdpsm_glcmb6 <- extract(sm_glcmb6,acdp_sp)
acdpsm_glcmb7 <- extract(sm_glcmb7,acdp_sp)
acdpsm_glcmb8A <- extract(sm_glcmb8A,acdp_sp)
acdpsm_glcmb11 <- extract(sm_glcmb11,acdp_sp)
acdpsm_glcmb12 <- extract(sm_glcmb12,acdp_sp)

# Model 5 bands +VI +texture
s2016_5 <- data.frame(cbind(acdp,acdpb2,acdpb3,acdpb4,acdpb5,acdpb6,acdpb7,acdpb8A,acdpb11,acdpb12,
                            acdpMCARI,acdpMARI,acdpNLI_RE1,acdpNLI_RE2,acdpNDVI705,acdpPSRI,
                            acdpmean_glcmb2,acdpmean_glcmb3,acdpmean_glcmb4,acdpmean_glcmb5,acdpmean_glcmb6,acdpmean_glcmb7,acdpmean_glcmb8A,acdpmean_glcmb11,acdpmean_glcmb12,
                            acdpvar_glcmb2,acdpvar_glcmb3,acdpvar_glcmb4,acdpvar_glcmb5,acdpvar_glcmb6,acdpvar_glcmb7,acdpvar_glcmb8A,acdpvar_glcmb11,acdpvar_glcmb12,
                            acdphomo_glcmb2,acdphomo_glcmb3,acdphomo_glcmb4,acdphomo_glcmb5,acdphomo_glcmb6,acdphomo_glcmb7,acdphomo_glcmb8A,acdphomo_glcmb11,acdphomo_glcmb12,
                            acdpcont_glcmb2,acdpcont_glcmb3,acdpcont_glcmb4,acdpcont_glcmb5,acdpcont_glcmb6,acdpcont_glcmb7,acdpcont_glcmb8A,acdpcont_glcmb11,acdpcont_glcmb12,
                            acdpdiss_glcmb2,acdpdiss_glcmb3,acdpdiss_glcmb4,acdpdiss_glcmb5,acdpdiss_glcmb6,acdpdiss_glcmb7,acdpdiss_glcmb8A,acdpdiss_glcmb11,acdpdiss_glcmb12,
                            acdpet_glcmb2,acdpet_glcmb3,acdpet_glcmb4,acdpet_glcmb5,acdpet_glcmb6,acdpet_glcmb7,acdpet_glcmb8A,acdpet_glcmb11,acdpet_glcmb12,
                            acdpsm_glcmb2,acdpsm_glcmb3,acdpsm_glcmb4,acdpsm_glcmb5,acdpsm_glcmb6,acdpsm_glcmb7,acdpsm_glcmb8A,acdpsm_glcmb11,acdpsm_glcmb12)) #,acdpcorr_glcmb2,acdpcorr_glcmb3,acdpcorr_glcmb4,acdpcorr_glcmb5,acdpcorr_glcmb6,acdpcorr_glcmb7,acdpcorr_glcmb8A,acdpcorr_glcmb11,acdpcorr_glcmb12))

names(s2016_5) <- c("acd","b2","b3","b4","b5","b6","b7","b8A","b11","b12",
                    "MCARI","MARI","NLI_RE1","NLI_RE2","NDVI705","PSRI",
                    "mean_glcmb2","mean_glcmb3","mean_glcmb4","mean_glcmb5","mean_glcmb6","mean_glcmb7","mean_glcmb8A","mean_glcmb11","mean_glcmb12",
                    "var_glcmb2","var_glcmb3","var_glcmb4","var_glcmb5","var_glcmb6","var_glcmb7","var_glcmb8A","var_glcmb11","var_glcmb12",
                    "homo_glcmb2","homo_glcmb3","homo_glcmb4","homo_glcmb5","homo_glcmb6","homo_glcmb7","homo_glcmb8A","homo_glcmb11","homo_glcmb12",
                    "cont_glcmb2","cont_glcmb3","cont_glcmb4","cont_glcmb5","cont_glcmb6","cont_glcmb7","cont_glcmb8A","cont_glcmb11","cont_glcmb12",
                    "diss_glcmb2","diss_glcmb3","diss_glcmb4","diss_glcmb5","diss_glcmb6","diss_glcmb7","diss_glcmb8A","diss_glcmb11","diss_glcmb12",
                    "et_glcmb2","et_glcmb3","et_glcmb4","et_glcmb5","et_glcmb6","et_glcmb7","et_glcmb8A","et_glcmb11","et_glcmb12",
                    "sm_glcmb2","sm_glcmb3","sm_glcmb4","sm_glcmb5","sm_glcmb6","sm_glcmb7","sm_glcmb8A","sm_glcmb11","sm_glcmb12") #,"corr_glcmb2","corr_glcmb3","corr_glcmb4","corr_glcmb5","corr_glcmb6","corr_glcmb7","corr_glcmb8A","corr_glcmb11","corr_glcmb12")

head(s2016_5)
write.csv(s2016_5,file="acdsentinel2016_5.csv")

s2016_5 = na.omit(s2016_5) # remove rows with NAN / NA values
s2016_5 <- s2016_5[is.finite(rowSums(s2016_5)),] # remove rows with Inf values

# tuning model
#t <- tuneRF(s2016_5[,-1], s2016_5[,1], ntreeTry=800, stepFactor=0.5, improve=0.05,trace=TRUE, plot=TRUE, doBest=FALSE)

rf_5 <- randomForest(acd ~ ., data = s2016_5,
                     ntree = 800,
                     mtry = 13,
                     importance = TRUE,
                     proximity = TRUE,
)#na.action = na.omit for missing value from cloud removal
print(rf_5)
varImpPlot(rf_5)


s2016 <- data.frame(cbind(acdp,acdpb2,acdpb3,acdpb4,acdpb5,acdpb6,acdpb7,acdpb8A,acdpb11,acdpb12,
                          acdpARVI,acdpCVI,acdpCI_G,acdpCI_RE1,acdpDVI,acdpEVI_RE1,acdpEVI_RE2,acdpEVI,acdpEVI_NIR2,acdpENDVI,acdpIRECI,acdpGNDVI,acdpGDVI,acdpGARI,acdpMSI,acdpMSR_RE1,acdpMSR_RE2,acdpMSR,acdpMSR_NIR2,acdpMCARI,acdpMARI,acdpNG,acdpNR,acdpNLI_RE1,acdpNLI_RE2,acdpNLI,acdpNLI_NIR2,acdpNNIR,acdpNDVI,acdpNDVI705,acdpNDII,acdpNBR,acdpNDWI,acdpPSSR,acdpPSRI, 
                          acdpRDVI,acdpSAVI,acdpTSAVI,acdpWDRVI_RE,acdpWDRVI,acdpVARIg,
                          acdpmean_glcmb2,acdpmean_glcmb3,acdpmean_glcmb4,acdpmean_glcmb5,acdpmean_glcmb6,acdpmean_glcmb7,acdpmean_glcmb8A,acdpmean_glcmb11,acdpmean_glcmb12,
                          acdpvar_glcmb2,acdpvar_glcmb3,acdpvar_glcmb4,acdpvar_glcmb5,acdpvar_glcmb6,acdpvar_glcmb7,acdpvar_glcmb8A,acdpvar_glcmb11,acdpvar_glcmb12,
                          acdphomo_glcmb2,acdphomo_glcmb3,acdphomo_glcmb4,acdphomo_glcmb5,acdphomo_glcmb6,acdphomo_glcmb7,acdphomo_glcmb8A,acdphomo_glcmb11,acdphomo_glcmb12,
                          acdpcont_glcmb2,acdpcont_glcmb3,acdpcont_glcmb4,acdpcont_glcmb5,acdpcont_glcmb6,acdpcont_glcmb7,acdpcont_glcmb8A,acdpcont_glcmb11,acdpcont_glcmb12,
                          acdpdiss_glcmb2,acdpdiss_glcmb3,acdpdiss_glcmb4,acdpdiss_glcmb5,acdpdiss_glcmb6,acdpdiss_glcmb7,acdpdiss_glcmb8A,acdpdiss_glcmb11,acdpdiss_glcmb12,
                          acdpet_glcmb2,acdpet_glcmb3,acdpet_glcmb4,acdpet_glcmb5,acdpet_glcmb6,acdpet_glcmb7,acdpet_glcmb8A,acdpet_glcmb11,acdpet_glcmb12,
                          acdpsm_glcmb2,acdpsm_glcmb3,acdpsm_glcmb4,acdpsm_glcmb5,acdpsm_glcmb6,acdpsm_glcmb7,acdpsm_glcmb8A,acdpsm_glcmb11,acdpsm_glcmb12)) #,acdpcorr_glcmb2,acdpcorr_glcmb3,acdpcorr_glcmb4,acdpcorr_glcmb5,acdpcorr_glcmb6,acdpcorr_glcmb7,acdpcorr_glcmb8A,acdpcorr_glcmb11,acdpcorr_glcmb12))

names(s2016) <- c("acd", "b2", "b3", "b4", "b5", "b6", "b7", "b8A", "b11", "b12",
                  "ARVI", "CVI","CI_G", "CI_RE1", "DVI","EVI_RE1", "EVI_RE2","EVI","EVI_NIR2","ENDVI","IRECI","GNDVI","GDVI","GARI","MSI","MSR_RE1","MSR_RE2","MSR","MSR_NIR2","MCARI","MARI","NG","NR","NLI_RE1","NLI_RE2","NLI","NLI_NIR2","NNIR","NDVI","NDVI705","NDII","NBR","NDWI","PSSR","PSRI","RDVI","SAVI","TSAVI","WDRVI_RE","WDRVI","VARIg",
                  "mean_glcmb2","mean_glcmb3","mean_glcmb4","mean_glcmb5","mean_glcmb6","mean_glcmb7","mean_glcmb8A","mean_glcmb11","mean_glcmb12",
                  "var_glcmb2","var_glcmb3","var_glcmb4","var_glcmb5","var_glcmb6","var_glcmb7","var_glcmb8A","var_glcmb11","var_glcmb12",
                  "homo_glcmb2","homo_glcmb3","homo_glcmb4","homo_glcmb5","homo_glcmb6","homo_glcmb7","homo_glcmb8A","homo_glcmb11","homo_glcmb12",
                  "cont_glcmb2","cont_glcmb3","cont_glcmb4","cont_glcmb5","cont_glcmb6","cont_glcmb7","cont_glcmb8A","cont_glcmb11","cont_glcmb12",
                  "diss_glcmb2","diss_glcmb3","diss_glcmb4","diss_glcmb5","diss_glcmb6","diss_glcmb7","diss_glcmb8A","diss_glcmb11","diss_glcmb12",
                  "et_glcmb2","et_glcmb3","et_glcmb4","et_glcmb5","et_glcmb6","et_glcmb7","et_glcmb8A","et_glcmb11","et_glcmb12",
                  "sm_glcmb2","sm_glcmb3","sm_glcmb4","sm_glcmb5","sm_glcmb6","sm_glcmb7","sm_glcmb8A","sm_glcmb11","sm_glcmb12") #,"corr_glcmb2","corr_glcmb3","corr_glcmb4","corr_glcmb5","corr_glcmb6","corr_glcmb7","corr_glcmb8A","corr_glcmb11","corr_glcmb12")
#str(s2016)
head(s2016)
write.csv(s2016,file="acdsentinel2016.csv")
#s2016$acd <- as.factor(s2016$acd)
#table(s$acd)
#sapply(s2016, function(x) sum(is.infinite(x))) #-- correlation have thousands of NA

s2016 = na.omit(s2016) # remove rows with NAN / NA values
s2016 <- s2016[is.finite(rowSums(s2016)),] # remove rows with Inf values

# tune model
t <- tuneRF(s2016[,-1], s2016[,1], ntreeTry=800, stepFactor=0.5, improve=0.05,
            trace=TRUE, plot=TRUE, doBest=FALSE)

rf <- randomForest(acd ~ ., data = s2016,
                   ntree = 800,
                   mtry = 18, 
                   importance = TRUE,
                   proximity = TRUE,
                   ) #na.action = na.omit for missing value from cloud removal
print(rf)
varImpPlot(rf)
plot(rf)

#Cross-validation with 5 fold

cval <- rfcv(s2016[,-1], s2016[,1], cv.fold=10, scale="log", step=0.5,
             mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE)

cval

# rename sequency of image stack with the same name than acd_2016
# stack Sentinel 2 images
image_stack = stack(sentinelb2,sentinelb3,sentinelb4,sentinelb5,sentinelb6,sentinelb7,sentinelb8A,sentinelb11,sentinelb12,
                    ARVI,CVI,CI_G,CI_RE1,DVI,EVI_RE1,EVI_RE2,EVI,EVI_NIR2,ENDVI,IRECI,GNDVI,GDVI,GARI,MSI,MSR_RE1,MSR_RE2,MSR,MSR_NIR2,MCARI,MARI,NG,NR,NLI_RE1,NLI_RE2,NLI,
                    NLI_NIR2,NNIR,NDVI,NDVI705,NDII,NBR,NDWI,PSSR,PSRI,RDVI,SAVI,TSAVI,WDRVI_RE,WDRVI,VARIg,
                    mean_glcmb2,mean_glcmb3,mean_glcmb4,mean_glcmb5,mean_glcmb6,mean_glcmb7,mean_glcmb8A,mean_glcmb11,mean_glcmb12,
                    var_glcmb2,var_glcmb3,var_glcmb4,var_glcmb5,var_glcmb6,var_glcmb7,var_glcmb8A,var_glcmb11,var_glcmb12,
                    homo_glcmb2,homo_glcmb3,homo_glcmb4,homo_glcmb5,homo_glcmb6,homo_glcmb7,homo_glcmb8A,homo_glcmb11,homo_glcmb12,
                    cont_glcmb2,cont_glcmb3,cont_glcmb4,cont_glcmb5,cont_glcmb6,cont_glcmb7,cont_glcmb8A,cont_glcmb11,cont_glcmb12,
                    diss_glcmb2,diss_glcmb3,diss_glcmb4,diss_glcmb5,diss_glcmb6,diss_glcmb7,diss_glcmb8A,diss_glcmb11,diss_glcmb12,
                    et_glcmb2,et_glcmb3,et_glcmb4,et_glcmb5,et_glcmb6,et_glcmb7,et_glcmb8A,et_glcmb11,et_glcmb12,
                    sm_glcmb2,sm_glcmb3,sm_glcmb4,sm_glcmb5,sm_glcmb6,sm_glcmb7,sm_glcmb8A,sm_glcmb11,sm_glcmb12) #,corr_glcmb2,corr_glcmb3,corr_glcmb4,corr_glcmb5,corr_glcmb6,corr_glcmb7,corr_glcmb8A,corr_glcmb11,corr_glcmb12)


image_stack[image_stack == 0] <- NaN

image_df = as.data.frame(image_stack)

names(image_df) <- c("b2", "b3", "b4", "b5", "b6", "b7", "b8A", "b11", "b12","ARVI", "CVI","CI_G", "CI_RE1", "DVI","EVI_RE1", "EVI_RE2","EVI","EVI_NIR2","ENDVI","IRECI","GNDVI","GDVI","GARI","MSI","MSR_RE1","MSR_RE2","MSR","MSR_NIR2","MCARI","MARI","NG","NR","NLI_RE1","NLI_RE2","NLI","NLI_NIR2","NNIR","NDVI","NDVI705","NDII","NBR","NDWI","PSSR","PSRI","RDVI","SAVI","TSAVI","WDRVI_RE","WDRVI","VARIg",
                         "mean_glcmb2","mean_glcmb3","mean_glcmb4","mean_glcmb5","mean_glcmb6","mean_glcmb7","mean_glcmb8A","mean_glcmb11","mean_glcmb12",
                         "var_glcmb2","var_glcmb3","var_glcmb4","var_glcmb5","var_glcmb6","var_glcmb7","var_glcmb8A","var_glcmb11","var_glcmb12",
                         "homo_glcmb2","homo_glcmb3","homo_glcmb4","homo_glcmb5","homo_glcmb6","homo_glcmb7","homo_glcmb8A","homo_glcmb11","homo_glcmb12",
                         "cont_glcmb2","cont_glcmb3","cont_glcmb4","cont_glcmb5","cont_glcmb6","cont_glcmb7","cont_glcmb8A","cont_glcmb11","cont_glcmb12",
                         "diss_glcmb2","diss_glcmb3","diss_glcmb4","diss_glcmb5","diss_glcmb6","diss_glcmb7","diss_glcmb8A","diss_glcmb11","diss_glcmb12",
                         "et_glcmb2","et_glcmb3","et_glcmb4","et_glcmb5","et_glcmb6","et_glcmb7","et_glcmb8A","et_glcmb11","et_glcmb12",
                         "sm_glcmb2","sm_glcmb3","sm_glcmb4","sm_glcmb5","sm_glcmb6","sm_glcmb7","sm_glcmb8A","sm_glcmb11","sm_glcmb12") #,"corr_glcmb2","corr_glcmb3","corr_glcmb4","corr_glcmb5","corr_glcmb6","corr_glcmb7","corr_glcmb8A","corr_glcmb11","corr_glcmb12")


#prediction

p1 <- predict(rf, image_df)

# use sentinel image to paste vales of prediction dataframe values into
image_out = sentinelb2
image_out[] = p1
writeRaster(image_out, "acd_sentinel2016.tif")

plot(image_out)