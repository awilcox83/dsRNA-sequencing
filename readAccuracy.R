
library(ShortRead)
library(ShadowRegression)
library(ggplot2)

# convert fastq files to rcnt

getReadCnts('.', 'phi6_15_40.fastq', 'Phi_15.rcnts')
getReadCnts('.', 'phi6_50_40.fastq', 'Phi_50.rcnts')
getReadCnts('.', 'phi6_90_40.fastq', 'Phi_90.rcnts')
getReadCnts('.', 'flu_0_40.fastq', 'Flu_0.rcnts')
getReadCnts('.', 'flu_15_40.fastq', 'Flu_15.rcnts')
getReadCnts('.', 'flu_50_40.fastq', 'Flu_50.rcnts')
getReadCnts('.', 'flu_90_40.fastq', 'Flu_90.rcnts')

# calculate error rates

PhiErrorRates15 = getErrorRates("Phi_15.rcnts", type="sub");
PhiErrorRates50 = getErrorRates("Phi_50.rcnts", type="sub");
PhiErrorRates90 = getErrorRates("Phi_90.rcnts", type="sub");

FluErrorRates0 = getErrorRates("flu_0.rcnts", type="sub");
FluErrorRates15 = getErrorRates("flu_15.rcnts", type="sub");
FluErrorRates50 = getErrorRates("flu_50.rcnts", type="sub");
FluErrorRates90 = getErrorRates("flu_90.rcnts", type="sub");

PhiErrorRates15$perReadER
PhiErrorRates50$perReadER
PhiErrorRates90$perReadER

FluErrorRates0$perReadER
FluErrorRates15$perReadER
FluErrorRates50$perReadER
FluErrorRates90$perReadER



# create a dataframe containing individual error rates and std errors.

dmsoTreatment <- c("15%", "50%", "90%")
errors <- c(PhiErrorRates15[["perReadER"]][["error rate"]], PhiErrorRates50[["perReadER"]][["error rate"]], PhiErrorRates90[["perReadER"]][["error rate"]])
stdErrors <-c(PhiErrorRates15[["perReadER"]][["standard error"]],PhiErrorRates50[["perReadER"]][["standard error"]],PhiErrorRates90[["perReadER"]][["standard error"]])
PhiDfErrorRate <- data.frame(dmsoTreatment, errors, stdErrors)
names(PhiDfErrorRate)<- c("DMSO", "Error Rate", "Std Error")

dmsoTreatment <- c("0%", "15%", "50%", "90%")
errors <- c(FluErrorRates0[["perReadER"]][["error rate"]], FluErrorRates15[["perReadER"]][["error rate"]], FluErrorRates50[["perReadER"]][["error rate"]], FluErrorRates90[["perReadER"]][["error rate"]])
stdErrors <-c(FluErrorRates0[["perReadER"]][["standard error"]],FluErrorRates15[["perReadER"]][["standard error"]],FluErrorRates50[["perReadER"]][["standard error"]],FluErrorRates90[["perReadER"]][["standard error"]])
FluDfErrorRate <- data.frame(dmsoTreatment, errors, stdErrors)
names(FluDfErrorRate)<- c("DMSO", "Error Rate", "Std Error")

# make bar charts to compare error rates

phiErrorRate = ggplot(data=PhiDfErrorRate,  aes(x=dmsoTreatment, y=errors, se=stdErrors)) + 
  geom_col(color="black",fill="white") +
  geom_errorbar(aes(ymin=errors-stdErrors, ymax=errors+stdErrors), width=.2,
                position=position_dodge(.9)) +
  xlab("DMSO treatment") + 
  ylab("error rate (per read)")
ggsave("phiErrorRate.pdf")

fluErrorRate = ggplot(data=FluDfErrorRate, aes(x=dmsoTreatment, y=errors, se=stdErrors)) + 
  geom_col(color="black",fill="white") +
  geom_errorbar(aes(ymin=errors-stdErrors, ymax=errors+stdErrors), width=.2,
                position=position_dodge(.9)) +
  xlab("DMSO treatment") + 
  ylab("error rate (per read)")
ggsave("fluErrorRate.pdf")


# make robust linear regression models for each

Philm15.rlm <- rlm(shadows ~ tags, data=PhiErrorRates15)
Philm50.rlm <- rlm(shadows ~ tags, data=PhiErrorRates50)
Philm90.rlm <- rlm(shadows ~ tags, data=PhiErrorRates90)

Flulm0.rlm <- rlm(shadows ~ tags, data=FluErrorRates0)
Flulm15.rlm <- rlm(shadows ~ tags, data=FluErrorRates15)
Flulm50.rlm <- rlm(shadows ~ tags, data=FluErrorRates50)
Flulm90.rlm <- rlm(shadows ~ tags, data=FluErrorRates90)


# plot all regressions on one chart and save:

PhiDf15 = data.frame(PhiErrorRates15$tags, PhiErrorRates15$shadows)
PhiDf50 = data.frame(PhiErrorRates50$tags, PhiErrorRates50$shadows)
PhiDf90 = data.frame(PhiErrorRates90$tags, PhiErrorRates90$shadows)
PhiDf15$DMSO <- "15"
PhiDf50$DMSO <- "50"
PhiDf90$DMSO <- "90"
names(PhiDf15)<- c("Tags", "Shadows", "DMSO")
names(PhiDf50)<- c("Tags", "Shadows", "DMSO")
names(PhiDf90)<- c("Tags", "Shadows", "DMSO")
phiData <- rbind(PhiDf15,PhiDf50, PhiDf90)
phiData.rlm <- rlm(Shadows ~ Tags*DMSO, data=phiData)

phiPlot = ggplot() +
  geom_smooth(data=PhiDf15, aes(x=Tags, y=Shadows, color="15% DMSO"), method="rlm", se=TRUE)+
  geom_smooth(data=PhiDf50, aes(x=Tags, y=Shadows, color="50% DMSO"), method="rlm", se=TRUE)+
  geom_smooth(data=PhiDf90, aes(x=Tags, y=Shadows, color="90% DMSO"), method="rlm", se=TRUE) +
  scale_colour_manual("", 
                      breaks = c("15% DMSO", "50% DMSO", "90% DMSO"),
                      values = c("green ", "blue", "red"))
ggsave("phiPlot.pdf")

# and for flu

FluDf0 = data.frame(FluErrorRates0$tags, FluErrorRates0$shadows)
FluDf15 = data.frame(FluErrorRates15$tags, FluErrorRates15$shadows)
FluDf50 = data.frame(FluErrorRates50$tags, FluErrorRates50$shadows)
FluDf90 = data.frame(FluErrorRates90$tags, FluErrorRates90$shadows)
FluDf0$DMSO <- "0"
FluDf15$DMSO <- "15"
FluDf50$DMSO <- "50"
FluDf90$DMSO <- "90"
names(FluDf0)<- c("Tags", "Shadows", "DMSO")
names(FluDf15)<- c("Tags", "Shadows", "DMSO")
names(FluDf50)<- c("Tags", "Shadows", "DMSO")
names(FluDf90)<- c("Tags", "Shadows", "DMSO")
FluData <- rbind(FluDf0, FluDf15,FluDf50, FluDf90)
FluData.rlm <- rlm(Shadows ~ Tags*DMSO, data=FluData)

FluPlot = ggplot() +
  geom_smooth(data=FluDf0, aes(x=Tags, y=Shadows, color="0% DMSO"), method="rlm", se=TRUE)+
  geom_smooth(data=FluDf15, aes(x=Tags, y=Shadows, color="15% DMSO"), method="rlm", se=TRUE)+
  geom_smooth(data=FluDf50, aes(x=Tags, y=Shadows, color="50% DMSO"), method="rlm", se=TRUE)+
  geom_smooth(data=FluDf90, aes(x=Tags, y=Shadows, color="90% DMSO"), method="rlm", se=TRUE) +
  scale_colour_manual("", 
                      breaks = c("0% DMSO", "15% DMSO", "50% DMSO", "90% DMSO"),
                      values = c("black", "green ", "blue", "red"))
ggsave("FluPlot.pdf")

