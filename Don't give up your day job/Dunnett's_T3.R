library(readxl)
library(PMCMRplus)

data <- read_excel("filtered_data.xlsx")
data$strain <- as.factor(data$strain)

fit_wet <- aov(wet ~ strain, data)
shapiro.test(residuals(fit_wet))
bartlett.test(wet ~ strain, data)
anova(fit_wet)

T3res_wet <- dunnettT3Test(wet ~ strain, data)
summary(T3res_wet)
summaryGroup(T3res_wet)

out <- capture.output(summary(T3res_wet))
cat("Dunnett's T3 wet summary", out, file = "summary_wet.txt", sep = "n", append = TRUE)

fit_dry <- aov(dry ~ strain, data)
shapiro.test(residuals(fit_dry))
bartlett.test(dry ~ strain, data)
anova(fit_dry)

T3res_dry <- dunnettT3Test(dry ~ strain, data)
summary(T3res_dry)
summaryGroup(T3res_dry)

out <- capture.output(summary(T3res_dry))
cat("Dunnett's T3 dry summary", out, file = "summary_dry.txt", sep = "n", append = TRUE)

nitrogen <- read_excel("raw_data.xlsx", sheet = "nitrogen")
nitrogen$strain <- as.factor(nitrogen$strain)

fit_nitrogen <- aov(value ~ strain, nitrogen)
shapiro.test(residuals(fit_nitrogen))
bartlett.test(value ~ strain, nitrogen)
anova(fit_nitrogen)

T3_nitrogen <- dunnettT3Test(value ~ strain, nitrogen)
summary(T3_nitrogen)
summaryGroup(T3_nitrogen)

out <- capture.output(summary(T3res_wet))
cat("Dunnett's T3 nitrogen summary", out, file = "summary_nitrogen.txt", sep = "n", append = TRUE)