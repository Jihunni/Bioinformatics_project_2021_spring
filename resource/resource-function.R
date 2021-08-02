
### Normalization ###
reference: https://medium.com/swlh/data-normalisation-with-r-6ef1d1947970
# normalized_age <- as.data.frame(lapply(normalized_clinic_process$age, norm_minmax))
# head(normalized_age)
# 
# normalized_age %>% summarise(Min = min(YearsExperience,na.rm = TRUE), Q1 = quantile(YearsExperience,probs = .25,na.rm = TRUE), Median = median(YearsExperience, na.rm = TRUE),Q3 = quantile(YearsExperience,probs = .75,na.rm = TRUE),Max = max(YearsExperience,na.rm = TRUE), Mean = mean(YearsExperience, na.rm = TRUE),SD = sd(YearsExperience, na.rm = TRUE), n = n(), Missing = sum(is.na(YearsExperience)))
# 
# normalise_salary %>% summarise(Min = min(Salary,na.rm = TRUE),Q1 = quantile(Salary,probs = .25,na.rm = TRUE),Median = median(Salary, na.rm = TRUE),Q3 = quantile(Salary,probs = .75,na.rm = TRUE),
#                                Max = max(Salary,na.rm = TRUE),
#                                Mean = mean(Salary, na.rm = TRUE),
#                                SD = sd(Salary, na.rm = TRUE),
#                                n = n(),
#                                Missing = sum(is.na(Salary)))