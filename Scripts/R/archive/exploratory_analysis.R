



# libraries ---------------------------------------------------------------


        library(tidyverse)
        library(nortest)
        library(FSA)
        library(tidymodels)




# data --------------------------------------------------------------------



        # set wd
        setwd("~")
        setwd("Documents/cerebral oxygenation/")

        # read in data
        df_hbox <- read.csv("data exports for analysis/98 first visit records with cerebral oxygenation stats_trimmed.csv") %>% as_tibble()


        # format
        df_hbox <- df_hbox %>% 
                filter(Status != "Other") %>% 
                mutate(Status = factor(Status, levels = c("HV", "UM", "CM")))
        
        df_hbox <- df_hbox %>% 
                mutate(Age = gsub(" Years, ", "_", Age)) %>% 
                mutate(Age = gsub(" Months ", "", Age)) %>% 
                separate(Age, into = c("Years", "Months"), sep = "_") %>% 
                mutate(Years = as.numeric(Years)) %>% 
                mutate(Months = as.numeric(Months)/12) %>% 
                mutate(Age = Years + Months) %>% 
                select(-c(Years, Months))
        
        df_hbox <- df_hbox %>% mutate(`%OxyHb` = cerebral.hbox/(cerebral.hbox + cerebral.hbdeox))
        

# test for normality ------------------------------------------------------

        
        # deoxyHb
        shapiro.test(df_hbox$cerebral.hbdeox)
        ad.test(df_hbox$cerebral.hbdeox)
        df_hbox %>% ggplot(aes(sample = cerebral.hbdeox)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # deoxyHb alpha
        shapiro.test(df_hbox$cerebral.hbdeox.alpha)
        ad.test(df_hbox$cerebral.hbdeox.alpha)
        df_hbox %>% ggplot(aes(sample = cerebral.hbdeox.alpha)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # oxyHb
        shapiro.test(df_hbox$cerebral.hbox)
        ad.test(df_hbox$cerebral.hbox)
        df_hbox %>% ggplot(aes(sample = cerebral.hbox)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # oxyHb alpha
        shapiro.test(df_hbox$cerebral.hbox.alpha)
        ad.test(df_hbox$cerebral.hbox.alpha)
        df_hbox %>% ggplot(aes(sample = cerebral.hbox.alpha)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # sat
        shapiro.test(df_hbox$cerebral.sat)
        ad.test(df_hbox$cerebral.sat)
        df_hbox %>% ggplot(aes(sample = cerebral.sat)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # normal
        
        # sat alpha
        shapiro.test(df_hbox$cerebral.sat.alpha)
        ad.test(df_hbox$cerebral.sat.alpha)
        df_hbox %>% ggplot(aes(sample = cerebral.sat.alpha)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # thc
        shapiro.test(df_hbox$cerebral.thc)
        ad.test(df_hbox$cerebral.thc)
        df_hbox %>% ggplot(aes(sample = cerebral.thc)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # thc alpha
        shapiro.test(df_hbox$cerebral.thc.alpha)
        ad.test(df_hbox$cerebral.thc.alpha)
        df_hbox %>% ggplot(aes(sample = cerebral.thc.alpha)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # thc norm std
        shapiro.test(df_hbox$cerebral.thc.norm.std)
        ad.test(df_hbox$cerebral.thc.norm.std)
        df_hbox %>% ggplot(aes(sample = cerebral.thc.norm.std)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # haptoglobin
        shapiro.test(df_hbox$Haptoglobin..mg.dl.)
        ad.test(df_hbox$Haptoglobin..mg.dl.)
        df_hbox %>% ggplot(aes(sample = Haptoglobin..mg.dl.)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        
        # Hemoglobin
        shapiro.test(df_hbox$Hemoglobin..uM.)
        ad.test(df_hbox$Hemoglobin..uM.)
        df_hbox %>% ggplot(aes(sample = Hemoglobin..uM.)) +
                stat_qq() +
                stat_qq_line(color = "red") +
                theme_bw() # not normal
        

        
        
        
# plots -----------------------------------------------------------------

        


        # Deoxy Hb
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.hbdeox)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.hbdeox") +
                theme(legend.position = "none")
        
        # Deoxy Hb alpha
        df_hbox %>% filter(cerebral.hbdeox.alpha < 3) %>% # removed two clear outliers just for visualization
                
                ggplot(aes(x = Status, y = cerebral.hbdeox.alpha)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.hbdeox.alpha") +
                theme(legend.position = "none")
        
        # Oxy Hb
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.hbox)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.hbox") +
                theme(legend.position = "none")
        
        # Oxy Hb alpha
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.hbox.alpha)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.hbox.alpha") +
                theme(legend.position = "none")
        
        # Oxy Hb/totalHb
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = `%OxyHb`)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("OxyHb/totalHb") +
                theme(legend.position = "none")
        
        # sat
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.sat)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.sat") +
                theme(legend.position = "none")
        
        # sat alpha
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.sat.alpha)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.sat.alpha") +
                theme(legend.position = "none")
        
        # THC
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.thc)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.thc") +
                theme(legend.position = "none")
        
        # THC alpha
        df_hbox %>% 
                
                ggplot(aes(x = Status, y = cerebral.thc.alpha)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.thc.alpha") +
                theme(legend.position = "none")
        
        # THC norm std
        df_hbox %>%
                
                ggplot(aes(x = Status, y = cerebral.thc.norm.std)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("cerebral.thc.norm.std") +
                theme(legend.position = "none")
        
        # Haptoglobin
        df_hbox %>%
                
                ggplot(aes(x = Status, y = log10(Haptoglobin..mg.dl.))) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("Haptoglobin (mg/dl)") +
                theme(legend.position = "none")
        
        # Hemoglobin
        df_hbox %>% filter(Hemoglobin..uM. < 80) %>% # remove outlier for visualization
                
                ggplot(aes(x = Status, y = Hemoglobin..uM.)) +
                geom_violin(aes(color = Status)) +
                geom_jitter(position = position_jitter(0.2), shape = 1) +
                stat_summary(fun = "median", geom = "crossbar", aes(color = Status), size = 0.2, width = 0.5) +
                
                theme_bw() +
                ggtitle("Hemoglobin (uM)") +
                theme(legend.position = "none")
        

        


# Kruskal-Wallis tests ----------------------------------------------------

       
        
        # Deoxy Hb
        kruskal.test(cerebral.hbdeox ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.hbdeox, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.hbdeox ~ Status, data = df_hbox, method = "bh")
        
        # Deoxy Hb alpha
        kruskal.test(cerebral.hbdeox.alpha ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.hbdeox.alpha, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.hbdeox.alpha ~ Status, data = df_hbox, method = "bh")
        
        # Oxy Hb 
        kruskal.test(cerebral.hbox ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.hbox, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.hbox ~ Status, data = df_hbox, method = "bh")
        
        # Oxy Hb alpha
        kruskal.test(cerebral.hbox.alpha ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.hbox.alpha, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.hbox.alpha ~ Status, data = df_hbox, method = "bh")
        
        # OxyHb/totalHb
        kruskal.test( cerebral.hbox/(cerebral.hbox + cerebral.hbdeox) ~ Status, data = df_hbox)

        # sat
        kruskal.test(cerebral.sat ~ Status, data = df_hbox)
        
        # sat alpha
        kruskal.test(cerebral.sat.alpha ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.sat.alpha, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.sat.alpha ~ Status, data = df_hbox, method = "bh")
        
        # THC
        kruskal.test(cerebral.thc ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.thc, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.thc ~ Status, data = df_hbox, method = "bh")
        
        # THC alpha
        kruskal.test(cerebral.thc.alpha ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.thc.alpha, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.thc.alpha ~ Status, data = df_hbox, method = "bh")
        
        # THC norm std
        kruskal.test(cerebral.thc.norm.std ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$cerebral.thc.norm.std, df_hbox$Status, p.adjust.method = "BH")
        dunnTest(cerebral.thc.norm.std ~ Status, data = df_hbox, method = "bh")
        
        # Haptoglobin
        kruskal.test(Haptoglobin..mg.dl. ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$Haptoglobin..mg.dl., df_hbox$Status, p.adjust.method = "BH")
        dunnTest(Haptoglobin..mg.dl. ~ Status, data = df_hbox, method = "bh")
        
        # Hemoglobin
        kruskal.test(Hemoglobin..uM. ~ Status, data = df_hbox)
        pairwise.wilcox.test(df_hbox$Hemoglobin..uM., df_hbox$Status, p.adjust.method = "BH")
        dunnTest(Hemoglobin..uM. ~ Status, data = df_hbox, method = "bh")
        

        


# spearman correlation with variables of interest -------------------------

        

        cor.test(df_hbox$Hematocrit, df_hbox$cerebral.thc.alpha, method = "spearman")
        cor.test(df_hbox$Lactate, df_hbox$cerebral.thc.alpha, method = "spearman")
        cor.test(df_hbox$Glucose, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$Sex, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$Temperature, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$Parasites, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$VO2, df_hbox$cerebral.thc.alpha, method = "spearman")
        cor.test(df_hbox$Whole.Blood.Nitrite, df_hbox$cerebral.thc.alpha, method = "spearman")
        cor.test(df_hbox$BP.Diastolic, df_hbox$cerebral.thc.alpha, method = "spearman")
        cor.test(df_hbox$BP.Systolic, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$Deoxygenation.Rate, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$HR, df_hbox$cerebral.thc.alpha, method = "spearman")
        # cor.test(df_hbox$Age, df_hbox$cerebral.thc.alpha, method = "spearman")
        
        # I don't know why Amit did that.... think I'd rather run kruskal wallis on everything
        
        # deoxyHb alpha
        aov(cerebral.hbdeox.alpha ~ Status + Age, data = df_hbox) %>% summary()
        
        # oxyHb alpha
        aov(cerebral.hbox.alpha ~ Status + Lactate, data = df_hbox) %>% summary()

        # sat
        aov(cerebral.sat ~ Status + Hematocrit + Lactate + Glucose + Age + Temperature +
                    Whole.Blood.Nitrite + BP.Diastolic + BP.Systolic, data = df_hbox) %>% summary() # nothing lol idk

        # sat alpha
        aov(cerebral.sat.alpha ~ Status + Hematocrit + Lactate + Glucose + Age + Temperature +
                    Whole.Blood.Nitrite + BP.Diastolic + BP.Systolic, data = df_hbox) %>% summary() # nothing here either
        
        # thc
        aov(cerebral.thc ~ Status + Hematocrit, data = df_hbox) %>% summary()

        # thc alpha
        aov(cerebral.thc.alpha ~ Status + Age + Temperature + BP.Diastolic, data = df_hbox) %>% summary()
        
        # thc norm std
        aov(cerebral.thc.norm.std ~ Status + Hematocrit + Lactate + Glucose + Age + Temperature +
                    Whole.Blood.Nitrite + BP.Diastolic + BP.Systolic, data = df_hbox) %>% summary() # nothing


        
        
        
        
        
        
        
        
                
