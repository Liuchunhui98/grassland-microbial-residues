# =========================================================
# =========================================================
Linear mixed-effects models (LMM)
# =========================================================
library(Matrix)
library(lme4)
library(sjPlot)
library(lmerTest)
library(readxl)
library(vegan)
library(ade4)
data<- read.csv('y.csv', row.names = 1)
print(colnames(data))
data.intercept= lmer (SOC~ Year + Layer+ Year * Layer+(1 | plot), data=data, REML = T)
summary (data.intercept)
plot(data.intercept)
tab_model(data.intercept)##fixed effects coefs
fixef(data.intercept)##random effects coefs
ranef(data.intercept)##all coef.
coef(data.intercept)
table(data$Year, data$Layer)

anova(data.intercept,type="III")
summary(data.intercept)
AIC(data.intercept, data.null)
anova(data.null, data.intercept, test = "Chisq") 
plot(data.intercept)
qqnorm(residuals(data.intercept)) 
tab_model(data.intercept)##fixed effects coefs
fixef(data.intercept)##random effects coefs
ranef(data.intercept)##all coef.
coef(data.intercept)
r.squaredGLMM(data.intercept) 
anova(data.intercept,type="I")
summary(data.intercept)

null_model = lmer(MNC~ (1 | plot), data = data, REML = T)  

qqnorm(residuals(data.intercept)) 
qqline(residuals(data.intercept), col = "red")

qqnorm(residuals(data.null)) 
qqline(residuals(data.null), col = "blue")

# Shapiro-Wilk
shapiro.test(residuals(data.intercept)) 
shapiro.test(residuals(data.null))  

models <- list(data.null, data.intercept)

model_selection <- model.sel(models)

print(model_selection)
head(model_selection)

model_selection_clean <- model_selection[!is.infinite(model_selection$delta) & !is.na(model_selection$delta), ]

best_models <- model_selection_clean[which(model_selection_clean$delta < 2), ]


print(best_models)
# =========================================================
# Random forest analysis
# =========================================================

# Random forest models
setwd('D:/R/yws')
#读取数据
df <- read.csv("D:/R/yws/y.csv", fileEncoding = "GBK")

print(colnames(df))
df1 <- na.omit(df)
set.seed(123)
df.rf <- randomForest(EPS.protein5~ aboveground.biomass + Root.biomass +
                        pH + EC + NO3 + NH4 + TN + TP +
                        SOC + MAOC + POC +passive.C + active.C + slow.C+
                        DOC + DON + DOP +
                        MBC + MBN + MBP + αG + βG + βX +
                        NAG + LAP + AP + Length + Angle + CUEest+ 
                        JS + Spore.density + Chao1 + Shannon,
                      data= df, 
                      ntree = 800,
                      importance=TRUE, proximity=TRUE)#一般都选择T
importance(df.rf, decreasing = F)
#检验每个变量的重要性并构建数据框
set.seed(123)
richness.rfP<-rfPermute(EPS.protein5~ aboveground.biomass + Root.biomass +
                          pH + EC + NO3 + NH4 + TN + TP +
                          SOC + MAOC + POC +passive.C + active.C + slow.C+
                          DOC + DON + DOP +
                          MBC + MBN + MBP + αG + βG + βX +
                          NAG + LAP + AP + Length + Angle + CUEest+ 
                        JS + Spore.density + Chao1 + Shannon, data = df, 
                        ntree = 800,
                        nrep = 800, 
                        num.cores = 2)#数据过多可以加大这个数字，多线程运行
richness.rfP<-rfPermute(MNC.MAOC~ aboveground.biomass + Root.biomass +
                          pH + EC + NO3 + NH4 + TN + TP +SOC+
                          MAOC + POC +passive.C + active.C + slow.C+
                          DOC + DON + DOP +
                          MBC + MBN + MBP + αG + βG + βX +
                          NAG + LAP + AP + Length + Angle + CUEest+ MNCMBC.M+
                          JS+ Spore.density + Chao1 + Shannon, data = df, 
                        ntree = 800,
                        nrep = 800, 
                        num.cores = 2)#数据过多可以加大这个数字，多线程运行
richness.data <- data.frame(importance(richness.rfP, decreasing = F))
# 创建颜色向量
colors <- c('#07326F', '#8BB1D4', '#D8E4F4', '#F3F9FE')
colors <- c('#4F6228', '#77933C', '#ACC17A', '#D7E4BD')
colors <- c('#3F014A', '#A797C7', '#DFD4E5', '#E6E6FA')

richness.data <- mutate(richness.data,
                        label = ifelse(X.IncMSE.pval < 0.001, '***',
                                       ifelse(X.IncMSE.pval < 0.01, '**',
                                              ifelse(X.IncMSE.pval < 0.05, '*', ''))))
#添加变量名
richness.data$name <- rownames(richness.data)
#变量名排序
richness.data$name <- factor(richness.data$name,
                             levels = richness.data$name)
# 创建颜色列
richness.data$color <- ifelse(richness.data$label == '***', colors[1],
                              ifelse(richness.data$label == '**', colors[2],
                                     ifelse(richness.data$label == '*', colors[3], colors[4])))
# 创建自定义排序列表
custom_order <- c("aboveground.biomass","Root.biomass",
                  "pH","EC","NO3","NH4","TN","TP",
                  "SOC","MAOC","POC","passive.C","active.C","slow.C",
                  "DOC","DON","DOP",
                  "MBC","MBN","MBP","αG","βG","βX",
                  "NAG","LAP","AP","Length","Angle","CUEest","MNCMBC.M",
                  "JS","Spore.density","Chao1","Shannon"
                  )
#倒序
custom_order_reversed <- rev(custom_order)
# 将'name'列转换为因子，并设置其级别为自定义排序
richness.data$name <- factor(richness.data$name, levels = custom_order_reversed)

# 处理X.IncMSE中小于0的值
richness.data$X.IncMSE[richness.data$X.IncMSE < 0] <- 0
p <- ggplot(richness.data, aes(name, X.IncMSE)) +
  geom_bar(aes(fill = label),
           stat = 'identity') +
  scale_fill_manual(values = brewer.pal(6, "Accent")) +
  geom_text(aes(y = X.IncMSE + 0.5,
                label = label)) +
  theme_classic() +
  labs(x = '',  # 隐藏x轴的标签
       y = 'Increase in MSE(%)') +
  theme(legend.position = '') +
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))#X从0开始,y垂直0点 

p

p1 <- ggplot(richness.data, aes(name, X.IncMSE, fill = color)) +
  geom_bar(stat = 'identity') +
  scale_fill_identity(guide = FALSE) +  # 使用自定义颜色，且不显示图例
  geom_text(aes(y = X.IncMSE + 1.5, label = label)) +
  theme_classic() +
  labs(x = '',
       y = 'IncMSE(%)   POM EPS.protein4 ') +
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank(),  # 隐藏y轴的文本标签
        axis.ticks.y = element_blank(), # 隐藏y轴的刻度线
        axis.title.y = element_blank()) + # 隐藏y轴的标题
  coord_flip() +  # 坐标轴翻转
  scale_y_continuous(expand = c(0, 0))#X从0开始,y垂直0点scale_y_continuous(breaks = seq(0, 25, by = 10), limits = c(0, 25), expand = c(0, 0)

p1

#读取环境变量和物种丰度矩阵
env <- df$MNC.MAOC
# 列名列表
cols_to_select <- c("Shannon","Chao1","Spore.density","JS","MNCMBC.M","CUEest",
                    "Angle","Length","AP","LAP","NAG",
                    "βX","βG","αG","MBP","MBN",
                    "MBC","DOP","DON","DOC","slow.C",
                    "active.C","passive.C","POC","MAOC","SOC",
                    "TP","TN","NH4","NO3","EC","pH",
                    "Root.biomass","aboveground.biomass")

# 选择特定的列
spe <- df[, cols_to_select]
#spe <- spe[rownames(env), ]
####环境变量和物种丰度的相关性分析
library(psych)
library(reshape2)
#可通过 psych 包函数 corr.test() 执行
#这里以 pearson 相关系数为例，暂且没对 p 值进行任何校正（可以通过 adjust 参数额外指定 p 值校正方法）
pearson <- corr.test(env, spe, method = 'pearson', adjust = 'none')
r <- data.frame(pearson$r)  #pearson 相关系数矩阵
p <- data.frame(pearson$p)  #p 值矩阵
#结果整理以便于作图
r$env <- rownames(r)
p$env <- rownames(p)
r <- melt(r, id = 'env')
p <- melt(p, id = 'env')
pearson <- cbind(r, p$value)
colnames(pearson) <- c('env', 'spe', 'pearson_correlation', 'p.value')
pearson$spe <- factor(pearson$spe, levels = colnames(spe))
head(pearson)  #整理好的环境变量和物种丰度的 pearson 相关性统计表

#ggplot2 作图，绘制环境变量和物种丰度的 pearson 相关性热图

# 自定义顺序列表
custom_order <- c("Shannon","Chao1","Spore.density","JS","MNCMBC.M","CUEest",
                  "Angle","Length","AP","LAP","NAG",
                  "βX","βG","αG","MBP","MBN",
                  "MBC","DOP","DON","DOC","slow.C",
                  "active.C","passive.C","POC","MAOC","SOC",
                  "TP","TN","NH4","NO3","EC","pH",
                  "Root.biomass","aboveground.biomass")
# 将 'env' 列转换为因子，并按照自定义顺序排序
pearson$env <- factor(pearson$env, levels =custom_order)

# 使用排序后的数据绘制热图
#liu文字
p2 <- ggplot() +
  geom_tile(data = pearson, aes(x = env, y = spe, fill = pearson_correlation)) +
  scale_fill_gradientn(colors = c('#A63446', 'white', '#747E95'), limit = c(-1, 1)) +
  theme(
    panel.background = element_blank(),   # 完全去掉背景
    panel.border = element_blank(),       # 去掉边框
    panel.grid = element_blank(),         # 去掉网格线
    axis.text.x = element_blank(),        # x轴文字可根据需要隐藏
    axis.ticks.x = element_blank(),       # x轴刻度线
    axis.ticks.y = element_blank(),       # y轴刻度线
    axis.text.y = element_text(color = 'black'),  # y轴文字保留
    legend.key = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(y = '', x = '', fill = '')

p2


p2 <- ggplot() +
  geom_tile(data = pearson, aes(x = env, y = spe, fill = pearson_correlation)) +
  scale_fill_gradientn(colors = c('#A63446', 'white', '#747E95'), limit = c(-1, 1)) +
  theme(panel.grid = element_line(), panel.background = element_rect(color = 'black'), 
        legend.key = element_blank(), legend.position = "bottom",
        #legend.margin = margin(t = -1, unit = "cm"),  
        #legend.box.margin = margin(t = 0, unit = "cm"),  
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio=1) +
  theme(axis.text.x = element_blank(),    
        axis.text.y = element_blank(),      
        axis.ticks = element_blank(),       
        axis.title.x = element_blank(),     
        axis.title.y = element_blank(),     
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y = '', x = '', fill = '')

p2

pearson[which(pearson$p.value<0.001),'sig'] <- '***'
pearson[which(pearson$p.value<0.01 & pearson$p.value>0.001),'sig'] <- '**'
pearson[which(pearson$p.value<0.05 & pearson$p.value>0.01),'sig'] <- '*'
head(pearson)  

p3 <- p2 +
  geom_text(data = pearson, aes(x = env, y = spe, label = sig), size = 4)

p3

p3+p1



xdata <- df[, c("aboveground.biomass","Root.biomass",
                      "pH","EC","NO3","NH4","TN","TP",
                      "SOC","MAOC","POC","passive.C","active.C","slow.C",
                      "DOC","DON","DOP",
                      "MBC","MBN","MBP","αG","βG","βX",
                      "NAG","LAP","AP","Length","Angle","CUEest",
                      "JS","Spore.density","Chao1","Shannon")]
ydata <- df$EPS.protein5

library(randomForest)
set.seed(123)
df.rf <- randomForest(x = xdata, y = ydata,
                      ntree = 800,
                      importance = TRUE,
                      proximity = TRUE)


# install.packages("permimp")  
library(permimp)
library(randomForest)

set.seed(123)
df.rf <- randomForest(x = xdata, y = ydata, ntree = 500, importance = TRUE)


set.seed(123)
#perm_results <- permimp(df.rf, conditional = TRUE, nperm = 99)
#print(perm_results)


pred <- predict(df.rf, xdata)


R2 <- 1 - sum((ydata - pred)^2) / sum((ydata - mean(ydata))^2)
R2

nperm <- 99
perm_R2 <- numeric(nperm)
set.seed(123)
for(i in 1:nperm){
  y_perm <- sample(ydata)
  rf_perm <- randomForest(x = xdata, y = y_perm, ntree = 500)
  pred_perm <- predict(rf_perm, xdata)
  perm_R2[i] <- 1 - sum((y_perm - pred_perm)^2) / sum((y_perm - mean(y_perm))^2)
}

P_value <- mean(perm_R2 >= R2)
R2; P_value

# =========================================================
#  Partial least squares path modeling (PLS-PM)
# =========================================================

# PLS-PM analyses
setwd("D:/R/yws")
library(vegan)
install.packages("devtools") 
library(devtools)
install_github("gastonstat/plspm")
install.packages("plspm")
library(plspm)
#AMs
data1<- read.csv("D:/R/yws/y.csv", fileEncoding = "GBK")

data1.scaled <- as.data.frame(scale(data1[sapply(data1, is.numeric)]))
colnames(data1.scaled)
plant=c(0,0,0,0,0,0,0,0)
soil= c(1,0,0,0,0,0,0,0)
AMF= c(1,1,0,0,0,0,0,0)
MB=c(1,1,1,0,0,0,0,0)
Mlife = c(1,1,1,1,0,0,0,0)
EPS= c(1,1,1,1,1,0,0,0)
MNC = c(1,1,1,1,1,0,0,0)
MAOC= c(1,1,1,1,1,1,1,0)
path_matrix = rbind(plant,soil,AMF,MB,Mlife,EPS,MNC,MAOC)
colnames(path_matrix) = rownames(path_matrix)
innerplot(path_matrix)
sat_blocks = list(c(5,6),c(8,9,11),c(40:41),c(22:24),c(25:26,33),c(58,59),c(62),14)
sat_mod <- rep("A", length(sat_blocks))
satpls = plspm(data1.scaled , path_matrix, sat_blocks, modes = sat_mod,  scaled = FALSE)
innerplot(satpls)
summary(satpls)
data1.scaled$CUEest = -1 * data1.scaled$CUEest

fit <- lm(MAOC ~ EPS.polysaccharide4 , data = data1.scaled)
library(car)
library(carData)
vif(lm(MB ~ soil + AMF + plant + Mlife, data = data))

colnames(data1.scaled)
plant=c(0,0,0,0,0,0,0,0,0)
soil= c(1,0,0,0,0,0,0,0,0)
AMF= c(1,1,0,0,0,0,0,0,0)
MB=c(1,1,1,0,0,0,0,0,0)
Mlife = c(1,1,1,1,0,0,0,0,0)
EPSb= c(1,1,1,1,1,0,0,0,0)
EPSt= c(1,1,1,1,1,0,0,0,0)
MNC = c(1,1,1,1,1,0,0,0,0)
MAOC= c(0,0,0,0,0,1,1,1,0)
path_matrix = rbind(plant,soil,AMF,MB,Mlife,EPSb,EPSt,MNC,MAOC)
colnames(path_matrix) = rownames(path_matrix)
innerplot(path_matrix)
sat_blocks = list(c(5:6),c(9,11,20),c(40:41),c(22:24),c(38:39),59,58,c(62),14)
sat_mod <- rep("A", length(sat_blocks))
satpls = plspm(data1.scaled , path_matrix, sat_blocks, modes = sat_mod,  scaled = FALSE)
innerplot(satpls)
summary(satpls)

writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
plant=c(0,0,0,0,0,0,0,0,0)
soil= c(1,0,0,0,0,0,0,0,0)
AMF= c(1,1,0,0,0,0,0,0,0)
MB=c(1,1,1,0,0,0,0,0,0)
Mlife = c(1,1,1,1,0,0,0,0,0)
EPSb= c(1,1,1,1,1,0,0,0,0)
EPSt= c(1,1,1,1,1,0,0,0,0)
MNC = c(1,1,1,1,1,0,0,0,0)
POC= c(0,0,0,0,0,1,1,1,0)
path_matrix = rbind(plant,soil,AMF,MB,Mlife,EPSb,EPSt,MNC,POC)
colnames(path_matrix) = rownames(path_matrix)
innerplot(path_matrix)
sat_blocks = list(c(5:6),c(9,20),c(40:41),c(22:24),c(38:39),72,71,c(75),15)
sat_mod <- rep("A", length(sat_blocks))
satpls = plspm(data1.scaled , path_matrix, sat_blocks, modes = sat_mod,  scaled = FALSE)
innerplot(satpls)
summary(satpls)
data1.scaled$CUEest = -1 * data1.scaled$CUEest

vif( lm(MAOC ~ MBC+MBN+MBP+EC+NO3+TN , data = data1.scaled) )
