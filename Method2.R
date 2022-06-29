install.packages("mice")
install.packages("Rlof")
install.packages("PerformanceAnalytics")
install.packages("CORElearn")
install.packages("caret")
install.packages("sampling")
library(mice)
library(PerformanceAnalytics)
library(Rlof)
library(tidyr)
library(CORElearn)
library(caret)
library(sampling)

# read data
df = as.data.frame(read.csv('original data.csv',stringsAsFactors = F))


# missing value
md.pattern(df)


# outlier
outlier = lof(df, k =c(5:10))
pch = rep(".",nrow(df))
pch[outlier]="+"
col = rep("black", nrow(df))
pairs(df,pch=pch, col= col)


# correlation
df_ma = as.matrix(df)
df_ma = apply(df_ma, 2, as.numeric)
chart.Correlation(df_ma, histogram=TRUE, pch=19)

# miss value
df[df$bare_nucleoli=='?', ]=NA
df = drop_na(df)


# feature:relief
df$class = as.factor(df$class)
relief = as.data.frame(attrEval(class~.,df,"Relief"))
relief$var = rownames(relief)
colnames(relief)[1]='value'
ggplot(data=relief, aes(x=var, y=value)) + 
  geom_bar(colour="black", fill="steelblue", width=.6, stat="identity") + 
  guides(fill=FALSE) +
  xlab("Variable") + ylab("Contribution value") +
  ggtitle("Relief of df")+theme(plot.title = element_text(hjust = 0.5))

# EXTRACTION
df = df[,-1]
df_data = df[,c(1,4,5,6,7,10)]


# Logistic Regression
# train and test set:Stratified Sampling
n=round(3/5*nrow(df_data)/3)
sub_train=strata(df_data,stratanames=("class"),size=rep(n,3),method="srswor")
data_train=df_data[sub_train$ID_unit,]
data_train = as.data.frame(data_train)
data_test=df_data[-sub_train$ID_unit,]
data_test = as.data.frame(data_test)
data_train$class = factor(data_train$class)
data_test$class = factor(data_test$class)


# train
tran_fit = trainControl(method = "cv", number = 5)
lr = train(class~.,data_train,method = "glm", preProcess = c('center','scale'), trControl = tran_fit)
lr

# prediction
pred_lr = predict(lr,newdata = data_test)
pred_lr
compare_test = data_test$class
confusionMatrix(pred_lr, compare_test)

