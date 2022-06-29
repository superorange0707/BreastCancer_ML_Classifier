# install and load required packages
install.packages("Amelia")
install.packages("psych")
install.packages("corrplot")
install.packages("gridExtra")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("scorecard")
install.packages("caret")
install.packages("pROC")
install.packages("scatterplot3d")
library(Amelia)
library(dplyr)
library(ggplot2)
library(psych)
library(corrplot)
library(gridExtra)
library(FactoMineR)
library(factoextra)
library(scorecard)
library(caret)
library(pROC)
library(scatterplot3d)

# load data
wbcd = read.csv('original data.csv',stringsAsFactors = F)
wbcd = as.data.frame(wbcd)


# Data Analysis And Data Pre processing
# check the missing data
missmap(wbcd, col = c("white", "lightblue"))


# check the Outliers
for (i in 1:ncol(wbcd))
  print(unique(wbcd[,i]))

# find value ? of bare_nucleoli
wbcd1 = wbcd[wbcd["bare_nucleoli"] =='?',]

# replace these missing value(using mode)
# get the mode of bare_nucleoli
mode = as.numeric(names(table(wbcd$bare_nucleoli))[which.max(table(wbcd$bare_nucleoli))])
wbcd$bare_nucleoli[wbcd$bare_nucleoli=="?"]= mode
wbcd$bare_nucleoli=as.numeric(wbcd$bare_nucleoli)
wbcd1$bare_nucleoli[wbcd1$bare_nucleoli=="?"]= mode

# boxplot
ggplot(stack(wbcd[,1:10]), aes(x = ind, y = values)) + geom_boxplot() + xlab("Attributes") + ylab("Vaule") + ggtitle("Boxplot of Attributes") + theme(plot.title = element_text(hjust = 0.5))

# remove duplicate data
wbcd = distinct(wbcd)

# drop the id column
wbcd = wbcd[,-1]


# statistics of data
describe(wbcd)


# distribution of data
clump_thickness=ggplot(wbcd,aes(x =clump_thickness))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
size_uniformity=ggplot(wbcd,aes(x =size_uniformity))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
shape_uniformity=ggplot(wbcd,aes(x =shape_uniformity))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
marginal_adhesion=ggplot(wbcd,aes(x =marginal_adhesion))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
epithelial_size=ggplot(wbcd,aes(x =epithelial_size))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
bare_nucleoli=ggplot(wbcd,aes(x =bare_nucleoli))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
bland_chromatin=ggplot(wbcd,aes(x =bland_chromatin))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
normal_nucleoli=ggplot(wbcd,aes(x =normal_nucleoli))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
mitoses=ggplot(wbcd,aes(x =mitoses))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 1, center = 0)+geom_density()
class=ggplot(wbcd,aes(x =class))+geom_histogram(aes(y=..density..),color="#88ada6", fill="#fffbf0",alpha=.25,binwidth = 2, center = 0)+geom_density()
grid.arrange(clump_thickness, size_uniformity, shape_uniformity, marginal_adhesion,epithelial_size, bare_nucleoli,bland_chromatin,normal_nucleoli,mitoses,class,top = "Distribution of Features",nrow = 2, ncol = 5)


# correlation
corrplot(cor(wbcd), method = "square",shade.col = NA, tl.col ="black", tl.srt = 45, order = "AOE")


# feature selection
pca = PCA(wbcd, graph = FALSE)
fviz_eig(pca, addlabels = TRUE)+theme(plot.title = element_text(hjust = 0.5))
fviz_contrib(pca, choice = "var", axes = 1:3, top = 10)+theme(plot.title = element_text(hjust = 0.5))

# get new data after feature selection
wbcd_choice = wbcd[,c(1,2,3,9,10)]
wbcd_choice$class=as.factor(wbcd_choice$class)


# data division
# train data 70% and test data 30%
data_split = split_df(wbcd_choice,ratios = c(0.7,0.3))
train_set = data_split$train
test_set = data_split$test



# Classifier Models

# K-fold cross validation:10
tran_fit = trainControl(method = "cv", number = 10)


# KNN_model
# # data standardization and train model
# try different k value
knn_model = train(class~.,train_set,method = "knn", preProcess = c('center','scale'), trControl = tran_fit, tuneLength = 10 )
# draw the plot of accuracy of different k
knn_k = as.data.frame(knn_model$results)
ggplot(data=knn_k, mapping = aes(x=factor(k), y=Accuracy, group = 1)) + geom_line(color="steelblue") + xlab("K Vaule") + ylab("Accuracy") +ggtitle("Accuracy of different K value")+theme(plot.title = element_text(hjust = 0.5))

# prediction
pred_knn = predict(knn_model,newdata = test_set)
true =test_set$class
confusionMatrix(pred_knn, true)




# SVM_Radial
# with optimization
svm_Radial_model = train(class~.,train_set,method = "svmRadial", preProcess = c('center','scale'), trControl = tran_fit, tuneGrid = expand.grid(C= c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1),sigma =c(0, 0.05:2)))
# draw the plot of different parameters of SVM
svm_pa = as.data.frame(svm_Radial_model$results)
scatterplot3d(svm_pa[,1:3], pch = 16, color="steelblue")

# prediction
pred_svm = predict(svm_Radial_model,newdata = test_set)
confusionMatrix(pred_svm, true)




# Random_Forest
rf_model = train(class~.,train_set,method = "rf", preProcess = c('center','scale'), trControl = tran_fit, tuneGrid = expand.grid(.mtry= c(1:20)), ntree = 800, nodesize = 2, maxnodes = 10)
rf_model
# draw the plot of accuracy of different mtry
rf_mtry = as.data.frame(rf_model$results)
ggplot(data=rf_mtry, mapping = aes(x=factor(mtry), y=Accuracy, group = 1)) + geom_line(color="steelblue") + xlab("mtry Vaule") + ylab("Accuracy") +ggtitle("Accuracy of different mtry value")+theme(plot.title = element_text(hjust = 0.5))

# prediction
pred_rf = predict(rf_model,newdata = test_set)
confusionMatrix(pred_rf, true)



# Naive bayes
nb_model = train(class~.,train_set,method = "nb", preProcess = c('center','scale'), trControl = tran_fit, tuneGrid = expand.grid(.fL=c(0), .usekernel=c(FALSE),.adjust=0.5))
nb_model

# prediction
pred_nb = predict(nb_model,newdata = test_set)
confusionMatrix(pred_nb, true)


# logistic regression
lr_model = train(class~.,train_set, method = "glm", preProcess = c('center','scale'), trControl = tran_fit)
lr_model

# prediction
pred_lr = predict(lr_model,newdata = test_set)
confusionMatrix(pred_lr, true)


# Evaluation
# Accuracy
model_name = c('KNN', 'SVM', 'RF', 'NB', 'LR')
model_accuracy=c(0.9545,0.9591,0.9636,0.95,0.9584)
model_kappa =c(0.8908,0.9004,0.9149,0.8815,0.9119)
model_info = data.frame(model_name, model_accuracy, model_kappa)
ggplot(data=model_info, mapping = aes(x=model_name, y=model_accuracy, group = 1)) + geom_line(color="steelblue") + geom_point()+ xlab("Model") + ylab("Accuracy") +ggtitle("Accuracy of different Model")+theme(plot.title = element_text(hjust = 0.5))


# Consistency
ggplot(data=model_info, mapping = aes(x=model_name, y=model_kappa, group = 1)) + geom_line(color="steelblue") + geom_point()+ xlab("Model") + ylab("Consistency") +ggtitle("Consistency of different Model")+theme(plot.title = element_text(hjust = 0.5))


# ROC: to evaluate the sensitivity and specificity
roc_knn = roc(test_set$class,factor(pred_knn,ordered = T))
roc_svm = roc(test_set$class,factor(pred_svm,ordered = T))
roc_rf = roc(test_set$class,factor(pred_rf,ordered = T))
roc_nb = roc(test_set$class,factor(pred_nb,ordered = T)) 
roc_lr = roc(test_set$class,factor(pred_lr,ordered = T)) 
# draw multiple-models' roc curve
plot(roc_knn,col="red",main="ROC curve of different Models",legacy.axes=T)
plot(roc_svm, add=TRUE, col="blue")
plot(roc_rf, add=TRUE, col="green")
plot(roc_nb, add=TRUE, col="black")
plot(roc_lr, add=TRUE, col="yellow")
legend("bottomright", legend=c("roc_knn","roc_svm","roc_rf","roc_nb","roc_lr"),col=c("red","blue","green","black","yellow"),lty=1, lwd=3, cex =0.6,bty="n",x.intersp=0.5,y.intersp=0.5,text.font=10,seg.len=1,inset=c(-0.1,0))


# AUC
model_AUC = c(roc_knn$auc, roc_svm$auc, roc_rf$auc, roc_nb$auc, roc_lr$auc)
model_info$model_AUC = model_AUC
ggplot(data=model_info, mapping = aes(x=model_name, y=model_AUC, group = 1)) + geom_line(color="steelblue") + geom_point()+ xlab("Model") + ylab("AUC") +ggtitle("AUC of different Model")+theme(plot.title = element_text(hjust = 0.5))
