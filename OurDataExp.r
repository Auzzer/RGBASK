# -*- encoding=utf-8 -*-

# 本脚本用来实现案例分析

library(readxl)


plot_clustering <- function(X, res, main){
    pch_set=c(1,2,3,5,7,8,9)
    plot(X[,1], X[,2],xlab = 'the cumulative confirmed number',ylab='the cumulative death toll',main=main, col='white')
    for(i in 1:length(res_grid_svm_knn)){
        if(res[i] == 0){
            points(x = X[i,1], y = X[i,2], col=res[i]+1, pch=15)
        }else{
            points(x = X[i,1], y = X[i,2], col=res[i]+1, pch=pch_set[res[i]])
        }
    }
}



# 设置工作路径
setwd('E:/projects/grid_svm_knn/data/')
setwd('D:/个人/坚果云同步/我的坚果云/论文/grid_svm_knn/data')
setwd('E:/坚果云同步/我的坚果云/论文/grid_svm_knn/data/')
# 读取数据
mydata <- readxl::read_excel('china_provincedata.xlsx')
mydata <- as.data.frame(mydata)
# 截取需要的数据
usedvar <- c('confirmedCount', 'confirmedIncr', 'dateId', 'provinceShortName')
mydata <- mydata[,usedvar]
mydata <- mydata[mydata$dateId == 20200208, ]
# 绘制图像
plot(mydata[-18,]$confirmedCount, mydata[-18,]$confirmedIncr)
X <- as.matrix(mydata[,c('confirmedCount', 'confirmedIncr')], rownames.force = FALSE)

res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=1, k=3, Index1=1)
print(res_grid_svm_knn)
plot_clustering(X, res_grid_svm_knn, 'Result')

# 由于相差太大，需要剔除噪声点再次聚类
X_new <- X[-which(res_grid_svm_knn == 2),]
res_grid_svm_knn <- grid_svm_knn(X_new, gamma=50,cost=1, k=3, Index1=6)
print(res_grid_svm_knn)
plot_clustering(X_new, res_grid_svm_knn, 'Result')


mydata_new <- mydata[-which(res_grid_svm_knn == 2),]
mydata_new$provinceShortName[which(res_grid_svm_knn == 5)]



















# 新案例


set.seed(123)
n1 <- 1000  # 上圆弧的样本点数量
n2 <- 200  # 下圆弧的样本点数量
n3 <- 300
X1 <- matrix(nrow = n1,ncol = 2)
mu <- c(0,0)
sigma <- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
for(i in 1:n1){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if(((x[1])^2 + x[2]^2) >= 1 & ((x[1])^2 + x[2]^2) <= 2 & x[2] < 0){
            X1[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
}

mu <- c(2,2)
sigma <- matrix(c(0.2,0,0,0.2),nrow = 2,ncol = 2)
X2 <- mvrnorm(n2,mu,sigma)

mu <- c(-2,2)
sigma <- matrix(c(0.2,0,0,0.2),nrow = 2,ncol = 2)
X3 <- mvrnorm(n3,mu,sigma)


# 合并数据
X <- rbind(X1,X2,X3)
g <- c(rep(1,n1),rep(2,n2),rep(3,n3))


# 噪声数据
n_noise <- 50
X_noise <- matrix(ncol=2,nrow = n_noise)
for(i in 1:n_noise){
    min_dist <- 0
    while(min_dist < 0.1){
        noise_x<- runif(1,-2,2)
        noise_y <- runif(1,-2,2)
        noise_one <- matrix(rep(c(noise_x, noise_y), times=dim(X)[1]), ncol = 2,byrow = TRUE)
        min_dist <- min(dist(noise_one, X))
    }
    X_noise[i,] <- noise_one[1,]
}


X <- rbind(X,X_noise)
g <- c(g,rep(4,n_noise))



# 无标签绘图
plot(X,pch=20,col='blue',main='DataSet2',xlab = expression(X[1]),ylab=expression(X[2]))
# 有标签绘图
plot(X,col=g,pch=20,main='组合形状数据',xlab = expression(X[1]),ylab=expression(X[2]))
# 保存数据
dataset4 <- cbind(X,g)
colnames(dataset4) <- c('X1','X2','g')
write.csv(dataset4, 'dataset4.csv', fileEncoding = 'utf-8', row.names = FALSE)

res_db <- dbscan(X,eps=0.2,minPts=8)
res_grid_svm_knn <- grid_svm_knn(X, gamma=10, cost=0.023, k=18, Index1=1213)









setwd('E:/坚果云同步/我的坚果云/论文/grid_svm_knn/data/')
# 读取数据
mydata <- read.csv('12-10-2020.csv')
# 截取需要的数据
usedvar <- c("Country_Region", "Confirmed", "Deaths")
mydata <- mydata[,usedvar]
# 对数据按照国家进行聚合处理
mydata_new <- rowsum(mydata[,c(2,3)], group=mydata$Country_Region)
data_name <- row.names(mydata_new)
# 绘制图像
plot(mydata_new)
# 整理一下数据
X <- as.matrix(mydata_new, rownames.force = FALSE)
# 对数据进行标准化处理
X_scale <- scale(X)



# 绘制图像
plot(mydata[-18,]$confirmedCount, mydata[-18,]$confirmedIncr)
X <- as.matrix(mydata[,c('confirmedCount', 'confirmedIncr')], rownames.force = FALSE)

res_grid_svm_knn <- grid_svm_knn(X_scale, gamma=1,cost=10, k=10, Index1=1)
print(res_grid_svm_knn)
plot_clustering(X_scale, res_grid_svm_knn, '聚类结果')

# 由于相差太大，需要剔除噪声点再次聚类
X_new <- X[-which(res_grid_svm_knn == 2),]
res_grid_svm_knn <- grid_svm_knn(X_new, gamma=50,cost=1, k=3, Index1=6)
print(res_grid_svm_knn)
plot_clustering(X_new, res_grid_svm_knn, '聚类结果')
legend('bottomright', c('other', 'Mexico', 'Brazil', 'India', 'US'), pch =c(1,2,3,5,7), col=c(2,3,4,5,6))

mydata_new <- mydata[-which(res_grid_svm_knn == 2),]
mydata_new$provinceShortName[which(res_grid_svm_knn == 2)]