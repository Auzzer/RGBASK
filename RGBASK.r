# -*- encoding: utf-8 -*-

# 目的: 实现基于网格的SVM近邻聚类算法(GCASK)，同时进行数值模拟以及复杂度分析


##########################################################
# 设置工作路径
##########################################################
setwd('E:/projects/grid_svm_knn/data')
# setwd('D:/个人/坚果云同步/我的坚果云/论文/grid_svm_knn/data')
###########################################################
# 基于网格的SVM近邻聚类算法相关代码
###########################################################

# 对数据空间进行网格化，返回各个维度上网格化的边界
Gridding <- function(X,l_grid=0){
    # 得到X的属性个数
    p <- ncol(X)
    # 得到X样本点个数
    n <- nrow(X)
    # 计算得到各个属性维度的最大值和最小值，确定数据空间的范围
    X_range <- matrix(ncol=2,nrow=p)
    # X_range的第一列为每个维度上的最大值
    X_range[,1] <- apply(X,2,max)
    # X_range的第二列为每个维度上的最小值
    X_range[,2] <- apply(X,2,min)
    # 若X_range的第一列和第二列相同，说明该属性取值都相同，应该删除该属性
    # Index <- (X_range[,1] == X_range[,2])
    # X <- X[,!Index]
    # X_range <- X_range[!Index,]
    # p <- ncol(X)
    # 计算整个空间的体积
    V <- cumprod(X_range[,1]-X_range[,2])[p]
    # 计算小网格的边长
    if (l_grid == 0){
        l_grid <- (V/n) ** (1/p)
    }
    # 计算每个维度上网格的个数
    n_grid <- ceiling((X_range[,1]-X_range[,2]) / l_grid)
    # 用DividingLine来存放各个维度上网格的划分边界线
    DivdingLine <- matrix(ncol=p,nrow=max(n_grid)+1)
    for(i in 1:p){
        DivdingLine[1:(n_grid[i]+1),i] <- X_range[i,2] + ((0:n_grid[i]) * l_grid)
    }
    Lst <- list(DivdingLine=DivdingLine,n_grid=n_grid)
    return(Lst)
}


# 判断格子里有没有点
box_point <- function(X,border){
    Index <- 1:nrow(X)
    p <- ncol(X)
    for(i in 1:p){
        new_X <- X[Index,]
        new_X <- matrix(new_X,ncol=p)
        Index <- which( (new_X[,i]>=border[i,1]) & (new_X[,i]<=border[i,2]) )
        if(length(Index) == 0){
            return(1)
        }
    }
    return(0)
}


# 生成小格子边界的索引
get_combinations <- function(n_grid,Index_matrix,Index = 1,row=1){
    for(i in 1:n_grid[Index]){
        if(Index == length(n_grid)){
            for(j in 1:n_grid[Index]){
                Index_matrix[which(is.na(Index_matrix[,Index]))[1],Index] <- j
                row <- row + 1
            }
            return(list(Index_matrix=Index_matrix, row=row))
        }else{
            Index_matrix[row:nrow(Index_matrix),Index] <- i
            Index <- Index + 1
            Lst <- get_combinations(n_grid,Index_matrix,Index,row)
            Index <- Index - 1
            row <- Lst$row
            Index_matrix <- Lst$Index_matrix
        }
    }
    return(Lst)
}


# 对稀疏区域进行插点
get_nothing_point <- function(X,l_grid=0){
    # 得到网格边界
    Lst <- Gridding(X,l_grid)
    DivdingLine <- Lst$DivdingLine
    n_grid <- Lst$n_grid
    # 得到网格索引矩阵
    Index_matrix <- matrix(ncol=length(n_grid),nrow=cumprod(n_grid)[length(n_grid)])
    Lst <- get_combinations(n_grid,Index_matrix,Index = 1,row=1)
    Index_matrix <- Lst$Index_matrix
    # 生成一个插点矩阵
    point_matrix <- matrix(nrow=nrow(Index_matrix),ncol=1)
    # 生成小格子border，并进行判别
    border <- matrix(ncol=2,nrow=ncol(X))
    for(i in 1:nrow(Index_matrix)){
        for(j in 1:ncol(Index_matrix)){
            border[j,1] <- DivdingLine[Index_matrix[i,j],j]
            border[j,2] <- DivdingLine[Index_matrix[i,j]+1,j]
        }
        point_matrix[i] <- box_point(X,border)
    }
    Index_matrix <- cbind(Index_matrix,point_matrix)
    # 得到稀疏区域点的坐标
    Index_matrix <- Index_matrix[Index_matrix[,ncol(Index_matrix)]==1 ,]
    nothing_point <- matrix(ncol=ncol(X),nrow=nrow(Index_matrix))
    for(i in 1:nrow(Index_matrix)){
        for(j in 1:ncol(X)){
            nothing_point[i,j] <- (DivdingLine[Index_matrix[i,j],j] + DivdingLine[Index_matrix[i,j]+1,j]) / 2
        }
    }
    return(nothing_point)
}

# 生成有监督学习的数据
get_newdata <- function(X,l_grid=0){
    nothing_point <- get_nothing_point(X,l_grid)
    label0 <- matrix(0,nrow=nrow(nothing_point),ncol=1)
    X1 <- cbind(nothing_point,label0)
    label1 <- matrix(1,nrow = nrow(X),ncol=1)
    X2 <- cbind(X,label1)
    newdata <- rbind(X1,X2)
    return(newdata)
}


# 找到最近的k个点
nearest_k_point <- function(x,X,C_S,k){
    x <- matrix(rep(x,nrow(X)), byrow = TRUE, ncol = ncol(X))
    dis_matrix <- diag((x-X) %*% t(x-X))
    dis_matrix[C_S] <- max(dis_matrix)+1
    Index <- c()
    kk <- sum(!((1:nrow(X)) %in% C_S))
    k <- min(k,kk)
    if(k==0){return(Index)}
    for(i in 1:k){
        Index <- c(Index,which.min(dis_matrix))
        dis_matrix[Index[i]] <- max(dis_matrix)+1
    }
    return(Index) # Index中最后的那个就是最远的那个
}


# 判断最近的k个点中有几个要生成新类
new_class_point <- function(x,X,Index,model){
    x <- matrix(rep(x,length(Index)), byrow = TRUE, ncol = ncol(X))
    med <- as.data.frame((x + X[Index,])/2)
    colnames(med) <- colnames(X)
    pred <- predict(model,med)
    # 需要生成新类的点
    newclass_point <- Index[!(as.logical(as.numeric(pred)-1))]
    # 不需要生成新类的点
    oldclass_point <- Index[as.logical(as.numeric(pred)-1)]
    point <- list(newclass_point=newclass_point, oldclass_point=oldclass_point)
    return(point)
}


NewClass <- function(class_list,k,Index){
    # k表示当前有几个类
    class_list[[k]] <- Index
    return(class_list)
}


# 函数：合并到上一类
CombindClass <- function(class_list,k,Index){
    # k表示当前类的个数
    class_list[[k]] <- c(class_list[[k]],Index)
    return(class_list)
}


# 融合形成新类
Combind <- function(X,C,C_S,model,class_list,k,k_class=0,Index1=1){
    # 若是第一个点，先生成一个类
    if(k_class == 0){
        k_class <- k_class+1
        class_list <- list()
        class_list[k_class] <- C[Index1]
        C_S <- c(C_S,C[Index1])
        #print(X[Index1,])
        #points(X[Index1,1],X[Index1,2],col=k_class,pch=18)
        #pause()
    }
    
    # 有了类之后可以进行判别聚类了
    # 只要还有点没有分类就执行代码
    while(length(C_S) <= length(C)){
        #readline('请按回车')
        # 找到聚类点
        #if(exists('new_center')){
        #    x <- new_center
        #}else{
        #    x <- X[C[Index1],]
        #}
        x <- X[C[Index1],]
        # 找到最近的k个点
        Index <- nearest_k_point(x,X,C_S,k)
        if(is.null(Index)){return(class_list)}
        # 对是否每个k近邻点需要生成新类进行分析
        point <- new_class_point(x,X,Index,model)
        # 判断是不是所有的点都需要生成新类
        if(length(point$newclass_point) == length(Index)){ # 所有点都要生成新类
            # 找到最近的那个点
            Index2 <- Index[1]
            C_S <- c(C_S,C[Index2])
            # 生成新类
            k_class <- k_class + 1
            class_list <- NewClass(class_list,k_class,C[Index2])
            #points(X[Index2,1],X[Index2,2],pch=20,col=k_class)
            #pause()
            Index1 <- Index2
        }else{
            Index2 <- point$oldclass_point
            C_S <- c(C_S,C[Index2])
            # 合并类
            class_list <- CombindClass(class_list,k_class,C[Index2])
            #points(X[Index2,1],X[Index2,2],pch=20,col=k_class)
            #pause()
            #Index1 <- Index2[length(Index2)] # 选择最远的同类点作为下一个迭代点
            #this_class_points <- X[class_list[[length(class_list)]],]
            #new_center <- apply(this_class_points, 2, mean)
            Index1 <- Index2[1] # 选择最近的同类点作为下一个迭代点
        }
    }
    return(class_list)
}


# 噪声点识别
get_noise_point<- function(dat,pred){
    index_1 <- sum(dat$g==0)
    noise_index <- dat$g == 1 & pred != dat$g
    noise_point <- dat[noise_index,]
    row.names(noise_point) <- as.numeric(row.names(noise_point)) - index_1
    not_noise_point_index <- dat$g == 1 & pred == dat$g
    not_noise_point <- dat[not_noise_point_index,]
    row.names(not_noise_point) <- as.numeric(row.names(not_noise_point)) - index_1
    res <- list(noise_point=noise_point, not_noise_point=not_noise_point)
    return(res)
}

grid_svm_knn <- function(newdata, gamma ,cost, k, Index1){
    # 若是没有加载e1071则加载该包
    if(!('e1071' %in% (.packages()))){library(e1071)}
    mydata <- get_newdata(newdata,l_grid = 0)
    dat <- as.data.frame(mydata)
    colnames(dat)[ncol(dat)] <- 'g'
    dat$g<- as.factor(dat$g)
    gamma <- 50
    model <- svm(g~.,data = dat,type='C',kernel='radial',gamma=gamma,cost=cost)
    #plot(model,dat)
    pred <- predict(model,dat[,-ncol(dat)])
    X <- as.matrix(newdata[,1:2])
    res <- get_noise_point(dat,pred)
    X <- as.matrix(res$not_noise_point[,1:2])
    index_X <- as.numeric(row.names(res$not_noise_point))
    C <- 1:nrow(X)
    C_S <- c()
    class_list <- list()
    class_list <- Combind(X,C,C_S,model,class_list,k = k, Index1=Index1)
    noise_point <- as.numeric(row.names(res$noise_point))
    #clusters_result <- list(class_list = class_list, noise_point = noise_point)
    # 对结果进行整理
    pred <- rep(0, dim(newdata)[1])
    for(i in 1:length(class_list)){
        pred[index_X[class_list[[i]]]] <- i
    }
    return(pred)
}

#########################################################################################
# 动态分位数Kmeans算法
#########################################################################################
Radius.KMeans <- function(x,p,q=0.5,t=1.5){
    lst <- FILTER(x,11,q,t)
    C <- lst$C
    R <- lst$R
    k <- COMBINE(C,R)
    km <- kmeans(x,k)
    plot(x,col=km$cluster,main = 'Radius.KMeans result')
    return(list(k = k,km = km))
}

# 中心筛选——对应2.3中的算法2
FILTER <- function(x,p,q=0.5,t=1.5){
    C <- matrix(nrow = 1,ncol = dim(x)[2])
    R <- c()
    # 生成中心点集合C以及对应的半径集合R
    for(k in 2:p){
        # 进行Kmeans聚类
        km <- kmeans(x,k)
        for (i in 1:k) {
            # 返回本次聚类属于第i类的样本
            samples <- x[km$cluster==i,]
            # 计算样本到对应聚类中心的距离
            d <- dist(km$centers[i,],samples)
            # 得到当前类聚类中心对应的半径
            r <- as.numeric(quantile(d,q))
            R <- c(R,r)
            # 聚类中心
            C <- rbind(C,km$centers[i,])
        }
    }
    # 绘制图像
    #plot(C,pch=19,col=2,xlim = c(-30,30),ylim = c(-30,30),main = 'Sizeable radiuses(中心筛选前)')
    library(plotrix)
    for(i in 1:length(R)){
        draw.circle(C[i,1],C[i,2],R[i])
    }
    C <- C[-1,]
    COL <- c()
    for (i in 1:length(R)) {
        for (j in 1:length(R)) {
            # 计算Condition1
            d_i_j <- sqrt(sum((C[i,] - C[j,])^2))
            MIN_i_j <- min(R[i],R[j])
            MAX_i_j <- max(R[i],R[j])
            Condition1 <- (d_i_j + MIN_i_j < (t * MAX_i_j))
            if ((i!=j) && Condition1) {
                col <- c(i,j)[which.max(R[c(i,j)])]
                COL <- c(COL,col)
            }
        }
    }
    COL <- unique(COL)
    C <- C[-COL,]
    R <- R[-COL]
    # 绘制图像
    #plot(C,pch=19,col=2,xlim = c(-5,20),ylim = c(-5,20),main = 'Sizeable radiuses(中心筛选后)')
    library(plotrix)
    for(i in 1:length(R)){
        draw.circle(C[i,1],C[i,2],R[i])
    }
    return(list(C=C,R=R))
}

# 中心融合——对应2.3中的算法3
COMBINE <- function(C,R){
    L_C <- length(R)
    # 生成pointer矩阵
    pointer <- matrix(0,nrow = L_C,ncol = L_C)
    s <- c()
    c_s <- 1:L_C
    for(i in 1:L_C){
        for(j in 1:L_C){
            d_i_j <- sqrt(sum((C[i,1:2] - C[j,1:2])^2))
            Condition2 <- d_i_j<(R[i]+R[j])
            if((i!=j) & Condition2){
                pointer[i,j] <- 1
                pointer[j,i] <- 1
            }
        }
    }
    k <- 0
    while (length(c_s) != 0) {
        lst <- RECURSION(pointer,s,c_s,c_s[1])
        s <- lst$s
        c_s <- lst$c_s
        k <- k + 1
    }
    return(k)
}

# 递归函数，完成类的寻找工作——对应2.3中的算法4
RECURSION <- function(pointer,s,c_s,ini_co){
    corrdin <- which(pointer[ini_co,] == 1)
    s <- c(s,ini_co)
    c_s <- c_s[which(c_s!=ini_co)]
    if(length(corrdin) != 0){
        for(i in 1:length(corrdin)){
            if(corrdin[i] %in% s){
                next
            }else{
                ini_co <- corrdin[i]
                lst <- RECURSION(pointer,s,c_s,ini_co)
                s <- lst$s
                c_s <- lst$c_s
            }
        }
    }
    return(list(s=s,c_s=c_s))
}

# 定义distance函数，计算类内所有点到聚类中心距离的排序
dist <- function(centers,X){
    d <- sqrt(diag((X - centers) %*% t(X - centers)))
    return(d)
}


######################################################################################
# 生成模拟数据
#####################################################################################
# 导入用到的包
library(MASS)
# 多图绘制
par(mfrow=c(2,2))
# =============================================案例1：三圆环数据=========================================
set.seed(123)
n1 <- 500  # 内部圆的样本点数量
n2 <- 500  # 中间圆环的样本点数量
n3 <- 500  # 外部圆环的样本点数量
X1 <- matrix(nrow = n1,ncol = 2)
mu <- c(0,0)
sigma <- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
for(i in 1:n1){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if((x[1]^2 + x[2]^2) >= 0 & (x[1]^2 + x[2]^2) <= 1){
            X1[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
    
}

X2 <- matrix(nrow = n2,ncol = 2)
mu <- c(0,0)
sigma <- matrix(c(9,0,0,9),nrow = 2,ncol = 2)
for(i in 1:n2){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if((x[1]^2 + x[2]^2) >= 2 & (x[1]^2 + x[2]^2) <= 3){
            X2[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
}
X3 <- matrix(nrow = n3,ncol = 2)
mu <- c(0,0)
sigma <- matrix(c(36,0,0,36),nrow = 2,ncol = 2)
for(i in 1:n2){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if((x[1]^2 + x[2]^2) >= 5 & (x[1]^2 + x[2]^2) <= 6){
            X3[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
}
# 合并数据
X <- rbind(X1,X2,X3)
g <- c(rep(1,n1),rep(2,n2),rep(3,n3))
# 无标签绘图
plot(X,pch=20,col='blue',main='DataSet1',xlab = expression(X[1]),ylab=expression(X[2]))
# 有标签绘图
plot(X,col=g,pch=20,main='圆环数据',xlab = expression(X[1]),ylab=expression(X[2]))
# 保存数据
dataset1 <- cbind(X,g)
colnames(dataset1) <- c('X1','X2','g')
write.csv(dataset1, 'dataset1.csv', fileEncoding = 'utf-8', row.names = FALSE)
# ============================================案例2：块状数据===============================================
n <- 300
mu1 <- c(0,0)
mu2 <- c(0,-3)
mu3 <- c(3,-3)
mu4 <- c(3,3)
mu5 <- c(-3,3)
sigma <- matrix(c(0.1,0,0,0.1),nrow = 2,ncol = 2)
X1 <- mvrnorm(n,mu1,sigma)
X2 <- mvrnorm(n,mu2,sigma)
X3 <- mvrnorm(n,mu3,sigma)
X4 <- mvrnorm(n,mu4,sigma)
X5 <- mvrnorm(n,mu5,sigma)
X <- rbind(X1,X2,X3,X4,X5)
g <- rep(1:5,each=n)
# 无标签绘图
plot(X,pch=20,col='blue',main='DataSet2',xlab = expression(X[1]),ylab=expression(X[2]))
# 有标签绘图
plot(X,col=g,pch=20,main='块状数据',xlab = expression(X[1]),ylab=expression(X[2]))
# 保存数据
dataset2 <- cbind(X,g)
colnames(dataset2) <- c('X1','X2','g')
write.csv(dataset2, 'dataset2.csv', fileEncoding = 'utf-8', row.names = FALSE)
# ======================================案例3：块状数据+噪声==================================================
n <- 300
mu1 <- c(0,0)
mu2 <- c(0,-3)
mu3 <- c(3,-2)
mu4 <- c(2,2)
mu5 <- c(-1.5,3)
sigma <- matrix(c(0.1,0,0,0.1),nrow = 2,ncol = 2)
X1 <- mvrnorm(n,mu1,sigma)
X2 <- mvrnorm(n,mu2,sigma)
X3 <- mvrnorm(n,mu3,sigma)
X4 <- mvrnorm(n,mu4,sigma)
X5 <- mvrnorm(n,mu5,sigma)
X <- rbind(X1,X2,X3,X4,X5)
g <- rep(1:5,each=n)
n_noise <- 100
X_noise <- matrix(ncol=2,nrow = n_noise)
for(i in 1:n_noise){
    min_dist <- 0
    while(min_dist < 0.5){
        noise_x<- runif(1,-4,4)
        noise_y <- runif(1,-4,4)
        noise_one <- matrix(rep(c(noise_x, noise_y), times=dim(X)[1]), ncol = 2,byrow = TRUE)
        min_dist <- min(dist(noise_one, X))
    }
    X_noise[i,] <- noise_one[1,]
}
X <- rbind(X,X_noise)
g <- c(g,rep(6,n_noise))
# 无标签绘图
plot(X,pch=20,col='blue',main='DataSet3',xlab = expression(X[1]),ylab=expression(X[2]))
# 有标签绘图
plot(X,col=g,pch=20,main='带噪声的块状数据',xlab = expression(X[1]),ylab=expression(X[2]))
# 保存数据
dataset3 <- cbind(X,g)
colnames(dataset3) <- c('X1','X2','g')
write.csv(dataset3, 'dataset3.csv', fileEncoding = 'utf-8', row.names = FALSE)
# ==========================================案例4：圆弧数据+噪声==============================================
set.seed(123)
n1 <- 1000  # 上圆弧的样本点数量
n2 <- 1000  # 下圆弧的样本点数量
X1 <- matrix(nrow = n1,ncol = 2)
mu <- c(-0.5,0)
sigma <- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
for(i in 1:n1){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if(((x[1]+0.5)^2 + x[2]^2) >= 1 & ((x[1]+0.5)^2 + x[2]^2) <= 2 & x[2] > -0.5){
            X1[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
}
X2 <- matrix(nrow = n2,ncol = 2)
mu <- c(0.5,0)
sigma <- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
for(i in 1:n1){
    x <- mvrnorm(1,mu,sigma)
    while(TRUE){
        if(((x[1]-0.5)^2 + x[2]^2) >= 1 & ((x[1]-0.5)^2 + x[2]^2) <= 2 & x[2] < 0.5){
            X2[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu,sigma)
        } 
    }
}
# 合并数据
X <- rbind(X1,X2)
g <- c(rep(1,n1),rep(2,n2))
# 噪声数据
n_noise <- 100
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
g <- c(g,rep(3,n_noise))
# 无标签绘图
plot(X,pch=20,col='blue',main='DataSet4',xlab = expression(X[1]),ylab=expression(X[2]))
# 有标签绘图
plot(X,col=g,pch=20,main='带噪声的弧形数据',xlab = expression(X[1]),ylab=expression(X[2]))
# 保存数据
dataset4 <- cbind(X,g)
colnames(dataset4) <- c('X1','X2','g')
write.csv(dataset4, 'dataset4.csv', fileEncoding = 'utf-8', row.names = FALSE)


# ==========================================案例3：混合数据+噪声==============================================
set.seed(123)
n1 <- 300  # 左边块的样本点数量
n2 <- 200  # 右边块的样本点数量
n3 <- 1000  # 圆弧的样本点数量
# 生成左边的块
X1 <- matrix(nrow = n1,ncol = 2)
mu1 <- c(-2,2)
sigma1 <- matrix(c(0.2,0,0,0.2),nrow = 2,ncol = 2)
X1 <- mvrnorm(n1,mu1,sigma1)
# 生成右边的块
X2 <- matrix(nrow = n2,ncol = 2)
mu2 <- c(2,2)
sigma2 <- matrix(c(0.2,0,0,0.2),nrow = 2,ncol = 2)
X2 <- mvrnorm(n2,mu2,sigma2)
# 生成圆弧
X3 <- matrix(nrow = n3,ncol = 2)
mu3 <- c(0,0)
sigma3 <- matrix(c(3,0,0,3),nrow = 2,ncol = 2)
for(i in 1:n3){
    x <- mvrnorm(1,mu3,sigma3)
    while(TRUE){
        if(((x[1])^2 + x[2]^2) >= 1 & ((x[1])^2 + x[2]^2) <= 1.5 & x[2] <= 0){
            X3[i,] <- x
            break()
        }else{
            x <- mvrnorm(1,mu3,sigma3)
        } 
    }
}
# 添加噪声
n_noise <- 40
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



X <- rbind(X1,X2, X3,X_noise)
plot(X,pch=20,col='blue',main='DataSet4',xlab = expression(X[1]),ylab=expression(X[2]))




#################################################################################
# 利用各种方法进行聚类
#################################################################################
# 绘制聚类结果图
plot_clustering <- function(X, res, main){
    pch_set=c(1,2,3,5,7,8,9)
    plot(X[,1], X[,2],xlab = expression(X[1]),ylab=expression(X[2]),main=main, col='white')
    for(i in 1:length(res_grid_svm_knn)){
        if(res[i] == 0){
            points(x = X[i,1], y = X[i,2], col=res[i]+1, pch=15)
        }else{
            points(x = X[i,1], y = X[i,2], col=res[i]+1, pch=pch_set[res[i]])
        }
    }
}

# 加载一些包
library(fpc)
library(dbscan)
library(ellipse)
X <- dataset1[,1:2]; g <- dataset1[,3]
# X <- dataset2[,1:2]; g <- dataset2[,3]
# X <- dataset3[,1:2]; g <- dataset3[,3]
# X <- dataset4[,1:2]; g <- dataset4[,3]
# grid_SVM_knn
res_grid_svm_knn <- grid_svm_knn(X, gamma=10 ,cost=10, k=18, Index1=1213)  # dataset1
#res_grid_svm_knn <- grid_svm_knn(X, gamma=5 ,cost=0.36, k=18, Index1=1)
res_grid_svm_knn <- grid_svm_knn(X, gamma=20,cost=0.01, k=18, Index1=1213)
res_grid_svm_knn <- grid_svm_knn(X, gamma=10, cost=0.02265, k=25, Index1=1213)
res_grid_svm_knn <- grid_svm_knn(X, gamma=10, cost=0.023, k=18, Index1=1213)
# 绘制聚类结果
plot_clustering(X, res_grid_svm_knn, 'DataSet1')
# dbscan
kNNdistplot(X,k = 20)  # 找最优eps
res_db <- dbscan(X,eps=0.2,minPts=8)
res_db <- res_db$cluster
plot_clustering(X, res_db, 'DataSet4聚类效果图')
# DP-Kmeans
res_DP_Kmeans <- Radius.KMeans(X,1213,q=0.5,t=1.2)
res_DP_Kmeans <- res_DP_Kmeans$km$cluster
plot_clustering(X, res_DP_Kmeans,'DataSet4聚类效果图')
# 汇总上面的结果并保存
res1 <- cbind(res_grid_svm_knn, res_db, res_DP_Kmeans)
write.csv(res1, 'res1.csv', row.names = FALSE, fileEncoding = 'utf-8')

res2 <- cbind(res_grid_svm_knn, res_db, res_DP_Kmeans)
write.csv(res2, 'res2.csv', row.names = FALSE, fileEncoding = 'utf-8')

res3 <- cbind(res_grid_svm_knn, res_db, res_DP_Kmeans)
write.csv(res3, 'res3.csv', row.names = FALSE, fileEncoding = 'utf-8')

res4 <- cbind(res_grid_svm_knn, res_db, res_DP_Kmeans)
write.csv(res4, 'res4.csv', row.names = FALSE, fileEncoding = 'utf-8')
# birch用python处理


# 解除多图绘制
par(mfrow=c(1,1))
###############################################################################
# 聚类效果评估
###############################################################################
library(NMI)
NMI(res_db$cluster, g)





###############################################################################
# 对算法的时间复杂度进行分析
###############################################################################
library(MASS)
create_data <- function(n){
    n <- n/2
    mu1 <- c(0,2)
    mu2 <- c(0,-2)
    sigma <- matrix(c(0.1,0,0,0.1),nrow = 2,ncol = 2)
    X1 <- mvrnorm(n,mu1,sigma)
    X2 <- mvrnorm(n,mu2,sigma)
    X <- rbind(X1,X2)
    return(X)
}
X <- create_data(500)
# 无标签绘图
plot(X,pch=20,col='blue',xlab = expression(X[1]),ylab=expression(X[2]))
res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=1, k=5, Index1=1213)

# 随n增大的时间消耗情况
n_set <- (1:20)*50
time_set <- matrix(nrow = length(n_set), ncol = 1)
for(i in 1:length(n_set)){
    print(i)
    n <- n_set[i]
    X <- create_data(n)
    t1 <- proc.time()
    res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=1, k=5, Index1=1213)
    t2 <- proc.time()
    t <- t2 - t1
    time_set[i] <- t[3][[1]]
}
plot(x = n_set, y = time_set, type = 'b', main='the influence of n to running time', ylab = 'Time', xlab='n')

# 随k增大的时间消耗情况
time_set <- matrix(nrow = length(k_set), ncol = 10)
for(j in 1:10){
    k_set <- 5:30
    for(i in 1:length(k_set)){
        print(paste0('正在处理 行',i,' 列',j))
        k <- k_set[i]
        X <- create_data(1000)
        t1 <- proc.time()
        res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=1, k=k, Index1=1213)
        t2 <- proc.time()
        t <- t2 - t1
        time_set[i,j] <- t[3][[1]]
    }
}
time_set_mean <- apply(time_set, 1, mean)
plot(x = k_set, y = time_set_mean, type = 'b', main='the influence of k to running time', ylab = 'Time', xlab='k')


k_set <- 1:30
time_set1 <- matrix(nrow = length(k_set), ncol = 1)
for(i in 1:length(k_set)){
    print(paste0('正在处理 行',i))
    k <- k_set[i]
    set.seed(1213)
    X <- create_data(2000)
    t1 <- proc.time()
    res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=1, k=k, Index1=1213)
    t2 <- proc.time()
    t <- t2 - t1
    time_set1[i] <- t[3][[1]]
}
plot(x = k_set, y = time_set1, type = 'b', main='the influence of k to running time', ylab = 'Time', xlab='k')


# 随gamma增大的时间消耗情况
gamma_set <- 10**(-6:6)
time_set <- matrix(nrow = length(gamma_set), ncol = 1)
for(i in 1:length(gamma_set )){
    print(i)
    gamma <- gamma_set[i]
    X <- create_data(1000)
    t1 <- proc.time()
    res_grid_svm_knn <- grid_svm_knn(X, gamma=gamma,cost=1, k=20, Index1=1213)
    t2 <- proc.time()
    t <- t2 - t1
    time_set[i] <- t[3][[1]]
}
plot(x = gamma_set, y = time_set, type = 'b', main='gamma对算法运行时间的影响', ylab = '消耗时间', xlab='gamma')


# 随cost增大的时间消耗情况
cost_set <- 10**(-4:4)
time_set <- matrix(nrow = length(cost_set), ncol = 1)
for(i in 1:length(cost_set )){
    print(i)
    cost <- cost_set[i]
    X <- create_data(3000)
    t1 <- proc.time()
    res_grid_svm_knn <- grid_svm_knn(X, gamma=1,cost=cost, k=20, Index1=1213)
    t2 <- proc.time()
    t <- t2 - t1
    time_set[i] <- t[3][[1]]
}
plot(x = cost_set, y = time_set, type = 'b', main='cost对算法运行时间的影响', ylab = '消耗时间', xlab='cost')


















