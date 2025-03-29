## (1)
### 代码
```
head(iris)
ncol(iris)
c(class(iris[[1]]),class(iris[[2]]),class(iris[[3]]),class(iris[[4]]),class(iris[[5]]))
```
### 结果
```
> head(iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
> ncol(iris)
[1] 5
> c(class(iris[[1]]),class(iris[[2]]),class(iris[[3]]),class(iris[[4]]),class(iris[[5]]))
[1] "numeric" "numeric" "numeric" "numeric" "factor"
```
### 答案
有五列  
从第一列到第五列分别为：数值型、数值型、数值型、数值型、因子型
## (2)
### 代码
```
Mean <- aggregate(Sepal.Length ~ Species, iris, mean)
SD <- aggregate(Sepal.Length ~ Species, iris, sd)
colnames(Mean) <- c("Species", "mean_value")
colnames(SD) <- c("Species", "sd")
write.csv(data.frame(Mean,SD$sd),"~/Desktop/大二下/生物信息学/practice 2/MSD.csv",quote = F)
```
### 结果
```
,Species,mean_value,SD.sd
1,setosa,5.006,0.352489687213451
2,versicolor,5.936,0.516171147063863
3,virginica,6.588,0.635879593274432
```
## (3)
### 代码
`summary(aov(Sepal.Width ~ Species,data=iris))`
### 结果
```
             Df Sum Sq Mean Sq F value Pr(>F)    
Species       2  11.35   5.672   49.16 <2e-16 ***
Residuals   147  16.96   0.115                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

