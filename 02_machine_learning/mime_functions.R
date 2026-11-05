source("config.R")  ## 统一路径配置，见 config.R`r`n# ============================
# 加载所需的R包
# ============================
library(tgp)             # 贝叶斯树高斯过程
library(openxlsx)        # Excel文件操作
library(seqinr)          # 序列分析
library(plyr)            # 数据处理工具
library(randomForestSRC) # 随机森林
library(glmnet)          # 广义线性模型和正则化
library(plsRglm)         # 偏最小二乘回归
library(gbm)             # 梯度提升机
library(caret)           # 机器学习算法训练和测试
library(mboost)          # 增强模型
library(e1071)           # SVM等算法
library(BART)            # 贝叶斯加法回归树
library(MASS)            # 广泛应用的统计方法
library(snowfall)        # 并行计算支持
library(xgboost)         # 极端梯度提升模型
library(ComplexHeatmap)  # 复杂热图绘制
library(RColorBrewer)    # 颜色选择
library(pROC)            # ROC曲线工具
library(circlize)        # 圆形可视化

# ============================
# 设置工作目录和参数
# ============================
min.selected.var = 1     # 设置变量选择的最小数目
max.selected.var = 20    # 设置变量选择的最大数目
top.models.for.roc = 72  # 设置绘制ROC曲线的模型数量

# 定义一个函数RunML，用于运行机器学习算法
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  # 清理和准备算法名称和参数
  method = gsub(" ", "", method) # 去除方法名称中的空格
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # 从方法名称中提取算法名称
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # 从方法名称中提取参
  
  # 根据提取的算法名称，准备相应的参数
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL  # 如果没有匹配到任何名称，返回NULL
  )
  
  # 输出正在运行的算法和使用的变量数
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")
  
  # 将传入的参数和提取的参数组合成一个新的参数列表
  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)
  
  # 使用do.call动态调用相应的算法实现函数
  obj <- do.call(what = paste0("Run", method_name),
                 args = args) 
  
  # 根据模式，输出不同的信息
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

# 定义用于运行Elastic Net正则化线性模型的函数
RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  # 使用交叉验证找到最优的正则化参数
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10)
  # 建立最终模型
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行Lasso正则化线性模型的函数
RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

# 定义用于运行Ridge正则化线性模型的函数
RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

# 定义用于运行逐步广义线性模型的函数
RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  
  if(ncol(Train_set)>50){
    vars <- apply(Train_set,2,var)
    Train_set <- Train_set[,order(vars,decreasing=TRUE)[1:50],drop=FALSE]
  }
  
  fit <- step(
    glm(
      Train_label[[classVar]] ~ .,
      family="binomial",
      data=as.data.frame(Train_set)
    ),
    direction=direction,
    trace=0
  )
  
  fit$subFeature <- names(coef(fit))
  fit$subFeature <- setdiff(fit$subFeature,
                            c("(Intercept)","Intercept"))
  
  if(mode=="Model") return(fit)
  if(mode=="Variable") return(ExtractVar(fit))
}


# 定义用于运行支持向量机的函数
RunSVM <- function(Train_set, Train_label, mode, classVar){
  # 将数据框转换为因子类型，适合SVM模型输入
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  # 建立SVM模型
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行线性判别分析的函数
RunLDA <- function(Train_set, Train_label, mode, classVar){
  # 准备数据，将类变量转换为因子类型
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  # 使用train函数建立LDA模型
  fit = train(eval(parse(text = paste(classVar, "~."))), 
              data = data, 
              method="lda",
              trControl = trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# ============================
# Naive Bayes
# ============================
RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  
  library(e1071)
  
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- naiveBayes(
    formula = as.formula(paste(classVar, "~ .")),
    data = data
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if(mode == "Model"){
    return(fit)
  }
  
  if(mode == "Variable"){
    return(ExtractVar(fit))
  }
}

# 定义用于运行梯度提升机的函数
RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  # 准备数据，将类变量和训练集绑定
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  # 建立GLMBoost模型
  # 使用较大的迭代次数以获得最佳性能
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = 100))
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行偏最小二乘回归和广义线性模型的函数
RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  # 使用交叉验证评估模型参数
  cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ ., 
                              data = as.data.frame(Train_set),
                              nt=10, verbose = FALSE)
  # 建立PLSRGLM模型
  fit <- plsRglm(Train_label[[classVar]], 
                 as.data.frame(Train_set), 
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行随机森林的函数
RunRF <- function(Train_set, Train_label, mode, classVar){
  # 设置随机森林参数，如树的最小节点大小
  rf_nodesize = 2 # 可根据需要调整
  # 准备数据，将类变量转换为因子
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  # 建立随机森林模型
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 2000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 定义用于运行梯度提升机的函数
RunGBM <- function(Train_set, Train_label, mode, classVar){
  # 获取样本数量
  n_samples <- nrow(Train_set)
  
  # 根据样本数量调整参数
  if (n_samples < 50) {
    # 样本数较少时，使用更小的参数
    n_minobsinnode <- 2
    cv_folds <- 3
    n_trees <- 100
    interaction_depth <- 2
  } else if (n_samples < 100) {
    # 样本数中等时
    n_minobsinnode <- 5
    cv_folds <- 5
    n_trees <- 500
    interaction_depth <- 3
  } else {
    # 样本数足够时，使用原来的参数
    n_minobsinnode <- 10
    cv_folds <- 10
    n_trees <- 1000
    interaction_depth <- 3
  }
  
  # 建立初步的GBM模型
  tryCatch({
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = n_trees,
               interaction.depth = interaction_depth,
               n.minobsinnode = n_minobsinnode,
               shrinkage = 0.001,
               cv.folds = cv_folds,
               n.cores = 1)  # 改为单核避免并行问题
    
    # 选择最优的迭代次数
    best <- which.min(fit$cv.error)
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = best,
               interaction.depth = interaction_depth,
               n.minobsinnode = n_minobsinnode,
               shrinkage = 0.001,
               n.cores = 1)
    
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }, error = function(e) {
    # 如果GBM失败，返回空结果
    warning(sprintf("GBM模型构建失败: %s", e$message))
    return(if(mode == "Model") NULL else c())
  })
}

# 定义用于运行XGBoost的函数
RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  # 创建交叉验证折叠
  indexes = createFolds(Train_label[[classVar]], k = 5, list=T)
  # 计算每折的最优模型参数
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(data = Train_set[-pt, ], 
                         label = Train_label[-pt, classVar])
    dtest = xgb.DMatrix(data = Train_set[pt, ], 
                        label = Train_label[pt, classVar])
    
    bst <- xgb.train(data=dtrain, 
                     params = list(max_depth=2, eta=1, objective="binary:logistic"),
                     nrounds=10, 
                     watchlist=list(train=dtrain, eval=dtest), 
                     verbose = 0)
    
    # 3.2.1 兼容: evaluation_log 可能为空
    if (!is.null(bst$evaluation_log) && nrow(bst$evaluation_log) > 0) {
      cols <- colnames(bst$evaluation_log)
      if ("eval_logloss" %in% cols) {
        which.min(bst$evaluation_log$eval_logloss)
      } else if ("test_logloss" %in% cols) {
        which.min(bst$evaluation_log$test_logloss)
      } else {
        5
      }
    } else {
      5
    }
  }))
  
  # 使用最常见的轮数建立最终模型
  nround <- as.numeric(names(which.max(table(CV))))
  
  # 3.2.1 修复: xgboost() 用 x/y 参数, y 必须是 factor
  # 创建的模型类是 "xgboost" (不是 "xgb.Booster")
  raw_fit <- xgboost(x = Train_set, 
                     y = as.factor(Train_label[[classVar]]), 
                     max.depth = 2, eta = 1, nrounds = nround, 
                     objective = "binary:logistic")
  
  # 3.2.1 关键修复: 用 list 包装模型, class 不要包含 "xgb.Booster"
  fit <- list(model = raw_fit, subFeature = colnames(Train_set))
  class(fit) <- "xgb.wrapper"
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# DRF模型的注释已被删除，因为此部分代码被注释掉了
# RunDRF <- function(Train_set, Train_label, mode, classVar){
#   Train_label <- data.frame(
#     "0" = as.numeric(Train_label == 0),
#     "1" = as.numeric(Train_label == 1)
#   )
#   fit <- drf(X = Train_set, 
#              Y = Train_label, 
#              compute.variable.importance = F)
#   fit$subFeature = colnames(Train_set)
#   
#   summary(predict(fit, functional = "mean", as.matrix(Train_set))$mean)
#   
#   if (mode == "Model") return(fit)
#   if (mode == "Variable") return(ExtractVar(fit))
# }

# 定义一个函数用于在执行过程中抑制输出
quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())  # 将输出重定向到临时文件
    on.exit(sink())  # 确保在函数退出时恢复正常输出
  }
  # 根据参数决定是否抑制消息
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

# 定义一个函数用于标准化数据
standarize.fun <- function(indata, centerFlag, scaleFlag) {  
  scale(indata, center=centerFlag, scale=scaleFlag)
}

# 定义一个函数用于批量处理数据的标准化
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)  # 保存原始样本名称
  # 如果没有指定队列，将所有数据视为一个队列
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)  # 根据队列分割数据
  }
  
  # 如果没有提供中心化标志，默认不进行中心化
  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  # 如果中心化标志是单一值，应用于所有数据
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  # 如果中心化标志没有命名，按顺序匹配
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }
  
  # 如果没有提供缩放标志，默认不进行缩放
  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  # 如果缩放标志是单一值，应用于所有数据
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  # 如果缩放标志没有命名，按顺序匹配
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  # 使用mapply函数对每个数据队列应用标准化函数
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  # 将处理后的数据按原始顺序重新组合
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

# 定义一个函数用于从模型中提取重要的变量
ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
    "glm" = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train" = fit$coefnames,
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0],
    "rfsrc" = {
      imp <- fit$importance[,1]
      
      imp <- imp[!is.na(imp)]
      
      if(length(imp)==0){
        character(0)
      }else{
        names(sort(imp,decreasing=TRUE)[1:min(50,length(imp))])
      }
    }, 
    "gbm" = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0],
    "xgb.Booster" = fit$subFeature,
    "xgboost" = fit$subFeature,        # <-- 3.2.1 新增: xgboost() 创建的模型
    "xgb.wrapper" = fit$subFeature,    # <-- 3.2.1 新增: wrapper 对象
    "naiveBayes" = fit$subFeature
  ))
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  Feature <- Feature[!is.na(Feature)]
  return(Feature)
}

# 定义一个函数用于计算预测得分
# ========== 3. 替换 CalPredictScore 函数 ==========
CalPredictScore <- function(fit, new_data, type = "lp"){
  # 3.2.1 兼容: 如果 fit 是 wrapper list, 取 model
  if (
    is.list(fit) &&
    !inherits(fit, c("glm", "lm")) &&
    "model" %in% names(fit)
  ) {
    
    actual_fit <- fit$model
    features <- fit$subFeature
    
  } else {
    
    actual_fit <- fit
    features <- fit$subFeature
    
  }
  
  new_data <- new_data[, features, drop = FALSE]
  original_samples <- rownames(new_data)
  
  if(class(actual_fit)[1] == "svm.formula"){
    na_rows <- apply(new_data, 1, function(x) any(is.na(x)))
    if(any(na_rows)){
      warning(sprintf("SVM预测：发现 %d 个含有NA的样本，将返回NA预测值", sum(na_rows)))
      valid_data <- new_data[!na_rows, , drop = FALSE]
    } else {
      valid_data <- new_data
      na_rows <- rep(FALSE, nrow(new_data))
    }
  } else {
    valid_data <- new_data
    na_rows <- rep(FALSE, nrow(new_data))
  }
  
  # 获取实际模型类
  fit_class <- class(actual_fit)[1]
  
  RS_valid <- quiet(switch(
    EXPR = fit_class,
    "lognet"      = predict(actual_fit, type = 'response', as.matrix(valid_data)),
    "glm"         = predict(actual_fit, type = 'response', as.data.frame(valid_data)),
    "svm.formula" = {
      pred <- predict(actual_fit, as.data.frame(valid_data), probability = TRUE)
      probs <- attr(pred, "probabilities")
      if("1" %in% colnames(probs)){
        probs[, "1"]
      } else if("0" %in% colnames(probs)){
        1 - probs[, "0"]
      } else {
        stop("SVM概率预测未包含预期的类别")
      }
    },
    "train"       = predict(actual_fit, valid_data, type = "prob")[[2]],
    "glmboost"    = predict(actual_fit, type = "response", as.data.frame(valid_data)),
    "plsRglmmodel" = {
      
      X <- as.matrix(
        valid_data[, actual_fit$subFeature, drop = FALSE]
      )
      
      score <- X %*% actual_fit$ww
      
      score <- as.data.frame(score)
      
      # 自动匹配FinalModel需要的变量名
      required_vars <- all.vars(
        formula(actual_fit$FinalModel)
      )
      
      required_vars <- setdiff(
        required_vars,
        "y"
      )
      
      colnames(score) <- required_vars
      
      predict(
        actual_fit$FinalModel,
        newdata = score,
        type = "response"
      )
    },
    "rfsrc"        = predict(actual_fit, as.data.frame(valid_data))$predicted[, "1"],
    "gbm"          = predict(actual_fit, type = 'response', as.data.frame(valid_data)),
    "xgboost"     = predict(actual_fit, newdata = valid_data, type = "response"),  # <-- 3.2.1 修复
    "xgb.Booster" = predict(actual_fit, as.matrix(valid_data)),
    "naiveBayes" = predict(object = actual_fit, type = "raw", newdata = valid_data)[, "1"]
  ))
  
  RS_valid = as.numeric(as.vector(RS_valid))
  RS = rep(NA, nrow(new_data))
  RS[!na_rows] = RS_valid
  names(RS) = original_samples
  return(RS)
}

# 定义一个函数用于预测类别
# ========== 4. 替换 PredictClass 函数 ==========
PredictClass <- function(fit, new_data){
  # 3.2.1 兼容: 如果 fit 是 wrapper list, 取 model
  # ===============================
  # 修复 glm/lm 自带 model 字段问题
  # ===============================
  
  if (
    is.list(fit) &&
    !inherits(fit, c("glm","lm","train")) &&
    "model" %in% names(fit)
  ) {
    
    actual_fit <- fit$model
    features <- fit$subFeature
    
  } else {
    
    actual_fit <- fit
    features <- fit$subFeature
    
  }
  
  new_data <- new_data[, features, drop = FALSE]
  original_samples <- rownames(new_data)
  
  if(class(actual_fit)[1] == "svm.formula"){
    na_rows <- apply(new_data, 1, function(x) any(is.na(x)))
    if(any(na_rows)){
      warning(sprintf("SVM分类预测：发现 %d 个含有NA的样本，将返回NA预测值", sum(na_rows)))
      valid_data <- new_data[!na_rows, , drop = FALSE]
    } else {
      valid_data <- new_data
      na_rows <- rep(FALSE, nrow(new_data))
    }
  } else {
    valid_data <- new_data
    na_rows <- rep(FALSE, nrow(new_data))
  }
  
  fit_class <- class(actual_fit)[1]
  
  label_valid <- quiet(switch(
    EXPR = fit_class,
    "lognet"      = predict(actual_fit, type = 'class', as.matrix(valid_data)),
    "glm"         = ifelse(test = predict(actual_fit, type = 'response', as.data.frame(valid_data))>0.5,
                           yes = "1", no = "0"),
    "svm.formula" = as.character(predict(actual_fit, as.data.frame(valid_data))),
    "train"       = predict(actual_fit, valid_data, type = "raw"),
    "glmboost"    = predict(actual_fit, type = "class", as.data.frame(valid_data)),
    "plsRglmmodel" = {
      
      X <- as.matrix(
        valid_data[, actual_fit$subFeature, drop = FALSE]
      )
      
      score <- X %*% actual_fit$ww
      
      score <- as.data.frame(score)
      
      required_vars <- all.vars(
        formula(actual_fit$FinalModel)
      )
      
      required_vars <- setdiff(
        required_vars,
        "y"
      )
      
      colnames(score) <- required_vars
      
      pred <- predict(
        actual_fit$FinalModel,
        newdata = score,
        type = "response"
      )
      
      ifelse(
        pred > 0.5,
        "1",
        "0"
      )
    },
    "rfsrc"        = predict(actual_fit, as.data.frame(valid_data))$class,
    "gbm"          = ifelse(test = predict(actual_fit, type = 'response', as.data.frame(valid_data))>0.5,
                            yes = "1", no = "0"),
    "xgboost"     = as.character(predict(actual_fit, newdata = valid_data, type = "class")),  # <-- 3.2.1 修复
    "xgb.Booster" = ifelse(test = predict(actual_fit, as.matrix(valid_data))>0.5,
                           yes = "1", no = "0"),
    "naiveBayes" = predict(object = actual_fit, type = "class", newdata = valid_data)
  ))
  
  label_valid = as.character(as.vector(label_valid))
  label = rep(NA_character_, nrow(new_data))
  label[!na_rows] = label_valid
  names(label) = original_samples
  return(label)
}

# 定义一个函数用于评估模型性能
RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){
  
  # 检查测试标签中是否存在队列指标
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }
  
  # 如果提供了训练集和训练标签，将它们与测试集合并
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])
    
    # 如果提供了训练名称，将其作为队列名称
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Train-GEO"
    }
    # 更新训练标签的列名，包括队列变量和类变量
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }
  
  # 计算预测得分
  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  
  # 准备输出数据，包括预测得分
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  # 按队列分组
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  
  # 计算每个队列的AUC值
  result <- unlist(lapply(Predict.out, function(data){
    # 检查该队列是否有两个类别
    unique_classes <- unique(data[[classVar]])
    if(length(unique_classes) < 2){
      warning(sprintf("队列 '%s' 只有一个类别 (%s)，无法计算AUC，返回NA",
                      unique(data[[cohortVar]]), paste(unique_classes, collapse=",")))
      return(NA)
    }
    # 检查是否有NA值
    valid_idx <- !is.na(data$RS) & !is.na(data[[classVar]])
    if(sum(valid_idx) < 2){
      warning(sprintf("队列 '%s' 有效样本数不足，返回NA",
                      unique(data[[cohortVar]])))
      return(NA)
    }
    # 计算AUC
    tryCatch({
      as.numeric(auc(suppressMessages(roc(data[[classVar]][valid_idx], data$RS[valid_idx]))))
    }, error = function(e){
      warning(sprintf("队列 '%s' 计算AUC时出错: %s",
                      unique(data[[cohortVar]]), e$message))
      return(NA)
    })
  }))
  
  return(result)
}

# 定义一个简单的热图绘制函数
SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows,
                          gene_counts = NULL){  # 添加基因数量参数
  Cindex_mat <- as.matrix(Cindex_mat)
  # 定义列注释
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  # 定义行注释，包括平均C指数的条形图
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  # 如果提供了基因数量信息，修改行名以包含基因数量
  row_labels <- rownames(Cindex_mat)
  if (!is.null(gene_counts)) {
    row_labels <- paste0(rownames(Cindex_mat), " (n=", gene_counts[rownames(Cindex_mat)], ")")
  }
  
  # 根据列数调整热图宽度
  n_cols <- ncol(Cindex_mat)
  if (n_cols > 50) {
    # 当列数太多时，使用较小的单位宽度
    cellwidth_adj <- 0.3
  } else {
    cellwidth_adj <- cellwidth
  }
  
  # 计算合理的热图尺寸
  heatmap_width <- min(n_cols * cellwidth_adj + 2, 50)  # 限制最大宽度为50厘米
  heatmap_height <- max(nrow(Cindex_mat) * cellheight, 10)  # 最小高度为10厘米
  
  # 创建热图（简化版本，移除 column_split 以避免问题）
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha,
          top_annotation = col_ha,
          col = circlize::colorRamp2(
            c(min(Cindex_mat, na.rm = TRUE),
              mean(Cindex_mat, na.rm = TRUE),
              max(Cindex_mat, na.rm = TRUE)),
            c("#4195C1",
              "#FFFFFF",
              "#FFBC90")
          ), # 定义颜色，从蓝到红通过白色渐变
          rect_gp = gpar(col = "black", lwd = 0.5), # 设置边框颜色为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 是否对行和列进行聚类
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_labels = row_labels,  # 使用修改后的行标签
          width = unit(heatmap_width, "cm"),
          height = unit(heatmap_height, "cm"),
          cell_fun = function(j, i, x, y, w, h, col) { # 在每个格子中添加文本
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 7))
          }
  )
}

