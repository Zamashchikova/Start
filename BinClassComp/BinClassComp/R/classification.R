#' Get data
#'
#' Get data of Sql Server
#'
#' @param name_object is name of oject
#'
#' @param query is query text
#'
#' @param connectionStringDataMine is a connection string
#'
#' @return list of data and y_name
#'
#' @export
bmi3<-function(x){
	b<-x^2
	return(b)
}



library(RODBC)
library(caret)
library(InformationValue)
library(glmnet)
library(pROC)
library(randomForest)
library(car)
library(RColorBrewer)
library(xgboost)
library(ROSE)

scriptPath <- function() {
	getSrcDirectory(scriptPath);
}

#' Get data
#'
#' Get data of Sql Server
#' @param y_name is name of outcome
#' @param query is query text
#' @param connectionStringDataMine is a connection string
#' @return list of data and y_name
#' @export
get_data<-function(y_name,query,connection_string)
{
	dbhandle <- odbcDriverConnect(connection_string)
	data<-sqlQuery(dbhandle,query)
	on.exit(odbcClose(dbhandle))
	return(list(y_name=y_name,data=data))
}


#' Z-test
#'
#' #' Z-test for binary features
#' @param x1 is first vecot
#' @param x2 is second vector
#' @return value of statisic
#' @export
z_test<-function(x1,x2)
{
	n1<-length(x1)
	n2<-length(x2)
	m1<-sum(x1)
	m2<-sum(x2)
	if((m1+m2)!=0 & (m1+m2)!=(n1+n2) & n1!=0 & n2!=0)
	{
		zd<-((m1+m2)/(n1+n2))*((n1+n2-m1-m2)/(n1+n2))*(1/n1+1/n2)
		z<-(m1/n1+1/(2*n1)-m2/n2-1/(2*n2))/sqrt(zd)
		p_value<-pnorm(z,mean = 0, sd = 1,lower.tail=FALSE)
	}
	else
	{
		z<-NA
		p_value<-1
	}

	return(p_value)
}

#' Wilcoxon test
#'
#' wilcoxon test for numeric vectors
#' @param x1 is first vecot
#' @param x2 is second vector
#' @return value of statisic
#' @export
w_test<-function(x1,x2)
{
	w<-wilcox.test(x1,x2)
	p_value<-w$p.value
	return(p_value)
}

#' Bartletta test
#'
#' Bartletta test for multicollinearity
#' @param x is matrix of features
#' @param alpha is significant level
#' @return value of statisic
#' @export
testBartlett<-function(x,alpha)
{
	if(ncol(x)<2)
	{
		return(FALSE)
	}
	m<-ncol(x)
	N<-nrow(x)
	R<-cor(x)
	detR<-(det(R))
	if( !is.na(detR) ){
		B<- -(N-1-(1/6)*(2*m+5))*log(abs(det(R)))
		pValue <- pchisq(B, df=m*(m-1)/2, lower.tail = FALSE)
		result <- ifelse(pValue > alpha, TRUE, FALSE)
		return(result)
	}
	return(FALSE)
}

#' Features selection
#'
#' Features selection for classification
#' @param env is data environment
#' @param alpha is significant level
#' @return significant columns
#' @export
significant_features<-function(env,alpha = 0.1)
{
	sig_colums<-c(env$y_number)
	data_study_0<-env$data_study[which(env$data_study[,env$y_number] == 0),]
	data_study_1<-env$data_study[which(env$data_study[,env$y_number] == 1),]
	for(i in (env$y_number + 1):ncol(env$data_study))
	{
		if(typeof(env$data_study[,i]) == "double")
		{
			if(length(unique(env$data_study[,i])) == 2)
			{
				p_value<-z_test(data_study_1[,i], data_study_0[,i])
			}
			else
			{
				#means
				p_value<-w_test(data_study_1[,i], data_study_0[,i])
				#variance
				if (p_value >= alpha/2 & p_value <= (1-alpha/2))
				{
					p_value<-mood.test(data_study_1[,i], data_study_0[,i])$p.value
				}
			}
			if(!is.nan(p_value)) {
				if (p_value <= alpha/2 | p_value >= (1-alpha/2)){
					sig_colums<-c(sig_colums,i)
				}
			}
		}
		else{
			values = factor(env$data_study[,i])
			pValue = chisq.test(values, env$data_study[,env$y_number])$p.value
			if(!is.nan(pValue)) {
				if (pValue <= alpha/2){
					sig_colums<-c(sig_colums,i)
				}
			}
		}
	}
	env$data_study<-env$data_study[,sig_colums]
	env$data_test<-env$data_test[,sig_colums]
	return(sig_colums)
}

#' Multicollinearity deletion
#'
#' Multicollinearity deletion
#' @param env is data environment
#' @param method is choice of vif or bartlett
#' @param alpha is significant level
#' @param vif_bound is bound for vif
#' @return remaining features
#' @export
multi_features<-function(env, method = 'bartlett', alpha = 0.05, vif_bound = 4)
{
	if(ncol(env$data_study)<=2)
	{
		return(seq(1:ncol(env$data_study)))
	}

	sig_colums<-c(env$y_number)
	is_numeric_columns<-apply(env$data_study,2,is.numeric)
	sig_colums<-c(sig_colums,which( colnames(env$data_study)== names(is_numeric_columns[which(is_numeric_columns == FALSE)])))

	if(method == 'vif')
	{
		ev<-names(env$data_study[,-env$y_number])
		vif_list<-NULL
		while((is.null(vif_list)|length(vif_list[vif_list>=vif_bound])>0) & length(ev)>1)
		{
			vif_list<-sapply(ev,function(e){
				r<-summary(lm(as.formula(paste(e,"~",paste(ev[ev!=e],collapse="+"))),data=env$data_study, family = binomial()))$r.squared
				return(ifelse(r==1,100000,sqrt(1/(1-r))))
			})
			ev<-ev[-(which.max(vif_list))]
		}
		sig_colums<-c(sig_colums, which(colnames(env$data_study) %in% ev))
	}
	out_colums<-c()
	if(method == 'bartlett')
	{
		while(ncol(env$data_study[,-c(sig_colums,out_colums),drop=FALSE]) >= 2 & testBartlett(env$data_study[,-c(sig_colums,out_colums),drop=FALSE],alpha = alpha)==FALSE  )
		{
			R<-cor(env$data_study[,-c(sig_colums,out_colums)])
			Robr<-solve(R,tol = 1e-40)
			Rm<-c()
			for( j in 1:ncol(env$data_study[,-c(sig_colums,out_colums)]))
			{
				Rm<-c(Rm,sqrt(1-1/Robr[j,j]))
			}
			out_colums<-c(out_colums,which(colnames(env$data_study) == names(env$data_study[,-c(sig_colums,out_colums)])[(which.max(Rm))]))
		}
		sig_colums<-seq(1,ncol(env$data_study))
		sig_colums<-c(sig_colums[!sig_colums %in% out_colums])
	}
	env$data_study<-env$data_study[,sig_colums]
	env$data_test<-env$data_test[,sig_colums]
	return(sig_colums)
}

#' Box-Cox Transformation
#'
#' Box-Cox Transformation
#' @param feature is feature for transformation
#' @param lamda is conversion parameter
#' @return transform feature
#' @export
box_cox<-function(feature,lamda)
{
	if(lamda!=0)
	{
		feature<-((feature)^lamda-1)/lamda
	}
	else
	{
		feature<-log(feature)
	}
	return(feature)
}

#' Features selection by AUC
#'
#' Features selection by AUC
#' @param env is data environment
#' @param alg is algoritm of forward/backwards selection
#' @return significant columns
#' @export
get_feature_step<-function(env,alg='b')
{
	model='glm'
	feature_number<-seq(env$y_number+1,ncol(env$data_study))
	if(alg=="f")
	{
		step<-TRUE
		auc_best<-0
		auc_best_number<-c()
		while(step)
		{
			auc<-sapply(feature_number,function(t){
				model_res<-get_model(env$data_study[,c(env$y_number,auc_best_number,t)],model)
				predict<-get_prediction(model_res,env$data_study[,c(auc_best_number,t),drop=FALSE],model)
				auc<-roc.curve(env$data_study[,env$y_number],predict,plotit = F)$auc
				return(auc)
			})
			auc_max<-max(auc)
			if (auc_best<=auc_max)
			{
				auc_best<-auc_max
				auc_best_number<-c(auc_best_number,feature_number[which.max(auc)])
				feature_number<-feature_number[-which.max(auc)]
			}
			else
			{
				step<-FALSE
			}
		}
		feature_number<-auc_best_number
	}
	if(alg=="b")
	{
		step<-TRUE
		model_res<-get_model(env$data_study,model)
		predict<-get_prediction(model_res,env$data_study,model)
		auc_all<-roc.curve(env$data_study[,env$y_number],predict,plotit = F)$auc
		feature_del<-c()
		while(step & length(feature_number)>1)
		{
			auc<-sapply(feature_number,function(t)
			{
				model_res<-get_model(env$data_study[,-c(t,feature_del)],model)
				predict<-get_prediction(model_res,env$data_study[,-c(t,feature_del)],model)
				auc<-roc.curve(env$data_study[,env$y_number],predict,plotit = F)$auc
				return(auc)
			})
			auc_max<-max(auc)
			if (auc_all<=auc_max)
			{
				feature_del<-c(feature_del,feature_number[which.max(auc)])
				feature_number<-feature_number[-which.max(auc)]
				auc_all<-auc_max
			}
			else
			{
				step<-FALSE
			}
		}
	}
	env$data_study<-env$data_study[,c(env$y_number,feature_number)]
	env$data_test<-env$data_test[,c(env$y_number,feature_number)]
}

#' Prepare data by box_cox
#'
#' Prepare data by box_cox transformation
#' @param env is data environment
#' @param lamda range for lamda
#' @return data frame with information about lamda
#' @export
prepare_box_cox<-function(env,lamda=seq(-4,4,0.01))
{
	model<-'glm'
	normalization(env)
	bc_info<-NULL
	tmp_study<-env$data_study
	for(i in (env$y_number+1):ncol(env$data_study))
	{
		auc<-sapply(lamda,function(t){
			tmp_study[,i]<-box_cox(env$data_study[,i],t)
			model_res<-get_model(tmp_study[,c(env$y_number,i)],model)
			predict<-get_prediction(model_res,tmp_study[,i,drop=FALSE],model)
			auc<- roc.curve(tmp_study[,env$y_number],predict,plotit = F)$auc
			return(auc)
		})
		lamda_best<-lamda[which.max(auc)]

		bc_info<-rbind(bc_info,c(names(env$data_study)[i],lamda_best))
		env$data_study[,i]<-box_cox(env$data_study[,i],lamda_best)
		if(is.null(env$data_test)==FALSE)
		{
			env$data_test[,i] <-box_cox(env$data_test[,i],lamda_best)
		}
	}
	bc_info<-as.data.frame(bc_info)
	names(bc_info)<-c("names","lamda")
	bc_info[,1]<-as.character(paste(bc_info[,1]))
	bc_info[,2]<-as.numeric(paste(bc_info[,2]))
	return(bc_info=bc_info)
}

#' Replace anomalies
#'
#' Replace anomalies
#' @param env is data environment
#' @param k is coefficient for box plot
#' @param max_part_outlier is max part outlier for replacement
#' @param step step for change k
#' @return data frame with change information
#' @export
replace_anomalies<-function(env, k=1.5, max_part_outlier=0.2, step=0.5)
{
	i<-env$y_number+1
	numbers_outlier<-c()
	max_count_outlier<-ceiling(max_part_outlier*nrow(env$data_study))

	repeat
	{
		xL<-log(env$data_study[,i]+abs(min(env$data_study[,i]))+1)
		Q1<-quantile(xL)[[2]]
		Q3<-quantile(xL)[[4]]
		if(Q1!=Q3)
		{
			x1<-Q1-k*(Q3-Q1)
			x2<-Q3+k*(Q3-Q1)

			outlier<-which( xL<x1 | xL>x2 )
			numbers_outlier<-unique(c(numbers_outlier,outlier))
			if(length(numbers_outlier)>max_count_outlier)
			{
				i<-env$y_number+1
				numbers_outlier<-c()
				k<-k+step
			}
		}
		if(i==ncol(env$data_study))
		{
			break
		}
		i<-i+1
	}
	max_min<-data.frame("names"=character(),"max"=numeric(),"min"=numeric())
	if(length(numbers_outlier)>0)
	{
		x_without_outlier<-env$data_study[-numbers_outlier,]
		for(i in (env$y_number + 1):ncol(env$data_study))
		{
			max<-max(x_without_outlier[,i])
			min<-min(x_without_outlier[,i])
			env$data_study[env$data_study[,i]>max,i]<-max
			env$data_study[env$data_study[,i]<min,i]<-min
			if(is.null(env$data_test)==FALSE)
			{
				env$data_test[env$data_test[,i]>max,i]<-max
				env$data_test[env$data_test[,i]<min,i]<-min
			}
			max_min<-rbind(max_min,data.frame("names"=names(env$data_study)[i],"max"=max,"min"=min))
		}
	}
	max_min$names<-as.character(paste(max_min$names))
	return(max_min=max_min)
}

#' Standartization of data
#'
#' Standartization of data
#' @param env is data environment
#' @export
standartization<-function(env)
{
	env$data_study[,-env$y_number] <- as.data.frame(scale(env$data_study[,-env$y_number,drop=FALSE]))
	env$data_study[is.na(env$data_study)]<-0
	if( !is.null(env$data_test))
	{
		env$data_test[,-env$y_number] <- as.data.frame(scale(env$data_test[,-env$y_number,drop=FALSE]))
		env$data_test[is.na(env$data_test)]<-0
	}
}

#' Normalization of data
#'
#' Normalization of data
#' @param env is data environment
#' @export
normalization<-function(env)
{
	env$data_study[,-env$y_number] <- do.call(cbind, apply(env$data_study[,-env$y_number,drop=FALSE],2,function(x){list(x+abs(min(x))+1)}))
	if( !is.null(env$data_test))
	{
		env$data_test[,-env$y_number] <- do.call(cbind, apply(env$data_test[,-env$y_number,drop=FALSE],2,function(x){list(x+abs(min(x))+1)}))
	}
}

#' Prepare data
#'
#' Prepare data study and test
#' @param data_study data study
#' @param data_test data test
#' @param y_name is a name of oucome colum
#' @param feature_select logical value is it necessary feature select
#' @param stand logical value is it necessary feature standardization
#' @param multi is logical value delete multicollinearity
#' @param anomalies is logical value using of method for optimal boarding
#' @param box_cox is logical value replace anomalies
#' @param feature_step  is logical value select features by AUC
#' @return list
#' @export
prepare_study<-function(data_study, data_test = NULL,y_name = 'EndCat',
												feature_select = FALSE,multi=TRUE,anomalies=TRUE,box_cox=FALSE,stand=FALSE,feature_step=FALSE)
{
	if(is.null(data_study) | nrow(data_study) == nrow(data_study[data_study$y == 1,]))
	{
		print('Incorect data study')
		return(NULL)
	}
	#reNaming
	env<-new.env()
	env$y<-"y"
	columns_numbers<-seq(which( colnames(data_study) == y_name),ncol(data_study))

	colnames(data_study)[which( colnames(data_study) == y_name)] <- env$y
	env$data_study <- data_study[,columns_numbers]

	if( !is.null(data_test))
	{
		colnames(data_test)[which( colnames(data_test) == y_name)] <- env$y
		env$data_test <- data_test[,columns_numbers]
	}
	else
	{
		env$data_test <- NULL
	}
	env$y_number<-which( colnames(env$data_study) == env$y)

	if(feature_select)
	{
		columns_numbers<-significant_features(env)
		if(length(columns_numbers)<2)
		{
			print(paste('Not enough feature = ', columns_numbers, 'feature_select = ', feature_select))
			return( NULL )
		}
	}
	if(multi)
	{
		columns_numbers<-multi_features(env)
		if(length(columns_numbers)<2)
		{
			print(paste('Not enough feature = ', columns_numbers, 'feature_select = ', feature_select))
			return( NULL )
		}
	}
	max_min_info<-NULL
	if(anomalies)
	{
		max_min_info<-replace_anomalies(env)
	}
	lamda_info<-NULL
	if(box_cox)
	{
		lamda_info<-prepare_box_cox(env)
		print(lamda_info)
	}
	if(feature_step)
	{
		columns_numbers<-get_feature_step(env)
		if(length(columns_numbers)<2)
		{
			print(paste('Not enough feature = ', columns_numbers, 'feature_select = ', feature_select))
			return( NULL )
		}
	}
	if(stand)
	{
		standartization(env)
	}
	return(list(data_study = env$data_study,
							data_test  = env$data_test,
							prepare_info = list(feature_select = feature_select,
																	multi=multi,
																	anomalies=anomalies,
																	box_cox=box_cox,
																	stand = stand,
																	feature_step=feature_step,
																	lamda_info=lamda_info,
																	max_min_info=max_min_info,
																	significant_features=names(env$data_study)[-1],
																	y_name = y_name,
																	y = env$y)
	))
}

#' Get model
#'
#' Comparison of several models
#' @param data_study data study
#' @param model is name of model
#' @param y_name is a name of oucome colum
#' @return object
#' @export
get_model<-function(data_study = NULL, model='glm', y_name = 'y')
{
	predictorsNames <- names(data_study)[!names(data_study) %in% c(y_name)]
	#get_feature_step
	if (model =='glmnetCV')
	{
		objModel<-cv.glmnet(x = as.matrix(data_study[,predictorsNames,drop=FALSE]), y = data_study[,y_name], alpha = 1
												, family = "binomial", type.measure = "auc")
	}
	if (model == 'glm')
	{
		objModel<-glm(as.formula(paste("y~",paste(predictorsNames,collapse="+"))),data = data_study, family = binomial())
	}
	if(model=='gbm')
	{
		data_study$y<-ifelse(data_study$y==1,'Return','No')

		objControl <- trainControl(method='repeatedcv', repeats = 5,number=5,returnResamp='none',
															 summaryFunction = twoClassSummary,
															 classProbs = TRUE,verboseIter = FALSE)

		objModel <- train(data_study[,predictorsNames,drop=FALSE], data_study[,y_name],
											method = model,
											trControl = objControl,
											metric = "ROC",
											preProc = c("center", "scale")
											#,tuneGrid = expand.grid(alpha = 1, lambda = seq(0,1,by = 0.001))
		)
	}
	if(model=='rf')
	{

		bestmtry <- tuneRF(data_study[,predictorsNames], as.factor(data_study[,y_name]),  stepFactor=2, improve=0.05, ntree=500)
		data_study[,y_name]<-as.factor(make.names(data_study[,y_name]))

		objControl <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary, number = 5)

		objModel<- train(as.formula(paste(y_name,"~",paste(predictorsNames,collapse="+"))), data = data_study,
										 method = model,mpty=bestmtry[which.min(bestmtry[,2]),1],# tuneGrid = rfGrid,
										 metric = "ROC", trControl = objControl, importance = TRUE)
	}
	if(model=='rpart')
	{
		data_study[,y_name]<-make.names(as.factor(data_study[,y_name]))
		objControl <- trainControl(method='repeatedcv', number=5,returnResamp='none',
															 repeats =20, summaryFunction = twoClassSummary,
															 classProbs = TRUE,verboseIter = FALSE)
		objModel<- train(data_study[,predictorsNames,drop=FALSE], data_study[,y_name],
										 method = model,
										 metric = "ROC", trControl = objControl)
	}
	return(objModel)
}

#' Get prediction
#'
#' Prediction function for all methods
#' @param objModel model object
#' @param data data for preidction
#' @param model is name of model
#' @param y_name is a name of oucome colum
#' @return vecotr of prediction values
#' @export
get_prediction<-function(objModel, data = NULL, model='glm', y_name = 'y')
{
	predictorsNames <- names(data)[!names(data) %in% y_name]
	if(model=='glmnetCV')
	{
		predictions <- predict(object = objModel, as.matrix(data[,predictorsNames,drop=FALSE]), s = "lambda.min", type = "response")[,1]
	}
	if(model=='glm')
	{
		predictions<-1/(1+exp(-predict.glm(objModel ,data[,predictorsNames,drop=FALSE])))
	}
	if(model=='gbm')
	{
		predictions <- predict(object=objModel, data[,predictorsNames,drop=FALSE], type='prob')$Return
	}
	if(model=='rf')
	{
		predictions <- predict(object=objModel, data[,predictorsNames,drop=FALSE], type='prob')$X1
	}
	if(model=='rpart')
	{
		predictions <- predict(object=objModel, data[,predictorsNames,drop=FALSE], type='prob')$X1
	}
	return(unname(predictions))
}

#' Get quality classification
#'
#' Get quality classification
#' @param y true values for response
#' @param est estimation values for response
#' @param optimal_border using of method for optimal boarding
#' @param plot_roc is logical value for plot roc
#' @return list of quality metrics
#' @export
get_quality<-function(y,est, optimal_border = 0.5, plot_roc = FALSE)
{
	est<-data.frame("est" = est,"y" = y)
	N<-nrow(est)
	TP<-nrow(est[est$est >= optimal_border & est$y==1,])
	TN<-nrow(est[est$est < optimal_border & est$y==0,])
	AP<-nrow(est[est$y==1,])
	AN<-nrow(est[est$y==0,])

	FP<-nrow(est[est$est >= optimal_border & est$y==0,])
	FN<-nrow(est[est$est < optimal_border & est$y==1,])

	accuracy<-(TP+TN)/N
	precision<-ifelse((TP+AN-TN)==0,0,TP/(TP+AN-TN))
	recall<-ifelse(AP==0,0,TP/AP)
	fmesure<-ifelse((precision+recall)==0,0,2*(precision*recall)/(precision+recall))
	roc_obj<-roc(est$y,est$est,quiet = TRUE)
	plot_res <- 0
	if(plot_roc){
		plot_res = plot(roc_obj)
	}
	est<-est[order(est$est,decreasing=TRUE),]
	partTop<-sum(est$y[1:AP])/AP

	return(list("Optimal_border" = optimal_border,
							"AP" = AP,
							"AN" = AN,
							"TP" = TP,
							"TN" = TN,
							"FP" = FP,
							"FN" = FN,
							'Auc'= roc.curve(est$y,est$est,plotit = F)$auc,
							'F_measure' = fmesure,
							'PartTop' = partTop,
							'Plot_auc' = plot_res))
}

#' Cross-Validation
#'
#' get quality cross-validation
#' @param xy data study
#' @param n_folds folds number
#' @param models is name of model
#' @param optimal_border using of method for optimal boarding
#' @param y_name is a name of oucome colum
#' @return result data.frame
#' @export
get_quality_cv<- function(xy, n_folds=10, model, y_name = 'y',optimal_border = FALSE)
{
	columns_numbers<-seq(which( colnames(xy) == y_name),ncol(xy))
	xy<-xy[sample(nrow(xy)),columns_numbers]
	if(n_folds>1)
	{
		xy_0<-xy[xy$y==0,]
		xy_1<-xy[xy$y==1,]
		folds_i_0 <- sample(rep(1:n_folds, length.out =nrow(xy_0)))
		folds_i_1 <- sample(rep(1:n_folds, length.out =nrow(xy_1)))
		quality<-data.frame()
		for(j in 1:n_folds)
		{
			study_xy <- rbind(xy_0[-which(folds_i_0 == j), ],xy_1[-which(folds_i_1 == j), ])
			test_xy <-  rbind(xy_0[ which(folds_i_0 == j), ],xy_1[ which(folds_i_1 == j), ])
			model_xy<-get_model(study_xy,model)
			predict_xy<-get_prediction(model_xy,test_xy,model)
			quality_xy<-get_quality(test_xy[,y_name],predict_xy,optimal_border,plot_roc = FALSE)
			model_info<-get_model_info(model_xy,model)
			quality_xy<-c(quality_xy,count_features_in_model=nrow(model_info[model_info[,2] < 0 | model_info[,2] > 0, ]))

			quality<-rbind(quality,quality_xy)
		}
		quality<-as.data.frame(as.list(apply(quality,2,"mean")))
	}
	if(n_folds==1)
	{
		model_xy<-get_model(xy,model)
		predict_xy<-get_prediction(model_xy,xy,model)
		quality_xy<-get_quality(xy[,y_name],predict_xy,optimal_border,plot_roc = FALSE)
		model_info<-get_model_info(model_xy,model)
		quality<-as.data.frame(c(quality_xy,count_features_in_model=nrow(model_info[model_info[,2] < 0 | model_info[,2] > 0, ])))
	}
	if(n_folds<1)
	{
		return(NULL)
	}


	return(quality)
}

#' model info
#'
#' Get model informatoin
#' @param obj_model model
#' @param model model name
#' @return result data.frame
#' @export
get_model_info<-function(obj_model,model)
{
	if(model =='glmnetCV')
	{
		coef<-as.data.frame(as.matrix(coef(obj_model)))
		model_info<-cbind(row.names(coef),coef)
		names(model_info)<-c("Evidence","Coef")
	}
	if (model == 'glm')
	{
		coef<-as.data.frame(coef(obj_model))
		model_info<-cbind(row.names(coef),coef)
		names(model_info)<-c("Evidence","Coef")
	}
	#carret
	if(model=='gbm')
	{
		model_info<-summary(obj_model)
		names(model_info)<-c("Evidence","Influence")

	}
	if(model=='rf')
	{
		inf<-varImp(object=obj_model)$importance
		model_info<-cbind(row.names(inf),inf)[-3]
		names(model_info)<-c("Evidence","Importance")
	}
	if(model=='rpart')
	{
		inf<-varImp(object=obj_model)$importance
		model_info<-cbind(row.names(inf),inf)[-3]
		names(model_info)<-c("Evidence","Importance")
	}
	row.names(model_info)<-NULL
	return(model_info)
}

#' Comparison classifications
#'
#' Comparison of several models
#' @param data_study data study
#' @param data_test data test
#' @param models is name of models
#' @param y_name is a name of oucome colum
#' @param cv_n_folds folds number
#' @param is_optimal_border  using of method for optimal boarding
#' @param stand logical vector is it necessary feature standardization
#' @param multi is logical vector delete multicollinearity
#' @param anomalies is logical vector using of method for optimal boarding
#' @param box_cox is logical vector replace anomalies
#' @param feature_step  is logical vector select features by AUC
#' @param feature_select is logical vector is it necessary feature select
#' @return list
#' @export
compare_classifications<-function(data_study, data_test = NULL, y_name = 'y', models=c("gbm","glm"),is_optimal_border=TRUE,
																	feature_select=c(FALSE),multi=c(TRUE,FALSE),anomalies=c(TRUE,FALSE),box_cox=c(FALSE),
																	feature_step=c(TRUE,FALSE),stand=TRUE,n_folds=10,plot_test=FALSE)
{
	options(digits=5)
	options(scipen = 999)

	name_list<-c()
	model_info_all<-data.frame("Evidence"=character())
	quality_test_all<-data.frame()
	quality_cv_all<-data.frame()
	quality_study_all<-data.frame()

	i<-1
	list_plot<-list()

	for (fsl in feature_select)
	{
		for(mu in multi)
		{
			for(an in anomalies)
			{
				for(bc in box_cox)
				{
					for(fst in feature_step)
					{
						start.time <- Sys.time()
						data<-prepare_study(data_study, data_test ,y_name ,fsl,mu,an,bc,stand,fst)
						end.time <- Sys.time()
						time_prepare <- end.time - start.time
						if(!is.null(data))
						{
							for (m in models)
							{
								#create name
								name_fsl<-ifelse(fsl,"fsl","nfsl")
								name_mu<-ifelse(mu,"mu","nmu")
								name_an<-ifelse(an,"an","nan")
								name_bc<-ifelse(bc,"bc","nbc")
								name_fst<-ifelse(fst,"fst","nfst")
								name_st<-ifelse(stand,"st","nst")
								name_full<-paste(m,name_fsl,name_mu,name_an,name_bc,name_fst,name_st,sep = "_")
								name_list<-c(name_list,name_full)
								print(name_full)

								#model info
								start.time <- Sys.time()
								model<-get_model(data$data_study,m)
								if(is_optimal_border)
								{
									optimal_border<-optimal_border<-get_optimal_border(model, data$data_study, m, 'y')
								}
								else
								{
									optimal_border<-0.5
								}

								end.time <- Sys.time()
								time_model <- end.time - start.time
								model_info<-get_model_info(model,m)
								count_features<-nrow(model_info[model_info[,2] < 0 | model_info[,2] > 0, ])

								model_info_all<-merge(model_info_all,model_info, by=c("Evidence"),all = TRUE)
								names(model_info_all)[ncol(model_info_all)]<-paste(names(model_info_all)[ncol(model_info_all)] , name_full ,sep = "_")

								model_info[model_info == 0] <- NA
								model_info[,-1]<-rank(-abs(model_info[,-1]),"max",na.last = "keep")
								model_info_all<-merge(model_info_all,model_info, by=c("Evidence"),all = TRUE)
								names(model_info_all)[ncol(model_info_all)]<-paste("rank",names(model_info_all)[ncol(model_info_all)-1] , name_full ,sep = "_")

								#quality study
								quality_study<-get_quality_cv(data$data_study, n_folds=1, m, "y",optimal_border)
								quality_study_all<-rbind(quality_study_all,cbind(Full_Name=name_full,
																																 Time_All=time_model+time_prepare,
																																 Count_features=count_features,
																																 quality_study[names(quality_study)!='Plot_auc']))
								#quality cv
								if(n_folds > 1)
								{
									quality_cv<-get_quality_cv(data$data_study, n_folds=n_folds, m, "y",optimal_border)
									quality_cv_all<-rbind(quality_cv_all,cbind(Full_Name=name_full,
																														 Time_All=time_model+time_prepare,
																														 Count_features=count_features,
																														 quality_cv[names(quality_cv)!='Plot_auc']))
								}
								#quality test
								if(!is.null(data_test))
								{
									start.time <- Sys.time()
									predict_test<-get_prediction(model,data$data_test,m)
									end.time <- Sys.time()
									time_predict <- end.time - start.time
									if(plot_test)
									{
										quality_test<-get_quality(data$data_test[,"y"],predict_test,optimal_border,plot_roc = TRUE)
										list_plot[[i]]<-quality_test$Plot_auc
										i=i+1
									}
									else
									{
										quality_test<-get_quality(data$data_test[,"y"],predict_test,optimal_border,plot_roc = FALSE)
									}
									quality_test_all<-rbind(quality_test_all,cbind(Full_Name=name_full,
																																 Time_All=time_predict,
																																 Count_features=count_features,
																																 as.data.frame(quality_test[names(quality_test)!='Plot_auc'])))
								}
							}
						}
					}
				}
			}
		}


	}

	if(!is.null(data_test))
	{
		if(plot_test)
		{
			qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
			color_list = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
			for(l in 1:length(list_plot))
			{
				#plot roc
				if(l==1)
				{
					plot(list_plot[[l]],col = color_list[l],lwd = 2, main = "ROC for test")
				}
				else
				{
					plot(list_plot[[l]],col = color_list[l],lwd = 2,add = TRUE)
				}
			}
			legend("bottom", lwd = 5, col=color_list, name_list,cex=0.60,ncol=3)
		}
	}
	return(list(quality_test=quality_test_all,
							quality_cv=quality_cv_all,
							quality_study=quality_study_all,
							model_info=model_info_all,
							list_plot=list_plot
	))
}


#' Get optimal border
#'
#' Get optimal border for get quality
#' @param model_obj is model derived from classification_go
#' @param data data study
#' @param model is name of model
#' @param y_name is a name of oucome colum
#' @return list
#' @export
get_optimal_border<-function(model_obj, data = NULL, model='glm', y_name = 'y')
{
	predictStudy<-get_prediction(model_obj,data,model,y_name = 'y')
	optimBorder<-optimalCutoff(actuals = data[,names(data) %in% y_name], predictedScores = predictStudy,optimiseFor ="Both" )
	return(optimBorder)
}


#' Classification
#'
#' Start classification
#' @param data_study is data study
#' @param data_test is data test
#' @param model is name of model
#' @param y_name is a name of oucome colum
#' @param feature_select is logical value is it necessary feature select
#' @param n_folds is folds number
#' @param is_optimal_border is logical value using of method for optimal boarding
#' @param stand is logical value is it necessary feature standardization
#' @param plot_roc is logical value is it necessary plot roc
#' @param multi is logical value delete multicollinearity
#' @param anomalies is logical value using of method for optimal boarding
#' @param box_cox is logical value replace anomalies
#' @param feature_step  is logical value select features by AUC
#' @return list
#' @export
classification_go<-function(data_study, data_test = NULL, model = 'glmnet', y_name = 'EndCat',is_optimal_border=FALSE,plot_roc=FALSE,
														feature_select=TRUE,multi=TRUE,anomalies=TRUE,box_cox=TRUE,feature_step=TRUE,stand=FALSE,n_folds=0)
{
	data<-prepare_study(data_study, data_test ,y_name,feature_select,multi,anomalies,box_cox,stand,feature_step)
	result_model<-get_model(data$data_study,model)
	print(result_model)
	if(is_optimal_border)
	{
		optimal_border<-get_optimal_border(result_model, data$data_study, model, 'y')
	}
	else
	{
		optimal_border<-0.5
	}


	quality_cv<-NULL
	quality_study<-get_quality_cv(data$data_study, n_folds=1, model, "y",optimal_border)
	if(n_folds > 1)
	{
		quality_cv<-get_quality_cv(data$data_study, n_folds=n_folds, model, "y",optimal_border)
	}
	quality_test<-NULL
	if(!is.null(data_test))
	{
		predict_test<-get_prediction(result_model,data$data_test,model)
		quality_test<-get_quality(data$data_test[,"y"],predict_test,optimal_border,plot_roc)
	}
	result<-list( "quality_cv"=quality_cv,
								"quality_test"=quality_test,
								"quality_study"=quality_study,
								"model"=list(obj_model=result_model,model_type=model),
								"prepare_info"=data$prepare_info
	)
	return(result)
}

#' Prediction
#'
#' Start prediction
#' @param data_predict is data for prediction
#' @param res_model is model derived from classification_go
#' @return result vector with predict
#' @export
prediction_go<-function(data_predict, res_model)
{
	columns_numbers<-which(colnames(data_predict) %in% res_model$prepare_info$significant_features)
	data_predict<-data_predict[,columns_numbers]

	if(res_model$prepare_info$anomalies)
	{
		for(i in 1:ncol(data_predict))
		{
			max<-res_model$prepare_info$max_min_info$max[res_model$prepare_info$max_min_info$name==names(data_predict)[i]]
			min<-res_model$prepare_info$max_min_info$min[res_model$prepare_info$max_min_info$name==names(data_predict)[i]]
			data_predict[data_predict[,i]>max,i]<-max
			data_predict[data_predict[,i]<min,i]<-min
		}
	}
	if(res_model$prepare_info$box_cox)
	{
		for(i in 1:ncol(data_predict))
		{
			data_predict[,i] <- box_cox(data_predict[,i],res_model$lamda_info$lamda[res_model$lamda_info$name==names(data_predict)[i]])
		}
	}
	if(res_model$prepare_info$stand)
	{
		data_predict <- as.data.frame(scale(data_predict))
	}
	predict<-get_prediction(res_model$model$obj_model,data_predict,res_model$model$model_type)
	return(predict)
}
