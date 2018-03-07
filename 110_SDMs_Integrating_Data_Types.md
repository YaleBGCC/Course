---
title: "Integrating different data types into SDMs"
author: "Cody Shank"
---



<!-- <div> -->
<!-- <iframe src="05_presentation/05_Spatial.html" width="100%" height="700px"> </iframe> -->
<!-- </div> -->

<!-- <div> -->
<!-- <object data="101_assets/SDM101_Intro.pdf" type="application/pdf" width="100%" height="600px">  -->
<!--   <p>It appears you don't have a PDF plugin for this browser. -->
<!--    No biggie... you can <a href="02_assets/02_DataWrangling.pdf">click here to -->
<!--   download the PDF file.</a></p>   -->
<!--  </object> -->
<!--  </div> -->
<!--  <p><a href="101_assets/SDM101_Intro.pdf">Download the PDF of the presentation</a></p>   -->


[<i class="fa fa-file-code-o fa-3x" aria-hidden="true"></i> The R Script associated with this page is available here](110_SDMs_Integrating_Data_Types.R).  Download this file and open it (or copy-paste into a new script) with RStudio so you can follow along.  

# Setup

Load the functions, and library (make sure 'dismo' is installed)

> Note that you don't need to wade through these functions to start, just run them and more to the next chunk where we begin using them.


```r
library(dismo)

logit = function(pp) { log(pp) - log(1-pp) }
expit = function(eta) {1/(1+exp(-eta))}

ObsInfo.po = function(param) {
  
  beta = param[1:dim(X.po)[2]]
  alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
  
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  p.po = expit(W.po %*% alpha)
  
  nxcovs = length(beta)
  nwcovs = length(alpha)
  
  nparams = nxcovs + nwcovs
  Hmat = matrix(nrow=nparams, ncol=nparams)
  
  #  beta partials
  for (i in 1:nxcovs) {
    for (j in 1:i) {
      Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
      Hmat[j,i] = Hmat[i,j]
    }
  }
  
  # alpha partials
  for (i in 1:nwcovs) {
    for (j in 1:i) {
      Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.po[,i] * W.po[,j] * p.po * (1-p.po))
      Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
    }
  }
  
  # alpha-beta partials
  for (i in 1:nwcovs) {
    for (j in 1:nxcovs) {
      Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
      Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
    }
  }
  
  Hmat
}

negLL.po = function(param) {
  
  beta = param[1:dim(X.po)[2]]
  alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
  
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  logL.po = sum(X.po %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
  
  (-1)*sum(logL.po)
}

negLL.pa = function(param) {
  beta = param[1:dim(X.pa)[2]]
  lambda.pa = exp(X.pa %*% beta)
  
  alpha = param[(dim(X.pa)[2]+1):(dim(X.pa)[2]+dim(W.pa)[2])]
  p.pa = expit(W.pa %*% alpha)
  
  logL.pa = rep(NA,n.pa)
  
  for (i in 1:n.pa) {
    yvec=y.pa[i,]
    navec=is.na(yvec)
    nd=sum(yvec[!navec])
    nj=sum(!navec)
    pvec=p.pa[i,]
    cp= (pvec^yvec)*((1-pvec)^(1-yvec))
    cp[navec]=1
    logL.pa[i]= log(prod(cp)*(1-exp(-lambda.pa[i]*area.pa[i])) + ifelse(nd==0,1,0)*exp(-lambda.pa[i]*area.pa[i]))
  }
  
  (-1)*sum(logL.pa)
}

negLL.poANDpa = function(param)  {
  
  param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
  param.pa = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pa)[2]))]
  negLL.po(param.po) + negLL.pa(param.pa)
}

rarefy = function(pointSample,r){
  pointSampleCopy = pointSample
  pointSampleTrain = pointSample[0,]
  pointSampleTest = pointSample[0,]
  nsamples=1
  while(dim(pointSampleCopy)[1] > 0){
    rowID = sample(nrow(pointSampleCopy), 1)
    s1 = pointSampleCopy[rowID, ]
    pointSampleTrain[nsamples,] = s1
    coordinates(s1) = c("X","Y")
    proj4string(s1) = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    coordinates(pointSampleCopy) = c("X","Y")
    proj4string(pointSampleCopy) = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    pointSampleCopy$distance = (spDistsN1(pointSampleCopy,s1)>r)
    pointSampleCopy = data.frame(pointSampleCopy)
    pointSampleCopy = pointSampleCopy[-rowID,]
    pointSampleTemp = pointSampleCopy[!pointSampleCopy$distance,]
    #rownames(pointSampleTemp) = NULL
    #since I reinstalled R, for some reason it adds an extra column during this function, I get rid of it below
    pointSampleTemp$distance = NULL
    pointSampleTemp = pointSampleTemp[,-dim(pointSampleTemp)[2]]
    pointSampleTest = rbind(pointSampleTest,pointSampleTemp)
    pointSampleCopy  = pointSampleCopy[pointSampleCopy$distance,]
    pointSampleCopy$distance = NULL
    #since I reinstalled R, for some reason it adds an extra column during this function, I get rid of it below
    pointSampleCopy = pointSampleCopy[,-dim(pointSampleCopy)[2]]
    nsamples=nsamples+1
    #rownames(pointSampleTrain) = NULL
  }
  return(list(pointSampleTrain,pointSampleTest,r))
}
```

Prepare the data (First try running the models on the thinned presence data)

Model settings

```r
minrecipCondNum = 1e-6
resolution = 4 # this is in km
site.area = resolution^2
n.samples = 6
thinning.radius = resolution*1000*2^0.5 # this is in m
```

Download data

> Read the the Notes in this section, or stuff will break, and we won't feel bad for you...


```r
#workspace = "/Users/Cody/Dropbox/YIBS exercise/"
#predictors.folder = paste0("/Users/Cody/Dropbox/YIBS exercise/","predictors")
# !!!NOTE: you must set predictors.folder to a path on your computer
predictors.folder='/Users/ctg/Dropbox/Projects/Workshops/YaleBGCCourses/110_assets'
isDownloaded=list.files(predictors.folder,pattern='predictors.zip',full.names=T)
if(!file.exists(isDownloaded)){
  download.file('https://cmerow.github.io/YaleBGCCourses/110_assets/predictors.zip',
                paste0(predictors.folder,'/predictors.zip'))
}
# !!! NOTE: you must go to the download folder and unzip the file before proceeding

#presenceData_all = read.csv(paste0(workspace,"presence_only_proj.csv"))
presenceData_all = read.csv('https://cmerow.github.io/YaleBGCCourses/110_assets/presence_only_proj.csv')
#paData_all = read.csv(paste0(workspace,"detection_histories_proj_10samples.csv"))
paData_all = read.csv('https://cmerow.github.io/YaleBGCCourses/110_assets/detection_histories_proj_10samples.csv')
```

Prep for modeling


```r
presenceData = rarefy(presenceData_all,thinning.radius)[[1]]
paData = rarefy(paData_all,thinning.radius)[[1]]
#paData = paData_all
#presenceData = presenceData_all

allStack.files = list.files(paste0(predictors.folder,'/predictors'),full.names=TRUE,pattern=".tif")
allStack = stack(allStack.files)

names.x.po.covs = c("chelsa_temp_seasonality","chelsa_precip_seasonality","chelsa_max_temp_warmest_month","chelsa_annual_precip",
                    "chelsa_temp_seasonality_sq","chelsa_precip_seasonality_sq","chelsa_max_temp_warmest_month_sq","chelsa_annual_precip_sq",
                    "treecover2000","distancePA","EVI","EVI_sq","forestloss","road_length","fire_density","mean_slope","treecover2000_EVI")
names.x.pa.covs = names.x.po.covs
names.w.po.covs = c("forest","protected","road_distance","road_distance_sq","mean_slope")
names.w.pa.covs = c("treecover2000","distancePA","road_distance","road_distance_sq","mean_slope",
                    "forestloss","road_length","fire_density","EVI")
names.w.pa.sampling.covs = c("Tapir","Cat","OnTrail")

names_allCovs = unique(c(names.x.po.covs,names.w.po.covs,names.x.pa.covs,names.w.pa.covs))

sgrid = stack(subset(allStack,names_allCovs))


xycov = na.omit(rasterToPoints(sgrid))
xy = xycov[,1:2]
des = as.data.frame(xycov[,3:dim(xycov)[2]])
names(des)=names_allCovs
X.back = as.matrix(cbind(rep(1, dim(des)[1]), des[,names.x.po.covs]))
W.back = as.matrix(cbind(rep(1, dim(des)[1]), des[,names.w.po.covs]))

po.des = data.frame(na.omit(extract(sgrid,presenceData[,c("X","Y")])))
n.po = dim(po.des)[1]
X.po = as.matrix(cbind(rep(1, dim(po.des)[1]), po.des[,names.x.po.covs]))
W.po = as.matrix(cbind(rep(1, dim(po.des)[1]), po.des[,names.w.po.covs]))

area.back = rep(site.area, dim(X.back)[1])

pa.extract = extract(sgrid,paData[,c("X","Y")])
filterNA = !is.na(pa.extract)
filterNA = apply(filterNA, 1, prod)
filterNA = filterNA==1
pa.des = as.data.frame(na.omit(pa.extract))
paDataSelect = paData[filterNA,]
X.pa = as.matrix(cbind(rep(1,dim(pa.des)[1]),pa.des[,names.x.pa.covs]))
W.pa.extract = as.matrix(pa.des[,names.w.pa.covs])
W.pa = as.matrix(cbind(rep(1, dim(W.pa.extract)[1]),W.pa.extract,paDataSelect[,names.w.pa.sampling.covs]))

# detection histories
y.pa = paDataSelect[,paste("sample",seq(1,n.samples),sep="_")]
y.pa = as.matrix(y.pa)

n.pa = dim(y.pa)[1]
J.pa =  dim(y.pa)[2]

area.pa = rep(site.area, n.pa)

betaGuess = rep(0, dim(X.po)[2])
alphaGuess.po = c(logit(n.po/nrow(xy)),rep(0, (dim(W.po)[2]-1))) # use naive detectability
#alphaGuess.po = rep(0, dim(W.po)[2]) 
alphaGuess.pa = c(logit(sum(y.pa,na.rm=T)/sum(y.pa==0,na.rm=T)),rep(0, (dim(W.pa)[2]-1))) # use naive detectability
#alphaGuess.pa = rep(0, dim(W.pa)[2])
```

Fit the Integrated SDM and predict intensity and occupancy surfaces. This may take 5-8 minutes.


```r
paramGuess = c(betaGuess, alphaGuess.po, alphaGuess.pa)
fit.poANDpa = NA

start.time = Sys.time()
possibleError.poANDpa = tryCatch(
  expr = (fit.poANDpa = optim(par=paramGuess, fn=negLL.poANDpa, method='BFGS', hessian=TRUE, control=list(maxit=200))),
  error=function(e) e
)
end.time = Sys.time()
elapsed.time = end.time - start.time 

if(!inherits(possibleError.poANDpa, "error")){
  recipCondNum.poANDpa = NA
  se.poANDpa = rep(NA, length(fit.poANDpa$par))
  if (fit.poANDpa$convergence==0) {
    hess = fit.poANDpa$hessian
    ev = eigen(hess)$values
    recipCondNum.poANDpa = ev[length(ev)]/ev[1]
    if (recipCondNum.poANDpa>minrecipCondNum) {
      vcv = chol2inv(chol(hess))
      se.poANDpa = sqrt(diag(vcv))
    }
  }	
  
  linearPredictor.poANDpa.sdm = X.back %*% fit.poANDpa$par[1:dim(X.back)[2]]
  predict.poANDpa.sdm= exp(linearPredictor.poANDpa.sdm)
  raster.poANDpa.sdm = raster(sgrid)
  cells = cellFromXY(raster.poANDpa.sdm, xy)
  raster.poANDpa.sdm[cells] = predict.poANDpa.sdm
  raster.poANDpa.sdm.psi = 1-exp(-1*raster.poANDpa.sdm*site.area)
  
  poANDpa.coefs = data.frame(t(unlist(fit.poANDpa$par)))
  poANDpa.se = data.frame(t(unlist(se.poANDpa)))
  poANDpa.population =sum(predict.poANDpa.sdm*site.area,na.rm=TRUE)
  
}
```

Fit Site-occupancy model only

```r
paramGuess = c(betaGuess, alphaGuess.pa)
fit.pa = NA

start.time = Sys.time()
possibleError.pa = tryCatch(
  expr = (fit.pa = optim(par=paramGuess, fn=negLL.pa, method='BFGS', hessian=TRUE,control=list(maxit=200))),
  error=function(e) e
)
end.time = Sys.time()
elapsed.time = end.time - start.time 

if(!inherits(possibleError.pa, "error")){
  recipCondNum.pa = NA
  se.pa = rep(NA, length(fit.pa$par))
  if (fit.pa$convergence==0) {
    hess = fit.pa$hessian
    ev = eigen(hess)$values
    recipCondNum.pa = ev[length(ev)]/ev[1]
    if (recipCondNum.pa>minrecipCondNum) {
      vcv = chol2inv(chol(hess))
      se.pa = sqrt(diag(vcv))
    }
  }
  
  linearPredictor.pa.sdm = X.back %*% fit.pa$par[1:dim(X.back)[2]]
  predict.pa.sdm = exp(linearPredictor.pa.sdm)
  raster.pa.sdm = raster(sgrid)
  cells = cellFromXY(raster.pa.sdm, xy)
  raster.pa.sdm[cells] = predict.pa.sdm
  raster.pa.sdm.psi = 1-exp(-1*raster.pa.sdm*resolution^2)
  
  pa.coefs = data.frame(t(unlist(fit.pa$par)))
  pa.se = data.frame(t(unlist(se.pa)))
  
  pa.population =sum(predict.pa.sdm*site.area,na.rm=TRUE)
  
}
```

Fit Presence-only model only


```r
paramGuess = c(betaGuess, alphaGuess.po)

start.time = Sys.time()
possibleError.po = tryCatch(
  expr = (fit.po = optim(par=paramGuess, fn=negLL.po, method='BFGS', hessian=FALSE, control=list(maxit=200))),
  error=function(e) e
)
end.time = Sys.time()
elapsed.time = end.time - start.time # this takes about 2-3 minutes to run

if(!inherits(possibleError.po, "error")){
  recipCondNum.po = NA
  se.po = rep(NA, length(fit.po$par))
  if (fit.po$convergence==0) {		
    hess = ObsInfo.po(fit.po$par)
    possibleError.hess = tryCatch(expr = (ev = eigen(hess)$values),error=function(e) e)
    if(!inherits(possibleError.hess, "error")){
      recipCondNum.po = ev[length(ev)]/ev[1]
      if (!is.na(recipCondNum.po)){
        if (recipCondNum.po>minrecipCondNum) {
          vcv = chol2inv(chol(hess))
          se.po = sqrt(diag(vcv))
        }
      }
      if (is.na(recipCondNum.po)) {se.po=rep(NA, length(fit.po$par))}
    }
  }else{se.po=rep(NA, length(fit.po$par))}
  
  linearPredictor.po.sdm = X.back %*% fit.po$par[1:dim(X.back)[2]]
  predict.po.sdm = exp(linearPredictor.po.sdm)
  raster.po.sdm = raster(sgrid)
  cells = cellFromXY(raster.po.sdm, xy)
  raster.po.sdm[cells] = predict.po.sdm
  raster.po.sdm.psi = 1-exp(-1*raster.po.sdm*resolution^2)
  
  po.coefs = data.frame(t(unlist(fit.po$par)))
  po.se = data.frame(t(unlist(se.po)))
  
  po.population =sum(predict.po.sdm*site.area,na.rm=TRUE)
}
```

Plot the intensity of the different models


```r
plot(raster.poANDpa.sdm)
```

![](110_SDMs_Integrating_Data_Types_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
plot(raster.pa.sdm)
```

![](110_SDMs_Integrating_Data_Types_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```r
plot(raster.po.sdm)
```

![](110_SDMs_Integrating_Data_Types_files/figure-html/unnamed-chunk-9-3.png)<!-- -->


Naming of coefficient and standard error outputs



Print out total population and look at coefficients and se's

```r
poANDpa.coefs
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["x.distancePA"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["x.EVI"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["x.EVI_sq"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["x.forestloss"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["x.road_length"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["x.fire_density"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["x.mean_slope"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000_EVI"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["alpha0.po"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["w.po.forest"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["w.po.protected"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance_sq"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["w.po.mean_slope"],"name":[24],"type":["dbl"],"align":["right"]},{"label":["alpha0.pa"],"name":[25],"type":["dbl"],"align":["right"]},{"label":["w.pa.treecover2000"],"name":[26],"type":["dbl"],"align":["right"]},{"label":["w.pa.distancePA"],"name":[27],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance"],"name":[28],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance_sq"],"name":[29],"type":["dbl"],"align":["right"]},{"label":["w.pa.mean_slope"],"name":[30],"type":["dbl"],"align":["right"]},{"label":["w.pa.forestloss"],"name":[31],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_length"],"name":[32],"type":["dbl"],"align":["right"]},{"label":["w.pa.fire_density"],"name":[33],"type":["dbl"],"align":["right"]},{"label":["w.pa.EVI"],"name":[34],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Tapir"],"name":[35],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Cat"],"name":[36],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.OnTrail"],"name":[37],"type":["dbl"],"align":["right"]}],"data":[{"1":"-0.5815728","2":"-0.4827651","3":"0.3224282","4":"0.1580765","5":"1.438296","6":"-0.323158","7":"-0.2579493","8":"0.0003289864","9":"-0.8878898","10":"1.264749","11":"-0.05291894","12":"-0.641176","13":"-0.1248524","14":"-0.08094696","15":"-0.2039516","16":"-0.3560615","17":"3.206404","18":"0.1167512","19":"-8.02003","20":"0.2398945","21":"1.48111","22":"-0.462504","23":"0.06467697","24":"-3.525587","25":"-3.127045","26":"0.1488322","27":"0.1273323","28":"-0.2816245","29":"-0.141174","30":"-0.8534065","31":"-0.2323479","32":"-0.2149084","33":"0.6722108","34":"0.4258688","35":"1.026512","36":"0.02506977","37":"0.5422595"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
poANDpa.se
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0.se"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality.se"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality.se"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month.se"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip.se"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq.se"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq.se"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq.se"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq.se"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000.se"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["x.distancePA.se"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["x.EVI.se"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["x.EVI_sq.se"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["x.forestloss.se"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["x.road_length.se"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["x.fire_density.se"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["x.mean_slope.se"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000_EVI.se"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["alpha0.po.se"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["w.po.forest.se"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["w.po.protected.se"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance.se"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance_sq.se"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["w.po.mean_slope.se"],"name":[24],"type":["dbl"],"align":["right"]},{"label":["alpha0.pa.se"],"name":[25],"type":["dbl"],"align":["right"]},{"label":["w.pa.treecover2000.se"],"name":[26],"type":["dbl"],"align":["right"]},{"label":["w.pa.distancePA.se"],"name":[27],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance.se"],"name":[28],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance_sq.se"],"name":[29],"type":["dbl"],"align":["right"]},{"label":["w.pa.mean_slope.se"],"name":[30],"type":["dbl"],"align":["right"]},{"label":["w.pa.forestloss.se"],"name":[31],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_length.se"],"name":[32],"type":["dbl"],"align":["right"]},{"label":["w.pa.fire_density.se"],"name":[33],"type":["dbl"],"align":["right"]},{"label":["w.pa.EVI.se"],"name":[34],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Tapir.se"],"name":[35],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Cat.se"],"name":[36],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.OnTrail.se"],"name":[37],"type":["dbl"],"align":["right"]}],"data":[{"1":"1.600649","2":"0.1128138","3":"0.1742078","4":"0.1714821","5":"0.3042578","6":"0.1065095","7":"0.1193622","8":"0.04132669","9":"0.2358779","10":"0.1779658","11":"0.1372942","12":"0.1412511","13":"0.1073277","14":"0.06351127","15":"0.1808518","16":"0.2144564","17":"1.936653","18":"0.1559999","19":"1.599514","20":"0.2239499","21":"0.2140094","22":"0.1314043","23":"0.04056721","24":"1.940785","25":"0.6644296","26":"0.6747717","27":"0.3868287","28":"0.479111","29":"0.1923115","30":"0.2830889","31":"0.1588521","32":"0.5618555","33":"1.059699","34":"0.4174205","35":"0.4663462","36":"0.5849214","37":"0.4744538"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
pa.coefs
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["x.distancePA"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["x.EVI"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["x.EVI_sq"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["x.forestloss"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["x.road_length"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["x.fire_density"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["x.mean_slope"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000_EVI"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["alpha0.pa"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["w.pa.treecover2000"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["w.pa.distancePA"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_distance_sq"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["w.pa.mean_slope"],"name":[24],"type":["dbl"],"align":["right"]},{"label":["w.pa.forestloss"],"name":[25],"type":["dbl"],"align":["right"]},{"label":["w.pa.road_length"],"name":[26],"type":["dbl"],"align":["right"]},{"label":["w.pa.fire_density"],"name":[27],"type":["dbl"],"align":["right"]},{"label":["w.pa.EVI"],"name":[28],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Tapir"],"name":[29],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.Cat"],"name":[30],"type":["dbl"],"align":["right"]},{"label":["w.pa.sampling.OnTrail"],"name":[31],"type":["dbl"],"align":["right"]}],"data":[{"1":"-0.0003147523","2":"-0.0001659311","3":"8.4881e-05","4":"-0.0002303436","5":"-8.370776e-05","6":"-0.0003278914","7":"-0.0001055936","8":"0.0004332846","9":"-0.0008800678","10":"-0.00033809","11":"-3.748323e-05","12":"-0.0002925802","13":"-0.0002895514","14":"0.0001275581","15":"6.540034e-05","16":"8.687111e-05","17":"0.0001532065","18":"-0.0003026035","19":"-4.038865","20":"0.3315309","21":"0.3346737","22":"-0.3705029","23":"-0.06583263","24":"-0.2911668","25":"-0.2036972","26":"0.0620944","27":"0.2459344","28":"0.6479038","29":"1.239141","30":"0.5133014","31":"0.2908128"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
pa.se
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0.se"],"name":[1],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality.se"],"name":[2],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality.se"],"name":[3],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month.se"],"name":[4],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_annual_precip.se"],"name":[5],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq.se"],"name":[6],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq.se"],"name":[7],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq.se"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq.se"],"name":[9],"type":["lgl"],"align":["right"]},{"label":["x.treecover2000.se"],"name":[10],"type":["lgl"],"align":["right"]},{"label":["x.distancePA.se"],"name":[11],"type":["lgl"],"align":["right"]},{"label":["x.EVI.se"],"name":[12],"type":["lgl"],"align":["right"]},{"label":["x.EVI_sq.se"],"name":[13],"type":["lgl"],"align":["right"]},{"label":["x.forestloss.se"],"name":[14],"type":["lgl"],"align":["right"]},{"label":["x.road_length.se"],"name":[15],"type":["lgl"],"align":["right"]},{"label":["x.fire_density.se"],"name":[16],"type":["lgl"],"align":["right"]},{"label":["x.mean_slope.se"],"name":[17],"type":["lgl"],"align":["right"]},{"label":["x.treecover2000_EVI.se"],"name":[18],"type":["lgl"],"align":["right"]},{"label":["alpha0.pa.se"],"name":[19],"type":["lgl"],"align":["right"]},{"label":["w.pa.treecover2000.se"],"name":[20],"type":["lgl"],"align":["right"]},{"label":["w.pa.distancePA.se"],"name":[21],"type":["lgl"],"align":["right"]},{"label":["w.pa.road_distance.se"],"name":[22],"type":["lgl"],"align":["right"]},{"label":["w.pa.road_distance_sq.se"],"name":[23],"type":["lgl"],"align":["right"]},{"label":["w.pa.mean_slope.se"],"name":[24],"type":["lgl"],"align":["right"]},{"label":["w.pa.forestloss.se"],"name":[25],"type":["lgl"],"align":["right"]},{"label":["w.pa.road_length.se"],"name":[26],"type":["lgl"],"align":["right"]},{"label":["w.pa.fire_density.se"],"name":[27],"type":["lgl"],"align":["right"]},{"label":["w.pa.EVI.se"],"name":[28],"type":["lgl"],"align":["right"]},{"label":["w.pa.sampling.Tapir.se"],"name":[29],"type":["lgl"],"align":["right"]},{"label":["w.pa.sampling.Cat.se"],"name":[30],"type":["lgl"],"align":["right"]},{"label":["w.pa.sampling.OnTrail.se"],"name":[31],"type":["lgl"],"align":["right"]}],"data":[{"1":"NA","2":"NA","3":"NA","4":"NA","5":"NA","6":"NA","7":"NA","8":"NA","9":"NA","10":"NA","11":"NA","12":"NA","13":"NA","14":"NA","15":"NA","16":"NA","17":"NA","18":"NA","19":"NA","20":"NA","21":"NA","22":"NA","23":"NA","24":"NA","25":"NA","26":"NA","27":"NA","28":"NA","29":"NA","30":"NA","31":"NA"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
po.coefs
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["x.distancePA"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["x.EVI"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["x.EVI_sq"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["x.forestloss"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["x.road_length"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["x.fire_density"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["x.mean_slope"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000_EVI"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["alpha0.po"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["w.po.forest"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["w.po.protected"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance_sq"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["w.po.mean_slope"],"name":[24],"type":["dbl"],"align":["right"]}],"data":[{"1":"3.934571","2":"-0.492576","3":"0.3548169","4":"0.1656607","5":"1.502789","6":"-0.2915151","7":"-0.2886869","8":"0.0009452964","9":"-0.9452048","10":"1.270242","11":"-0.1052937","12":"-0.7153668","13":"-0.1469811","14":"-0.09987559","15":"-0.317552","16":"-0.3387835","17":"9.912818","18":"0.1553599","19":"-12.48636","20":"0.280741","21":"1.456716","22":"-0.5252029","23":"0.07684464","24":"-10.29325"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
po.se
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["beta0.se"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality.se"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality.se"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month.se"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip.se"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_temp_seasonality_sq.se"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_precip_seasonality_sq.se"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_max_temp_warmest_month_sq.se"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["x.chelsa_annual_precip_sq.se"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000.se"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["x.distancePA.se"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["x.EVI.se"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["x.EVI_sq.se"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["x.forestloss.se"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["x.road_length.se"],"name":[15],"type":["dbl"],"align":["right"]},{"label":["x.fire_density.se"],"name":[16],"type":["dbl"],"align":["right"]},{"label":["x.mean_slope.se"],"name":[17],"type":["dbl"],"align":["right"]},{"label":["x.treecover2000_EVI.se"],"name":[18],"type":["dbl"],"align":["right"]},{"label":["alpha0.po.se"],"name":[19],"type":["dbl"],"align":["right"]},{"label":["w.po.forest.se"],"name":[20],"type":["dbl"],"align":["right"]},{"label":["w.po.protected.se"],"name":[21],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance.se"],"name":[22],"type":["dbl"],"align":["right"]},{"label":["w.po.road_distance_sq.se"],"name":[23],"type":["dbl"],"align":["right"]},{"label":["w.po.mean_slope.se"],"name":[24],"type":["dbl"],"align":["right"]}],"data":[{"1":"6.248012","2":"0.1128055","3":"0.1788559","4":"0.1705399","5":"0.3078827","6":"0.1079978","7":"0.122022","8":"0.04112033","9":"0.2491692","10":"0.1817512","11":"0.14383","12":"0.1488298","13":"0.1106662","14":"0.06616428","15":"0.2048053","16":"0.2146877","17":"6.312799","18":"0.1591219","19":"6.25752","20":"0.2384649","21":"0.2293832","22":"0.1378663","23":"0.04114941","24":"6.301258"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
