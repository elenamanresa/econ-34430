
The goal of the following simulation is to develop our understanding of the two-way fixed effect model. See the original paper by [Abowd Kramartz and Margolis](http://onlinelibrary.wiley.com/doi/10.1111/1468-0262.00020/full).


# Constructing Employer-Employee matched data

## Simulating a network

One central piece is to have a network of workers and firms over time. We then start by simulating such an object. The rest of the simulation will focus on adding wages to this model. As we know, a central issue of the network will be the number of movers.

We are going to model the mobility between workers and firms. Given a transition matrix we can solve for a stationary distribution, and then construct our panel from there.


```r
nk = 30
nl = 10


## Calibration to CHK data
# psi_sd = 0.16
# alpha_sd = 0.29
# w_sigma = 0.37

alpha_sd = 1
psi_sd   = 1

# let's draw some FE
psi   = qnorm(1:nk/(nk+1)) * alpha_sd
alpha = qnorm(1:nl/(nl+1)) * psi_sd

# let's assume moving PR is fixed
lambda = 0.05

csort = 0.5 # sorting effect
cnetw = 0.2 # network effect
csig  = 0.5 

# lets create type specific transition matrices
# we are going to use joint normal centered on different values
G = array(0,c(nl,nk,nk)) # def: conditional distribution k,k' given l
for (l in 1:nl) for (k in 1:nk) {
  #G[l,k,] = dnorm(csort*(psi[k] - alpha[l])) * dnorm(cnetw*(psi[k] - psi))
  G[l,k,] = dnorm( psi - cnetw *psi[k] - csort * alpha[l],sd = csig)
  #G[l,k,] = (dnorm(csort*(psi[k] - alpha[l])) + dnorm(cnetw*(psi[k] - psi)))/2
  G[l,k,] = G[l,k,]/sum(G[l,k,])
}


# we then solve for the stationary distribution over psis for each alpha value
H = array(1/nk,c(nl,nk))
for (l in 1:nl) {
 M = G[l,,]
 for (i in 1:100) {
   H[l,] = t(G[l,,]) %*% H[l,]
 }
}

Plot1=wireframe(G[1,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
Plot2=wireframe(G[nl,,],aspect = c(1,1),xlab = "previous firm",ylab="next firm")
grid.arrange(Plot1, Plot2,nrow=1)
```

![](index_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

And we can plot the joint distribution of matches


```r
wireframe(H,aspect = c(1,1),xlab = "worker",ylab="firm")
```

![](index_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The next step is to simulate our network given our transitions rules. 


```r
nt = 5
ni = 100000

# we simulate a panel
network    = array(0,c(ni,nt))
spellcount = array(0,c(ni,nt))
A = rep(0,ni)

for (i in 1:ni) {
  # we draw the worker type
  l = sample.int(nl,1)
  A[i]=l 
  # at time 1, we draw from H
  network[i,1] = sample.int(nk,1,prob = H[l,])
  for (t in 2:nt) {
    if (runif(1)<lambda) { #movers
      network[i,t] = sample.int(nk,1,prob = G[l,network[i,t-1],])
      spellcount[i,t] = spellcount[i,t-1] +1
    } else { # stayers
      network[i,t]    = network[i,t-1]
      spellcount[i,t] = spellcount[i,t-1]
    }
  }
}

data  = data.table(melt(network,c('i','t'))) 
data2 = data.table(melt(spellcount,c('i','t')))
setnames(data,"value","k") # firm type
data <- data[,spell := data2$value] # number of spells
data <- data[,l := A[i],i] # worker type
data <- data[,alpha := alpha[l],l] # worker FE
data <- data[,psi := psi[k],k] # firm FE
```


The final step is assigning identities to firms. We are going to do this by randomly assigning firm ids to spells.


```r
firm_size = 10
within_firm_count = ni/(firm_size*nk*nt)

dspell <- data[,list(len=.N),list(i,spell,k)] #length of the spell in each firm by worker
dspell <- dspell[,fid := sample( 1: pmax(1,sum(len)/(firm_size*nt) ) ,.N,replace=TRUE) , k]
dspell <- dspell[,fid := .GRP, list(k,fid)]

setkey(data,i,spell)
setkey(dspell,i,spell)

data <- data[, fid:= dspell[data,fid]]
```
<span class="label label-success">Task 1</span> Compute:
 - mean firm size, in the crossection, expect something like 15.
 - mean number of movers per firm in total in our panel (a worker that moved from firm i to j is counted as mover in firm i as well as in firm j).


```r
#firm size
firm_size_t = data[,.N,list(fid,t)]
setkey(firm_size_t,fid)
firm_size = firm_size_t[,mean(N),fid]
data <- data[, firm_size:= firm_size[fid,V1]]
Nfirm = firm_size[.N][[1]]

#movers per firm
firm_movers_t = data[,max(spell),i] # we identify the movers
firm_movers_t = firm_movers_t[,mover:=(V1>0),i]
data_aux <- data
data_aux <- data_aux[, mover:= firm_movers_t[i,mover]]
firm_mover <- data_aux[mover==TRUE, unique(i),fid]
firm_mover <- firm_mover[,ones:=1]
firm_mover <- firm_mover[,sum(ones),fid]

data <- data[, movers_size:= firm_mover[fid,V1]]

# PLOTS
hist(data[t==1]$firm_size) # at any t, the firm size distribution is the same
```

![](index_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
hist(data$movers_size)
```

![](index_files/figure-html/unnamed-chunk-6-2.png)<!-- -->


## Simulating AKM wages

We start with just AKM wages, which is log additive with some noise.


```r
w_sigma = 0.8
data <- data[, lw := alpha + psi + w_sigma * rnorm(.N) ]
```

<span class="label label-success">Task 2</span> We use this generated data to create the event study plot from Card-Heining-Kline:

 1. Compute the mean wage within firm
 2. group firms into quartiles
 3. select workers around a move (2 periods pre, 2 periods post)
 4. compute wages before/after the mover for each transitions (from each quartile to each quartile)


```r
data[,mean_w:=mean(lw),fid]
q1 = quantile(data$mean_w,.25)
q2 = quantile(data$mean_w,.5)
q3 = quantile(data$mean_w,.75)

data[,flag_q:=0]
data[mean_w <= q1,flag_q:=1]
data[mean_w > q3,flag_q:=4]
data[(mean_w < q2) & (mean_w>q1) ,flag_q:=2]
data[(mean_w < q3) & (mean_w >q2) ,flag_q:=3]

# select one time movers at time 3:
data[,max_spell:=max(spell),i]

data_aux <- data[max_spell == 1]
data_aux[,flag_ckh1:=0]
data_aux[,flag_ckh2:=0]

data_aux[spell == 1 & t == 3,flag_ckh1:=1]
data_aux[spell == 0 & t == 2,flag_ckh2:=1]

data_aux[,flag_ckh1:=max(flag_ckh1),i]
data_aux[,flag_ckh2:=max(flag_ckh2),i]

data_aux[,flag_ckh3:=flag_ckh1*flag_ckh2]

data_ckh <- data_aux[flag_ckh3==1,]
data_aux <- NULL

# identify the different type of movers q1->q4, q2->q4, q3->q4, q4->q1, etc.
data_ckh[, flag_ini:=0]
data_ckh[, flag_end:=0]

data_ckh[(t<3), flag_ini:=flag_q]
data_ckh[, flag_ini:=max(flag_ini),i]

data_ckh[(t>4), flag_end:=flag_q]
data_ckh[, flag_end:=max(flag_end),i]

# mean of log-wages by movement
data_ckh <- data_ckh[,mean(lw),list(t,flag_ini,flag_end)]
```


# Estimating AKM

This requires to first extract a large connected set, and then to estimate the linear problem with many dummies.

## Extracting the connected set

Because we are not going to deal with extremely large data-sets, we can use of the shelf algorithms to extract the connected set. Use the function `conComp` from the package `ggm` to extract the connected set from our data. To do so you will need to construct first an adgency matrix between the firms. Here is how I would proceed to construct the adjency matrix:



```r
setkey(data,i,t)

adjMat <- function(data){ 
  setkey(data,i,t)
  data[ ,fid.l1 := data[J(i,t-1),fid]]
  data <- data[complete.cases(data)]
  nfirms = data[,max(fid)]
  jdata = data[fid.l1!=fid][,list(fid,fid.l1)]
  jdata = unique(jdata)
  setkey(jdata,fid,fid.l1)
  # jdata holds the indices where a sparse matrix should be 1.
  amat = Matrix::sparseMatrix(i=jdata[,c(fid,fid.l1)],j=jdata[,c(fid.l1,fid)],dims=c(nfirms,nfirms))
 amat
}

concomp <- function(data){
  setkey(data,i,t)
  data[ ,fid.l1 := data[J(i,t-1),fid]]
  data <- data[complete.cases(data)]
  cf = lfe::compfactor(list(f1=data[,factor(fid)],f2=data[,factor(fid.l1)]))
  fr = data.frame(f1=data[,factor(fid)],f2=data[,factor(fid.l1)],cf)
  data = data[fr$cf==1]
  return(data)
}


data <- concomp(data)
# compute the number of firms:
aa = unique(data$fid)
Nfirm_con = length(aa)

print("proportion of firms in the connected set is")
```

```
## [1] "proportion of firms in the connected set is"
```

```r
print(Nfirm_con/Nfirm)
```

```
## [1] 0.9637456
```


## Estimating worker and firms FE

Finally, we are going to implement the AKM estimator. This can be done simply by updating in turn the worker FE and the firms FE. 

Start by appending 2 new columns `alpha_hat` and `psi_hat` to your data. Then loop over the following:

 1. update `alpha_hat` by just taking the mean within `i` net of firm FE, use data.table to do that very efficiently
 2. update `psi_hat`   by just taking the mean within `fid` net of worker FE, use data.table to do that very efficiently
 
 
Estimating this model is non-trivial. We are going to take the approach of [Guimareas and Portugal](https://s3.amazonaws.com/academia.edu.documents/42461471/Variable_selection_in_linear_regression20160209-26324-14r87y1.pdf?AWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1525336798&Signature=I2SdtwWgl0U6e28114dmaiLjlA4%3D&response-content-disposition=inline%3B%20filename%3DVariable_selection_in_linear_regression.pdf#page=130), who propose a *ZigZag* estimator in the Stata Journal. [Simen Gaure](https://www.sciencedirect.com/science/article/pii/S0167947313001266) proposes an almost equivalent approach and develops the [`lfe`](https://cran.r-project.org/web/packages/lfe/lfe.pdf) package for R. The idea in both approaches starts with the formulation

$$
\mathbf{Y} = \mathbf{Z}\beta + \mathbf{D}\alpha + \epsilon
$$

where $\mathbf{Z}$ is $(N,k)$, and $\mathbf{D}$ is a  $(N,G_1)$ matrix of dummy variables: element $(i,j)$ of $\mathbf{D}$ is $1$ if $i$ is associated to $j$.

Both papers show the recursive relationship

$$
\left[ \begin{array}{c}
\beta &=& (\mathbf{Z}'\mathbf{Z})^{-1} \mathbf{Z}'(\mathbf{Y}-\mathbf{D}\alpha) \\
\alpha &=& (\mathbf{D}'\mathbf{D})^{-1} \mathbf{D}'(\mathbf{Y}-\mathbf{Z}\beta)
\end{array}
\right]
$$

* line 2 is just `data[,mean(resid(lm(y ~ z))),by=i]`
* Dimensionality of $\mathbf{D}$ becomes irrelevant.
* Same principle for more than 1 FE:

$$
\mathbf{Y} = \mathbf{Z}\beta + \mathbf{D}_1\alpha + \mathbf{D}_2\gamma + \epsilon
$$

and recursive structure: 

$$
\left[ \begin{array}{c}
\beta &=& (\mathbf{Z}'\mathbf{Z})^{-1} \mathbf{Z}'(\mathbf{Y}-\mathbf{D}_1\alpha-\mathbf{D}_2\gamma) \\
\alpha &=& (\mathbf{D}_1'\mathbf{D}_1)^{-1} \mathbf{D}_1'(\mathbf{Y}-\mathbf{Z}\beta-\mathbf{D}_2\gamma) \\
\gamma &=& (\mathbf{D}_2'\mathbf{D}_2)^{-1} \mathbf{D}_2'(\mathbf{Y}-\mathbf{Z}\beta-\mathbf{D}_1\alpha)
\end{array}
\right]
$$

* First, store a `mod = data[,lm(y ~ z + D1 + D2)]`
* line 2 is just `data[,mean(resid(mod) + coef(mod)[3]*D2),by=i]`
* line 3 is just `data[,mean(resid(mod) + coef(mod)[2]*D1),by=j]`

Guimaraes+Portugal do this in an iterative fashion until the difference in mean squared error (MSE) of successive estimates of line 0 becomes small.

 


```r
guimaraesPortugal <- function(data,tol=1e-6){
  setkey(data,i,t)
  data[,fid.l1 := data[J(i,t-1)][["fid"]]]
  data[,move:=fid!=fid.l1]
  mvid <- data[(move),list(id=unique(i))]
  data[,alpha_hat := 0.0]
  data[,psi_hat := 0.0]
  delta = Inf
  MSE1 = Inf
  i = 0
  while (delta>tol){
    i = i + 1
    mod = data[,lm(lw ~ psi_hat + alpha_hat)]
    coefs = coef(mod)
    if (i==1){
      coefs[2:3] <- 0
    }
    data[,res := resid(mod)]
    MSE0 = data[.(mvid$id),mean(res^2)] #You can increase speed by focusing on movers only first to recover the `psi`.
    delta = abs(MSE0 - MSE1)
    MSE1 = MSE0
    data[,alpha_hat := mean(res + coefs[3]*alpha_hat),by=i ]
    data[,psi_hat := mean(res + coefs[2]*psi_hat),by=fid ]
    if (i==1 || (i %% 10 == 0)){
     # flog.info("iter=%d, MSE=%f, delta(MSE)=%f",i,MSE0,delta)  
    }
  }
  return(data)
}
dd=guimaraesPortugal(data)

# estimated coefs on a simple lm must be one
summary(lm(lw~alpha_hat+psi_hat,data))
```

```
## 
## Call:
## lm(formula = lw ~ alpha_hat + psi_hat, data = data)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -3.2183 -0.4546  0.0001  0.4552  3.2072 
## 
## Coefficients:
##              Estimate Std. Error  t value Pr(>|t|)    
## (Intercept) 0.0061970  0.0010926    5.672 1.42e-08 ***
## alpha_hat   1.0000516  0.0009803 1020.179  < 2e-16 ***
## psi_hat     1.0001096  0.0012066  828.862  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.6815 on 388974 degrees of freedom
## Multiple R-squared:  0.7897,	Adjusted R-squared:  0.7897 
## F-statistic: 7.305e+05 on 2 and 388974 DF,  p-value: < 2.2e-16
```

## Computing the variance decomposition

We now compute the variance decomposition and compare it with its infeasible counterpart, using alpha and psi, instead of $\hat{\alpha}$ and $\hat{\psi}$.

The variance decomposition is: $(var(\hat{\alpha}),var(\hat{\psi}),cov(\hat{\alpha},\hat{\psi}))$:

```r
c(var(data$alpha_hat),var(data$psi_hat),cov(data$alpha_hat,data$psi_hat))
```

```
## [1]  1.2881528  0.8502136 -0.1972321
```


The true variance decomposition is: $(var(\alpha),var(\psi),cov(\alpha,\psi))$

```r
c(var(data$alpha),var(data$psi),cov(data$alpha,data$psi))
```

```
## [1] 0.6202044 0.3470182 0.2961474
```
