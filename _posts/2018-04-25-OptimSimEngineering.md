---
title: "Optimization for (simulation) engineering"
author: "Y. Richet"
layout: post
date:   2018-04-25 00:00:00 +0100
categories: r jekyll
---


_"Efficient Global Optimization of Expensive Black-Box Functions"_
Jones, Schonlau, Welch, Journal of Global Optimization, December 1998

<hr/>


Engineering objective: optimize $f$ function/simulator, with lowest $f$ evaluations as possible.


![plot of chunk unnamed-chunk-1](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-1-1.png)




## Basic idea

create some $\color{blue}{models\ of\ f}$ based on few evaluations $x={X}$


![plot of chunk unnamed-chunk-2](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-2-1.png)



## (simple) Kriging

$$\color{blue}M(x) = \mathcal{N}(m(x),s^2(x))$$


* $m(x) = C(x)^T C(X)^{-1} f(X)$
* $s^2(x) = c(x) - C(x)^T C(X)^{-1} C(x)$
* $C$ is the covariance kernel  $C(.) = C(X,.)$, $c(.) = C(x,.)$



![plot of chunk unnamed-chunk-3](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-3-1.png)



## Efficient Global Optimization



![plot of chunk unnamed-chunk-4](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-4-1.png)


</div>

Let's define the _Expected Improvement_:

$$\color{violet}{EI}(x) = E[(min \{ f(X) \} - M(x))^+]$$

_which is (also) analytical thanks to $M$ properties..._



$$u(x) = { { min \{ f(X) \} - m(x) } \over { s(x) } }$$



![plot of chunk unnamed-chunk-5](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-5-1.png)



$$\color{violet}{EI}(x) = s(x)\ \big(\ u(x)p_{\mathcal{N}}(u(x)) + d_{\mathcal{N}}(u(x))\ \big)$$

<hr/>

EGO: Maximize $EI(x)$ (*), compute $f$ there, add to $X$, ... Repeat until ...

<hr/>

* __+__ good trade-off between exploration and exploitation
* __+__ requires few evaluations of $f$
* __-__ often lead to add close points to each others ...  
Which is not very comfortable for kriging numerical stability
* __-__ "one step lookahead" (myopic) strategy
* __-__ rely on model suitability to $f$

(*) using standard optimization algorithm: BFGS, PSO, DiRect, ...


## EGO - step 0


![plot of chunk unnamed-chunk-6](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-6-1.png)
![plot of chunk unnamed-chunk-7](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-7-1.png)



## EGO - step 1


![plot of chunk unnamed-chunk-8](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-8-1.png)
![plot of chunk unnamed-chunk-9](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-9-1.png)



## EGO - step 2


![plot of chunk unnamed-chunk-10](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-10-1.png)
![plot of chunk unnamed-chunk-11](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-11-1.png)




## EGO - step 3


![plot of chunk unnamed-chunk-12](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-12-1.png)
![plot of chunk unnamed-chunk-13](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-13-1.png)




## EGO - step 4


![plot of chunk unnamed-chunk-14](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-14-1.png)
![plot of chunk unnamed-chunk-15](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-15-1.png)




## EGO - step 5


![plot of chunk unnamed-chunk-16](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-16-1.png)
![plot of chunk unnamed-chunk-17](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-17-1.png)



## EGO - step 6


![plot of chunk unnamed-chunk-18](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-18-1.png)
![plot of chunk unnamed-chunk-19](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-19-1.png)



## EGO - step 7


![plot of chunk unnamed-chunk-20](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-20-1.png)
![plot of chunk unnamed-chunk-21](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-21-1.png)



## EGO - step 8


![plot of chunk unnamed-chunk-22](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-22-1.png)
![plot of chunk unnamed-chunk-23](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-23-1.png)



## EGO - step 9


![plot of chunk unnamed-chunk-24](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-24-1.png)
![plot of chunk unnamed-chunk-25](/www/figure/source/2018-04-25-OptimSimEngineering/unnamed-chunk-25-1.png)




