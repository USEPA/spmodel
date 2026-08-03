# splm print output is unchanged (point-referenced, default)

    Code
      print(spmod)
    Output
      
      Call:
      splm(formula = z ~ water + tarp, data = caribou, spcov_type = "exponential", 
          xcoord = x, ycoord = y)
      
      
      Coefficients (fixed):
      (Intercept)       waterY     tarpnone    tarpshade  
          2.05021     -0.08336      0.08006      0.28663  
      
      
      Coefficients (exponential spatial covariance):
            de        ie     range  
       0.10672   0.02244  18.01186  
      

---

    Code
      print(summary(spmod))
    Output
      
      Call:
      splm(formula = z ~ water + tarp, data = caribou, spcov_type = "exponential", 
          xcoord = x, ycoord = y)
      
      Residuals:
           Min       1Q   Median       3Q      Max 
      -0.41321 -0.20784 -0.11238  0.02915  0.45415 
      
      Coefficients (fixed):
                  Estimate Std. Error z value Pr(>|z|)    
      (Intercept)  2.05021    0.30373   6.750 1.48e-11 ***
      waterY      -0.08336    0.06443  -1.294 0.195745    
      tarpnone     0.08006    0.07750   1.033 0.301564    
      tarpshade    0.28663    0.07657   3.743 0.000181 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Pseudo R-squared: 0.3972
      
      Coefficients (exponential spatial covariance):
            de       ie    range 
       0.10672  0.02244 18.01186 

# splm print output is unchanged (anisotropy shown)

    Code
      print(spmod)
    Output
      
      Call:
      splm(formula = z ~ water + tarp, data = caribou, spcov_type = "exponential", 
          xcoord = x, ycoord = y, anisotropy = TRUE)
      
      
      Coefficients (fixed):
      (Intercept)       waterY     tarpnone    tarpshade  
          2.05728     -0.12276      0.06734      0.26739  
      
      
      Coefficients (exponential spatial covariance):
             de         ie      range     rotate      scale  
      6.041e-02  6.041e-06  4.487e+00  1.081e+00  2.425e-01  
      

# splm print output is unchanged (none/ie collapse)

    Code
      print(spmod)
    Output
      
      Call:
      splm(formula = z ~ water + tarp, data = caribou, spcov_type = "none")
      
      
      Coefficients (fixed):
      (Intercept)       waterY     tarpnone    tarpshade  
           1.9600      -0.0352       0.0497       0.2509  
      
      
      Coefficients (none spatial covariance):
           ie  
      0.03635  
      

# spautor print output is unchanged (de/range only)

    Code
      print(spmod)
    Output
      
      Call:
      spautor(formula = y ~ x, data = exdata_poly, spcov_type = "car", 
          estmethod = "reml")
      
      
      Coefficients (fixed):
      (Intercept)            x  
          -0.1579      -0.1314  
      
      
      Coefficients (car spatial covariance):
          de   range  
      3.8458  0.1481  
      

---

    Code
      print(summary(spmod))
    Output
      
      Call:
      spautor(formula = y ~ x, data = exdata_poly, spcov_type = "car", 
          estmethod = "reml")
      
      Residuals:
          Min      1Q  Median      3Q     Max 
      -1.5263 -0.7825 -0.0775  0.6407  2.2922 
      
      Coefficients (fixed):
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept)  -0.1579     0.1459  -1.082    0.279
      x            -0.1314     0.1237  -1.062    0.288
      
      Pseudo R-squared: 0.02343
      
      Coefficients (car spatial covariance):
          de  range 
      3.8458 0.1481 

# spautor print output is unchanged (de/ie/range/extra all present)

    Code
      print(spmod)
    Output
      
      Call:
      spautor(formula = y ~ x, data = exdata_Upoly, spcov_type = "car", 
          estmethod = "reml")
      
      
      Coefficients (fixed):
      (Intercept)            x  
          -0.0900      -0.1035  
      
      
      Coefficients (car spatial covariance):
          de   range   extra  
      4.2151  0.3173  0.5812  
      

---

    Code
      print(summary(spmod))
    Output
      
      Call:
      spautor(formula = y ~ x, data = exdata_Upoly, spcov_type = "car", 
          estmethod = "reml")
      
      Residuals:
           Min       1Q   Median       3Q      Max 
      -1.57759 -0.85906  0.02955  0.63673  2.21781 
      
      Coefficients (fixed):
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept)  -0.0900     0.1746  -0.515    0.606
      x            -0.1035     0.1375  -0.752    0.452
      
      Pseudo R-squared: 0.01316
      
      Coefficients (car spatial covariance):
          de  range  extra 
      4.2151 0.3173 0.5812 

# spglm print output is unchanged

    Code
      print(spmod)
    Output
      
      Call:
      spglm(formula = y_pois ~ x, family = "poisson", data = exdata, 
          spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
      
      
      Coefficients (fixed):
      (Intercept)            x  
           0.7371       0.0630  
      
      
      Coefficients (exponential spatial covariance):
          de      ie   range  
      0.2708  0.0001  0.2109  
      
      
      Coefficients (Dispersion for poisson family):
      dispersion  
               1  
      

---

    Code
      print(summary(spmod))
    Output
      
      Call:
      spglm(formula = y_pois ~ x, family = "poisson", data = exdata, 
          spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord)
      
      Deviance Residuals:
           Min       1Q   Median       3Q      Max 
      -1.83138 -0.67823 -0.09114  0.42672  1.34076 
      
      Coefficients (fixed):
                  Estimate Std. Error z value Pr(>|z|)    
      (Intercept)  0.73708    0.10883   6.773 1.26e-11 ***
      x            0.06300    0.08243   0.764    0.445    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Pseudo R-squared: 0.007079
      
      Coefficients (exponential spatial covariance):
          de     ie  range 
      0.2708 0.0001 0.2109 
      
      Coefficients (Dispersion for poisson family):
      dispersion 
               1 

# spgautor print output is unchanged (de/ie/range/extra all present)

    Code
      print(spmod)
    Output
      
      Call:
      spgautor(formula = abs(y) ~ x, family = "Gamma", data = exdata_Upoly, 
          spcov_type = "car", estmethod = "reml")
      
      
      Coefficients (fixed):
      (Intercept)            x  
         -0.26235      0.01073  
      
      
      Coefficients (car spatial covariance):
             de      range      extra  
      9.663e-03  8.738e-01  9.202e-06  
      
      
      Coefficients (Dispersion for Gamma family):
      dispersion  
           1.637  
      

---

    Code
      print(summary(spmod))
    Output
      
      Call:
      spgautor(formula = abs(y) ~ x, family = "Gamma", data = exdata_Upoly, 
          spcov_type = "car", estmethod = "reml")
      
      Deviance Residuals:
          Min      1Q  Median      3Q     Max 
      -1.9559 -0.7188 -0.1259  0.4029  1.1949 
      
      Coefficients (fixed):
                  Estimate Std. Error z value Pr(>|z|)  
      (Intercept) -0.26235    0.11952  -2.195   0.0282 *
      x            0.01073    0.10831   0.099   0.9211  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Pseudo R-squared: 0.0001923
      
      Coefficients (car spatial covariance):
             de     range     extra 
      9.663e-03 8.738e-01 9.202e-06 
      
      Coefficients (Dispersion for Gamma family):
      dispersion 
           1.637 

