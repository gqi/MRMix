## Individual functions for variance calculation
f1 = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi*(sy2+theta^2*sx2))*exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

f2 = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi*sigma2)*exp(-(betay-theta*betax)^2/(2*sigma2))
}

f = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    pi0*f1(betax, betay, sx2, sy2, theta) + (1-pi0)*f2(betax, betay, theta, sigma2)
}

## Derivative functions to compute the estimating equation
pf2_psigma2 = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*((sigma2)^(-5/2)*(betay-theta*betax)^2-(sigma2)^(-3/2))*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f2_p2sigma2 = function(betax, betay, theta, sigma2){
    1/4/sqrt(2*pi)*(3*sigma2^(-5/2)-6*(sigma2)^(-7/2)*(betay-theta*betax)^2+sigma2^(-9/2)*(betay-theta*betax)^4)*exp(-(betay-theta*betax)^2/2/sigma2)
}

pf1_ptheta = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi)*(-theta*sx2*(sy2+theta^2*sx2)^(-3/2)+betax*(betay-theta*betax)*(sy2+theta^2*sx2)^(-3/2)+(betay-theta*betax)^2*theta*sx2*(sy2+theta^2*sx2)^(-5/2))*exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

pf2_ptheta = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi)*sigma2^(-3/2)*betax*(betay-theta*betax)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f2_ptheta_psigma2 = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*betax*(betay-theta*betax)*(sigma2^(-7/2)*(betay-theta*betax)^2-3*sigma2^(-5/2))*exp(-(betay-theta*betax)^2/2/sigma2)
}

### Further derivatives over theta
p3f2_p2sigma2_ptheta = function(betax, betay, theta, sigma2){
    1/4/sqrt(2*pi)*betax*(betay-theta*betax)*
        (15*sigma2^(-7/2)-10*sigma2^(-9/2)*(betay-theta*betax)^2+sigma2^(-11/2)*(betay-theta*betax)^4)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f1_p2theta = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi)*(-(betax^2+sx2)*(sy2+theta^2*sx2)^(-3/2)-(5*betax*(betay-theta*betax)-3*theta*sx2)*theta*sx2*(sy2+theta^2*sx2)^(-5/2)+
                      sx2*(betay-theta*betax)^2*((sy2+theta^2*sx2)^(-5/2)-5*theta^2*sx2*(sy2+theta^2*sx2)^(-7/2))+
                      (-theta*sx2*(sy2+theta^2*sx2)^(-3/2)+betax*(betay-theta*betax)*(sy2+theta^2*sx2)^(-3/2)+(betay-theta*betax)^2*theta*sx2*(sy2+theta^2*sx2)^(-5/2))*
                      (betax*(betay-theta*betax)/(sy2+theta^2*sx2)+(betay-theta*betax)^2*theta*sx2/(sy2+theta^2*sx2)^2))*exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

p2f2_p2theta = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi)*(-betax^2*sigma2^(-3/2)+betax^2*(betay-theta*betax)^2*sigma2^(-5/2))*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

p3f2_p2theta_psigma2 = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*betax^2*(3*sigma2^(-5/2)-6*sigma2^(-7/2)*(betay-theta*betax)^2+sigma2^(-9/2)*(betay-theta*betax)^4)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

### Derivatives over betax
pf1_pbetax = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi)*(sy2+theta^2*sx2)^(-3/2)*theta*(betay-theta*betax)*
        exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

pf2_pbetax = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi)*sigma2^(-3/2)*theta*(betay-theta*betax)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f2_psigma2_pbetax = function(betax, betay, theta, sigma2){
    theta/2/sqrt(2*pi)*(-3*sigma2^(-5/2)*(betay-theta*betax)+sigma2^(-7/2)*(betay-theta*betax)^3)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

p3f2_p2sigma2_pbetax = function(betax, betay, theta, sigma2){
    theta/4/sqrt(2*pi)*(15*sigma2^(-7/2)*(betay-theta*betax)-10*sigma2^(-9/2)*(betay-theta*betax)^3+
                            sigma2^(-11/2)*(betay-theta*betax)^5)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f1_ptheta_pbetax = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi)*((betay-2*theta*betax)*(sy2+theta^2*sx2)^(-3/2)-
                      2*theta^2*(betay-theta*betax)*sx2*(sy2+theta^2*sx2)^(-5/2)+
                      (-theta*sx2*(sy2+theta^2*sx2)^(-3/2)+betax*(betay-theta*betax)*(sy2+theta^2*sx2)^(-3/2)+(betay-theta*betax)^2*theta*sx2*(sy2+theta^2*sx2)^(-5/2))*
                      theta*(betay-theta*betax)/(sy2+theta^2*sx2))*exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

p2f2_ptheta_pbetax = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi)*sigma2^(-3/2)*(betay-2*theta*betax+
                                    theta/sigma2*betax*(betay-theta*betax)^2)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p3f2_ptheta_psigma2_pbetax = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*((betay-2*theta*betax)*(sigma2^(-7/2)*(betay-theta*betax)^2-3*sigma2^(-5/2))-
                        5*theta*sigma2^(-7/2)*betax*(betay-theta*betax)^2+theta*sigma2^(-9/2)*betax*(betay-theta*betax)^4)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

### Derivatives over betay
pf1_pbetay = function(betax, betay, sx2, sy2, theta){
    -1/sqrt(2*pi)*(sy2+theta^2*sx2)^(-3/2)*(betay-theta*betax)*
        exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

pf2_pbetay = function(betax, betay, theta, sigma2){
    -1/sqrt(2*pi)*sigma2^(-3/2)*(betay-theta*betax)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f2_psigma2_pbetay = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*(3*sigma2^(-5/2)*(betay-theta*betax)-sigma2^(-7/2)*(betay-theta*betax)^3)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

p3f2_p2sigma2_pbetay = function(betax, betay, theta, sigma2){
    1/4/sqrt(2*pi)*(-15*sigma2^(-7/2)*(betay-theta*betax)+10*sigma2^(-9/2)*(betay-theta*betax)^3-
                        sigma2^(-11/2)*(betay-theta*betax)^5)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p2f1_ptheta_pbetay = function(betax, betay, sx2, sy2, theta){
    1/sqrt(2*pi)*(betax*(sy2+theta^2*sx2)^(-3/2)+
                      2*theta*(betay-theta*betax)*sx2*(sy2+theta^2*sx2)^(-5/2)+
                      (-theta*sx2*(sy2+theta^2*sx2)^(-3/2)+betax*(betay-theta*betax)*(sy2+theta^2*sx2)^(-3/2)+(betay-theta*betax)^2*theta*sx2*(sy2+theta^2*sx2)^(-5/2))*
                      (-1)*(betay-theta*betax)/(sy2+theta^2*sx2))*exp(-(betay-theta*betax)^2/2/(sy2+theta^2*sx2))
}

p2f2_ptheta_pbetay = function(betax, betay, theta, sigma2){
    1/sqrt(2*pi)*sigma2^(-3/2)*(betax-
                                    1/sigma2*betax*(betay-theta*betax)^2)*exp(-(betay-theta*betax)^2/2/sigma2)
}

p3f2_ptheta_psigma2_pbetay = function(betax, betay, theta, sigma2){
    1/2/sqrt(2*pi)*(betax*(sigma2^(-7/2)*(betay-theta*betax)^2-3*sigma2^(-5/2))+
                        5*sigma2^(-7/2)*betax*(betay-theta*betax)^2-sigma2^(-9/2)*betax*(betay-theta*betax)^4)*
        exp(-(betay-theta*betax)^2/2/sigma2)
}

## First derivatives of the likelihood function
l = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    log(f(betax, betay, sx2, sy2, theta, pi0, sigma2))
}

pl_ppi0 = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2)) / f(betax, betay, sx2, sy2, theta, pi0, sigma2)
}

pl_psigma2 = function(sigma2){
    (1-pi0)*pf2_psigma2(betax, betay, theta, sigma2) / f(betax, betay, sx2, sy2, theta, pi0, sigma2)
}

## Second derivatives of the likelihood
p2l_p2sigma2 = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)*(p2f2_p2sigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                 pf2_psigma2(betax, betay, theta, sigma2)^2*(1-pi0))/
        f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2
}

p2l_ppi0_ptheta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    ((pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
         (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))/
        f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2
}

p2l_ppi0_psigma2 = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (-pf2_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
         (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2))/
        f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2
}

p2l_psigma2_ptheta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)*(p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                 pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))/
        f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2
}

## Third derivatives of the likelihood - over betax
p3l_p2sigma2_pbetax = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ((p3f2_p2sigma2_pbetax(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+p2f2_p2sigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2))-
              2*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2)*p2f2_psigma2_pbetax(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_p2sigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-pf2_psigma2(betax, betay, theta, sigma2)^2*(1-pi0))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2)))
}

p3l_ppi0_ptheta_pbetax = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    1/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        (((p2f1_ptheta_pbetax(betax, betay, sx2, sy2, theta)-p2f2_ptheta_pbetax(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              (pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2))-
              (pf1_pbetax(betax, betay, sx2, sy2, theta)-pf2_pbetax(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))-
              (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*p2f1_ptheta_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_ptheta_pbetax(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             ((pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2)))
}

p3l_ppi0_psigma2_pbetax = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    -1/(f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4)*
        ((p2f2_psigma2_pbetax(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2))+
              (1-pi0)*(pf1_pbetax(betax, betay, sx2, sy2, theta)-pf2_pbetax(betax, betay, theta, sigma2))*pf2_psigma2(betax, betay, theta, sigma2)+
              (1-pi0)*(f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*p2f2_psigma2_pbetax(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (pf2_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2)))
}

p3l_psigma2_ptheta_pbetax = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ####################
    #########################
        ((p3f2_ptheta_psigma2_pbetax(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2))-
              p2f2_psigma2_pbetax(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))-
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*p2f1_ptheta_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_ptheta_pbetax(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetax(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetax(betax, betay, theta, sigma2)))
}

## Third derivative of the likelihood - over betay
p3l_p2sigma2_pbetay = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ((p3f2_p2sigma2_pbetay(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+p2f2_p2sigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2))-
              2*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2)*p2f2_psigma2_pbetay(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_p2sigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-pf2_psigma2(betax, betay, theta, sigma2)^2*(1-pi0))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2)))
}

p3l_ppi0_ptheta_pbetay = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    1/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        (((p2f1_ptheta_pbetay(betax, betay, sx2, sy2, theta)-p2f2_ptheta_pbetay(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              (pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2))-
              (pf1_pbetay(betax, betay, sx2, sy2, theta)-pf2_pbetay(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))-
              (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*p2f1_ptheta_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_ptheta_pbetay(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             ((pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2)))
}

p3l_ppi0_psigma2_pbetay = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    -1/(f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4)*
        ((p2f2_psigma2_pbetay(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2))+
              (1-pi0)*(pf1_pbetay(betax, betay, sx2, sy2, theta)-pf2_pbetay(betax, betay, theta, sigma2))*pf2_psigma2(betax, betay, theta, sigma2)+
              (1-pi0)*(f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*p2f2_psigma2_pbetay(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (pf2_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2)))
}

p3l_psigma2_ptheta_pbetay = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ((p3f2_ptheta_psigma2_pbetay(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2))-
              p2f2_psigma2_pbetay(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))-
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*p2f1_ptheta_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_ptheta_pbetay(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_pbetay(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_pbetay(betax, betay, theta, sigma2)))
}

## Third derivative of the likelihood - over theta
p3l_p2sigma2_ptheta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ((p3f2_p2sigma2_ptheta(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+p2f2_p2sigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))-
              2*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2)*p2f2_ptheta_psigma2(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_p2sigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-pf2_psigma2(betax, betay, theta, sigma2)^2*(1-pi0))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))
}

p3l_ppi0_p2theta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    1/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        (((p2f1_p2theta(betax, betay, sx2, sy2, theta)-p2f2_p2theta(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
              (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*p2f1_p2theta(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_p2theta(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             ((pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))
}

p3l_ppi0_psigma2_ptheta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    -1/(f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4)*
        ((p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2))+
              (1-pi0)*(pf1_ptheta(betax, betay, sx2, sy2, theta)-pf2_ptheta(betax, betay, theta, sigma2))*pf2_psigma2(betax, betay, theta, sigma2)+
              (1-pi0)*(f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*p2f2_ptheta_psigma2(betax, betay, theta, sigma2))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (pf2_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)+
                  (f1(betax, betay, sx2, sy2, theta)-f2(betax, betay, theta, sigma2))*(1-pi0)*pf2_psigma2(betax, betay, theta, sigma2))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))
}

p3l_psigma2_p2theta = function(betax, betay, sx2, sy2, theta, pi0, sigma2){
    (1-pi0)/f(betax, betay, sx2, sy2, theta, pi0, sigma2)^4*
        ((p3f2_p2theta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
              pf2_psigma2(betax, betay, theta, sigma2)*(pi0*p2f1_p2theta(betax, betay, sx2, sy2, theta)+(1-pi0)*p2f2_p2theta(betax, betay, theta, sigma2)))*f(betax, betay, sx2, sy2, theta, pi0, sigma2)^2-
             (p2f2_ptheta_psigma2(betax, betay, theta, sigma2)*f(betax, betay, sx2, sy2, theta, pi0, sigma2)-
                  pf2_psigma2(betax, betay, theta, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))*
             2*f(betax, betay, sx2, sy2, theta, pi0, sigma2)*(pi0*pf1_ptheta(betax, betay, sx2, sy2, theta)+(1-pi0)*pf2_ptheta(betax, betay, theta, sigma2)))
}
