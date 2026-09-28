library(nleqslv)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)

# alpha0: type I error rate, 1 - alpha1: power
get_n1_UIP <- function(mu1_star, sigma1, alpha0 = 0.025, alpha1 = 0.2){
  z1 = qnorm(1-alpha0, lower.tail = T)
  z2 = qnorm(1-alpha1, lower.tail = T)
  return( sigma1^2 * ( (z1 + z2)/mu1_star )^2 )
}

target_strapp <- function(n1, a0, n0, sigma0,
                          mu1_star, sigma1, z1, z2){
  mu1_star*(n1 + a0*n0) + sigma1*sqrt(n1 + a0*n0)*z1 + sigma1*sqrt(n1)*z2
}

get_n1_strapp <- function(a0, n0, sigma0, mu1_star, sigma1,
                          xstart = 500,
                          alpha0 = 0.025, alpha1 = 0.2){
  z1 = qnorm(1-alpha0, lower.tail = T)
  z2 = qnorm(1-alpha1, lower.tail = T)
  n1 = nleqslv(xstart, target_strapp, a0 = a0, n0 = n0, sigma0 = sigma0,
          mu1_star = mu1_star, sigma1 = sigma1, z1 = z1, z2 = z2,
          control = list(ftol = 10^(-14)),
          method = "Newton")$x
  
  return(n1)
}

target_pp <- function(n1, a0, n0, sigma0,
                      mu1_star, sigma1, z1, z2){
  mu1_star*(n1 + a0*n0*sigma1/sigma0) + sigma1*sqrt(n1 + a0*n0*(sigma1^2)/(sigma0^2))*z1 + sigma1*sqrt(n1)*z2
}

get_n1_pp <- function(a0, n0, sigma0, mu1_star, sigma1,
                          xstart = 500,
                          alpha0 = 0.025, alpha1 = 0.2){
  z1 = qnorm(1-alpha0, lower.tail = T)
  z2 = qnorm(1-alpha1, lower.tail = T)
  n1 = nleqslv(xstart, target_pp, a0 = a0, n0 = n0, sigma0 = sigma0,
               mu1_star = mu1_star, sigma1 = sigma1, z1 = z1, z2 = z2,
               control = list(ftol = 10^(-14)),
               method = "Newton")$x
  
  return(n1)
}

get_df <- function(a0_val, n0, sigma0,
                    mu1_star, sigma1, 
                    alpha0 = 0.025, alpha1 = 0.2){
  n1_UIP    <- rep(get_n1_UIP(mu1_star = mu1_star, sigma1 = sigma1), length(a0_val))
  n1_strapp <- sapply(a0_val, get_n1_strapp, n0 = n0, sigma0 = sigma0,
                      mu1_star = mu1_star, sigma1 = sigma1)
  n1_pp     <- sapply(a0_val, get_n1_pp, n0 = n0, sigma0 = sigma0,
                      mu1_star = mu1_star, sigma1 = sigma1)
  df <- data.frame(n1 = c(n1_UIP, n1_strapp, n1_pp),
                   a0 = rep(a0_val, 3),
                   prior = rep(c("UIP", "straPP", "PP"), each = length(a0_val)))
  return(df)
  
}

a0_val <- seq(0, 1, length = 1000)

# case (a)
n0 = 85; sigma0 = 3
mu1_star = -0.2; sigma1 = 1
dfa <- get_df(a0_val = a0_val, n0 = n0, sigma0 = sigma0,
              mu1_star = mu1_star, sigma1 = sigma1)
plt_a = dfa %>%
  ggplot( aes(x=a0, y=n1, linetype = prior, shape = prior, color = prior)) +
  geom_line() +
  labs(y = TeX("$n$"), x = TeX("$a_0$"), 
       #subtitle = TeX("$n_0 = 85, \\bar{y}_0=0.6, \\sigma_0 = 3,
      #              \\mu_1^* = -0.2, \\sigma_1=1$"),
       tag = TeX(r"(\textbf{(a)} \overset{\normalsize{$n_0 = 85, \bar{y}_0=-0.6, \sigma_0 = 3$}}{\normalsize{$\mu_1^* = -0.2, \sigma_1=1$}$})"))+
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.key.width= unit(2, 'cm'),
        plot.tag.position = "bottom")

# case (b)
n0 = 50; sigma0 = 1
mu1_star = -0.75; sigma1 = 3
dfb <- get_df(a0_val = a0_val, n0 = n0, sigma0 = sigma0,
              mu1_star = mu1_star, sigma1 = sigma1)
plt_b = dfb %>%
  ggplot( aes(x=a0, y=n1, linetype = prior, shape = prior, color = prior)) +
  geom_line() +
  labs(y = TeX("$n$"), x = TeX("$a_0$"), 
       tag = TeX(r"(\textbf{(b)} \overset{\normalsize{$n_0 = 50, \bar{y}_0=-0.25, \sigma_0 = 1$}}{\normalsize{$\mu_1^* = -0.75, \sigma_1=3$}$})"),
       )+
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.key.width= unit(2, 'cm'),
        #plot.caption = element_text(hjust = 0.1, size =12),
        #plot.caption.position = "panel",
        plot.tag.position = "bottom"
        )

# case (c)
n0 = 125; sigma0 = 1
mu1_star = -0.6; sigma1 = 4
dfc <- get_df(a0_val = a0_val, n0 = n0, sigma0 = sigma0,
              mu1_star = mu1_star, sigma1 = sigma1)
plt_c = dfc %>%
  ggplot( aes(x=a0, y=n1, linetype = prior, shape = prior, color = prior)) +
  geom_line() +
  labs(y = TeX("$n$"), x = TeX("$a_0$"), 
       tag = TeX(r"(\textbf{(c)} \overset{\normalsize{$n_0 = 125, \bar{y}_0=-0.15, \sigma_0 = 1$}}{\normalsize{$\mu_1^* = -0.6, \sigma_1=4$}$})"))+
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.key.width= unit(2, 'cm'),
        plot.tag.position = "bottom")

plts = ggarrange(plt_a, plt_b, plt_c, nrow = 1, 
                 #labels = c("(a)", "(b)", "(c)"),
          common.legend = T, legend = "bottom")
annotate_figure(plts, top = text_grob("Required Sample Sizes for Different Priors", 
                                               face = "bold", size = 14)
)
#width: 1200, height: 480
                