library(ggplot2) # for graph
require(deSolve) # for the "ode" function
library(dplyr)
#library(bbmle) # for ML method  install.packages("bbmle", repos="http://R-Forge.R-project.org")



#Equation#######################
#
# Susceptible
# Infectious_gamma
# Hospitalized_hospL
# Critical (severe)_deltaL (severe length)
# Dead_zeta (prop)
# Recovered
# Vaccinated
#
# adult/young1(20-30)/young2(40-50)/child-nonvac/insufficientvac/fullvac-(original/variant)-(hoMe-hospitalized)
#
# vaccine---reduce infection and hospitalization (eventually severity and death)
# variant---increase transmission and severity (eventually death, but not hospitalization)
#

sir_1 <- function(rt_ini, rt_rebLength, rt_reb, rt_var, soe_eff, rt_age, rt_same, 
                  soe_time1, soe_time2, soe_time3, soe_time4, soe_time5, soe_time6,
                  gamma, deltaC_temp, deltaY1_temp, deltaY2_temp, deltaA_temp, deltaL, del_var,
                  hosp_rateC, hosp_rateY, hosp_rateA, hospL_temp,
                  zetaC_temp, zetaY1_temp, zetaY2_temp, zetaA_temp,
                  vaceff_inf, vaceff_death_temp, vac_brk,
                  
                  popA, popY1, popY2, popC,
                  
                  
                  San0, Sai0, Saf0,
                  Sy1n0, Sy1i0, Sy1f0,
                  Sy2n0, Sy2i0, Sy2f0,
                  Scn0, Sci0, Scf0,
                  
                  Ianom0, Ianoh0, Ianvm0, Ianvh0, Iaiom0, Iaioh0, Iaivm0, Iaivh0, Iafom0, Iafoh0, Iafvm0, Iafvh0,
                  Iy1nom0, Iy1noh0, Iy1nvm0, Iy1nvh0, Iy1iom0, Iy1ioh0, Iy1ivm0, Iy1ivh0, Iy1fom0, Iy1foh0, Iy1fvm0, Iy1fvh0,
                  Iy2nom0, Iy2noh0, Iy2nvm0, Iy2nvh0, Iy2iom0, Iy2ioh0, Iy2ivm0, Iy2ivh0, Iy2fom0, Iy2foh0, Iy2fvm0, Iy2fvh0, 
                  Icnom0, Icnoh0, Icnvm0, Icnvh0, Iciom0, Icioh0, Icivm0, Icivh0, Icfom0, Icfoh0, Icfvm0, Icfvh0,
                  
                  Hano0, Hanv0, Haio0, Haiv0, Hafo0, Hafv0,
                  Hy1no0, Hy1nv0, Hy1io0, Hy1iv0, Hy1fo0, Hy1fv0,
                  Hy2no0, Hy2nv0, Hy2io0, Hy2iv0, Hy2fo0, Hy2fv0,
                  Hcno0, Hcnv0, Hcio0, Hciv0, Hcfo0, Hcfv0,
                  
                  Cano0, Canv0, Caio0, Caiv0, Cafo0, Cafv0,
                  Cy1no0, Cy1nv0, Cy1io0, Cy1iv0, Cy1fo0, Cy1fv0,
                  Cy2no0, Cy2nv0, Cy2io0, Cy2iv0, Cy2fo0, Cy2fv0,
                  Ccno0, Ccnv0, Ccio0, Cciv0, Ccfo0, Ccfv0,
                  
                  Dano0, Danv0, Daio0, Daiv0, Dafo0, Dafv0,
                  Dy1no0, Dy1nv0, Dy1io0, Dy1iv0, Dy1fo0, Dy1fv0,
                  Dy2no0, Dy2nv0, Dy2io0, Dy2iv0, Dy2fo0, Dy2fv0,
                  Dcno0, Dcnv0, Dcio0, Dciv0, Dcfo0, Dcfv0,
                  
                  Rano0, Ranv0, Raio0, Raiv0, Rafo0, Rafv0,
                  Ry1no0, Ry1nv0, Ry1io0, Ry1iv0, Ry1fo0, Ry1fv0,
                  Ry2no0, Ry2nv0, Ry2io0, Ry2iv0, Ry2fo0, Ry2fv0,
                  Rcno0, Rcnv0, Rcio0, Rciv0, Rcfo0, Rcfv0,
                  
                  Vaf0,
                  Vy1f0,
                  Vy2f0,
                  Vcf0,
                  
                  times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      jinryu_to_rt <- 0.025

      rt <- rt_ini
      

      #soe_time
      if (soe_time1 < time) {
        if (time < soe_time1+40) {
          rt <- rt * soe_eff
        }
      }
        
      if (soe_time2 < time) {
        if (time < soe_time2+40) {
          rt <- rt * soe_eff
        }
      }
      
      if (soe_time3 < time) {
        if (time < soe_time3+40) {
          rt <- rt * soe_eff
        }
      }

      if (soe_time4 < time) {
        if (time < soe_time4+40) {
          rt <- rt * soe_eff
        }
      }
      
      if (soe_time5 < time) {
        if (time < soe_time5+40) {
          rt <- rt * soe_eff
        }
      }
      
      if (soe_time6 < time) {
        if (time < soe_time6+40) {
          rt <- rt * soe_eff
        }
      }
      
      
      
            

      #var
      rt_varnow <- rt * rt_var
      
      
      #vaccines

      adult_vaccination <- 0
      young_vaccination <- 0
              

      vaccine_delay <- 1/28
      
      
      
      #parameter adjustment
      
      vaceff_death <- (1-vaceff_death_temp) / (1-vaceff_inf)
      
      hospL <- 1 / (hospL_temp - 1/gamma)
      
      deltaA <- deltaA_temp / hosp_rateA
      deltaY1 <- deltaY1_temp / hosp_rateY
      deltaY2 <- deltaY2_temp / hosp_rateY
      deltaC <- deltaC_temp / hosp_rateC
      
      zetaA <- zetaA_temp / deltaA_temp
      zetaY1 <- zetaY1_temp / deltaY1_temp
      zetaY2 <- zetaY2_temp / deltaY2_temp
      zetaC <- zetaC_temp / deltaC_temp
      

      

            

      
      
      infa     <- Ianom+Ianoh+Iaiom+Iaioh+Iafom*(1-vac_brk)+Iafoh*(1-vac_brk)
      infa_var <- Ianvm+Ianvh+Iaivm+Iaivh+Iafvm*(1-vac_brk)+Iafvh*(1-vac_brk)
      
      infc      <- Icnom+Icnoh+Iciom+Icioh+Icfom*(1-vac_brk)+Icfoh*(1-vac_brk)
      infc_var  <- Icnvm+Icnvh+Icivm+Icivh+Icfvm*(1-vac_brk)+Icfvh*(1-vac_brk)
      
      infy1     <- (Iy1nom+Iy1noh+Iy1iom+Iy1ioh+Iy1fom*(1-vac_brk)+Iy1foh*(1-vac_brk))*rt_age
      infy1_var <- (Iy1nvm+Iy1nvh+Iy1ivm+Iy1ivh+Iy1fvm*(1-vac_brk)+Iy1fvh*(1-vac_brk))*rt_age
      
      infy2     <- (Iy2nom+Iy2noh+Iy2iom+Iy2ioh+Iy2fom*(1-vac_brk)+Iy2foh*(1-vac_brk))*rt_age
      infy2_var <- (Iy2nvm+Iy2nvh+Iy2ivm+Iy2ivh+Iy2fvm*(1-vac_brk)+Iy2fvh*(1-vac_brk))*rt_age
      

      
      
      
      
      
      
      
      
      dSan <- -rt / popA * gamma * (infy1+infy2) * San -rt*rt_same / popA * gamma * infa * San -rt / popA * gamma * infc * San -rt_varnow / popA * gamma * (infy1_var+infy2_var) * San -rt_varnow*rt_same / popA * gamma * infa_var * San  -rt_varnow / popA * gamma * infc_var * San - adult_vaccination
      dSai <- -rt / popA * gamma * (infy1+infy2) * Sai -rt*rt_same / popA * gamma * infa * Sai -rt / popA * gamma * infc * Sai -rt_varnow / popA * gamma * (infy1_var+infy2_var) * Sai -rt_varnow*rt_same / popA * gamma * infa_var * Sai  -rt_varnow / popA * gamma * infc_var * Sai + adult_vaccination - Sai * vaccine_delay
      dSaf <- -rt / popA * gamma * (infy1+infy2) * Saf -rt*rt_same / popA * gamma * infa * Saf -rt / popA * gamma * infc * Saf -rt_varnow / popA * gamma * (infy1_var+infy2_var) * Saf -rt_varnow*rt_same / popA * gamma * infa_var * Saf  -rt_varnow / popA * gamma * infc_var * Saf                     + Sai * vaccine_delay * (1-vaceff_inf)

      dSy1n <- -rt*rt_same / popY1 * gamma * infy1 * Sy1n -rt / popY1 * gamma * infy2 * Sy1n -rt / popY1 * gamma * (infa+infc) * Sy1n -rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1n -rt_varnow / popY1 * gamma * infy2_var * Sy1n -rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1n - young_vaccination
      dSy1i <- -rt*rt_same / popY1 * gamma * infy1 * Sy1i -rt / popY1 * gamma * infy2 * Sy1i -rt / popY1 * gamma * (infa+infc) * Sy1i -rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1i -rt_varnow / popY1 * gamma * infy2_var * Sy1i -rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1i + young_vaccination - Sy1i * vaccine_delay
      dSy1f <- -rt*rt_same / popY1 * gamma * infy1 * Sy1f -rt / popY1 * gamma * infy2 * Sy1f -rt / popY1 * gamma * (infa+infc) * Sy1f -rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1f -rt_varnow / popY1 * gamma * infy2_var * Sy1f -rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1f                     + Sy1i * vaccine_delay * (1-vaceff_inf)
      
      dSy2n <- -rt / popY2 * gamma * infy1 * Sy2n -rt*rt_same / popY2 * gamma * infy2 * Sy2n -rt / popY2 * gamma * (infa+infc) * Sy2n -rt_varnow / popY2 * gamma * infy1_var * Sy2n -rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2n -rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2n - young_vaccination
      dSy2i <- -rt / popY2 * gamma * infy1 * Sy2i -rt*rt_same / popY2 * gamma * infy2 * Sy2i -rt / popY2 * gamma * (infa+infc) * Sy2i -rt_varnow / popY2 * gamma * infy1_var * Sy2i -rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2i -rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2i + young_vaccination - Sy2i * vaccine_delay
      dSy2f <- -rt / popY2 * gamma * infy1 * Sy2f -rt*rt_same / popY2 * gamma * infy2 * Sy2f -rt / popY2 * gamma * (infa+infc) * Sy2f -rt_varnow / popY2 * gamma * infy1_var * Sy2f -rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2f -rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2f                     + Sy2i * vaccine_delay * (1-vaceff_inf)
      
      dScn <- -rt / popC * gamma * (infy1+infy2) * Scn -rt / popC * gamma * infa * Scn -rt*rt_same / popC * gamma * infc * Scn -rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scn -rt_varnow / popC * gamma * infa_var * Scn -rt_varnow*rt_same / popC * gamma * infc_var * Scn 
      dSci <- -rt / popC * gamma * (infy1+infy2) * Sci -rt / popC * gamma * infa * Sci -rt*rt_same / popC * gamma * infc * Sci -rt_varnow / popC * gamma * (infy1_var+infy2_var) * Sci -rt_varnow / popC * gamma * infa_var * Sci -rt_varnow*rt_same / popC * gamma * infc_var * Sci - Sci * vaccine_delay
      dScf <- -rt / popC * gamma * (infy1+infy2) * Scf -rt / popC * gamma * infa * Scf -rt*rt_same / popC * gamma * infc * Scf -rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scf -rt_varnow / popC * gamma * infa_var * Scf -rt_varnow*rt_same / popC * gamma * infc_var * Scf + Sci * vaccine_delay * (1-vaceff_inf)

      
      
      
      
      dIanom <- (rt / popA * gamma * (infy1+infy2) * San                +rt*rt_same / popA * gamma * infa * San             +rt / popA * gamma * infc * San)            * (1 - hosp_rateA)  - gamma * Ianom
      dIanoh <- (rt / popA * gamma * (infy1+infy2) * San                +rt*rt_same / popA * gamma * infa * San             +rt / popA * gamma * infc * San)            * hosp_rateA        - gamma * Ianoh
      dIanvm <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * San +rt_varnow*rt_same / popA * gamma * infa_var * San  +rt_varnow / popA * gamma * infc_var * San) * (1 - hosp_rateA)  - gamma * Ianvm
      dIanvh <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * San +rt_varnow*rt_same / popA * gamma * infa_var * San  +rt_varnow / popA * gamma * infc_var * San) * hosp_rateA        - gamma * Ianvh
      
      dIaiom <- (rt / popA * gamma * (infy1+infy2) * Sai                +rt*rt_same / popA * gamma * infa * Sai             +rt / popA * gamma * infc * Sai)            * (1 - hosp_rateA)  - gamma * Iaiom
      dIaioh <- (rt / popA * gamma * (infy1+infy2) * Sai                +rt*rt_same / popA * gamma * infa * Sai             +rt / popA * gamma * infc * Sai)            * hosp_rateA        - gamma * Iaioh
      dIaivm <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * Sai +rt_varnow*rt_same / popA * gamma * infa_var * Sai  +rt_varnow / popA * gamma * infc_var * Sai) * (1 - hosp_rateA)  - gamma * Iaivm
      dIaivh <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * Sai +rt_varnow*rt_same / popA * gamma * infa_var * Sai  +rt_varnow / popA * gamma * infc_var * Sai) * hosp_rateA        - gamma * Iaivh
      
      dIafom <- (rt / popA * gamma * (infy1+infy2) * Saf                +rt*rt_same / popA * gamma * infa * Saf             +rt / popA * gamma * infc * Saf)            * (1 - hosp_rateA*vaceff_death)  - gamma * Iafom
      dIafoh <- (rt / popA * gamma * (infy1+infy2) * Saf                +rt*rt_same / popA * gamma * infa * Saf             +rt / popA * gamma * infc * Saf)            * hosp_rateA*vaceff_death        - gamma * Iafoh
      dIafvm <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * Saf +rt_varnow*rt_same / popA * gamma * infa_var * Saf  +rt_varnow / popA * gamma * infc_var * Saf) * (1 - hosp_rateA*vaceff_death)  - gamma * Iafvm
      dIafvh <- (rt_varnow / popA * gamma * (infy1_var+infy2_var) * Saf +rt_varnow*rt_same / popA * gamma * infa_var * Saf  +rt_varnow / popA * gamma * infc_var * Saf) * hosp_rateA*vaceff_death        - gamma * Iafvh
      
      
      
      dIy1nom <- (rt*rt_same / popY1 * gamma * infy1 * Sy1n            +rt / popY1 * gamma * infy2 * Sy1n            +rt / popY1 * gamma * (infa+infc) * Sy1n)                * (1 - hosp_rateY)  - gamma * Iy1nom
      dIy1noh <- (rt*rt_same / popY1 * gamma * infy1 * Sy1n            +rt / popY1 * gamma * infy2 * Sy1n            +rt / popY1 * gamma * (infa+infc) * Sy1n)                * hosp_rateY        - gamma * Iy1noh
      dIy1nvm <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1n +rt_varnow / popY1 * gamma * infy2_var * Sy1n +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1n) * (1 - hosp_rateY)  - gamma * Iy1nvm
      dIy1nvh <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1n +rt_varnow / popY1 * gamma * infy2_var * Sy1n +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1n) * hosp_rateY        - gamma * Iy1nvh
      
      dIy1iom <- (rt*rt_same / popY1 * gamma * infy1 * Sy1i            +rt / popY1 * gamma * infy2 * Sy1i            +rt / popY1 * gamma * (infa+infc) * Sy1i)                * (1 - hosp_rateY)  - gamma * Iy1iom
      dIy1ioh <- (rt*rt_same / popY1 * gamma * infy1 * Sy1i            +rt / popY1 * gamma * infy2 * Sy1i            +rt / popY1 * gamma * (infa+infc) * Sy1i)                * hosp_rateY        - gamma * Iy1ioh
      dIy1ivm <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1i +rt_varnow / popY1 * gamma * infy2_var * Sy1i +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1i) * (1 - hosp_rateY)  - gamma * Iy1ivm
      dIy1ivh <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1i +rt_varnow / popY1 * gamma * infy2_var * Sy1i +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1i) * hosp_rateY        - gamma * Iy1ivh
      
      dIy1fom <- (rt*rt_same / popY1 * gamma * infy1 * Sy1f            +rt / popY1 * gamma * infy2 * Sy1f            +rt / popY1 * gamma * (infa+infc) * Sy1f)                * (1 - hosp_rateY*vaceff_death)  - gamma * Iy1fom
      dIy1foh <- (rt*rt_same / popY1 * gamma * infy1 * Sy1f            +rt / popY1 * gamma * infy2 * Sy1f            +rt / popY1 * gamma * (infa+infc) * Sy1f)                * hosp_rateY*vaceff_death        - gamma * Iy1foh
      dIy1fvm <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1f +rt_varnow / popY1 * gamma * infy2_var * Sy1f +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1f) * (1 - hosp_rateY*vaceff_death)  - gamma * Iy1fvm
      dIy1fvh <- (rt_varnow*rt_same / popY1 * gamma * infy1_var * Sy1f +rt_varnow / popY1 * gamma * infy2_var * Sy1f +rt_varnow / popY1 * gamma * (infa_var+infc_var) * Sy1f) * hosp_rateY*vaceff_death        - gamma * Iy1fvh
      
      
      
      dIy2nom <- (rt / popY2 * gamma * infy1 * Sy2n            +rt*rt_same / popY2 * gamma * infy2 * Sy2n            +rt / popY2 * gamma * (infa+infc) * Sy2n)                * (1 - hosp_rateY)  - gamma * Iy2nom
      dIy2noh <- (rt / popY2 * gamma * infy1 * Sy2n            +rt*rt_same / popY2 * gamma * infy2 * Sy2n            +rt / popY2 * gamma * (infa+infc) * Sy2n)                * hosp_rateY        - gamma * Iy2noh
      dIy2nvm <- (rt_varnow / popY2 * gamma * infy1_var * Sy2n +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2n +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2n) * (1 - hosp_rateY)  - gamma * Iy2nvm
      dIy2nvh <- (rt_varnow / popY2 * gamma * infy1_var * Sy2n +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2n +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2n) * hosp_rateY        - gamma * Iy2nvh
      
      dIy2iom <- (rt / popY2 * gamma * infy1 * Sy2i            +rt*rt_same / popY2 * gamma * infy2 * Sy2i            +rt / popY2 * gamma * (infa+infc) * Sy2i)                       * (1 - hosp_rateY)  - gamma * Iy2iom
      dIy2ioh <- (rt / popY2 * gamma * infy1 * Sy2i            +rt*rt_same / popY2 * gamma * infy2 * Sy2i            +rt / popY2 * gamma * (infa+infc) * Sy2i)                       * hosp_rateY        - gamma * Iy2ioh
      dIy2ivm <- (rt_varnow / popY2 * gamma * infy1_var * Sy2i +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2i +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2i) * (1 - hosp_rateY)  - gamma * Iy2ivm
      dIy2ivh <- (rt_varnow / popY2 * gamma * infy1_var * Sy2i +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2i +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2i) * hosp_rateY        - gamma * Iy2ivh
      
      dIy2fom <- (rt / popY2 * gamma * infy1 * Sy2f            +rt*rt_same / popY2 * gamma * infy2 * Sy2f            +rt / popY2 * gamma * (infa+infc) * Sy2f)                * (1 - hosp_rateY*vaceff_death)  - gamma * Iy2fom
      dIy2foh <- (rt / popY2 * gamma * infy1 * Sy2f            +rt*rt_same / popY2 * gamma * infy2 * Sy2f            +rt / popY2 * gamma * (infa+infc) * Sy2f)                * hosp_rateY*vaceff_death        - gamma * Iy2foh
      dIy2fvm <- (rt_varnow / popY2 * gamma * infy1_var * Sy2f +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2f +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2f) * (1 - hosp_rateY*vaceff_death)  - gamma * Iy2fvm
      dIy2fvh <- (rt_varnow / popY2 * gamma * infy1_var * Sy2f +rt_varnow*rt_same / popY2 * gamma * infy2_var * Sy2f +rt_varnow / popY2 * gamma * (infa_var+infc_var) * Sy2f) * hosp_rateY*vaceff_death        - gamma * Iy2fvh

      
      
      dIcnom <- (rt / popC * gamma * (infy1+infy2) * Scn                +rt / popC * gamma * infa * Scn            +rt*rt_same / popC * gamma * infc * Scn)            * (1 - hosp_rateC)  - gamma * Icnom
      dIcnoh <- (rt / popC * gamma * (infy1+infy2) * Scn                +rt / popC * gamma * infa * Scn            +rt*rt_same / popC * gamma * infc * Scn)            * hosp_rateC        - gamma * Icnoh
      dIcnvm <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scn +rt_varnow / popC * gamma * infa_var * Scn +rt_varnow*rt_same / popC * gamma * infc_var * Scn) * (1 - hosp_rateC)  - gamma * Icnvm
      dIcnvh <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scn +rt_varnow / popC * gamma * infa_var * Scn +rt_varnow*rt_same / popC * gamma * infc_var * Scn) * hosp_rateC        - gamma * Icnvh
      
      dIciom <- (rt / popC * gamma * (infy1+infy2) * Sci                +rt / popC * gamma * infa * Sci            +rt*rt_same / popC * gamma * infc * Sci)                       * (1 - hosp_rateC)  - gamma * Iciom
      dIcioh <- (rt / popC * gamma * (infy1+infy2) * Sci                +rt / popC * gamma * infa * Sci            +rt*rt_same / popC * gamma * infc * Sci)                       * hosp_rateC        - gamma * Icioh
      dIcivm <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Sci +rt_varnow / popC * gamma * infa_var * Sci +rt_varnow*rt_same / popC * gamma * infc_var * Sci) * (1 - hosp_rateC)  - gamma * Icivm
      dIcivh <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Sci +rt_varnow / popC * gamma * infa_var * Sci +rt_varnow*rt_same / popC * gamma * infc_var * Sci) * hosp_rateC        - gamma * Icivh
      
      dIcfom <- (rt / popC * gamma * (infy1+infy2) * Scf                +rt / popC * gamma * infa * Scf            +rt*rt_same / popC * gamma * infc * Scf)                       * (1 - hosp_rateC*vaceff_death)  - gamma * Icfom
      dIcfoh <- (rt / popC * gamma * (infy1+infy2) * Scf                +rt / popC * gamma * infa * Scf            +rt*rt_same / popC * gamma * infc * Scf)                       * hosp_rateC*vaceff_death        - gamma * Icfoh
      dIcfvm <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scf +rt_varnow / popC * gamma * infa_var * Scf +rt_varnow*rt_same / popC * gamma * infc_var * Scf) * (1 - hosp_rateC*vaceff_death)  - gamma * Icfvm
      dIcfvh <- (rt_varnow / popC * gamma * (infy1_var+infy2_var) * Scf +rt_varnow / popC * gamma * infa_var * Scf +rt_varnow*rt_same / popC * gamma * infc_var * Scf) * hosp_rateC*vaceff_death        - gamma * Icfvh
      
      
      
      
      
      dHano <- gamma * Ianoh * (1 - deltaA)         - hospL * Hano
      dHanv <- gamma * Ianvh * (1 - deltaA*del_var) - hospL * Hanv

      dHaio <- gamma * Iaioh * (1 - deltaA)         - hospL * Haio
      dHaiv <- gamma * Iaivh * (1 - deltaA*del_var) - hospL * Haiv
      
      dHafo <- gamma * Iafoh * (1 - deltaA)         - hospL * Hafo
      dHafv <- gamma * Iafvh * (1 - deltaA*del_var) - hospL * Hafv
      
      
      
      dHy1no <- gamma * Iy1noh * (1 - deltaY1)         - hospL * Hy1no
      dHy1nv <- gamma * Iy1nvh * (1 - deltaY1*del_var) - hospL * Hy1nv
      
      dHy1io <- gamma * Iy1ioh * (1 - deltaY1)         - hospL * Hy1io
      dHy1iv <- gamma * Iy1ivh * (1 - deltaY1*del_var) - hospL * Hy1iv
      
      dHy1fo <- gamma * Iy1foh * (1 - deltaY1)         - hospL * Hy1fo
      dHy1fv <- gamma * Iy1fvh * (1 - deltaY1*del_var) - hospL * Hy1fv
      
      
      
      dHy2no <- gamma * Iy2noh * (1 - deltaY2)         - hospL * Hy2no
      dHy2nv <- gamma * Iy2nvh * (1 - deltaY2*del_var) - hospL * Hy2nv
      
      dHy2io <- gamma * Iy2ioh * (1 - deltaY2)         - hospL * Hy2io
      dHy2iv <- gamma * Iy2ivh * (1 - deltaY2*del_var) - hospL * Hy2iv
      
      dHy2fo <- gamma * Iy2foh * (1 - deltaY2)         - hospL * Hy2fo
      dHy2fv <- gamma * Iy2fvh * (1 - deltaY2*del_var) - hospL * Hy2fv
      
      
      
      dHcno <- gamma * Icnoh * (1 - deltaC)         - hospL * Hcno
      dHcnv <- gamma * Icnvh * (1 - deltaC*del_var) - hospL * Hcnv
      
      dHcio <- gamma * Icioh * (1 - deltaC)         - hospL * Hcio
      dHciv <- gamma * Icivh * (1 - deltaC*del_var) - hospL * Hciv
      
      dHcfo <- gamma * Icfoh * (1 - deltaC)         - hospL * Hcfo
      dHcfv <- gamma * Icfvh * (1 - deltaC*del_var) - hospL * Hcfv
      
      
      
      
      
      dCano <- gamma * Ianoh * deltaA         - deltaL * Cano
      dCanv <- gamma * Ianvh * deltaA*del_var - deltaL * Canv
      
      dCaio <- gamma * Iaioh * deltaA         - deltaL * Caio
      dCaiv <- gamma * Iaivh * deltaA*del_var - deltaL * Caiv
      
      dCafo <- gamma * Iafoh * deltaA         - deltaL * Cafo
      dCafv <- gamma * Iafvh * deltaA*del_var - deltaL * Cafv

      
      
      dCy1no <- gamma * Iy1noh * deltaY1         - deltaL * Cy1no
      dCy1nv <- gamma * Iy1nvh * deltaY1*del_var - deltaL * Cy1nv
      
      dCy1io <- gamma * Iy1ioh * deltaY1         - deltaL * Cy1io
      dCy1iv <- gamma * Iy1ivh * deltaY1*del_var - deltaL * Cy1iv
      
      dCy1fo <- gamma * Iy1foh * deltaY1         - deltaL * Cy1fo
      dCy1fv <- gamma * Iy1fvh * deltaY1*del_var - deltaL * Cy1fv
      
      
      
      dCy2no <- gamma * Iy2noh * deltaY2         - deltaL * Cy2no
      dCy2nv <- gamma * Iy2nvh * deltaY2*del_var - deltaL * Cy2nv
      
      dCy2io <- gamma * Iy2ioh * deltaY2         - deltaL * Cy2io
      dCy2iv <- gamma * Iy2ivh * deltaY2*del_var - deltaL * Cy2iv
      
      dCy2fo <- gamma * Iy2foh * deltaY2         - deltaL * Cy2fo
      dCy2fv <- gamma * Iy2fvh * deltaY2*del_var - deltaL * Cy2fv
      
      
      
      dCcno <- gamma * Icnoh * deltaC         - deltaL * Ccno
      dCcnv <- gamma * Icnvh * deltaC*del_var - deltaL * Ccnv
      
      dCcio <- gamma * Icioh * deltaC         - deltaL * Ccio
      dCciv <- gamma * Icivh * deltaC*del_var - deltaL * Cciv
      
      dCcfo <- gamma * Icfoh * deltaC         - deltaL * Ccfo
      dCcfv <- gamma * Icfvh * deltaC*del_var - deltaL * Ccfv
      
      

      
      
      dRano <- gamma * Ianom + hospL * Hano + deltaL * Cano * (1 - zetaA)
      dRanv <- gamma * Ianvm + hospL * Hanv + deltaL * Canv * (1 - zetaA)
      
      dRaio <- gamma * Iaiom + hospL * Haio + deltaL * Caio * (1 - zetaA)
      dRaiv <- gamma * Iaivm + hospL * Haiv + deltaL * Caiv * (1 - zetaA)
      
      dRafo <- gamma * Iafom + hospL * Hafo + deltaL * Cafo * (1 - zetaA)
      dRafv <- gamma * Iafvm + hospL * Hafv + deltaL * Cafv * (1 - zetaA)
      
      
      
      dRy1no <- gamma * Iy1nom + hospL * Hy1no + deltaL * Cy1no * (1 - zetaY1)
      dRy1nv <- gamma * Iy1nvm + hospL * Hy1nv + deltaL * Cy1nv * (1 - zetaY1)
      
      dRy1io <- gamma * Iy1iom + hospL * Hy1io + deltaL * Cy1io * (1 - zetaY1)
      dRy1iv <- gamma * Iy1ivm + hospL * Hy1iv + deltaL * Cy1iv * (1 - zetaY1)
      
      dRy1fo <- gamma * Iy1fom + hospL * Hy1fo + deltaL * Cy1fo * (1 - zetaY1)
      dRy1fv <- gamma * Iy1fvm + hospL * Hy1fv + deltaL * Cy1fv * (1 - zetaY1)
      
      
      
      dRy2no <- gamma * Iy2nom + hospL * Hy2no + deltaL * Cy2no * (1 - zetaY2)
      dRy2nv <- gamma * Iy2nvm + hospL * Hy2nv + deltaL * Cy2nv * (1 - zetaY2)
      
      dRy2io <- gamma * Iy2iom + hospL * Hy2io + deltaL * Cy2io * (1 - zetaY2)
      dRy2iv <- gamma * Iy2ivm + hospL * Hy2iv + deltaL * Cy2iv * (1 - zetaY2)
      
      dRy2fo <- gamma * Iy2fom + hospL * Hy2fo + deltaL * Cy2fo * (1 - zetaY2)
      dRy2fv <- gamma * Iy2fvm + hospL * Hy2fv + deltaL * Cy2fv * (1 - zetaY2)
      
      
      
      dRcno <- gamma * Icnom + hospL * Hcno + deltaL * Ccno * (1 - zetaC)
      dRcnv <- gamma * Icnvm + hospL * Hcnv + deltaL * Ccnv * (1 - zetaC)
      
      dRcio <- gamma * Iciom + hospL * Hcio + deltaL * Ccio * (1 - zetaC)
      dRciv <- gamma * Icivm + hospL * Hciv + deltaL * Cciv * (1 - zetaC)
      
      dRcfo <- gamma * Icfom + hospL * Hcfo + deltaL * Ccfo * (1 - zetaC)
      dRcfv <- gamma * Icfvm + hospL * Hcfv + deltaL * Ccfv * (1 - zetaC)
      
      
      
      
      
      dDano <- deltaL * Cano * zetaA
      dDanv <- deltaL * Canv * zetaA
      
      dDaio <- deltaL * Caio * zetaA
      dDaiv <- deltaL * Caiv * zetaA
      
      dDafo <- deltaL * Cafo * zetaA
      dDafv <- deltaL * Cafv * zetaA
      
      
      
      dDy1no <- deltaL * Cy1no * zetaY1
      dDy1nv <- deltaL * Cy1nv * zetaY1
      
      dDy1io <- deltaL * Cy1io * zetaY1
      dDy1iv <- deltaL * Cy1iv * zetaY1
      
      dDy1fo <- deltaL * Cy1fo * zetaY1
      dDy1fv <- deltaL * Cy1fv * zetaY1
      
      
      
      dDy2no <- deltaL * Cy2no * zetaY2
      dDy2nv <- deltaL * Cy2nv * zetaY2
      
      dDy2io <- deltaL * Cy2io * zetaY2
      dDy2iv <- deltaL * Cy2iv * zetaY2
      
      dDy2fo <- deltaL * Cy2fo * zetaY2
      dDy2fv <- deltaL * Cy2fv * zetaY2
      
      
      
      dDcno <- deltaL * Ccno * zetaC
      dDcnv <- deltaL * Ccnv * zetaC
      
      dDcio <- deltaL * Ccio * zetaC
      dDciv <- deltaL * Cciv * zetaC
      
      dDcfo <- deltaL * Ccfo * zetaC
      dDcfv <- deltaL * Ccfv * zetaC
      
      
      
      
      dVaf <-  Sai  * vaccine_delay * vaceff_inf
      dVy1f <- Sy1i * vaccine_delay * vaceff_inf
      dVy2f <- Sy2i * vaccine_delay * vaceff_inf
      dVcf <-  Sci  * vaccine_delay * vaceff_inf
      
      

      
      
      
      
      
      
            
      return(list(c(dSan,dSai,dSaf,dSy1n,dSy1i,dSy1f,dSy2n,dSy2i,dSy2f,dScn,dSci,dScf,dIanom,dIanoh,dIanvm,dIanvh,dIaiom,dIaioh,dIaivm,dIaivh,dIafom,dIafoh,dIafvm,dIafvh,dIy1nom,dIy1noh,dIy1nvm,dIy1nvh,dIy1iom,dIy1ioh,dIy1ivm,dIy1ivh,dIy1fom,dIy1foh,dIy1fvm,dIy1fvh,dIy2nom,dIy2noh,dIy2nvm,dIy2nvh,dIy2iom,dIy2ioh,dIy2ivm,dIy2ivh,dIy2fom,dIy2foh,dIy2fvm,dIy2fvh,dIcnom,dIcnoh,dIcnvm,dIcnvh,dIciom,dIcioh,dIcivm,dIcivh,dIcfom,dIcfoh,dIcfvm,dIcfvh,dHano,dHanv,dHaio,dHaiv,dHafo,dHafv,dHy1no,dHy1nv,dHy1io,dHy1iv,dHy1fo,dHy1fv,dHy2no,dHy2nv,dHy2io,dHy2iv,dHy2fo,dHy2fv,dHcno,dHcnv,dHcio,dHciv,dHcfo,dHcfv,dCano,dCanv,dCaio,dCaiv,dCafo,dCafv,dCy1no,dCy1nv,dCy1io,dCy1iv,dCy1fo,dCy1fv,dCy2no,dCy2nv,dCy2io,dCy2iv,dCy2fo,dCy2fv,dCcno,dCcnv,dCcio,dCciv,dCcfo,dCcfv,dDano,dDanv,dDaio,dDaiv,dDafo,dDafv,dDy1no,dDy1nv,dDy1io,dDy1iv,dDy1fo,dDy1fv,dDy2no,dDy2nv,dDy2io,dDy2iv,dDy2fo,dDy2fv,dDcno,dDcnv,dDcio,dDciv,dDcfo,dDcfv,dRano,dRanv,dRaio,dRaiv,dRafo,dRafv,dRy1no,dRy1nv,dRy1io,dRy1iv,dRy1fo,dRy1fv,dRy2no,dRy2nv,dRy2io,dRy2iv,dRy2fo,dRy2fv,dRcno,dRcnv,dRcio,dRciv,dRcfo,dRcfv,dVaf,dVy1f,dVy2f,dVcf)))
      })
  }
  
  # the parameters values:
  parameters_values <- c(rt_ini=rt_ini, rt_rebLength=rt_rebLength, rt_reb=rt_reb, rt_var=rt_var, soe_eff=soe_eff, rt_age=rt_age, rt_same=rt_same, 
                         soe_time1=soe_time1, soe_time2=soe_time2, soe_time3=soe_time3, soe_time4=soe_time4, soe_time5=soe_time5, soe_time6=soe_time6,
                         gamma=gamma, deltaC_temp=deltaC_temp, deltaY1_temp=deltaY1_temp, deltaY2_temp=deltaY2_temp, deltaA_temp=deltaA_temp, deltaL=deltaL, del_var=del_var,
                         hosp_rateC=hosp_rateC, hosp_rateY=hosp_rateY, hosp_rateA=hosp_rateA, hospL_temp=hospL_temp,
                         zetaC_temp=zetaC_temp, zetaY1_temp=zetaY1_temp, zetaY2_temp=zetaY2_temp, zetaA_temp=zetaA_temp,
                         vaceff_inf=vaceff_inf, vaceff_death_temp=vaceff_death_temp, vac_brk=vac_brk,
                         popA=popA, popY1=popY1, popY2=popY2, popC=popC)
  
  # the initial values of variables:
  initial_values <- c(San=San0, Sai=Sai0, Saf=Saf0,
                      Sy1n=Sy1n0, Sy1i=Sy1i0, Sy1f=Sy1f0,
                      Sy2n=Sy2n0, Sy2i=Sy2i0, Sy2f=Sy2f0,
                      Scn=Scn0, Sci=Sci0, Scf=Scf0,
                      
                      Ianom=Ianom0, Ianoh=Ianoh0, Ianvm=Ianvm0, Ianvh=Ianvh0, Iaiom=Iaiom0, Iaioh=Iaioh0, Iaivm=Iaivm0, Iaivh=Iaivh0, Iafom=Iafom0, Iafoh=Iafoh0, Iafvm=Iafvm0, Iafvh=Iafvh0,
                      Iy1nom=Iy1nom0, Iy1noh=Iy1noh0, Iy1nvm=Iy1nvm0, Iy1nvh=Iy1nvh0, Iy1iom=Iy1iom0, Iy1ioh=Iy1ioh0, Iy1ivm=Iy1ivm0, Iy1ivh=Iy1ivh0, Iy1fom=Iy1fom0, Iy1foh=Iy1foh0, Iy1fvm=Iy1fvm0, Iy1fvh=Iy1fvh0,
                      Iy2nom=Iy2nom0, Iy2noh=Iy2noh0, Iy2nvm=Iy2nvm0, Iy2nvh=Iy2nvh0, Iy2iom=Iy2iom0, Iy2ioh=Iy2ioh0, Iy2ivm=Iy2ivm0, Iy2ivh=Iy2ivh0, Iy2fom=Iy2fom0, Iy2foh=Iy2foh0, Iy2fvm=Iy2fvm0, Iy2fvh=Iy2fvh0,
                      Icnom=Icnom0, Icnoh=Icnoh0, Icnvm=Icnvm0, Icnvh=Icnvh0, Iciom=Iciom0, Icioh=Icioh0, Icivm=Icivm0, Icivh=Icivh0, Icfom=Icfom0, Icfoh=Icfoh0, Icfvm=Icfvm0, Icfvh=Icfvh0,
                      
                      Hano=Hano0, Hanv=Hanv0, Haio=Haio0, Haiv=Haiv0, Hafo=Hafo0, Hafv=Hafv0,
                      Hy1no=Hy1no0, Hy1nv=Hy1nv0, Hy1io=Hy1io0, Hy1iv=Hy1iv0, Hy1fo=Hy1fo0, Hy1fv=Hy1fv0,
                      Hy2no=Hy2no0, Hy2nv=Hy2nv0, Hy2io=Hy2io0, Hy2iv=Hy2iv0, Hy2fo=Hy2fo0, Hy2fv=Hy2fv0,
                      Hcno=Hcno0, Hcnv=Hcnv0, Hcio=Hcio0, Hciv=Hciv0, Hcfo=Hcfo0, Hcfv=Hcfv0,
                      
                      Cano=Cano0, Canv=Canv0, Caio=Caio0, Caiv=Caiv0, Cafo=Cafo0, Cafv=Cafv0,
                      Cy1no=Cy1no0, Cy1nv=Cy1nv0, Cy1io=Cy1io0, Cy1iv=Cy1iv0, Cy1fo=Cy1fo0, Cy1fv=Cy1fv0,
                      Cy2no=Cy2no0, Cy2nv=Cy2nv0, Cy2io=Cy2io0, Cy2iv=Cy2iv0, Cy2fo=Cy2fo0, Cy2fv=Cy2fv0,
                      Ccno=Ccno0, Ccnv=Ccnv0, Ccio=Ccio0, Cciv=Cciv0, Ccfo=Ccfo0, Ccfv=Ccfv0,
                      
                      Dano=Dano0, Danv=Danv0, Daio=Daio0, Daiv=Daiv0, Dafo=Dafo0, Dafv=Dafv0,
                      Dy1no=Dy1no0, Dy1nv=Dy1nv0, Dy1io=Dy1io0, Dy1iv=Dy1iv0, Dy1fo=Dy1fo0, Dy1fv=Dy1fv0,
                      Dy2no=Dy2no0, Dy2nv=Dy2nv0, Dy2io=Dy2io0, Dy2iv=Dy2iv0, Dy2fo=Dy2fo0, Dy2fv=Dy2fv0,
                      Dcno=Dcno0, Dcnv=Dcnv0, Dcio=Dcio0, Dciv=Dciv0, Dcfo=Dcfo0, Dcfv=Dcfv0,
                      
                      Rano=Rano0, Ranv=Ranv0, Raio=Raio0, Raiv=Raiv0, Rafo=Rafo0, Rafv=Rafv0,
                      Ry1no=Ry1no0, Ry1nv=Ry1nv0, Ry1io=Ry1io0, Ry1iv=Ry1iv0, Ry1fo=Ry1fo0, Ry1fv=Ry1fv0,
                      Ry2no=Ry2no0, Ry2nv=Ry2nv0, Ry2io=Ry2io0, Ry2iv=Ry2iv0, Ry2fo=Ry2fo0, Ry2fv=Ry2fv0,
                      Rcno=Rcno0, Rcnv=Rcnv0, Rcio=Rcio0, Rciv=Rciv0, Rcfo=Rcfo0, Rcfv=Rcfv0,
                      
                      Vaf=Vaf0,
                      Vy1f=Vy1f0,
                      Vy2f=Vy2f0,
                      Vcf=Vcf0)
  
  # solving
  out <- ode(method="rk4",initial_values, times, sir_equations, parameters_values)
  #out <- lsode(initial_values, times, sir_equations, parameters_values,atol=1e-20,maxsteps=10000)
  
  #method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
  #           "euler", "rk4", "ode23", "ode45", "radau", 
  #           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration")

    
  # returning the output:
  simulationresult <- as.data.frame(out)
  return(simulationresult)
}


#
#Simulations
#
# 



list_rt_ini <- seq(0.01, 2.00, by = 0.01)
                 
list_del_var <- c(1)
list_vacC12 <- c(0, 0.45, 0.60, 0.75, 0.9)
list_vacY1 <-  c(0, 0.45, 0.60, 0.75, 0.9)
list_vacY2 <-  c(0, 0.60, 0.70, 0.80, 0.9)
list_vacA <-   c(0, 0.80, 0.85, 0.90, 0.9)
list_rt_reb <- c(0)
list_reb_len <- c(0)
list_soe1 <- c(999, 999, 999, 999, 999)
list_soe2 <- c(999, 999, 999, 999, 999)
list_soe3 <- c(999, 999, 999, 999, 999)
list_soe4 <- c(999, 999, 999, 999, 999)
list_soe5 <- c(999, 999, 999, 999, 999)
list_soe6 <- c(999, 999, 999, 999, 999)
list_age <- c(1,4)
list_same <- c(1,4)
list_vacInf <- c(0.7,0.9)
list_vacDis <- c(0.9,0.95)
list_vacBrk <- c(0.25,0.5)

i <- 0

for (set_rt_ini in list_rt_ini) {
  for (set_del_var in list_del_var) {
    for (select_vac in 1:5) {
      for (select_vacEff in 1:2) {
          for (set_rt_reb in list_rt_reb) {
            for (set_reb_len in list_reb_len) {
                    for (select_AgeSame in 1:2) {
              

      i = i+1
      
      vacC12 <- list_vacC12[select_vac]
      vacY1 <- list_vacY1[select_vac]
      vacY2 <- list_vacY2[select_vac]
      vacA <- list_vacA[select_vac]
      
      soe_time1 <- list_soe1[select_vac]
      soe_time2 <- list_soe2[select_vac]
      soe_time3 <- list_soe3[select_vac]
      soe_time4 <- list_soe4[select_vac]
      soe_time5 <- list_soe5[select_vac]
      soe_time6 <- list_soe6[select_vac]
      
      hosp_rateC = 0.1
      hosp_rateY = 0.2
      hosp_rateA = 0.8
      
      vaceff_inf <- list_vacInf[select_vacEff]
      vaceff_death_temp <- list_vacDis[select_vacEff]
      vac_brk <- list_vacBrk[select_vacEff]
      
      rt_age <-  list_age[select_AgeSame]
      rt_same <- list_same[select_AgeSame]
      
      timeL <- 150
      
      # adult/young1(20-30)/young2(40-50)/child

      result <- sir_1(rt_ini=set_rt_ini, rt_rebLength=set_reb_len, rt_reb=set_rt_reb, rt_var=1, soe_eff=0.5, rt_age = rt_age, rt_same = rt_same,
                      
                      soe_time1=soe_time1, soe_time2=soe_time2, soe_time3=soe_time3, soe_time4=soe_time4, soe_time5=soe_time5, soe_time6=soe_time6,
                      
                      gamma = 1/5, deltaC_temp = 0.03/100, deltaY1_temp = 0.06/100, deltaY2_temp = 1.0/100, deltaA_temp = 8.5/100, deltaL = 1/14, del_var=set_del_var,
                      
                      hosp_rateC = hosp_rateC, hosp_rateY = hosp_rateY, hosp_rateA = hosp_rateA, hospL_temp = 10, 
                      
                      zetaC_temp = 0.005/100, zetaY1_temp = 0.01/100, zetaY2_temp = 0.20/100, zetaA_temp = 5.7/100,
                      
                      vaceff_inf = vaceff_inf, vaceff_death_temp = vaceff_death_temp, vac_brk = vac_brk,
                      
                      popA=3700*(10^4), popY1=3250*(10^4), popY2=3350*(10^4), popC=1000*(10^4)+1300*(10^4),

                      
                      
                      San0 = 3700*0.99*(10^4)*(1-vacA),
                      Sai0 = 0, Saf0 = 3700*0.99*(10^4)*vacA*(1-vaceff_inf),
                      
                      Sy1n0 = 3250*0.99*(10^4)*(1-vacY1),
                      Sy1i0 = 0, Sy1f0 = 3250*0.99*(10^4)*vacY1*(1-vaceff_inf),
                      
                      Sy2n0 = 3350*0.99*(10^4)*(1-vacY2),
                      Sy2i0 = 0, Sy2f0 = 3350*0.99*(10^4)*vacY2*(1-vaceff_inf),
                      
                      Scn0 = 1000*0.99*(10^4)*(1-vacC12) + 1300*0.99*(10^4),
                      Sci0 = 0, Scf0 = 1000*0.99*(10^4)*vacC12*(1-vaceff_inf),
                      
                      Ianom0 = 100*(1-hosp_rateA), Ianoh0 = 100*hosp_rateA,
                      Ianvm0 = 0, Ianvh0 = 0,
                      Iaiom0=0, Iaioh0=0, Iaivm0=0, Iaivh0=0, Iafom0=0, Iafoh0=0, Iafvm0=0, Iafvh0=0,
                      
                      Iy1nom0 = 100*(1-hosp_rateY), Iy1noh0 = 100*hosp_rateY,
                      Iy1nvm0 = 0, Iy1nvh0 = 0,
                      Iy1iom0=0, Iy1ioh0=0, Iy1ivm0=0, Iy1ivh0=0, Iy1fom0=0, Iy1foh0=0, Iy1fvm0=0, Iy1fvh0=0,
                      
                      Iy2nom0 = 100*(1-hosp_rateY), Iy2noh0 = 100*hosp_rateY,
                      Iy2nvm0 = 0, Iy2nvh0 = 0,
                      Iy2iom0=0, Iy2ioh0=0, Iy2ivm0=0, Iy2ivh0=0, Iy2fom0=0, Iy2foh0=0, Iy2fvm0=0, Iy2fvh0=0, 
                      
                      Icnom0 = 100*(1-hosp_rateC), Icnoh0 = 100*hosp_rateC,
                      Icnvm0 = 0, Icnvh0 = 0,
                      Iciom0=0, Icioh0=0, Icivm0=0, Icivh0=0, Icfom0=0, Icfoh0=0, Icfvm0=0, Icfvh0=0,
                      
                      Hano0 = 0, Hanv0=0, Haio0=0, Haiv0=0, Hafo0=0, Hafv0=0,
                      Hy1no0 = 0, Hy1nv0=0, Hy1io0=0, Hy1iv0=0, Hy1fo0=0, Hy1fv0=0,
                      Hy2no0 = 0, Hy2nv0=0, Hy2io0=0, Hy2iv0=0, Hy2fo0=0, Hy2fv0=0,
                      Hcno0 = 0, Hcnv0=0, Hcio0=0, Hciv0=0, Hcfo0=0, Hcfv0=0,
                      
                      Cano0 = 0, Canv0=0, Caio0=0, Caiv0=0, Cafo0=0, Cafv0=0,
                      Cy1no0 = 0, Cy1nv0=0, Cy1io0=0, Cy1iv0=0, Cy1fo0=0, Cy1fv0=0,
                      Cy2no0 = 0, Cy2nv0=0, Cy2io0=0, Cy2iv0=0, Cy2fo0=0, Cy2fv0=0,
                      Ccno0 = 0, Ccnv0=0, Ccio0=0, Cciv0=0, Ccfo0=0, Ccfv0=0,
                      
                      Dano0=0, Danv0=0, Daio0=0, Daiv0=0, Dafo0=0, Dafv0=0,
                      Dy1no0=0, Dy1nv0=0, Dy1io0=0, Dy1iv0=0, Dy1fo0=0, Dy1fv0=0,
                      Dy2no0=0, Dy2nv0=0, Dy2io0=0, Dy2iv0=0, Dy2fo0=0, Dy2fv0=0,
                      Dcno0=0, Dcnv0=0, Dcio0=0, Dciv0=0, Dcfo0=0, Dcfv0=0,
                      
                      Rano0 = 3700*0.01*(10^4), Ranv0=0, Raio0=0, Raiv0=0, Rafo0=0, Rafv0=0,
                      Ry1no0 = 3250*0.01*(10^4), Ry1nv0=0, Ry1io0=0, Ry1iv0=0, Ry1fo0=0, Ry1fv0=0,
                      Ry2no0 = 3350*0.01*(10^4), Ry2nv0=0, Ry2io0=0, Ry2iv0=0, Ry2fo0=0, Ry2fv0=0,
                      Rcno0 = 2300*0.01*(10^4), Rcnv0=0, Rcio0=0, Rciv0=0, Rcfo0=0, Rcfv0=0,
                      
                      Vaf0 = 3700*0.99*(10^4)*vacA*vaceff_inf,
                      Vy1f0 = 3250*0.99*(10^4)*vacY1*vaceff_inf,
                      Vy2f0 = 3350*0.99*(10^4)*vacY2*vaceff_inf,
                      Vcf0 = 1000*0.99*(10^4)*vacC12*vaceff_inf,
                      
                      times = seq(0, timeL, by=1))

      
      


  var_rt_ini <- rep(set_rt_ini, timeL+1)
  var_del_var <- rep(set_del_var, timeL+1)
  var_vacC12 <- rep(vacC12, timeL+1)
  var_vacY1 <- rep(vacY1, timeL+1)
  var_vacY2 <- rep(vacY2, timeL+1)
  var_vacA <- rep(vacA, timeL+1)
  var_rt_reb <- rep(set_rt_reb, timeL+1)
  var_reb_len <- rep(set_reb_len, timeL+1)
  var_age <- rep(rt_age, timeL+1)
  var_same <- rep(rt_same, timeL+1)
  var_vacEffInf <- rep(vaceff_inf, timeL+1)
  var_vacEffDis <- rep(vaceff_death_temp, timeL+1)
  var_vacEffBrk <- rep(vac_brk, timeL+1)

  result <- cbind(result,var_rt_ini)
  result <- cbind(result,var_del_var)
  result <- cbind(result,var_vacC12)
  result <- cbind(result,var_vacY1)
  result <- cbind(result,var_vacY2)
  result <- cbind(result,var_vacA)
  result <- cbind(result,var_rt_reb)
  result <- cbind(result,var_reb_len)
  result <- cbind(result,var_age)
  result <- cbind(result,var_same)
  result <- cbind(result,var_vacEffInf)
  result <- cbind(result,var_vacEffDis)
  result <- cbind(result,var_vacEffBrk)
  
  
  if (i == 1) {
    comb_result <- result
  } else {
    comb_result <- rbind(comb_result,result)
  }
    
  }}}}}}}

hosp <- comb_result$Ianoh+comb_result$Ianvh+comb_result$Iaioh+comb_result$Iaivh+comb_result$Iafoh+comb_result$Iafvh+comb_result$Iy1noh+comb_result$Iy1nvh+comb_result$Iy1ioh+comb_result$Iy1ivh+comb_result$Iy1foh+comb_result$Iy1fvh+comb_result$Iy2noh+comb_result$Iy2nvh+comb_result$Iy2ioh+comb_result$Iy2ivh+comb_result$Iy2foh+comb_result$Iy2fvh+comb_result$Icnoh+comb_result$Icnvh+comb_result$Icioh+comb_result$Icivh+comb_result$Icfoh+comb_result$Icfvh+comb_result$Hano+comb_result$Hanv+comb_result$Haio+comb_result$Haiv+comb_result$Hafo+comb_result$Hafv+comb_result$Hy1no+comb_result$Hy1nv+comb_result$Hy1io+comb_result$Hy1iv+comb_result$Hy1fo+comb_result$Hy1fv+comb_result$Hy2no+comb_result$Hy2nv+comb_result$Hy2io+comb_result$Hy2iv+comb_result$Hy2fo+comb_result$Hy2fv+comb_result$Hcno+comb_result$Hcnv+comb_result$Hcio+comb_result$Hciv+comb_result$Hcfo+comb_result$Hcfv+comb_result$Cano+comb_result$Canv+comb_result$Caio+comb_result$Caiv+comb_result$Cafo+comb_result$Cafv+comb_result$Cy1no+comb_result$Cy1nv+comb_result$Cy1io+comb_result$Cy1iv+comb_result$Cy1fo+comb_result$Cy1fv+comb_result$Cy2no+comb_result$Cy2nv+comb_result$Cy2io+comb_result$Cy2iv+comb_result$Cy2fo+comb_result$Cy2fv+comb_result$Ccno+comb_result$Ccnv+comb_result$Ccio+comb_result$Cciv+comb_result$Ccfo+comb_result$Ccfv
severe <- comb_result$Cano+comb_result$Canv+comb_result$Caio+comb_result$Caiv+comb_result$Cafo+comb_result$Cafv+comb_result$Cy1no+comb_result$Cy1nv+comb_result$Cy1io+comb_result$Cy1iv+comb_result$Cy1fo+comb_result$Cy1fv+comb_result$Cy2no+comb_result$Cy2nv+comb_result$Cy2io+comb_result$Cy2iv+comb_result$Cy2fo+comb_result$Cy2fv+comb_result$Ccno+comb_result$Ccnv+comb_result$Ccio+comb_result$Cciv+comb_result$Ccfo+comb_result$Ccfv
death <- comb_result$Dano+comb_result$Danv+comb_result$Daio+comb_result$Daiv+comb_result$Dafo+comb_result$Dafv+comb_result$Dy1no+comb_result$Dy1nv+comb_result$Dy1io+comb_result$Dy1iv+comb_result$Dy1fo+comb_result$Dy1fv+comb_result$Dy2no+comb_result$Dy2nv+comb_result$Dy2io+comb_result$Dy2iv+comb_result$Dy2fo+comb_result$Dy2fv+comb_result$Dcno+comb_result$Dcnv+comb_result$Dcio+comb_result$Dciv+comb_result$Dcfo+comb_result$Dcfv
new_inf <- (comb_result$Ianom+comb_result$Ianoh+comb_result$Iaiom+comb_result$Iaioh+comb_result$Iafom+comb_result$Iafoh +comb_result$ Icnom+comb_result$Icnoh+comb_result$Iciom+comb_result$Icioh+comb_result$Icfom+comb_result$Icfoh +comb_result$ Ianvm+comb_result$Ianvh+comb_result$Iaivm+comb_result$Iaivh+comb_result$Iafvm+comb_result$Iafvh +comb_result$ Icnvm+comb_result$Icnvh+comb_result$Icivm+comb_result$Icivh+comb_result$Icfvm+comb_result$Icfvh +comb_result$ Iy1nom+comb_result$Iy1noh+comb_result$Iy1iom+comb_result$Iy1ioh+comb_result$Iy1fom+comb_result$Iy1foh +comb_result$ Iy2nom+comb_result$Iy2noh+comb_result$Iy2iom+comb_result$Iy2ioh+comb_result$Iy2fom+comb_result$Iy2foh +comb_result$ Iy1nvm+comb_result$Iy1nvh+comb_result$Iy1ivm+comb_result$Iy1ivh+comb_result$Iy1fvm+comb_result$Iy1fvh +comb_result$ Iy2nvm+comb_result$Iy2nvh+comb_result$Iy2ivm+comb_result$Iy2ivh+comb_result$Iy2fvm+comb_result$Iy2fvh)/5

comb_result <- cbind(comb_result, hosp)
comb_result <- cbind(comb_result, severe)
comb_result <- cbind(comb_result, death)
comb_result <- cbind(comb_result, new_inf)

