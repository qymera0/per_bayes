
# 0 LOAD PACKAGES ---------------------------------------------------------

library(lme4)
library(extraoperators)
library(JWileymisc)
library(multilevelTools)

# 1 LOAD DATA -------------------------------------------------------------

data("aces_daily", package = 'JWileymisc')

## 2 MIXED EFFECTS --------------------------------------------------------

strictControl <-
        lmerControl(
                optCtrl = list(
                        algorithm = 'NLOPT_LN_NELDERMEAD',
                        xtol_abs = 1e-12,
                        ftol_abs = 1e-12
                )
        )


m1 <-
        lmer(
                NegAff ~ STRESS + (1 + STRESS | UserID),
                data = aces_daily,
                control = strictControl
        )

mdlDiagnosis <-
        modelDiagnostics(
                m1,
                ev.perc = .001
        )

plot(
        mdlDiagnosis,
        ask = FALSE,
        ncol = 2,
        nrow = 3
)

modelPerformance(m1)

summary(m1)

# Model tests

mt3 <- modelTest(m1)

names(mt3)

APAStyler(mt3)
