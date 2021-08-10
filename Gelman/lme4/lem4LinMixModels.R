
# 0 LOAD PACKAGES ---------------------------------------------------------

library(lme4)

# 1 INTRODUCTION ----------------------------------------------------------

## 1.1 Linear Mixed Models ------------------------------------------------


## 1.2 Example ------------------------------------------------------------

data("sleepstudy")

fm1 <-
        lmer(
                Reaction ~ Days + (Days | Subject),
                sleepstudy
        )

summary(fm1)

RShowDoc("lmerperf", package = "lme4")


## 1.3 High-level modular structure ---------------------------------------

# Recriate fm1 as step by step


# Parse formula

parsedFormula <-
        lFormula(
                formula = Reaction ~ Days + (Days | Subject),
                data = sleepstudy
        )

# Create the -log-likelihood functions

devianceFunction <-
        do.call(
                mkLmerDevfun, parsedFormula
        )

# Run optimization

optimizerOutput <-
        optimizeLmer(
                devianceFunction
        )

# Output

mkMerMod(
        rho = environment(devianceFunction),
        opt = optimizerOutput,
        reTrms = parsedFormula$reTrms,
        fr = parsedFormula$fr
)


# 2 FORMULA MODULE --------------------------------------------------------


# # 2.2 Understanding mixed model formulas --------------------------------

# formula without correlation parameter

fm2 <-
        lmer(
                Reaction ~ Days + (Days || Subject),
                sleepstudy
        )

summary(fm2)
