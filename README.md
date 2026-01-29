# Generalized Linear Models — Applied Statistical Modeling Project

## Overview

This repository contains an applied statistical modeling project focused on **Generalized Linear Models (GLMs)**. The project aims to connect statistical theory with practical implementation by developing, fitting, and critically evaluating GLMs on real-world datasets. The emphasis is on interpretability, model diagnostics, and robustness rather than purely predictive performance.

## Project Objectives

The main objectives of the project are:

* Understand the theoretical foundations of GLMs, including link functions, exponential family distributions, and likelihood-based inference.
* Implement core estimation algorithms manually, with particular focus on the **Iteratively Reweighted Least Squares (IRLS)** algorithm.
* Fit and compare Poisson and logistic regression models in applied settings.
* Evaluate models using likelihood-based criteria and nested model comparisons.
* Address common real-world issues such as overdispersion and parameter uncertainty.

## Methods and Implementation

The project combines theoretical analysis with hands-on implementation. A custom IRLS algorithm was developed from scratch to estimate GLM parameters, enabling direct inspection of convergence behavior and numerical stability. The estimates obtained from the manual implementation were compared against standard `glm` functions, showing very close agreement and validating the implementation.

Model performance was assessed using deviance, log-likelihood, and information criteria. Nested models were systematically compared to evaluate the contribution of individual predictors. Model coefficients were interpreted in terms of their practical meaning (e.g., rate ratios in Poisson regression), ensuring that statistical results remained interpretable in real-world contexts.

## Uncertainty and Robustness Analysis

To assess uncertainty beyond asymptotic approximations, bootstrap resampling methods were applied to estimate confidence intervals for selected parameters. This allowed direct comparison between bootstrap-based and asymptotic intervals and provided insight into model stability, particularly for weak or borderline effects.

The project also examines overdispersion in count data, discussing its implications for inference and the limitations of standard Poisson assumptions.

## Key Takeaways

* Manual implementation of GLM estimation algorithms provides deeper insight into model behavior and assumptions.
* Likelihood-based model comparison helps balance model complexity and explanatory power.
* Bootstrap methods are valuable for uncertainty quantification when classical assumptions may be unreliable.
* Interpretability and diagnostic analysis are essential components of responsible statistical modeling.

## Repository Structure

The repository is organized to reflect the progression of the project, from algorithm implementation to analysis and interpretation. Code, results, and report sections are structured to promote clarity, reproducibility, and ease of review.

## Intended Audience

This project is intended for students, data analysts, and researchers interested in statistical modeling, particularly those seeking a deeper understanding of GLMs beyond black-box use of existing libraries.
