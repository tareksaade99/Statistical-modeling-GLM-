\section*{Generalized Linear Models -- Applied Statistical Modeling Project}

\subsection*{Overview}
This repository contains an applied statistical modeling project focused on \textbf{Generalized Linear Models (GLMs)}. The objective of the project is to bridge theoretical understanding and practical implementation by developing, fitting, and critically evaluating GLMs on real-world datasets. Emphasis is placed on interpretability, model diagnostics, and robustness rather than purely predictive performance.

\subsection*{Project Objectives}
The main goals of the project are:
\begin{itemize}
\item To understand the theoretical foundations of GLMs, including link functions, exponential family distributions, and likelihood-based inference.
\item To implement core estimation algorithms manually, with particular focus on the \textbf{Iteratively Reweighted Least Squares (IRLS)} procedure.
\item To fit and compare Poisson and logistic regression models in applied settings.
\item To assess model quality using likelihood-based criteria and nested model comparisons.
\item To address common real-world issues such as overdispersion and parameter uncertainty.
\end{itemize}

\subsection*{Methods and Implementation}
The project includes both theoretical analysis and hands-on implementation. A custom IRLS algorithm was developed from scratch to estimate GLM parameters, allowing direct inspection of convergence behavior and numerical stability. Results obtained from the manual implementation were systematically compared to those produced by standard \texttt{glm} functions, showing near-identical estimates and validating the correctness of the approach.

Model evaluation relied on deviance, log-likelihood, and information criteria, with explicit comparisons between nested models to assess the contribution of additional predictors. Coefficients were interpreted in terms of their practical meaning (e.g., rate ratios in Poisson models), ensuring statistical results were connected to real-world implications.

\subsection*{Uncertainty and Robustness Analysis}
To quantify uncertainty beyond asymptotic approximations, bootstrap resampling techniques were applied to estimate confidence intervals for selected parameters. This allowed comparison between bootstrap-based and asymptotic intervals and provided insight into model stability, particularly in borderline or weak-effect scenarios.

Special attention was given to overdispersion in count data, discussing its impact on inference and the limitations of standard Poisson assumptions when applied to complex datasets.

\subsection*{Key Takeaways}
\begin{itemize}
\item Manual implementation of GLM estimation algorithms provides valuable insight into model behavior and assumptions.
\item Likelihood-based model comparison is essential for understanding the trade-off between model complexity and explanatory power.
\item Bootstrap methods offer a practical tool for assessing uncertainty when classical assumptions may be unreliable.
\item Interpretability and diagnostic analysis are central to responsible statistical modeling.
\end{itemize}

\subsection*{Repository Structure}
The repository is organized to reflect the progression of the project, from model implementation to analysis and interpretation. Code, results, and written report sections are structured to support clarity, reproducibility, and ease of review.

\subsection*{Intended Audience}
This project is intended for students, data analysts, and researchers with an interest in statistical modeling, particularly those seeking a deeper understanding of GLMs beyond black-box usage of existing libraries.
