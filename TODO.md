rMEA TODO
================
Johann R. Kleinbub

1.0.0.9010
----------

*date: "3/21/2018"*

-   Identify outliers in MEAlagplot()
-   'GET' functions, e.g. getCCF, getMEA
-   add a $rand object to *MEA* class
-   all functions accepting *MEAlist* objects should accept a list instead and try to parse as *MEAlist* internally
-   Calculate % of simultaneous zeroes (when zeros are during at least 5s)
-   Implement peak picking algorithms results in $ccfRes
-   clean code: debug comments, italian comments, old commented code
-   export long vector plots
-   diagnostic plots with 2 rows, on the second showing random section in great detail, e.g. 10 s of nonempty signal
-   add outlier removal to the readme example
-   **smart** confidence interval analysis for bootstrap
