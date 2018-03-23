rMEA TODO
================
Johann R. Kleinbub

1.0.0.9011
----------

*date: "3/21/2018"*

### easy stuff for next release

-   Identify outliers in MEAlagplot()
-   'GET' functions, e.g. ~~getCCF~~, getMEA
-   add a $rand object to *MEA* class
-   ~~all functions accepting *MEAlist* objects now should accept a list instead and try to parse as *MEAlist* internally~~
-   Calculate % of simultaneous zeroes (when zeros are during at least 5s)
-   clean code: debug comments, italian comments, old commented code
-   diagnostic plots with 2 rows, on the second showing random section in great detail, e.g. 10 s of nonempty signal
-   add outlier removal to the readme example
-   export nowness % information

### advanced wishlist

-   export long vector plots
-   Implement peak picking algorithms results in $ccfRes
-   **smart** confidence interval analysis for bootstrap
-   dynamic window size!!!
-   cross-validation to define a confidence band for each lag value, instead of mean (for MEAlagplots and group/bootstrap comparisons)
