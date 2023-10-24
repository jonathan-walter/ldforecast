# ldforecast
Derived datasets and code accompanying the manuscript (in review) Forecasting the spread of an invasive forest-defoliating insect

Data files included are: 
a) gridded occupancy files derived from Slow the Spread program trap catch data. Names follow the convention GridTrapsYEAR_full_presenceN.tif, where N is a number indicating the trap catch threshold corresponding to occupancy (1, 3, or 10 moths).
b) Annual total area of defoliation attributed to Lymantria dispar; Defol.csv.
3) Annual climatic suitability layers. Names follow the convention csuit_diff_YEAR.tif. These are annual deviations from average climatic suitability by grid cell.
c) Long-term average environmental suitability; csuit_mean_2000_2009.tif.
d) Forest canopy cover fraction; forest_canopy_5km.tif.
e) Omniscape outputs further organized into subdirectories. Subdirectory names give the year and search radius.

Code included are R scripts that fit each model.
