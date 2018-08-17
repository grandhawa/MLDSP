%  
% Path ->  /Users/djones/work/mWork/Fathom
% 
%   cl                  - same as clc
%   clz                 - same as clear + clc + close all
%   f_AIC               - compute AIC (or BIC) for least-squares based methods or classifiers
%   f_CCorA             - canonical correlation analysis (CCorA)
%   f_CCorA_PL          - CCorA based on P. Legendre's vegan code
%   f_CCorAplot         - plot results from a CCorA analysis
%   f_NSIDCimport       - import a National Snow & Ice Data Center (NSIDC) sea ice coverage binary file
%   f_NSIDCinterp       - interpolate sea ice coverage values for a specific location
%   f_adjustP           - adjust p-values for multiple comparison tests
%   f_ancova            - nonparametric (permutation-based) 1-way ANCOVA
%   f_ancovaPW          - a posteriori, multiple-comparison tests for Homogeneous Slopes in ANCOVA
%   f_anosim            - 1-way Analysis of Similarity (ANOSIM)
%   f_anosim2           - 2-way crossed ANOSIM with no replication
%   f_arrow             - draw an arrow
%   f_author            - display author and license information for the 'Fathom Toolbox for Matlab'
%   f_balance           - get indices of elements to create a balanced design
%   f_bartlett          - Bartlett's test for Homogeneity of Variances
%   f_beep              - play a custom sound file
%   f_betainc           - incomplete beta function ('non-regularized')
%   f_bioenv            - correlation between distance matrix and all subsets of X (BEST or BIOENV)
%   f_biplot            - eigenvector-based 2-d distance biplot
%   f_biplotEnv2        - create environmental vectors for 2-d nMDS ordination distance biplot
%   f_biplotEnv3        - create environmental vectors for 3-d nMDS ordination distance biplot
%   f_biplotPca2        - create 2-d PCA distance biplot
%   f_biplotSpecies     - create species vectors for ordination distance biplot
%   f_boot              - bootstrap resampling with replacement
%   f_bootCI            - bootstrapped confidence interval of the mean for univariate data
%   f_braycurtis        - legacy Bray-Curtis symmetric distance matrix function
%   f_brokenstick       - determine # of significant ordination dimensions via "Broken-Stick model"
%   f_brush2idx         - create index to rows of data currently selected by 'brushing'
%   f_bubble            - create B&W bubble plots
%   f_butter            - smooth columns of time series data via Butterworth lowpass filter
%   f_cap               - Canonical Analysis of Principal Coordinates (CAP) using ANY distance matrix
%   f_capMLE            - maximum likelihood estimator of mixture proportions based on CAP
%   f_capOptimal        - get optimal value of m for f_cap
%   f_capPlot           - plot results from a CAP analysis
%   f_cda               - Canonical Discriminant Analysis
%   f_cda632            - estimate CDA error rate using bootstrap .632+ rule
%   f_cdaBCV            - bootstrap cross-validation for Canonical Discriminant Analysis
%   f_cdaCV             - leave-one-out cross validation for Canonical Discriminant Analysis
%   f_cdaClass          - classify unknown observations for f_cda
%   f_cdaPlot           - plot results from a Canonical Discriminant Analysis
%   f_cellstr2cat       - convert a cell array of strings to a single categorical variable
%   f_cellstr2num       - convert a cell array of strings to a single categorical variable
%   f_center            - center data on column mean
%   f_central           - central tendency (mean, median, or mode)
%   f_centroid          - returns coordinates of the centroid of X, optionally partitioned into groups
%   f_chanceClass       - classification success expected by chance (= proportional chance criterion)
%   f_chisqIndep        - Chi-square test of independence (= homogeneity of proportions)
%   f_clock             - current date and time as hyphen ('-') delimited character string
%   f_cluster           - UPGMA-based cluster analysis of a symmetric dissimilarity matrix
%   f_colormap          - create a custom colormap
%   f_compet            - get index of column containing the maximum value of each row
%   f_confEllipse       - parametric or bootstrapped confidence ellipse for bivariate data
%   f_convHull          - convex hull for bivariate data
%   f_copy              - copies a numerical-array to the clipboard
%   f_corr              - Pearson's, Spearman's, or Kendall's correlation between 2 vectors
%   f_corrSign          - determine if a correlation b/n 2 vectors is positive or negative
%   f_correlogram       - multivariate Mantel correlogram
%   f_cov               - returns covariance (= dispersion) matrix from X
%   f_covPool           - covariance matrix pooled across groups
%   f_cv                - coefficient of variation
%   f_deg2utm           - convert lat/lon vectors into UTM coordinates (WGS84)
%   f_delaunay          - create a Delaunay triangulation from 2-D spatial coordinates
%   f_depthCM           - 'Depth of Center of Mass' for plankton survey data
%   f_depthMean         - weighted 'Mean Depth of Occurrence' for plankton survey data
%   f_designMatrix      - create ANOVA design matrix using dummy variables
%   f_diag              - replace diagonals of a square matrix
%   f_dis               - create symmetric dissimilarity (or distance) matrix
%   f_dis2sim           - convert square symmetric distance matrix to a similarity matrix
%   f_disprof           - dissimilarity profile analysis (DISPROF, SIMPROF)
%   f_disprof_clust     - dissimilarity profile analysis (DISPROF) of a cluster analysis dendrogram
%   f_disprof_clustPlot - plot dendrogram of a DISPROF-based custer analysis
%   f_disprof_clust_bin - compare 2 DISPROF binary connectivity matrices
%   f_dnn               - distance-based Nearest Neighbor Graph from 2-D spatial coordinates
%   f_dummy             - dummy coding of categorical variables
%   f_dummy2cat         - convert dummy codes to a single categorical variable
%   f_eig               - eigenanalysis of a square matrix
%   f_eigenMaps         - create Moran's eigenvector maps from spatial coordinates
%   f_eigenMapsStepwise - stepwise selection of Moran's Eigenvector Maps (MEM's)
%   f_ekmanDepth        - depth of the Ekman layer
%   f_empPDF            - empirical probability density function
%   f_errRate           - error rate for a classifier
%   f_euclid            - legacy Euclidean distance matrix function
%   f_export            - write ASCII delimited file.
%   f_exportDods        - export DODS bathymetry for import into Surfer
%   f_exportR           - export data for import into R
%   f_extractFields     - extract structure fields & combine into a single matrix
%   f_figArea           - finds the area of figure objects
%   f_filterMA          - smooth columns of a matrix using a Moving Average (or Median)
%   f_filterSinclair    - filter/smooth columns of a matrix via 11-point Moving Median + Average
%   f_findCell          - index to rows of Y that match X
%   f_firstOccur        - returns indices of the first occurrence of unique elements of input vector
%   f_gabriel           - create a Gabriel graph from 2-D spatial coordinates
%   f_getNLOM           - download Navy Layered Ocean Model SSH Nowcast imagery
%   f_gower             - Gower's centered matrix
%   f_graphviz_mds      - create a Graphviz DOT file to perform multidimensional scaling (MDS)
%   f_graphviz_mst      - export a minimal spanning tree to Graphviz DOT format
%   f_graphviz_neato    - create an undirected graph using Graphviz
%   f_greenwood         - Greenwood's statistic
%   f_greenwood_cdf     - cumulative probability density function for Greenwood's statistic
%   f_greenwood_par     - get parameters of the distribution of Greenwood's statistic for n>2
%   f_greenwood_pdf     - probability density function for Greenwood's statistic
%   f_greenwood_plt     - plot the PDF or CDF for Greenwood's statistic
%   f_greenwood_rnd     - randomized distribution of Greenwood's statistic
%   f_gregorian         - convert Julian date to Gregorian
%   f_grpBoot           - within-group bootstrap sampling (fixed size)
%   f_grpBootMix        - within-group bootstrap sampling (variable size), with mixing among columns
%   f_grpMean           - returns the mean/stdv/se of X (column-wise) separately for groups
%   f_grpOutlier        - logical index identifying outliers, separately for each group
%   f_grpPlot           - group plotting function
%   f_grpRel            - returns the relative proportion of X (column-wise) separately for groups
%   f_grpResample       - within-group resampling or bootstrapping (variable size)
%   f_grpSize           - get the number of observations in each category for a grouping variable
%   f_gshhs2Shp         - export a shapefile from the GSHHS database for use in Surfer/Arcview
%   f_halfchange        - scale nMDS configuration to half-change
%   f_hellinger         - Hellinger transform data or create a symmetric dissimilarity matrix
%   f_helmert           - Helmert orthogonal contrast codes for an ANOVA design matrix
%   f_hist              - create a normalized histogram plot of a column vector
%   f_ibc               - indicator-species based cluster analysis
%   f_importCSV         - import numeric data from a comma-separated-values (*.csv) file
%   f_importShapefile   - import ArcView shapefile for M_Map Mapping Toolbox
%   f_importSurfer      - import Surfer *.bln blanking file for M_Map
%   f_importSurferBln   - import Surfer *.bln blanking file for M_Map
%   f_importSurferGrd   - import Surfer GRID (*.grd) file
%   f_indVal            - species indicator values
%   f_inv               - matrix inversion via "\" (left division)
%   f_isAbsent          - determine which rows of A are absent from B
%   f_isDST             - determine if dates occur during Daylight Savings Time
%   f_isOdd             - determine if integer is odd or even
%   f_isPresent         - determine which rows of A are present in B
%   f_isScalar          - determine if input is a scalar
%   f_issymSim          - determine if input is square symmetric similarity matrix
%   f_issymdis          - determine if input is square symmetric distance matrix
%   f_julian            - convert date vector to Julian date
%   f_julianSub         - get indices to subsample a series of Julian dates
%   f_kde               - Botev's nonparametric kernel density estimator for univariate data
%   f_labelplot         - create 2-d label plot
%   f_latlong           - compute terrestrial distance matrix
%   f_leapYear          - determine if a year is a leap year
%   f_lenfreq           - length-frequency plot (adjusted for sampling effort)
%   f_ll2utm            - legacy function (use f_deg2utm instead)
%   f_lowpass           - smooth columns of a matrix via a lowpass filter
%   f_lowpass40         - smooth columns of a matrix via a 6th order, 40-hr lowpass filter
%   f_mad               - median absolute deviation
%   f_mantel            - standardized Mantel statistic for 2 symmetric distance matrices
%   f_mle               - maximum likelihood estimator of mixture proportions via an EM algorithm
%   f_mleCook           - maximum likelihood estimation via Cook's constrained corrected classification
%   f_mocnessAb         - returns abundance (# under 10 m^2 of sea surface)
%   f_modelMatrix       - create model matrix for a Mantel Test
%   f_month2num         - convert month string to number
%   f_moran             - Moran's coefficient of spatial autocorrelation
%   f_mregress          - Multiple Linear Regression via Least Squares Estimation
%   f_mst               - Minimum Spanning Tree from a symmetric distance matrix
%   f_mst_mex           - create a minimal spanning tree (MEX version)
%   f_mst_old           - legacy version of f_mst
%   f_multicomb         - generate permutation distribution of grouping labels
%   f_nan2ave           - replace missing values (NaN's) with average value (columnwise) by group
%   f_ncap              - Nonlinear Canonical Analysis of Principal Coordinates (NCAP)
%   f_ncapOptimal       - evaluate a range of values of m for f_ncap
%   f_nmds              - Nonmetric Multidimensional Scaling (nMDS)
%   f_nmdsPlot          - plot results from a nonmetric Multidimensional Scaling (nMDS)
%   f_nnMLP             - feedforward multilayer perceptron neural network for classification
%   f_nnMLP632          - estimate MLP error rate using bootstrap .632+ rule
%   f_nnMLPcv           - leave-one-out cross validation for f_nnMLP
%   f_normal            - normalize values of X
%   f_npDisp            - distance-based measures of multivariate dispersion
%   f_npDispAll         - similar to f_npDisp, but returns distance to all centroids
%   f_npDispPlot        - plot distance-based measures of Homogeneity in Multivariate Dispersion
%   f_npMLE             - distance-based, maximum likelihood estimator
%   f_npManova          - nonparametric (permutation-based) MANOVA for ANY distance matrix
%   f_npManovaPW        - a posteriori, multiple-comparison tests
%   f_num2cell          - convert a vector of numbers to a cell array of strings
%   f_numAtLength       - estimate number of fish within length classes for visual survey data
%   f_origin            - mark the origin of a graph with horizontal/vertical line(s)
%   f_outlier           - calculate outlier measures for a symmetric proximity (similarity) matrix
%   f_outlyingness      - normalized measure of outlyingness
%   f_outlyingness_test - permutation test of maximum outlyingness
%   f_parse             - parse matrix to structure
%   f_pca               - Principal Component Analysis of a data matrix
%   f_pcaPlot           - plot results from a Principal Components Analysis
%   f_pcnm              - Principal Coordinates of Neighbor Matrices (PCNM) from spatial coordinates
%   f_pcoa              - Principal Coordinates Analysis (PCoA)
%   f_pcoaPlot          - plot results from a Principal Coordinates Analysis
%   f_pdf               - export current figure to PDF file
%   f_pdfCrop           - crop a PDF file to its bounding box
%   f_pdfDistill        - distill Postscript files to PDF using Acrobat Distiller
%   f_pdfLatex          - typeset a LaTeX file
%   f_pdfMerge          - merge a series of PS (or PDF) files into a single PDF
%   f_pdfSplit          - split a multipage PDF file into separate single-page PDF's
%   f_perError          - percent error between measured and accepted values
%   f_permtest          - two sample permutation test of means
%   f_perturb           - randomly perturb the values of a matrix, separately for each column
%   f_plotBarsH         - horizontal bar graphs
%   f_plotBarsV         - vertical bar graphs
%   f_plotError         - add error bars to current figure
%   f_plotNeigh         - plotting function for f_delaunay, f_dnn, f_gabriel, f_mst, & f_relNeigh
%   f_plotUSGS          - plot USGS Coastwatch SST image with M_Map Toolbox
%   f_pnn               - Probabilistic Neural Network
%   f_pnn632            - estimate PNN error using bootstrap .632+ rule
%   f_pnnAIC            - stepwise selection of explanatory variables in PNN using AIC (or BIC)
%   f_pnnCV             - leave-one-out cross validation for PNN
%   f_pnnSm             - determine optimal smoothing factor for a PNN
%   f_poly              - create a matrix of polynomials of degree d
%   f_prd               - percent relative difference between 2 sets of values
%   f_procrustes        - Procrustes rotation of Y to X
%   f_qSort             - quick sort elements i:j of a vector, ascending
%   f_randDir           - draw random numbers from a Dirichlet distribution
%   f_randRange         - returns n random integers ranging from min to max
%   f_randSub           - extract a random subset of N rows from matrix X
%   f_randWH            - Wichmann & Hill's (2006) good pseudo-random number generator
%   f_range             - return the min and max values of a vector
%   f_ranging           - scale columns of x to range [0 to 1] or [-1 to 1]
%   f_ranks             - ranks data in X (column-wise) with averaging of ties
%   f_rda               - Redundancy Analysis (RDA)
%   f_rdaAIC            - stepwise forward selection of explanatory variables in RDA using AIC (or BIC)
%   f_rdaAnova          - permutation-based (M)ANOVA via RDA
%   f_rdaDB             - distance-based Redundancy Analysis (db-RDA)
%   f_rdaDB_AIC         - stepwise forward selection of explanatory variables in db-RDA using AIC (or BIC)
%   f_rdaDB_Stepwise    - stepwise selection of explanatory variables in db-RDA based on F-stat
%   f_rdaDB_manova      - nonparametric (permutation-based) MANOVA via db-RDA
%   f_rdaPlot           - ordination distance biplot for a Redundancy Analysis (RDA)
%   f_rdaStepwise       - stepwise selection of explanatory variables in RDA based on F-stat
%   f_readcwf           - import a Coastwatch Satellite SST file
%   f_recode            - recode elements of vector as consecutive integers
%   f_relNeigh          - create a Relative Neighbor graph from 2-D spatial coordinates
%   f_rename            - rename variable without memory reallocation
%   f_renameBatch       - rename variables across multiple files
%   f_renameField       - rename a structure field
%   f_rewrap            - wraps vector into symmetric distance matrix (reverses f_unwrap)
%   f_rgb               - utility program for selecting color of plot symbols
%   f_round             - round data to specified number of decimal places
%   f_rsd               - percent relative standard deviation
%   f_selectW           - AIC-based selection of optimal spatial weighting matrices
%   f_shadeBox          - shade subsets of a time series plot
%   f_shuffle           - randomly sorts vector, matrix, or symmetric distance matrix
%   f_simper            - similarity percentages (SIMPER) & species contributions among 2 groups
%   f_smooth            - smooth columns of a matrix using a Moving Average filter
%   f_sort              - sort elements of X (ascending) and create sorting/unsorting indices
%   f_stdErr            - returns the standard error
%   f_stnd              - standardize values of a matrix, column-wise (= z-scores)
%   f_struct2flat       - turn a structure into a flat table & optionally export
%   f_style             - utility program for selecting line styles
%   f_sub               - subsample data every 6 hrs
%   f_subsetDisPW       - extract subsets of distance matrix based on all pairs of a grouping factor
%   f_subsetPW          - extract subsets of column vector based on all pairs of a grouping factor
%   f_svd               - singular value decomposition of a matrix
%   f_swapInt           - swap 2 integer values
%   f_symb              - utility program for selecting plot symbols
%   f_table             - display a matrix X in the command window
%   f_trajectory        - get distance along a trajectory
%   f_transform         - several methods for data transformation
%   f_unique            - returns unsorted list of unique values
%   f_unwrap            - unwrap lower tri-diagonal (w/o diag) of symmetric distance matrix into a column vector
%   f_utm2deg           - convert vectors of UTM coordinates into Lat/Lon vectors (WGS84)
%   f_variogram         - multivariate empirical variogram
%   f_vecAngle          - counter-clockwise angle between 2 points
%   f_vecDiagram        - plot progressive vector diagrams
%   f_vecMagDir         - get magnitude & direction from u,v vector components
%   f_vecPlot           - plot time series of velocity vectors
%   f_vecRot            - rotate vectors (U,V) by angle THETA
%   f_vecRotAngle       - determine the angle of an isobath
%   f_vecTrans          - transform 2d vector coordinates
%   f_vecTrans3d        - transform 3d vector coordinates
%   f_vecUV             - returns U,V components of a vector given its magnitude & direction
%   f_vectorfit         - plot environmental ordination vectors via multiple linear regression
%   f_vonBert           - von Bertalanffy growth function
%   f_vonBertAge        - age a fish of known length with the von Bertalanffy growth equation
%   f_vonBertModel      - fit a von Bertalanffy growth model using nonlinear regression 
%   f_vonBertModelInit  - initial parameter estimate of a von Bertalanffy growth model
%   f_vonBert_La        - calculate La parameter for the von Bertalanffy growth function
%   f_vonBert_t0        - calculate t0 parameter for the von Bertalanffy growth function
%   f_wascores          - weighted-averages scores of species for a site ordination
%   f_windCMAN          - process CMAN (or NDBC) historical wind data
%   f_windstress        - wind stress in dynes/cm^2
%   f_writeTxt          - write character array to a text file
%   f_wtMean            - weighted mean
%   f_xMatrix           - design matrix of contrast codes for ANOVA linear models
%   f_xdiss             - calculate Extended Dissimilarities from a symmetric distance matrix
%   f_xmlExtract        - extract node, name, attribute, or tag from XML tree
%   finder              - open present working directory in a MacOS X's 'FINDER' window
%   m_2D_surf           - draw a 2D surface on an M_Map map
%   m_bubble            - M_Map compatible version of f_bubble
%   m_extContours       - extract contours from an M_Map map
%   m_ginput            - get lon/lat coordinates from an M_Map figure
%   m_scaleBar          - create horizontal scale bar for M_Map plots
%   m_subset            - create a subset of an M_Map usercoast file
%   saver               - save workspace to file tagged in 'fname' variable to current directory
