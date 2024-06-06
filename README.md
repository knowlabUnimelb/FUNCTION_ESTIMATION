# FUNCTION_ESTIMATION
Function Estimation data from Little, Shiffrin &amp; Laham (2022) https://psyarxiv.com/5c9na/

Project contains four files:

function_estimation_stimuli.dat 
-- scatterplot data shown to observers

-- Columns:
  - subject code
  - item code
  - x value
  - y value
  - number of points in scatter plot
  - scale of scatter plot (1 = zoomed in, 2 = zoomed out)
  - generating function (1 = linear, 2 = quadratic, 3 = cubic)

function_estimation_responses.dat
-- subject drawn function downsampled to 40 points

-- Columns
  - subject code (matches stimulus file)
  - item code  (matches stimulus file)
  - x value
  - y value

Supplement - Function Questionairre.pdf 
-- Function estimation stimuli

# ARCHIVED CODE AND DATA
Folders and files associated with: 
    Little, D. R., Shiffrin, R., & Laham, S. (2024). Function estimation: Quantifying individual differences in hand-drawn functions. Memory & Cognition. [Accepted 26-May-2024]

/analysis_code/
- /gaussian_process/
-- contains the code to estimate parameters for the gaussian process 

/expert_classification/
- Contains expert classification of responses into linear, quadratic, cubic, and data-tracking functions

/questionairre/
- original questionairre and code to generate the graphs used in the questionairre

/raw_data_scans/
- original scans of questionairres 

/summary_data/
- summary responses and stimuli

=========================================================================

Analysis steps 
1) Participants completed the function estimation questionairre using pencil 
  and paper. This resulted in some idiosyncracies that needed to be corrected
  (e.g., non-function drawings; overlapping x-values; multiple lines)

2) Questionairres were scanned by a Research Assistant to a pdf and then .JPG 
  images were created from the pdf for each graph and response

3) Datathief (https://datathief.org/) B. Tummers, DataThief III. 2006 <https://datathief.org/>
   was used to transform the drawings into useable numbers for further 
   pre-processing
    - Some useful functions for this step:
        a) Rename data (using renameData.m) if necessary to give new subject number
        b) Run organizeData.m - combines the individual scanned files into a single mat
        c) Run compileDataFile.m - writes the data in the form used by the GP analysis (in the GP folder)

4) Pre-preprocessing steps are combined in preprocessing_pipeline_script.m
    a) preprocess_1.m will display all of the graphs and drawing data for a given subject
        - Usage: preprocess_1(subject, ploton, saveon)
        - if saveon is true, then it will additionally create a mat file which 
          identifies the trials that have multiple responses. This mat file can 
          then be passed to preprocess_2_correct_multiple_responses.m

    b) preprocess_2_correct_multiple_responses.m
        - Usage: preprocess_2_correct_multiple_responses(subject)
        -- This function will read the outputfile of preprocess_1_out_s%_data.mat 
           or a multiple_resp_corrected_*.mat file (if exists)
        -- If there are multiple responses to any graph, the function will 
            present the graph and prompt the user to select one to keep
    c) preprocess_3.m
        - Usage: preprocess_3(subject, ploton, saveon)
        -- reads the preprocess_1_out_s%_data.mat or multiple_resp_corrected_*.mat file (if exists)
        -- The output is a mat file containing finalData matrices of sorted data, 
           containing only one response with the noise data adjusted to the 
           appropriate size and the draw data replicated so that there are only 
           nan's where there are no observed or response data - all other trial
           information is replicated to make indexing easier
        -- output files are final_*_data.mat

5) Next step is to fit the GP regression gp_analysis_gp1.m using only the covSEiso kernel
    - Output is saved in fits/gpAnalysis_Full_s*.mat

6) Other analyses then conducted on estimated parameters of covSEiso 
   a) /analysis_code/plot_expert_classifications
        - .fig files show plot of responses color coded by expert classification
        - handParmClassifications.m - plot all parms color coded by expert category

   b) /analysis_code/dp_classified_clusters
        - contains figures to generate graphs and responses color coded by cluster
        - contains code to edit figures for publication (reordering figures based on 
          stimulus features so that they can be presented in an ordered fashion)

   c) /analysis_code/linear_discriminant_anlaysis
        - use GRT toolbox (Alfonso-Reese, 2006) to fit boundaries to estimated parameters

        A. Alfonso-Reese, L. (2006). General recognition theory of categorization: A MATLAB toolbox. Behavior research methods, 38, 579-583.

   d) /analysis_code/parameter_cluster_analysis
        - use https://github.com/jacobeisenstein/DPMM toolbox to cluster parameters
          using a Dirichlet Process clustering model

Other useful files in /analysis_code/gaussian_process/:
- showAllSubjectsData.m
-- plot a single graph with all subject responses

- showSubjectData
-- plot all graphs for a single subject showing downsampled drawing locations

- socfun_subject_numbers.txt
-- list of subject numbers
