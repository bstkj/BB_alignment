** README **

---------------------------------------------------------------------------
Short summaries of the top-level functions. These are the files that are 
located within the main BB_alignment folder, but not included in any subfolder
(e.g. visualize_raw, visualize_MIP, ..., setup, model2report, generated_marked_tif).
---------------------------------------------------------------------------
1) visualize_*: 
Each of the visualize_* functions is meant to display a certain aspect or 

2) train_model
Takes an .nd2 raw image file and trains a model of the BB distribution within
the given cell. 

3) setup
Run this directly in the BB_alignment folder, and before doing anything else.
It sets up all the paths so that MATLAB knows where to search for your functions.

4) model2report
Produces various graphs/numerical figures summarizing certain properties of
the trained model.

5) generate_marked_tif
Takes an .nd2 raw image, picks out BBs using our pipeline and produces a tif
copy of the raw image with colored markings indicating the BB positions in 
each z-layer of the image. 

---------------------------------------------------------------------------
Short summaries of the subfolders within the main BB_alignment folder
---------------------------------------------------------------------------
1) nd2_files:
Contains our raw image files.

2) models:
Contains the models we train from our raw image files.

3) generated_tifs:
Contains the marked tifs that we generate from processing our raw image files.

4) utilities:
Contains all the helper functions to our top-level functions.

5) demos: 
Contain scripts that 'show' you how to use the top level functions.

6) development:
Contains functions/scripts that are being used for testing new extensions or 
improvements to the current pipeline.

7) bfmatlab:
Contains bioformats for MATLAB package that allows us to process our raw 
.nd2 image files.

8) getting_rid_of:
Scripts/functions that we are almost certainly going to remove or replace.