This Directory provides the supporting scripts to calibrate the BRICK model (add link), process model results and produce figures...

## Folders

**calib_driver**
- R Scripts for model calibration for high temperature scenario
  1.  CARL_calib_driver.R - generalized version of the calibration script
  2.  CARL_calib_driver_standard.R - calibration of model to data and setting the standard priors
  3.  CARL_calib_driver_expert.R - calibration of model to expert predictions and adjusting the standard priors accordingly
  4.  CARL_calib_driver_priors.R - calibration of model BRICK data and priors provided by past liturature
  5.  CARL_calib_driver_complete.R- calibration of model to combine standard data, expert predictions and prior data

**processingPipeline** 
- RScripts for processing BRICK model from calibrations discribed above for high temperature scenario
  1.  CARL_processingPipeline_standard.R - process model calibrated to data and setting the standard priors
  2.  CARL_processingPipeline_expert.R - process model calibrated to expert predictions and adjusting the standard priors accordingly
  3.  CARL_processingPipeline_priors.R - process model calibrated BRICK data and priors provided by past liturature
  4.  CARL_cprocessingPipeline_complete.R- cprocess model calibrated to combine standard data, expert predictions and prior data
 
**analysis_and_plots**
- RScripts for data extraction, analysis and figure generation
