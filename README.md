# gfd-optical-flow
Interface for viewing experiment data and processing OpticalFlow calculations within MATLAB.
---
# INFORMATION
This package contains three key components:
  1. A dynamic class object *exptdata.m*.
  - This class (to be treated as a repository) executes in-built functions and stores 
      information about an experiment. Refer to the documentation for more information
      regarding its full functionality.
      
   **Usage (THROUGH COMMAND-LINE)**
   (TO LOAD)
   ```javascript
   >> filepath_handle = 'cd\filename.nc'     % Example path
   >> varname = exptdata(filepath_handle);
   ```
    
   **(TO OPERATE)**
   ```javascript
   >> vel_struct = varname.getVelocity(1)    % Example OpticalFlow velocity calculation.
   >> [x,z,t,ds,dt] = varname.getDimensions; % Example dimensional vector output
   ```
      
  2. A single-running-instance (SRI) application *app_selector.mlapp*.
    - This GUI program behaves as an explorer for navigating and filtering experiments within
      a specified directory. Experiments may be 'checked' for loading and 'selected' for 
      viewing. 
      Experiments that are 'checked' will be assigned to the workspace in a cellular list
      containing only their path on the directory.
      
      - Usage (THROUGH COMMAND-LINE)
      ```javascript 
      >> app_selector
      ```
      
      
    
  3. A dynamic SRI application 'app_velocityfield.mlapp'
    - This GUI program is intended for previewing OpticalFlow calculations on a specified 
      experiment. It offers two algorithms for calculating optical flow velocities:
      
      +  [Horn-Schunck's algorithm](https://doi.org/10.1016/0004-3702(81)90024-2)
         - Spatiotemporal derivative estimation
         - Weighted cross-correlation averaging kernel
         - Iterative velocity computation

      +  [Lie-Shun's algorithm](https://doi.org/10.1017/S0022112008003273)
            - Coarse field estimatation using Horn-Schunck algorithm
            - Refinement with inverted matrix and tolerance-threshold loop 
      
      
      - Usage (THROUGH APP_SELECTOR)
      ```javascript
      > Select an experiment
      > Press Optical-Flow pushbutton
      ```
      
      (THROUGH COMMAND-LINE)
      ```javascript
      >> app_velocityfield
      ```
---
# GETTING STARTED
