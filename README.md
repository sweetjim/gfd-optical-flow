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
### **You must have a version of MATLAB onwards of the 2020 release and Image Processing Toolbox installed.**

1. Clone this package to an accessible directory.
3. Open Matlab (version 2021 onwards with Image Processing Toolbox installed) and navigate to the package directory.
4. Execute _app_selector_ in the command line and select the root folder which contains experiment data. The program should then display a tree-node file explorer with folders and .nc files.

    - **Experiment data is assumed to be contained in NetCDF files.**
    - **By default, the image data is assumed to be stored in a variable titled _'SST'_. This can be changed by opening _exptdata.loadImage_ and changing _'SST'_ to the desired variable**
    ```javascript
    [line 659] im = ncread(obj.ncinfo.Filename,'SST',[1 1 imageindex],[Inf Inf 1]);
    ```
5. Using the checkmark boxes beside the file/folder nodes, _mark_ experiments for preview. Experiments that are marked will be available for preview in a figure window.
  - Use the _image slider_ to preview the output of the .nc file.
  - Use the _experiment slider_ to preview the output of different .nc files.
  - **Marked experiments will be imported to the 'base' workspace under the variable _expts_ as a cellular list of filepaths.**
6. With an experiment selected (and marked), press _Optical Flow_ to open the Optical-Flow interface application.

---
# Saving OpticalFlow data
If the Optical-Flow app is open
The script _sequencer.m_ is an example of processing

---
# MAKING CHANGES TO DEFAULT SETTINGS
### Directory of experiments
  - Modify _set_startup.m_
### NetCDF variable
  - Modify _exptdata.loadImage_
### Reference frame (cartesian or polar)
  - Modify _exptdata.refframe_
