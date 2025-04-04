# hdps-tmle-sim
Code and materials for a simulation study evaluating high-dimensional proxy adjustment using hdPS, machine learning (LASSO, SL), and TMLE. Includes plasmode simulations based on NHANES data to assess bias, coverage, and the impact of model complexity in high-dimensional settings.

## Run Simulations for Each Method in Each Simulation Scenario

### 1. Navigate to the `simulation` Folder

### 2. Navigate to the Corresponding Folder for Each Scenario
- For each simulation scenario (i.e., frequent exposure & outcome, rare exposure, and rare outcome), go to the appropriate folder (i.e., `frequent`, `rare_exposure`, and `rare_outcome`, respectively).

### 3. Running Simulations for Each Method:

**3.1 For Each Method Except SL Methods:**

3.1.1 Navigate to the Corresponding Folder for Each Method
- Each folder contains 1 to 3 files. For folders containing more than one file, the file name will indicate the SL library used, with filenames ending in `_1` for the 1-learner library, `_3` for the 3-learner library, and `_4` for the 4-learner library.

3.1.2 Run the File to Generate the Results
- Run the file to generate the simulation results. To adjust the number of simulation iterations, modify the values of `iterStart` and `iterStop`. The total number of simulations generated will be: `iterStop - iterStart + 1`.

3.1.3 View the Results
- The estimated results will be provided in the data frame `resultDF`.

**3.2 For the SL Methods:**

3.2.1 Run the Corresponding TMLE Method
- Follow the instructions in section `3.1` to run the corresponding TMLE method. For example, to generate simulation results for the `hdPS.SL` method using the 1-learner library, first run the simulation for the `hdPS.TMLE` method with the 1-learner library. The simulation code for SL methods with the 3- and 4-learner libraries is in the same file (denoted by `_3_4` at the end of the file name), and the corresponding TMLE method for both the 3- and 4-learner libraries needs to be run.

3.2.2 Copy the Weight Folder Generated by the TMLE Method
- After running the TMLE method, copy the generated weight folder named `weight` followed by the corresponding `method name(s)` and paste it into the SL folder. Do **not** change the folder name.


3.2.3 Run the SL File to Generate the Results
- Run the SL file to generate the results. As with other methods, the number of simulation iterations can be adjusted by modifying the `iterStart` and `iterStop` values. The total number of simulations generated will be: `iterStop - iterStart + 1`.


3.2.4 View the Results
- The estimated results will be provided in the data frame `resultDF`.

