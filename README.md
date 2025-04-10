# **Local Interpolation Lab**

### **Overview**
This repository contains the activity for Week 12 of AGR 333: Local Interpolation Lab. The lab focuses on comparing different methods of spatial interpolation using R. Interpolation techniques are used to "gap fill" spatial data, predicting values at unsampled locations based on sparse point observations.

The dataset used in this lab is `soja98`, which contains gridded measurements of soil fertility and soybean yield. This activity will focus specifically on potassium (`K`) data.

---

### **Learning Objectives**
By completing this lab, students will:
1. Understand the concepts behind spatial interpolation methods such as **Inverse Distance Weighting (IDW)** and **Ordinary Kriging**.
2. Learn how to create grids for interpolation and visualize results.
3. Evaluate interpolation accuracy using **Root Mean Squared Error (RMSE)**.
4. Explore the effect of sample size on interpolation accuracy.
5. Compare geostatistical methods (e.g., Kriging) with non-geostatistical methods (e.g., IDW).

---

### **Prerequisites**
Before starting this lab, ensure you have:
- Basic knowledge of R programming.
- Familiarity with spatial data concepts.
- R installed on your computer along with RStudio.

---

### **Required R Packages**
The following R packages are required for this lab:
- `geoR`: For geostatistical analysis.
- `gstat`: For spatial interpolation methods.
- `sp`: For spatial data manipulation.
- `raster`: For raster data handling.
- `RColorBrewer`: For color palettes used in visualization.
- `yardstick`: For calculating RMSE as an alternative method.

You can install these packages using:
```R
install.packages(c("geoR", "gstat", "sp", "raster", "RColorBrewer", "yardstick"))
```

---

### **Dataset**
The dataset used in this lab is `soja98`, which is included in the `geoR` package. You can load it directly into R using:
```R
library(geoR)
data(soja98)
```

---

### **Instructions**
Follow the steps below to complete the lab:

#### 1. Load and Prepare Data
- Load the `soja98` dataset and inspect its structure using `help()` and `head()`.
- Select only relevant columns (`X`, `Y`, and potassium (`K`)) and convert the dataframe into a spatial points object.

#### 2. Inverse Distance Weighting (IDW)
- Create a grid for interpolation using `extent()` and `expand.grid()`.
- Perform IDW interpolation using calibration data (`cal`) and visualize results with `spplot()`.
- Evaluate accuracy using RMSE for both calibration and validation datasets.

#### 3. Effect of Sample Size on RMSE
- Simulate different sample sizes using a loop to observe their effect on RMSE.
- Plot RMSE against sample size to analyze trends.

#### 4. Ordinary Kriging
- Fit spherical and exponential variogram models to calibration data using `variogram()` and `fit.variogram()`.
- Perform kriging using both models and visualize predictions with `spplot()`.
- Evaluate kriging accuracy using RMSE for calibration and validation datasets.

---

### **Files in This Repository**
1. **`README.md`**: This file provides an overview of the lab activity.
2. **`local_interpolation_lab.R`**: The main R script containing all tasks for the lab.
3. **`AGR_333_Interpolation-Lab-Guide-Week-12.docx`**: Detailed instructions and background information about the lab.

---

### **Usage**
1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/your_username/local-interpolation-lab.git
   ```
2. Open the R script (`local_interpolation_lab.R`) in RStudio.
3. Follow the step-by-step instructions provided in the script to complete the lab tasks.

---

### **Key Notes**
1. **Grid Resolution**: Adjusting grid cell size (`cellSize`) affects both computational demands and interpolation accuracy. Experiment with different values (e.g., 10m vs 5m).
2. **RMSE Interpretation**: Validation RMSE is typically larger than calibration RMSE due to overfitting when predicting calibration data with itself.
3. **Kriging vs IDW**: Kriging accounts for spatial structure through semivariance, making it more robust under sparse sampling conditions compared to IDW.

---

### **Reflection Questions**
After completing the lab, consider:
1. How does sample size impact interpolation accuracy?
2. Why does validation RMSE differ from calibration RMSE?
3. Under what conditions does kriging outperform IDW?

