# **NE 204 Lab 2 Materials**
## Chris Lamb
### November 21st, 2022
___
### **Make sure to download data from link below before running any code**

***Important to note that `Lab-1/Data/Trapezoid_Heights/` has been modified so download the whole folder and replace it from before***

The directories should be structured such that: `Data/Combined_Data/` and `Data/Trapezoid_Heights/` are both directories that exist within the `NE-204/Lab-1/` directory. Outside the Lab 1 directory, there should be a python file containing all the relevant imports and created functions for this course `NE204_Functions.py`. Similarly, there should be a structure to the `Data` directory within `NE204/Lab-2/` such that `Data/Combined_Data/` and `Data/Clusters-Optimization/` exist with `Clusters-Optimization/` containing 20 files: `means_cluster_n.npy` and `stds_cluster_n.npy` with n ranging from 0 to 9.
#### [Lab 1 Data Link](https://drive.google.com/drive/folders/1z9oz9HB6gjCqFKw-JAUBIuuRUfO1F-cn?usp=sharing)
#### [Lab 2 Data Link](https://drive.google.com/drive/folders/1RqpK3-Ft9CDd36E8rBJOeZe8C1V5dYDv?usp=sharing)
___
#### Operation of the code to generate figures
> **Usage:**
> <br />Open the `NE-204-Lab-2-Final.ipynb` file in a Jupyter notebook environment with the Anaconda Navigator app or as below in a terminal:
>> `jupyter notebook NE-204-Lab-2-Final.ipynb`
>
> This will open the notebook in a local webpage and allow for operation of the code.
>
> Make sure to run the first cell as it imports all the functions from the `NE204_Functions.py` file as well as some local variables defined within this file.
>
> From there operation should be as straightforward as running cells within sections to reproduce the plots shown within the cells.

As an important note, optimization of the shaping parameters for each individual cluster takes quite a while. Doing so for all 10 clusters took on the order of 9 hours to run on my laptop, results may vary.
