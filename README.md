# QALYs-Lost-due-to-Genital-Herpes
### Lifetime Quality-Adjusted Life Years Lost Due to Genital Herpes Acquired in the United States in 2018
This repository includes codes to run model simulations and to visualize simulation outcomes for the manuscript **"Lifetime Quality-Adjusted Life Years Lost Due to Genital Herpes Acquired in the United States in 2018: mathematical modeling stud"** by Shiying You, Reza Yaesoubi, Kyueun Lee, Yunfei Li, Samuel T. Eppink, Katherine K. Hsu, Harrell W. Chesson, Thomas L. Gift, Andrés A. Berruti, Joshua A. Salomon, Minttu M. Rönn.


### Packages
All analyses were conducted in Python 3.8 on MacOS using packages including [NumPy](https://numpy.org/), [Pandas](https://pandas.pydata.org/), [SciPy](https://scipy.org/), [statsmodles](https://www.statsmodels.org/stable/index.html), [Matplotlib](https://matplotlib.org/), and [SimPy](https://github.com/yaesoubilab/SimPy) ([the version](https://github.com/yaesoubilab/SimPy/commit/f8a0249c1c384fdfb2e4520a79cab3e985652b10) committed on 2021 July)


### Organization of The Repository
- **analyses**: 
- **classes**:
- **source**:
- **supports**:
- **tree_outputs**:
- **visualization**:


### Template Code for Running Simulation for genital HSV-1 using MacOS Terminal
git clone https://github.com/yaesoubilab/SimPy.git

git clone https://github.com/shiying88/QALYs-Lost-due-to-Genital-Herpes.git

pip install numpy pandas matplotlib scipy statsmodels

export PYTHONPATH="${PYTHONPATH}:/{path_when_you_clone_SimPy}/SimPy"

cd {path_to_root_folder_of_QALYs-Lost-due-to-Genital-Herpes}

python3 analyses/RunHSV1ProbTree.py
