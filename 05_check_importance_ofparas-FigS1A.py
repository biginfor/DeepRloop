import matplotlib
matplotlib.use('TkAgg')
import optuna
import os
import random
SEED = 1
random.seed(SEED)
import numpy as np
np.random.seed(SEED)
#
study_name = "test"
#
working_dir = './data/'
os.chdir(working_dir)
#
loaded_study = optuna.load_study(study_name=study_name, storage="sqlite:///db.sqlite3")
fig_data = optuna.visualization.plot_param_importances(loaded_study)
fig_data.write_image(working_dir+"paras_importance.pdf")

