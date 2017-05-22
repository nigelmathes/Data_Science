import tensorflow as tf
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import numpy as np

# Columns:
# 0    1    2    3    4    5    6    7    8         9    10    11
# ID  FUV  Wav  SN   Flux eFlux EW  FWHM Line_FUV  Mass  z    label (0 = no, 1 = yes)

COLUMNS = ["id", "fuv", "wav", "sn", "flux", "eflux", "ew", "fwhm", "line_fuv", "mass", "z", "label"]
LABEL_COLUMN = 'label'
CONTINUOUS_COLUMNS = ["fuv", "mass", "z"]

# Define data set filenames
LYA_TRAINING = pd.read_csv("zcatalog_training.csv",names=COLUMNS)
LYA_TEST = pd.read_csv("zcatalog_test.csv",names=COLUMNS)

def input_fn(df):
    # Create a dictionary mapping from continuous feature column k to
    # the values in the column stored in a constant Tensor
    continuous_cols = {k: tf.constant(df[k].values) for k in CONTINUOUS_COLUMNS}
    feature_cols = dict(continuous_cols.items())
    label = tf.constant(df[LABEL_COLUMN].values)
    return feature_cols, label

def train_input_fn():
    return input_fn(LYA_TRAINING)

def eval_input_fn():
    return input_fn(LYA_TEST)

def build_estimator(model_dir, model_type):
    """Build an estimator."""

    # Continuous base columns.
    fuv = tf.contrib.layers.real_valued_column("fuv")
    mass = tf.contrib.layers.real_valued_column("mass")
    z = tf.contrib.layers.real_valued_column("z")

    deep_columns = [fuv, mass, z]

    # Deep Neural Net
    if model_type == "dnn":
        m = tf.contrib.learn.DNNClassifier(model_dir=model_dir,
                                           feature_columns=deep_columns,
                                           hidden_units=[200, 100])
    else:
        m = 0
    return m

def main(_):
    # Define model directory or create it
    model_dir = "MODEL"
    print("Model directory = %s" % model_dir)

    # ======================== Deep Neural Net ==================================

    # Define number of steps to train
    train_steps = 2000

    # Build the model for DNN
    m = build_estimator(model_dir, 'dnn')

    # Fit the model for DNN
    m.fit(input_fn=lambda: input_fn(LYA_TRAINING), steps=train_steps)

    # Evaluate the accuracy of the model for DNN
    results = m.evaluate(input_fn=lambda: input_fn(LYA_TEST), steps=1)
    for key in sorted(results):
        print("%s: %s" % (key, results[key]))
    print "\n====== Done with DNN Classifier ======\n"

    # ======================== RANDOM FOREST ===================================
    clf = RandomForestClassifier(n_jobs=2)
    y, _ = pd.factorize(LYA_TRAINING[LABEL_COLUMN])
    clf.fit(LYA_TRAINING[CONTINUOUS_COLUMNS],y)

    predictions = clf.predict(LYA_TEST[CONTINUOUS_COLUMNS])
    pd.crosstab(LYA_TEST[LABEL_COLUMN], predictions, rownames=['actual'], colnames=['predictions'])
    print "Random Forest Classification Accuracy: {}".format(np.sum(LYA_TEST[LABEL_COLUMN]==predictions) / float(len(LYA_TEST)))

tf.app.run(main=main)
