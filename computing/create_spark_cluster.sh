#!/bin/bash

DEFAULT_NUM_WORKERS=2

CLUSTER_NAME=$1
NUM_WORKERS=${2:-$DEFAULT_NUM_WORKERS}

# create the cluster
echo "Creating cluster $CLUSTER_NAME"
gcloud dataproc clusters create "$CLUSTER_NAME" \
    --metadata "JUPYTER_PORT=8785,JUPYTER_CONDA_PACKAGES=seaborn:boto3:s3fs:py-xgboost:numpy:pandas:scikit-learn:scipy=0.19.1:plotly:matplotlib:pysal,JUPYTER_CONDA_FORGE_PACKAGES=tornado=4.5.3:fuzzywuzzy:geopy:google-auth:google-cloud-storage:google-cloud-core" \
    --initialization-actions \
		gs://conrad-project/gcloud/jupyter.sh,gs://conrad-project/gcloud/bootstrap_gcloud.sh,gs://dataproc-initialization-actions/ganglia/ganglia.sh \
    --initialization-action-timeout "60m" \
    --subnet default --zone northamerica-northeast1-b \
    --tags "spark,dataproc-cluster-$CLUSTER_NAME" \
    --project conrad-197521 \
    --num-workers $NUM_WORKERS \
    --master-machine-type "n1-standard-1" \
    --worker-machine-type "n1-standard-1" \
    --properties "spark-env:PYTHONPATH=/repos/main/" \
    --scopes "https://www.googleapis.com/auth/cloud-platform,https://www.googleapis.com/auth/source.read_write"
