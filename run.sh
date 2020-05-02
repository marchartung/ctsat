#!/bin/bash
/usr/sbin/sshd -D &

BASENAME="${0##*/}"
log () {
  echo "${BASENAME} - ${1}"
}
HOST_FILE_PATH="/tmp/hostfile"
#aws s3 cp $S3_INPUT $SCRATCH_DIR
#tar -xvf $SCRATCH_DIR/*.tar.gz -C $SCRATCH_DIR

sleep 2
echo main node: ${AWS_BATCH_JOB_MAIN_NODE_INDEX}
echo this node: ${AWS_BATCH_JOB_NODE_INDEX}
echo Downloading problem from S3: ${S3_BKT}/${COMP_S3_PROBLEM_PATH}

if [[ "${COMP_S3_PROBLEM_PATH}" == *".xz" ]];
then
  aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/test.cnf.xz
  unxz supervised-scripts/test.cnf.xz
else
  aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/test.cnf
fi

echo "Starting $${S3_BKT}/${COMP_S3_PROBLEM_PATH}"

  # REPLACE THE FOLLOWING LINE WITH YOUR PARTICULAR SOLVER
time /ctsat/parallel/ctsat_parallel_release -nthreads=${NUM_PROCESSES}  supervised-scripts/test.cnf
