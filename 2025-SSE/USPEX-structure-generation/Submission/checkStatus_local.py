from __future__ import absolute_import
import argparse
import glob
import os

import subprocess
from subprocess import check_output

_author_ = 'etikhonov'


def checkStatus_local(jobID):
    u"""
    This function is to check if the submitted job is done or not
    One needs to do a little edit based on your own case.
    1   : whichCluster (0: no-job-script, 1: local submission, 2: remote submission)
    Step1: the command to check job by ID. 
    Step2: to find the keywords from screen message to determine if the job is done
    Below is just a sample:
    -------------------------------------------------------------------------------
    Job id                    Name             User            Time Use S Queue
    ------------------------- ---------------- --------------- -------- - -----
    2455453.nano              USPEX            qzhu            02:28:42 R cfn_gen04 
    -------------------------------------------------------------------------------
    If the job is still running, it will show as above.
    
    If there is no key words like 'R/Q Cfn_gen04', it indicates the job is done.
    :param jobID: 
    :return: doneOr
    """

    # Step 1
    try:
        output = check_output('squeue -j {}'.format(jobID), shell=True, universal_newlines=True)
    except:
        output = "JOBDONE"

    # process = subprocess.Popen(['qstat', str(jobID)], stdout=subprocess.PIPE)
    # output, err = process.communicate()

    # Step 2
    doneOr = True
    if ' R ' in output or ' PD ' in output:
        doneOr = False
    if doneOr:
        for file in glob.glob('USPEX*'):
            os.remove(file)  # to remove the log file
    print(str(doneOr))
    return doneOr

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', dest='jobID', type=int)
    args = parser.parse_args()

    isDone = checkStatus_local(jobID=args.jobID)
    print('<CALLRESULT>')
    print(int(isDone))
