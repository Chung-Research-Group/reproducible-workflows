import argparse
import glob
import os
import time

import subprocess
from subprocess import check_output


SLEEP_TIME = 20
SUBMISSION_SCRIPT = 'sleep {}'.format(SLEEP_TIME)
INDEX = 1
WORKING_DIR = '~'



def test_submission_scripts(whichCluster : int):
    if whichCluster == 1:
        from checkStatus_local import checkStatus_local
        from submitJob_local import submitJob_local
        jobID = submitJob_local(index=INDEX, commandExecutable=SUBMISSION_SCRIPT)
        print(str(jobID))
        assert(not checkStatus_local(jobID=jobID))
    elif whichCluster == 2:
        from checkStatus_remote import checkStatus_remote
        from submitJob_remote import submitJob_remote
        jobID = submitJob_remote(workingDir=WORKING_DIR, index=INDEX, commandExecutable=SUBMISSION_SCRIPT)
        assert(checkStatus_remote(jobID=jobID, workingDir=WORKING_DIR, index=INDEX))
        time.sleep(SLEEP_TIME)
        assert(not checkStatus_remote(jobID=jobID))
    else:
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest='whichCluster', type=int)
    args = parser.parse_args()

    test_submission_scripts(whichCluster=args.whichCluster)
