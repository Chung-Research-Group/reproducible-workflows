
import argparse
import os

from subprocess import check_output

def checkStatus_remote(jobID, workingDir, index):
    """
    This routine is to check if the submitted job is done or not
    One needs to do a little edit based on your own case.
    Step1: Specify the PATH to put your calculation folder
    Step2: Check JobID, the exact command to check job by jobID
    :param jobID:
    :param index:
    :param workingDir:
    :return:
    """
    # Step 1
    Home = '/home/etikhonov'  # 'pwd' of your home directory of your remote machine
    Address = 'rurik'  # Your target supercomputer: username@address or ssh alias
    # example of address: user@somedomain.edu -p 2222
    Path = Home + '/' + workingDir + '/CalcFold' + str(index)  # just keep it

    # Step 2
    output = str(check_output('ssh ' + Address + ' qstat ' + str(jobID), shell=True))
    # If you using full adress without ssh alias, you must provide valid ssh private key like there:
    # output = str(check_output('ssh -i ~/.ssh/id_rsa ' + Address + ' /usr/bin/qstat ' + str(jobID), shell=True))

    if not ' R ' in output or not ' Q ' in output:
        doneOr = True
        # [nothing, nothing] = unix(['scp -i ~/.ssh/id_rsa ' Address ':' Path '/OUTCAR ./']) %OUTCAR is not necessary by default
        os.system('scp ' + Address + ':' + Path + '/OSZICAR ./')  # For reading enthalpy/energy
        os.system('scp ' + Address + ':' + Path + '/CONTCAR ./')  # For reading structural info
        # Edit ssh command as above!
    else:
        doneOr = False
    return doneOr


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', dest='jobID', type=int)
    parser.add_argument('-i', dest='index', type=int)
    parser.add_argument('-f', dest='workingDir', type=str)
    args = parser.parse_args()

    isDone = checkStatus_remote(jobID=args.jobID, workingDir=args.workingDir, index=args.index)
    print('<CALLRESULT>')
    print((int(isDone)))
