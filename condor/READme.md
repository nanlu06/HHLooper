============================================

CONDOR scripts

===========================================

1. Before running condor, make sure to change the following items:
   ```
   https://github.com/LPC-HH/HHLooper/blob/master/condor/makeCondorSubmit.py#L24-L36
   https://github.com/LPC-HH/HHLooper/blob/master/condor/makeCondorSubmit.py#L67
   https://github.com/LPC-HH/HHLooper/blob/master/condor/proto_condor_submit#L1
   https://github.com/LPC-HH/HHLooper/blob/master/condor/run_myprog.sh#L2
   https://github.com/LPC-HH/HHLooper/blob/master/condor/run_myprog.sh#L6-L8
   https://github.com/LPC-HH/HHLooper/blob/master/condor/haddOutput.py#L24-L34
   https://github.com/LPC-HH/HHLooper/blob/master/condor/haddOutput.py#L36 
   https://github.com/LPC-HH/HHLooper/blob/master/condor/haddOutput.py#L39-L40
   https://github.com/LPC-HH/HHLooper/blob/master/condor/haddOutput.py#L44-L47

and have appropriate addresses in the `run_myprog.sh`, `proto_condor_submit`, `makeCondorsubmit.py` and `haddOutput.py`

2. Also create the following directories before submitting jobs

   a. `condor/condor_output/condor_logs`

   b. `condor/condor_submit`
3. To initiate condor jobs, do `python makeCondoSubmit.py`
4. You can resubmit failed jobs by doing `condor_submit condor_submit/submit_condor_(appropriate job name)`. You can check the log files, out files and error files in `condor/condor_output/condor_logs`
5. Once all jobs complete successfully, do `python haddOutput.py`
