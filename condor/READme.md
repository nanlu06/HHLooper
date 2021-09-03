============================================

CONDOR scripts

===========================================

1. Before running condor, make sure to change the following items:

2 . give appropriate addresses in the run_myprog.sh, proto_condor_submit and makecondorsubmit.py

4. to initiate condor jobs, do python makeCondoSubmit.py

5. Also create the following directories before submitting jobs

   a. condor/condor_output/condor_logs

   b. condor/condor_submit

6. You can resubmit failed jobs by doing condor_submit condor_submit/submit_condor_job (appropriate job name)
7. Once all jobs complete successfully, do `python haddOutput.py`