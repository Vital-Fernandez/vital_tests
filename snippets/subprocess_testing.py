import subprocess

# source ~/anaconda3/etc/profile.d/conda.sh
run_list_command = ['source /home/vital/anaconda3/bin/activate',
                    'conda activate iraf27',
                    'python /home/vital/PycharmProjects/vital_tests/snippets/matplotlib_plots.py',
                    'conda deactivate']
run_command = ' && '.join(run_list_command)
bash_command = f"bash -c '{run_command}'"

p1 = subprocess.run(bash_command, shell=True, capture_output=True)

# Print the output
print('- Stdout\n', p1.stdout.decode())

# Check for subprocess run
run_check = True if p1.returncode == 0 else False
print('-- Success subprocess' if run_check else '-- Failed subprocess')
if not run_check:
    print('- Stderr\n', p1.stderr.decode())


# import subprocess
#
# # source ~/anaconda3/etc/profile.d/conda.sh
# run_command = 'pwd'
# run_command += ' && ls'
# run_command += ' && which python'
# run_command += ' && source /home/vital/anaconda3/bin/activate'
# run_command += ' && conda activate iraf27'
# run_command += ' && which python'
# run_command += ' && conda deactivate'
# bash_command = f"bash -c '{run_command}'"
# p1 = subprocess.run(bash_command, shell=True, capture_output=True)
#
# # Print the output
# print('- Stdout\n', p1.stdout.decode())
#
# # Check for subprocess run
# run_check = True if p1.returncode == 0 else False
# print('-- Success subprocess' if run_check else '-- Failed subprocess')
# if not run_check:
#     print('- Stderr\n', p1.stderr.decode())

# https://stackoverflow.com/questions/48433478/how-to-activate-an-anaconda-environment-within-a-python-script-on-a-remote-machi

# https://stackoverflow.com/questions/7040592/calling-the-source-command-from-subprocess-popen