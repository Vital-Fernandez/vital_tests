from os import environ, chdir
from subprocess import Popen, PIPE, STDOUT

def run_script(ScriptName, ScriptFolder, cloudy_address='/home/vital/Cloudy/source/cloudy.exe',
               bins_folder="/usr/sbin:/sbin:/home/vital/.my_bin:"):

    # Move to the script folder
    chdir(ScriptFolder)

    # Adding the cloudy path to the environment
    my_env = environ
    my_env["PATH"] = bins_folder + my_env["PATH"]  # This variable should be adapted to the computer

    # Preparing the command
    Command = 'cloudy {}'.format(ScriptName[0:ScriptName.rfind('.')]) # Script name without extension

    # Run the command
    print("--Launching command:")
    print(Command)
    p = Popen(Command, shell=True, stdout=PIPE, stderr=STDOUT, env=my_env)

    # Terminal output in terminal
    if len(p.stdout.readlines()) > 0:
        print('-- Code output wording\n')
        for line in p.stdout.readlines():
            print(line)

    return