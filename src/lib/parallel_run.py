import os
import subprocess
import multiprocessing
from contextlib import redirect_stdout, redirect_stderr

import pathos.multiprocessing


def redirect_command(command, output_file):
    return f"{{ {command} ; }} >> '{output_file}' 2>&1"


def make_command_sequence(commands, log_file):
    commands.append("echo Finished successfully")
    command_sequence = " && ".join(commands)
    return redirect_command(command_sequence, log_file)


class SimpleRunner:
    def __init__(self, threads=120) -> None:
        self.threads = threads

    def add(self, command):
        print(f"Running command w/{self.threads} threads: {command}")
        try:
            p = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, executable="/bin/bash")
        except Exception as e:
            out = e.output
            assert False, f"Failed command: {command}\nError: {out}"


class LocalRunner:
    def __init__(self, log_file, threads=120):
        log_dir = os.path.dirname(os.path.abspath(log_file))
        assert os.path.exists(log_dir), f"{log_dir} folder doesn't exist"
        self.commands = []
        self.log_file = log_file
        self.threads = threads

    def add(self, command):
        self.commands.append(command)

    def run(self, dry=False, wait=True):
        if not self.commands:
            print("No commands to run")
            return

        p = None
        if dry:
            commands = '\n\t'.join(self.commands)
            print(f"Would run the following command sequence:\n\t{commands}\n")
        elif wait:
            command_sequence = make_command_sequence(self.commands, self.log_file)
            print(f"Running the following command sequence:\n{command_sequence}")
            try:
                p = subprocess.check_output(command_sequence, shell=True, stderr=subprocess.STDOUT,
                                            executable="/bin/bash")
            except Exception as e:
                out = e.output
                assert False, f"Failed command: {command_sequence}\nError: {out}"

            #return_code = os.system(command_sequence)
            #assert return_code == 0
        else:
            command_sequence = make_command_sequence(self.commands, self.log_file)            
            print(f"Running the following command sequence in parallel:\n{command_sequence}")
            p = subprocess.Popen(command_sequence, shell=True)

        self.commands = []
        return p


class PMACSRunner:
    def __init__(self, log_file, name="pybsub", threads=16, memory_gb=None, **kwargs):
        assert(os.path.exists(os.path.dirname(log_file)))
        
        self.log_file = log_file
        from bsub import bsub
        
        r = ""         
        if memory_gb is not None and memory_gb > 10:
          r += f"rusage [mem={memory_gb * 1024}]"
          
        if threads > 1:
          r += " span[hosts=1]"
        
        if r != "":
          kwargs["R"] = r
          
        if memory_gb is not None:
          kwargs["M"] = str(memory_gb * 1024)
        
        self.submit = bsub(name, n=str(threads), **kwargs)
        self.commands = []
        self.threads = threads

    def add(self, command):
        self.commands.append(command)

    def run(self, dry=False, wait=False):
        assert not wait
        command_sequence = make_command_sequence(self.commands, self.log_file)
        if dry:
            print(f"Would run the following command sequence:\n{command_sequence}\n")
        else:
            print(f"Submitting a job with the following command sequence:\n{command_sequence}")
            job = self.submit(command_sequence)
            print("Job id:", job.job_id, "\n")
            
        self.commands = []


def redirect_outputs(action, log_file):
    with open(log_file, "a") as f:
        with redirect_stderr(f):
            with redirect_stdout(f):
                return action()


def capture_stdout(command, suppress_stderr=False):
    res = subprocess.run(command, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL if suppress_stderr else None)

    return res.stdout.strip()


def start(action, log_file=None):
    if log_file is None:
        p = multiprocessing.Process(target=action)
    else:
        p = multiprocessing.Process(target=redirect_outputs, args=(action, log_file))

    p.start()
    return p


def smap(action):
    action()


def parallel_run(*actions, max_processes=20):
    with pathos.multiprocessing.ProcessPool(processes=max_processes) as pool:
        pool.map(smap, actions)


def wait_all(processes):
    for p in processes:
        p.join()
        