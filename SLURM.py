class SLURM:
    def __init__(self, command, memory, time):
        self.command = command
        self.memory  = memory
        self.time    = time


    #Make a command
    def make_slurm_command:
        cmd = 'sbatch  ' + self.return_script + ' --memory=' +  self.memory + \
              ' --time=' + self.time 
