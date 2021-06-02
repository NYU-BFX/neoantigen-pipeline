import os
from Config import Config
from datetime import datetime

class arcasHLA:
    def __init__(self, BAM, config, log, prefix, outDir):
        self.BAM    =  BAM
        self.config = config
        self.log    =  log
        self.prefix =  prefix
        self.outDir = outDir



    def return_BAM(self):
        return str(self.BAM)


    def return_prefix(self):
        return str(self.prefix)


    def return_outDir(self):
        return str(self.outDir)


    def return_log(self):
        return  str(self.log)
 
    #If directory to run arcas does not exists, create it! 
    def exists_directory(self, directory):
        try:
            os.mkdir(directory)
        except:
            print ('Directory already exists ' +  directory)



    def arcasHLA_log(self, cmd):
        #Check if log directory exists
        logfile = self.return_log()
        try : 
            os.mkdir(logfile)
        except :
            print ("LOG exists")

       #Make arcaslog file
        try:
            with open(logfile + "/arcasHLA.log", "w") as outfile:
              
                ts = datetime.now()
                outfile.write(str(ts) + '>> Running ' +  cmd )
        except IOError:
            print ('Cannot open file ' +  outfile)


    #Return values from config file
    def return_config(self):
        config = Config(self.config)
        #Find where is arcasHLA location 
        arcas = config.find_in_config('arcasHLA',  'arcasHLA')
        return str(arcas)





    #First step of arcasHLA
    def test_arcas(self):
         
        cmd = self.return_config() + ' '   + self.return_BAM() +  ' -o ' + self.return_prefix() + ' --paired -t 8'
        print (cmd)
        arcasHLA_dir =  self.return_outDir() + '/' +  self.return_prefix()
        self.exists_directory(arcasHLA_dir)
        return (cmd)


    #Second step of arcasHLA
    def typing_arcas(self):
        cmd = self.return_config + '  genotype  ' +  self.return_fastq1()  + '  ' +  self.return_fastq2() + '  -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o '  + self.return_prefix() + '  -t 8 -v'
        print (cmd)
        

