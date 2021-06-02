import configparser

class Config:
    def __init__(self, config_file):
        self.config_file = config_file

    def return_config_file(self):
        return str(self.config_file)

    #Find values of the root
    def find_in_config(self, root, branch):
       config = configparser.ConfigParser()
       config.read(self.return_config_file())
       return (config[root][branch])

