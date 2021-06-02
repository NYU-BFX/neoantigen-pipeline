from Config import Config

class Align:
    def __init__(self, fastq2, genomeDir, ngs, path, prefix, reference, reference_gtf,  reference_path):
        self.fastq2         = fastq2
        self.genomeDir      = genomeDir
        self.ngs            = ngs
        self.path           = path
        self.prefix         = prefix
        self.reference      = reference
        self.reference_gtf = reference_gtf
        self.reference_path  = reference_path

    #Return FASTQ2
    def return_fastq2(self):
        return  str(self.fastq2)


    #Return which ngs: rna-seq or exome
    def return_ngs(self):
        return str(self.ngs)

    #Return path where STAR exists
    def return_path(self):
	    return str(self.path)


    #Return which human genome
    def return_reference(self):
        return str(self.reference)


    #Return path to human genome
    def return_reference_path(self):
        return str(self.reference_path)

    #Return genomeDir location
    def return_genomeDir(self):
        return str(self.genomeDir)


    #Return reference_gtf
    def return_reference_gtf(self):
       return str(self.reference_gtf)


    #Return prefix
    def return_prefix(self):
       return str(self.prefix)
    

    # Make the STAR command. Also the length of read is hardcoded 
    def make_command(self):
        cmd = 'STAR --runThreadN 16 --genomeDir ' + self.return_genomeDir() + \
        '  --sjdbGTFfile  '  +  self.return_reference_gtf() + \
        '  --sjdbOverhang 49 --readFilesIn  '  + self.return_fastq2() + \
        ' --quantMode TranscriptomeSAM --outFileNamePrefix '  +  self.return_prefix() 
        print ("Hello" +  cmd )



r1=Align("fastq2" , "./", "Helleoe2 ","RNA-seq", "hg10", "", "fhfhgd", "HEHED" )

r1.make_command()
