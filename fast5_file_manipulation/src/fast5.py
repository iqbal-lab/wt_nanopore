import logging
import h5py

class fast5File(object):
    """
    Handler for ONT's fast5 files. 

    ...
    Parameters
    ----------
    filepath : str
        Path to the fast5file.


    Attributes
    ----------
    fq : str
        fasta string

    Methods
    -------

    """
    def __init__(self,filepath):
        self.filepath = filepath
        self.hdf = h5py.File(self.filepath, 'r')

    def get_fq(self,analysis_name="Basecall_2D_000",basecalled_name="BaseCalled_2D"):
        datasetname = '/Analyses/%s/%s/Fastq' % (analysis_name, basecalled_name)
        try:
            # logging.warning("Found Fasta File in %s : %s" % (self.filepath,datasetname ))
            fq = self.hdf[datasetname][()]
            return fq[:-1] + "\n"
        except KeyError,e:
            # logging.warning("No Fasta File in %s : %s" % (self.filepath,datasetname ))
            return None
    def close(self):
        self.hdf.close()

    @property 
    def twoD(self):
        """
        Returns the 2D basecalled FQ
        """
        return self.get_fq(basecalled_name='BaseCalled_2D')
    @property 
    def complement(self):
        """
        Returns the basecalled complement FQ
        """
        return self.get_fq(basecalled_name='BaseCalled_complement')
    @property 
    def template(self):
        """
        Returns the basecalled template FQ
        """
        return self.get_fq(basecalled_name='BaseCalled_template')

    @property 
    def analyses(self):
        return self.hdf['Analyses']
    @property 
    def has2d(self):
        try:
            return "BaseCalled_2D" in self.hdf['Analyses']['Basecall_2D_000']
        except KeyError:
            return False
    @property 
    def hasComplement(self):
        try:
            return "BaseCalled_complement" in self.hdf['Analyses']['Basecall_2D_000']
        except KeyError:
            return False

    @property 
    def hasTemplate(self):
        try:
            return "BaseCalled_template" in self.hdf['Analyses']['Basecall_2D_000']
        except KeyError:
            return False

    # @property 
    # def HairpinAlign(self):
    #     """"""
    #     return self.get_fq(basecalled_name='HairpinAlign')
