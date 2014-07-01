## Module for reporting error assesment
import logging
import subprocess
import shlex
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROJECT_PATH = os.path.join(BASE_DIR,os.pardir)
PROJECT_PATH = os.path.abspath(PROJECT_PATH)
class Reporter(object):
	"""Generate an error report for ONT reads

	...
    Parameters
    ----------
    outfileDir : str
    	Output file directory 
    counter : object
    	Couter object
    firstRead : object
    	arbirtrary aligned read from (s/b)amfile
    latexTemplate : str
    	Path to LaTeX report template

	"""
	def __init__(self,ID,statsDir,outfileDir,stats={},firstRead=None,latexTemplate="../src/template.tex",):
		self.template = open(latexTemplate,'r')
		self.docString =  TexString(string=self.template.read())
		self.outfileDir = outfileDir
		self.read = firstRead
		self.latexWriten = False
		self.statsDir = statsDir
		self.imgDir = statsDir + 'img/'
		self.ID = ID
		self.stats = stats
		print self.stats

	def renderTemplate(self):
		## Metadata
		self.docString.replace(k='ID',v=self.ID)
		self.docString.replace(k='numReads',v=str(self.stats['numReads']),escape=False)
		self.docString.replace(k='numAlignedReads',v=str(self.stats['numAlignedReads']),escape=False)
		self.docString.replace(k='numBases',v=str(self.stats['numBases']),escape=False)
		self.docString.replace(k='numAlignedBases',v=str(self.stats['numAlignedBases']),escape=False)
		self.docString.replace(k='M',v=str(self.stats['M']),escape=False)
		self.docString.replace(k='D',v=str(self.stats['D']),escape=False)
		self.docString.replace(k='I',v=str(self.stats['I']),escape=False)
		self.docString.replace(k='S',v=str(self.stats['S']),escape=False)
		self.docString.replace(k='H',v=str(self.stats['H']),escape=False)
		self.docString.replace(k='longestRead',v=str(self.stats['longestRead']),escape=False)
		self.docString.replace(k='longestAlignedRead',v=str(self.stats['longestAlignedRead']),escape=False)
		## Figures

		self.docString.replace(k='histogram_of_read_length',v=self.imgDir + "histogram_of_read_length.png" ,escape=False)
		self.docString.replace(k='position_versus_coverage',v=self.imgDir + "position_versus_coverage.png" ,escape=False)
		self.docString.replace(k='read_length_vs_aligned_read_length',v=self.imgDir + "read_length_vs_aligned_read_length.png" ,escape=False)

	def writeLatex(self):
		self.renderTemplate()
		logging.info("Writing LaTeX to %s " % (self.outfileDir+'alignment_report.tex'))
		with open(self.outfileDir+'alignment_report.tex','w') as outfile:
			outfile.write(self.docString.text)
		self.latexWriten = True

	def generatePdfReport(self):
		if not self.latexWriten:
			self.writeLatex()
		cmd =  'pdflatex -output-directory=%s %alignment_report.tex  ' % (self.outfileDir,self.outfileDir)
		print cmd
		logging.info("Running %s" % cmd)
		proc=subprocess.Popen(shlex.split(cmd))
		proc.communicate()




class TexString(object):
	"""docstring for TexString"""
	def __init__(self, string):
		self.text = str(string)
		self._latex_special_chars = {
								    '&':  r'\&',
								    '%':  r'\%',
								    '$':  r'\$',
								    '#':  r'\#',
								    '_':  r'\_',
								    '{':  r'\{',
								    '}':  r'\}',
								    '~':  r'\lettertilde{}',
								    '^':  r'\letterhat{}',
								    '\\': r'\letterbackslash{}',
								    '\n': r'\\\\',
									}

	def replace(self,k, v,escape=True):
		if escape:
			v = self.escape(v)
		self.text = self.text.replace("{{%s}}"%k,v)

	def escape(self,s):
		return ''.join(self._latex_special_chars.get(c, c) for c in s)

		


		