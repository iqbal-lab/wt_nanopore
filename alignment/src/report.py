## Module for reporting error assesment
import logging
import subprocess
import shlex
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROJECT_PATH = os.path.join(BASE_DIR,os.pardir)
PROJECT_PATH = os.path.abspath(PROJECT_PATH)
print PROJECT_PATH
class Reporter(object):
	"""Generate an error report for ONT reads

	...
    Parameters
    ----------
    opt : Optparse.opt
        Options from Optparse
    outfileDir : str
    	Output file directory 
    counter : object
    	Couter object
    firstRead : object
    	arbirtrary aligned read from (s/b)amfile
    latexTemplate : str
    	Path to LaTeX report template

	"""
	def __init__(self,opt,outfileDir,counter,firstRead=None,latexTemplate="../data/template.tex",):
		self.template = open(latexTemplate,'r')
		self.docString =  TexString(string=self.template.read())
		self.outfileDir = outfileDir
		self.opt = opt
		self.counter = counter
		self.read = firstRead

		self.latexWriten = False
		self.imgDir = PROJECT_PATH + '/proto_err/results/%s/img/' % (self.opt.runID  )
		

	def renderOptions(self):
		outstr  = ""
		d =  vars(self.opt)
		# for k,v in d.iteritems():
		# 	tmpStr = "".join([str(k),' : ',str(v),'\\'])
		# 	outstr = outstr.join(tmpStr)
		# print outstr
		outstr = str(d)
		return outstr

	# def addCaptionToTable(self,tableString,captionString):
	# 	caption = "\caption{%s}"
		


	def renderTemplate(self):
		## Metadata
		pass

	def writeLatex(self):
		self.renderTemplate()
		logging.info("Writing LaTeX to %s " % (self.outfileDir+'report.tex'))
		with open(self.outfileDir+'report.tex','w') as outfile:
			outfile.write(self.docString.text)
		self.latexWriten = True

	def generatePdfReport(self):
		if not self.latexWriten:
			self.writeLatex()
		cmd =  'pdflatex -output-directory=%s %sreport.tex  ' % (self.outfileDir,self.outfileDir)
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

		


		