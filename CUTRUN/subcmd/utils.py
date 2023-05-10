#!/hpcf/apps/python/install/2.7.13/bin/python
import os
import paramiko
import pickle
import getpass
import sys
import uuid
import glob
import subprocess
import re
import datetime
import time
from difflib import SequenceMatcher
import pandas as pd
# useful sub functions
p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
username = getpass.getuser()
def dos2unix(file):
	os.system(p_dir+"../bin/dos2unix "+file)
def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
def multireplace(myString, myDict):
	## keywords format is {{string}}
	for k in myDict:
		myString = myString.replace("{{"+str(k)+"}}",str(myDict[k]))
	return myString
def set_label(command,jid):
	command2 = 'attr -s command -V "%s" %s > /dev/null 2>&1'%(command,jid)
	# print (command2)
	os.system(command2)
class NGS_pipeline:

	# private variable, not accessible by users
	outputs_dict = {}
	outputs_dict['log_files'] = []
	outputs_dict['peak_files'] = []
	outputs_dict['bam_files'] = []
	outputs_dict['figures_tables'] = []
	outputs_dict['bdg_files'] = []
	outputs_dict['bw_files'] = []
	outputs_dict['QC_files'] = []
	submission_id_dict = {}
	parameter_dict = {}
	attachments = []
	urlDict = {}
	fastq_input_list=[]				
	
	def __init__(self,args,logger,dep=1):
		self.args = args
		self.logger = logger
		self.parameter_dict['project_name']=args.subcmd
		self.dep = dep
		if args.short:
			self.Num_cores = "2"
		else:
			self.Num_cores = "6"
		# parameter_dict has 10 additional parameters that need to be specified for every run
		self.effectiveGenomeSize={}
		self.effectiveGenomeSize['hg19'] = '2451960000'
		self.effectiveGenomeSize['hg19_20copy'] = '2451960000'
		self.effectiveGenomeSize['t2t'] = '2451960000'
		self.effectiveGenomeSize['hg19_lenti'] = '2451960000'
		self.effectiveGenomeSize['hg19_dctcf'] = '2451960000'
		self.effectiveGenomeSize['hg38'] = '2782746084'
		self.effectiveGenomeSize['mm9'] = '1865500000'
		self.effectiveGenomeSize['mm10'] = '2650000000'
		# http://seqanswers.com/forums/showthread.php?t=4167
		# https://www.biostars.org/p/55271/
		set_label(" ".join(sys.argv),self.args.jid)
	def submit_array_job(self):
		# most jobs should be split and then submit
		header="""

#BSUB -P {{project_name}}
#BSUB -o {{output_message}}_%J_%I.out -e {{output_message}}_%J_%I.err
#BSUB -n {{number_cores}}
#BSUB -q {{queue}}
#BSUB -R "span[hosts=1] rusage[mem={{memory_request}}]"
#BSUB -J "{{job_id}}[1-{{number_lines}}]"
module purge
{{dependencies}}

id=$LSB_JOBINDEX
COL1=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $1}'`
COL2=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $2}'`
COL3=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $3}'`
COL4=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $4}'`
COL5=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $5}'`
COL6=`head -n $id {{sample_list}}|tail -n1|awk -F "\t" '{print $6}'`
LINE=`head -n $id {{sample_list}}|tail -n1`

{{commands}}

		"""
		if hasattr(self.args, 'debug'):
			if self.args.debug:
				self.submission_id_dict[self.parameter_dict['job_id']] = "123"
				return 1
		if self.args.short:
			self.parameter_dict['queue']="short"	
			if self.parameter_dict['job_id'] == "BWA":
				self.parameter_dict['memory_request']=4000
			else:
				self.parameter_dict['memory_request']=1500
		if ".bw" in self.parameter_dict['job_id']:
			self.parameter_dict['queue']="standard"	
		# prepare job lsf
		if self.args.priority:
			self.parameter_dict['queue'] = "priority"
		header = multireplace(header,self.parameter_dict)
		write_file(self.parameter_dict['job_script_file'],header)
		dos2unix(self.parameter_dict['job_script_file'])
		# submit job and extract submission id
		job_id_regex = b"Job <([0-9]+)>"
		pipe = subprocess.Popen('bsub < '+self.parameter_dict['job_script_file'], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
		output = pipe.communicate()[0]
		match = re.match(job_id_regex,output)
		try:
			submission_id = match.group(1).decode('utf-8')
		except:
			self.logger.error(self.parameter_dict['job_script_file']+' Job is failed to submit. ')
			self.logger.error('This failure is likely to be caused by wrong input format. Program continues.')
			self.logger.error('You will still receive an email when jobs are done.')
			# exit()

		# organize outputs
		# every lsf job will have:
		# - .lsf
		# - .err
		# - .out
		self.outputs_dict['log_files'].append(self.parameter_dict['job_script_file'])
		for i in range(1,int(self.parameter_dict['number_lines'])+1):
			self.outputs_dict['log_files'].append(self.parameter_dict['output_message']+"_"+submission_id+"_"+str(i)+".out")
			self.outputs_dict['log_files'].append(self.parameter_dict['output_message']+"_"+submission_id+"_"+str(i)+".err")
		self.logger.info(self.parameter_dict['job_id']+" has been submitted. Job ID: "+submission_id)
		self.submission_id_dict[self.parameter_dict['job_id']] = submission_id
		
		pass


	def run_fastqc(self):
		commands = "fastqc ${COL1}"

		# prepare input list, the LSF job is to go through every line
		fastqc_input = self.args.jid + ".fastqc.input"
		write_file(fastqc_input,"\n".join(self.fastq_input_list))
		dos2unix(fastqc_input)
		self.outputs_dict['log_files'].append(fastqc_input)

		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".fastqc.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=4000
		self.parameter_dict['job_id']="QC1"
		self.parameter_dict['sample_list']=fastqc_input
		self.parameter_dict['number_lines']=len(self.fastq_input_list)
		self.parameter_dict['dependencies']='module load fastqc/0.11.5'
		self.parameter_dict['commands']=commands
		self.parameter_dict['job_script_file']=self.args.jid +".fastqc.lsf"

		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_list:
			if "fastq.gz" in i:
				self.outputs_dict['QC_files'].append(i.replace(".fastq.gz","_fastqc.html"))
				self.outputs_dict['QC_files'].append(i.replace(".fastq.gz","_fastqc.zip"))
			else:
				self.outputs_dict['QC_files'].append(i.replace(".fastq","_fastqc.html"))
				self.outputs_dict['QC_files'].append(i.replace(".fastq","_fastqc.zip"))

		pass

	def submit_single_command(self,command,output_dict={},attachments=[],queue="standard",number_cores=1,memory_request=4000):

		# most jobs should be split and then submit
		header="""

#BSUB -P {{project_name}}
#BSUB -o {{output_message}}_%J_%I.out -e {{output_message}}_%J_%I.err
#BSUB -n {{number_cores}}
#BSUB -q {{queue}}
#BSUB -R "span[hosts=1] rusage[mem={{memory_request}}]"
#BSUB -J "{{job_id}}[1-{{number_lines}}]"

{{commands}}

		"""

		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +"."+self.args.subcmd+".message"
		self.parameter_dict['number_cores']=number_cores
		self.parameter_dict['queue']=queue
		self.parameter_dict['memory_request']=memory_request
		self.parameter_dict['job_id']=self.args.subcmd
		self.parameter_dict['number_lines']=1
		self.parameter_dict['commands']=command
		self.parameter_dict['job_script_file']=self.args.jid +"."+self.args.subcmd+".lsf"
		if self.args.short:
			self.parameter_dict['queue']="short"	
		# prepare job lsf
		header = multireplace(header,self.parameter_dict)
		write_file(self.parameter_dict['job_script_file'],header)
		dos2unix(self.parameter_dict['job_script_file'])
		# submit job and extract submission id
		job_id_regex = b"Job <([0-9]+)>"
		pipe = subprocess.Popen('bsub < '+self.parameter_dict['job_script_file'], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
		output = pipe.communicate()[0]
		match = re.match(job_id_regex,output)
		try:
			submission_id = match.group(1).decode('utf-8')
		except:
			self.logger.error(self.parameter_dict['job_script_file']+' Job is failed to submit. Exit program...')
			exit()

		# lsf outputs
		# every lsf job will have:
		# - .lsf
		# - .err
		# - .out
		self.outputs_dict['log_files'].append(self.parameter_dict['job_script_file'])
		for i in range(1,int(self.parameter_dict['number_lines'])+1):
			self.outputs_dict['log_files'].append(self.parameter_dict['output_message']+"_"+submission_id+"_"+str(i)+".out")
			self.outputs_dict['log_files'].append(self.parameter_dict['output_message']+"_"+submission_id+"_"+str(i)+".err")
		self.logger.info(self.parameter_dict['job_id']+" has been submitted. Job ID: "+submission_id)
		self.submission_id_dict[self.parameter_dict['job_id']] = submission_id

		# commmand-specific outputs
		for k in output_dict:
			try:
				self.outputs_dict[k] += output_dict[k]
			except:
				self.outputs_dict[k] = output_dict[k]
		self.attachments+=attachments

		pass			
	def run_skewer(self):
		Num_cores = "2"
		if self.args.single_end:
			commands =  p_dir+"../bin/skewer-0.2.2-linux-x86_64 -t {{Num_cores}} -z -x {{adaptor_x}} ${COL1} -o ${COL2}"
		else:
			commands = p_dir+"../bin/skewer-0.2.2-linux-x86_64 -t {{Num_cores}} -x {{adaptor_x}} -y {{adaptor_y}} ${COL1} ${COL2} -z -o ${COL3}"

		# If trimming adapters from Nextera runs should cut the reads at CTGTCTCTTATACACATCT instead of the usual AGATCGGAAGAGC. Use of cutadapt, trim_galore or similar program is recommended with custom adapter specified.
		# prepare input list, the LSF job is to go through every line
		commands = multireplace(commands, {'Num_cores':Num_cores,\
			'adaptor_x':self.args.adaptor_x,\
			'adaptor_y':self.args.adaptor_y})


		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".skewer.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=8000
		self.parameter_dict['job_id']="trim"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=len(self.fastq_input_lines)
		self.parameter_dict['dependencies']=''
		self.parameter_dict['commands']=commands
		self.parameter_dict['job_script_file']=self.args.jid +".skewer.lsf"

		# submit job
		self.submit_array_job()
		
		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['log_files'].append(i[-1]+"-trimmed.log")
	

		pass

	def run_BWA(self):

		Num_cores = self.Num_cores
		read1 = "${COL1}"
		read2 = "${COL2}"
		uid = "${COL3}"
		commands = []
		if self.args.index_file != p_dir+'../hg19/bwa_16a_index/hg19.fa':
			genome_index = self.args.index_file
		else:
			genome_index = p_dir+"../hg19/bwa_16a_index/hg19.fa".replace("hg19",str(self.args.genome).lower())
		
		command_read1 = "bwa mem -t "+Num_cores+" "+genome_index+" " + read1 + " |samtools view -@ "+Num_cores+" -bS - > " + uid + ".read1.bam"
		command_single = "bwa mem -t "+Num_cores+" "+genome_index+" " + read1 + " |samtools view -@ "+Num_cores+" -bS - > " + uid + ".bam"
		command_pair = "bwa mem -t "+Num_cores+" "+genome_index+" " + read1 + " " + read2 + " |samtools view -@ "+Num_cores+" -bS - > " + uid + ".bam"
		# command_pair = "bwa mem -C -t "+Num_cores+" "+genome_index+" " + read1 + " " + read2 + " | sed 's/1:N:0/BC:Z/g; s/2:N:0/BC:Z/g' |samtools view -@ "+Num_cores+" -bS - > " + uid + ".bam"
		bwa_input = self.args.jid + ".bwa.input"
		outlines = []
		if self.args.subcmd == "atac_seq":				
			for i in self.fastq_input_lines:
				if self.args.single_end:
					outlines.append("\t".join(map(lambda x:i[-1]+x,['-trimmed.fastq.gz',''])))
				else:
					outlines.append("\t".join(map(lambda x:i[-1]+x,['-trimmed-pair1.fastq.gz','-trimmed-pair2.fastq.gz',''])))
		else:				
			for i in self.fastq_input_lines:
				outlines.append("\t".join(i))
		write_file(bwa_input,"\n".join(outlines))
		dos2unix(bwa_input)
		self.outputs_dict['log_files'].append(bwa_input)


		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "chip_seq_pair":
			commands.append(command_pair)
			commands.append(command_read1)
		elif self.args.subcmd == "chip_seq_single":
			commands.append(command_single.replace("{COL3}","{COL2}"))
		elif self.args.subcmd == "atac_seq":
			if self.args.single_end:
				commands.append(command_single.replace("{COL3}","{COL2}"))
				commands.append("rm ${COL1}")
			else:
				commands.append(command_pair)
				commands.append("rm ${COL1}")
				commands.append("rm ${COL2}")
		else:
			commands.append(command_pair)
		# if self.args.subcmd == "atac_seq":	
			# commands.append("rm ${COL1}")
			# if not self.args.single_end:
				# commands.append("rm ${COL2}")


		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".BWA.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=6000
		self.parameter_dict['job_id']="BWA"
		self.parameter_dict['sample_list']=bwa_input
		self.parameter_dict['number_lines']=len(self.fastq_input_lines)
		# depends on trimming
		dependencies = ['module load bwa/0.7.16a','module load samtools/1.7']
		if self.args.subcmd == "atac_seq":
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['trim']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)


		# self.parameter_dict['dependencies']="\n".join(['module load bwa/0.7.16a','module load samtools/1.7'])
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".BWA.lsf"

		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['bam_files'].append(i[-1]+".bam")
			if self.args.subcmd == "chip_seq_pair":
				self.outputs_dict['bam_files'].append(i[-1]+".read1.bam")

		pass

	def run_samtools_markdup_rmdup(self):
		commands = []
		Num_cores = self.Num_cores
		command_samtools = "samtools command -@ "+str(Num_cores)
		command_stat = "samtools stats UID.bam > UID.bam.stat"
		command_sort_by_name = command_samtools.replace("command","sort") + " -n -o UID.name.st.bam UID.bam"
		command_fixmate = command_samtools.replace("command","fixmate") + " -m UID.name.st.bam UID.fixmate.bam"
		command_sort2 = command_samtools.replace("command","sort") + " -o UID.fixmate.st.bam UID.fixmate.bam"
		command_markdup = command_samtools.replace("command","markdup") + " UID.fixmate.st.bam UID.markdup.bam"
		command_markdup_r = command_samtools.replace("command","markdup") + " -r UID.markdup.bam UID.rmdup.bam"
		command_delete = "rm UID.name.st.bam UID.fixmate.bam UID.fixmate.st.bam"
		command_flagstat = command_samtools.replace("command","flagstat") + " UID.markdup.bam > UID.markdup.report"
		command_index_markdup = command_samtools.replace("command","index") + " UID.markdup.bam"
		command_unique = command_samtools.replace("command","view") + " -q 1 -b UID.markdup.bam > UID.markdup.uq.bam"
		command_index_unique = command_samtools.replace("command","index") + " UID.markdup.uq.bam"
		command_index_rmdup = command_samtools.replace("command","index") + " UID.rmdup.bam"
		command_unique_rmdup = command_samtools.replace("command","view") + " -q 1 -b UID.rmdup.bam > UID.rmdup.uq.bam"
		command_index_unique_rmdup = command_samtools.replace("command","index") + " UID.rmdup.uq.bam"
		command_chrm = command_samtools.replace("command","view") + " UID.markdup.bam chrM -b > UID.markdup.chrM.bam"
		command_flagstat_chrM = command_samtools.replace("command","flagstat") + " UID.markdup.chrM.bam > UID.markdup.chrM.report"
		command_rmchrM = "samtools idxstats UID.rmdup.uq.bam | cut -f 1 | grep -v  chrM| xargs samtools view -@ "+str(Num_cores)+\
		" -b UID.rmdup.uq.bam > UID.rmdup.uq.rmchrM.bam;samtools index UID.rmdup.uq.rmchrM.bam"
		command_rmchrM_raw = "samtools idxstats UID.markdup.bam | cut -f 1 | grep -v  chrM| xargs samtools view -@ "+str(Num_cores)+\
			" -b UID.markdup.bam > UID.markdup.rmchrM.bam;samtools index UID.markdup.rmchrM.bam"
		command_rmchrM_raw_flagstat = command_samtools.replace("command","flagstat") + " UID.markdup.rmchrM.bam > UID.markdup.rmchrM.report"
		
		commands.append(command_stat)
		commands.append(command_sort_by_name)
		commands.append(command_fixmate)
		commands.append(command_sort2)
		commands.append(command_markdup)
		commands.append(command_markdup_r)
		commands.append(command_delete)
		commands.append(command_flagstat)
		commands.append(command_index_markdup)
		commands.append(command_unique)
		commands.append(command_index_unique)
		
		if self.args.subcmd in ['cut_run','cut_run_histone','chip_seq_pair','atac_seq']:
			## add keep only properly paired:
			command_view_f3 = command_samtools.replace("command","view") + " -b -h -f 3 -F 4 -F 8 -o UID.rmdup.f3.bam UID.rmdup.bam; rm UID.rmdup.bam; mv UID.rmdup.f3.bam UID.rmdup.bam"
			commands.append(command_view_f3)
		commands.append(command_index_rmdup)
		commands.append(command_unique_rmdup)
		commands.append(command_index_unique_rmdup)
		commands.append(command_chrm)
		commands.append(command_flagstat_chrM)
		commands.append(command_rmchrM)
		commands.append(command_rmchrM_raw)
		commands.append(command_rmchrM_raw_flagstat)
		# for cut_run
		bam_to_bed_markdup = "bedtools bamtobed -bedpe -i UID.markdup.bam | awk '$1==$4 && $6-$2 < 1000 {print $0}' | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > UID.markdup.fragments.bed; bedtools genomecov -bg -i UID.markdup.fragments.bed -g %s > UID.markdup.fragments.bdg"%(self.args.chrom_size)
		bam_to_bed_rmdup = "bedtools bamtobed -bedpe -i UID.rmdup.bam | awk '$1==$4 && $6-$2 < 1000 {print $0}' | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > UID.rmdup.fragments.bed; bedtools genomecov -bg -i UID.rmdup.fragments.bed -g %s > UID.rmdup.fragments.bdg"%(self.args.chrom_size)
		bam_to_bed_rmdup_uq = "bedtools bamtobed -bedpe -i UID.rmdup.uq.bam | awk '$1==$4 && $6-$2 < 1000 {print $0}' | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > UID.rmdup.uq.fragments.bed; bedtools genomecov -bg -i UID.rmdup.uq.fragments.bed -g %s > UID.rmdup.uq.fragments.bdg"%(self.args.chrom_size)
		if "cut_run" in self.args.subcmd:
			commands.append(bam_to_bed_markdup)
			commands.append(bam_to_bed_rmdup)
			commands.append(bam_to_bed_rmdup_uq)
		if "single" in self.args.subcmd:
			uid = "${COL2}"
		else:
			uid = "${COL3}"
		try:
			if self.args.single_end:
				uid = "${COL2}"
		except:
			print ("uid",uid)
		myList1 = map(lambda x:x.replace("UID",uid),commands)
		myList2 = map(lambda x:x.replace("UID",uid+".read1"),commands)
		commands = myList1
		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "chip_seq_pair":
			commands = myList1 + myList2
		read_mapping_percent_load1 = "module load conda3"
		read_mapping_percent_load2 = "source activate /home/yli11/.conda/envs/rseqc"
		read_mapping_percent_markdup = "read_distribution.py -i UID.markdup.bam -r /home/yli11/Data/RSEQC_bed/%s_RefSeq.bed > UID.markdup.read_distribution.tsv"%(self.args.genome)
		read_mapping_percent_markdup_uq = "read_distribution.py -i UID.markdup.uq.bam -r /home/yli11/Data/RSEQC_bed/%s_RefSeq.bed > UID.markdup.uq.read_distribution.tsv"%(self.args.genome)
		read_mapping_percent_rmdup_uq = "read_distribution.py -i UID.rmdup.uq.bam -r /home/yli11/Data/RSEQC_bed/%s_RefSeq.bed > UID.rmdup.uq.read_distribution.tsv"%(self.args.genome)
		read_mapping_percent_rmdup = "read_distribution.py -i UID.rmdup.bam -r /home/yli11/Data/RSEQC_bed/%s_RefSeq.bed > UID.rmdup.read_distribution.tsv"%(self.args.genome)
		commands.append(read_mapping_percent_load1)
		commands.append(read_mapping_percent_load2)
		commands+= [x.replace("UID",uid) for x in [read_mapping_percent_markdup,read_mapping_percent_markdup_uq,read_mapping_percent_rmdup_uq,read_mapping_percent_rmdup]]
		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".markdup.rmdup.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['memory_request']=12000
		self.parameter_dict['job_id']="rmdup"
		self.parameter_dict['queue']="standard"
		self.parameter_dict['number_lines']=str(len(self.fastq_input_lines))
		self.parameter_dict['sample_list']=self.args.input
		# depends on BWA
		dependencies = ['module load samtools/1.7','module load bedtools']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['BWA']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".mk.rm.lsf"


		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.fragments.bdg")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.fragments.bed")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.bam.bai")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.chrM.bam*")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.fragments.bdg")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.fragments.bed")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.bam.bai")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.fragments.bdg")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.fragments.bed")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.uq.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.bam.bai")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.uq.bam.bai")
			self.outputs_dict['QC_files'].append(i[-1]+".markdup.report")
			self.outputs_dict['QC_files'].append(i[-1]+".markdup.chrM.report")
			self.outputs_dict['QC_files'].append(i[-1]+".bam.stat")
			self.outputs_dict['QC_files'].append(i[-1]+"*tsv")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.rmchrM.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.rmchrM.bam")
			self.outputs_dict['bam_files'].append(i[-1]+".rmdup.uq.rmchrM.bam.bai")
			self.outputs_dict['bam_files'].append(i[-1]+".markdup.rmchrM.bam.bai")
			self.outputs_dict['QC_files'].append(i[-1]+".markdup.rmchrM.report")

			if self.args.subcmd == "chip_seq_pair":
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.bam.bai")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.uq.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.uq.bam.bai")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.chrM.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.bam.bai")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.uq.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.uq.bam.bai")
				self.outputs_dict['QC_files'].append(i[-1]+".read1.markdup.report")
				self.outputs_dict['QC_files'].append(i[-1]+".read1.markdup.chrM.report")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.uq.rmchrM.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.rmdup.uq.rmchrM.bam.bai")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.rmchrM.bam")
				self.outputs_dict['bam_files'].append(i[-1]+".read1.markdup.rmchrM.bam.bai")
				self.outputs_dict['QC_files'].append(i[-1]+".read1.markdup.rmchrM.report")				
				self.outputs_dict['QC_files'].append(i[-1]+".read1.bam.stat")				
		pass

	def run_spp(self):
		Num_cores = 2
		command = "Rscript "+p_dir+"../bin/run_spp.R -c=UID.markdup.uq.bam -p="+str(Num_cores)+" -savp -odir=. -rf > UID.spp.log\nparse_spp_log.py UID.spp.log" 

		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "chip_seq_pair":
			command = command.replace("UID","${COL3}.read1")
		elif self.args.subcmd == "chip_seq_single":
			command = command.replace("UID","${COL2}")
		else:
			command = command.replace("UID","${COL3}")


		# define LSF parameter
		self.parameter_dict['output_message']=self.args.jid +".spp.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['queue']="standard"
		# self.parameter_dict['queue']="normal"
		self.parameter_dict['memory_request']="6000"
		self.parameter_dict['job_id']="spp"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=str(len(self.fastq_input_lines))
		# depends on rmdup
		dependencies = ['module load R/3.4.0','module load samtools/1.7']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['rmdup']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)

		self.parameter_dict['commands']=command
		self.parameter_dict['job_script_file']=self.args.jid +".spp.lsf"


		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_lines:
			if self.args.subcmd == "chip_seq_pair":
				self.outputs_dict['peak_files'].append(i[-1]+".read1.markdup.uq.pdf")
				self.attachments.append(i[-1]+".read1.markdup.uq.pdf")
				self.outputs_dict['peak_files'].append(i[-1]+".read1.spp.log")
				self.outputs_dict['peak_files'].append(i[-1]+".read1.spp.log_mqc.tsv")
			else:
				self.outputs_dict['peak_files'].append(i[-1]+".markdup.uq.pdf")
				self.attachments.append(i[-1]+".markdup.uq.pdf")
				self.outputs_dict['peak_files'].append(i[-1]+".spp.log")
				self.outputs_dict['peak_files'].append(i[-1]+".spp.log_mqc.tsv")
		pass

	def run_bamCoverage(self,input,output,type):
		Num_cores = self.Num_cores
		if self.args.effectiveGenomeSize != "2451960000":
			effectiveGenomeSize = self.args.effectiveGenomeSize
		else:
			effectiveGenomeSize = self.effectiveGenomeSize[str(self.args.genome).lower()]	
		
		try:
			command = "bamCoverage -b "+input+" -o "+output+" --smoothLength=200 --ignoreForNormalization chrX chrM   --effectiveGenomeSize "+effectiveGenomeSize+" --numberOfProcessors "+Num_cores
		except:
			command = "bamCoverage -b "+input+" -o "+output+" --normalizeUsing CPM --numberOfProcessors "+Num_cores									 
		# paired-end data will add --centerReads parameter
		if not "single" in self.args.subcmd:
			command += " --centerReads"
		else:
			if self.args.design_matrix != None:
				# sometimes we use chip-seq pipeline to just mapping reads and generate bw, this -e will affect output
				command += " -e 300"
		if "cut_run_histone" == self.args.subcmd:
			command = "bamCoverage -b "+input+" -o "+output+" --binSize 20 --smoothLength=50 -e --ignoreForNormalization chrX chrM   --effectiveGenomeSize "+effectiveGenomeSize+" --numberOfProcessors "+Num_cores
		commands=[command]
		# for normalization
		command = "bamCoverage -b "+input+" -o "+output.replace(".bw",".RPGC.norm.bw")+" --smoothLength=200 --binSize 20 --ignoreForNormalization chrX chrM   --effectiveGenomeSize "+effectiveGenomeSize+" --numberOfProcessors "+Num_cores+" --normalizeUsing RPGC"
		if not "single" in self.args.subcmd:
			command += " --centerReads"
		commands.append(command)
		# for stranded
		try:
			if self.args.stranded:
				command = "bamCoverage -of bedgraph --filterRNAstrand reverse -b "+input+" -o "+output.replace(".bw",".pos.RPGC.norm.bdg")+" --smoothLength=200 --binSize 20 --ignoreForNormalization chrX chrM   --effectiveGenomeSize "+effectiveGenomeSize+" --numberOfProcessors "+Num_cores+" --normalizeUsing RPGC"
				if not "single" in self.args.subcmd:
					command += " --centerReads"
				commands.append(command)
				command = "bamCoverage -of bedgraph --filterRNAstrand forward -b "+input+" -o "+output.replace(".bw",".neg.RPGC.norm.bdg")+" --smoothLength=200 --binSize 20 --ignoreForNormalization chrX chrM   --effectiveGenomeSize "+effectiveGenomeSize+" --numberOfProcessors "+Num_cores+" --normalizeUsing RPGC"
				if not "single" in self.args.subcmd:
					command += " --centerReads"
				commands.append(command)
		except:
			pass
		# define LSF parameter	
		self.parameter_dict['output_message']=self.args.jid +"."+type+".bw.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=4500
		self.parameter_dict['job_id']=type+".bw"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=str(len(self.fastq_input_lines))
		# depends on rmdup
		# dependencies = ['#BSUB -R "select[rhel7]"','module load python/2.7.15-rhel7'] # not working anymore, 11/14/2021
		dependencies = ['#BSUB -R "select[rhel7]"','module load conda3/202011','source activate /home/yli11/.conda/envs/captureC']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['rmdup']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +"."+type+".bw.lsf"

		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['bw_files'].append(i[-1]+"."+type+".bw")
			self.outputs_dict['bw_files'].append(i[-1]+".rmdup.neg.RPGC.norm.bw")
			self.outputs_dict['bw_files'].append(i[-1]+".rmdup.pos.RPGC.norm.bw")
			self.outputs_dict['bw_files'].append(i[-1]+"."+type+".*.bw")
			self.outputs_dict['bw_files'].append(i[-1]+"."+type+".*.bdg")
		pass

	def run_lib_complexity(self):
		if self.args.subcmd == "atac_seq":
			command_lib_complexity = """bedtools bamtobed -i UID.bam | awk 'BEGIN{OFS="\\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > UID.lib.complexity; echo -e "# plot_type: 'table'\n# section_name: 'Library Complexity'\nTotal Reads\\tDistinct Reads\\tOne Read\\tTwo Reads\\tNRF_ATAC\\tPBC1_ATAC\\tPBC2_ATAC" > UID.lib.complexity.tsv; cat UID.lib.complexity >> UID.lib.complexity.tsv """
		else:
			command_lib_complexity = """bedtools bamtobed -i UID.bam | awk 'BEGIN{OFS="\\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > UID.lib.complexity; echo -e "# plot_type: 'table'\n# section_name: 'Library Complexity'\nTotal Reads\\tDistinct Reads\\tOne Read\\tTwo Reads\\tNRF\\tPBC1\\tPBC2" > UID.lib.complexity.tsv; cat UID.lib.complexity >> UID.lib.complexity.tsv """

		self.parameter_dict['output_message']=self.args.jid +".lib.complexity.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=10000
		self.parameter_dict['job_id']="libx"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=str(len(self.fastq_input_lines))
		# depends on rmdup
		dependencies = ['module load bedtools/2.25.0']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['BWA']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		if "single" in self.args.subcmd:
			uid = "${COL2}"
		else:
			uid = "${COL3}"
		try:
			if self.args.single_end:
				uid = "${COL2}"
		except:
			print ("uid",uid)
		self.parameter_dict['commands']=command_lib_complexity.replace("UID",uid)
		self.parameter_dict['job_script_file']=self.args.jid +".libx.lsf"
		# submit job
		self.submit_array_job()

		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['QC_files'].append(i[-1]+".lib.complexity.tsv")
			self.outputs_dict['QC_files'].append(i[-1]+".lib.complexity")
		pass

	def run_macs2_narrowPeak(self):
		
		try:
			if self.args.design_matrix == None:
				no_control_flag = True # only true when design matrix is None
			else:
				no_control_flag = False
		except:
			no_control_flag = False
		commands = []

		commands.append("module load conda3")
		# can't resolve this problem
		# all(elementNROWS(gal) < 3)
		tmp_name = str(uuid.uuid4()).split("-")[-1]
		commands.append("source activate /home/yli11/.conda/envs/r_env")
		commands.append("samtools view -b treatment_uid.markdup.bam chr1 > treatment_uid.chr1.%s.bam;samtools index treatment_uid.chr1.%s.bam"%(tmp_name,tmp_name))
		commands.append("run_ATACseqQC.R treatment_uid.chr1.%s.bam %s"%(tmp_name,self.args.genome))
		commands.append("rm treatment_uid.chr1.%s.bam*"%(tmp_name))
		# for efficiency, only using chr1
		# commands.append("run_ATACseqQC.R treatment_uid.bam %s"%(self.args.genome))
		commands.append("module load macs2/2.1.1")
		# atac-seq parameter is different
		# paired-end and single-end are different
		# paired-end data will call twice, one for raw bam (atac-seq will be raw.rmchrM.bam), one for uq.rmdup.bam
		# single-end data will always use raw bam
		if self.args.subcmd == "atac_seq":
			# ref: https://www.biostars.org/p/209592/
			# ref: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
			macs2_pealcall = "macs2 callpeak --nomodel --shift -100 --extsize 200 -t treatment_uid.bam -B -n NAME"
		elif "single" in self.args.subcmd:
			macs2_pealcall = "macs2 callpeak -t treatment_uid.bam -c control_uid.bam --keep-dup all -n NAME -B"
		else:
			macs2_pealcall = "macs2 callpeak -t treatment_uid.bam -c control_uid.bam -f BAMPE -B --keep-dup all -n NAME"
		if self.args.genome == "mm9" or self.args.genome == "mm10":
			macs2_pealcall += " -g mm"
		try:
			if self.args.nomodel:
				macs2_pealcall += " --nomodel"
		except:
			nomodel = False
		# local bias track and call peak
		# ref: https://github.com/taoliu/MACS/wiki/Build-Signal-Track#Run_MACS2_bdgcmp_to_generate_foldenrichment_and_logLR_track
		# ref: https://github.com/taoliu/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands#Step_2_Decide_the_fragment_length_d
		macs2_bdgcmp_qpois = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m qpois -o NAME_treat_pvalue.bdg"
		macs2_bdgcmp_FE = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m FE -o NAME_FE.bdg"
		macs2_bdgcmp_logLR = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m logLR -o NAME_logLR.bdg -p 0.01"
		macs2_bdgpeakcall = "macs2 bdgpeakcall -i NAME_treat_pvalue.bdg -c 1.301 -l 100 -g 75 -o NAME_bdgpeaks.bed"
		in_files = ['NAME_treat_pvalue.bdg','NAME_FE.bdg','NAME_logLR.bdg']
		out_files = ['NAME_treat_pvalue.bw','NAME_FE.bw','NAME_logLR.bw']
		
		# have to manually adjust , because the short queue is rhel6, and by default, the node is rhel7.
		if self.args.chrom_size != p_dir+"../hg19/hg19.chrom.sizes":
			chrom_size = self.args.chrom_size
		else:
			chrom_size = p_dir+"../hg19/hg19.chrom.sizes".replace("hg19",str(self.args.genome).lower())
		if self.args.short:
			command_wigToBigWig = p_dir+"../bin/wigToBigWig INFILE " + chrom_size + " OUTFILE"
		else:
			command_wigToBigWig = p_dir+"../bin/wigToBigWig INFILE " + chrom_size + " OUTFILE"
		if self.args.Blacklist != p_dir+'../hg19/hg19.blacklist.bed':
			Blacklist = self.args.Blacklist
		else:
			Blacklist = p_dir+"../hg19/hg19.blacklist.bed".replace("hg19",str(self.args.genome).lower())			
			
		remove_backlist1 = "bedtools intersect -a NAME_peaks.narrowPeak -b "+Blacklist+" -v > NAME_peaks.rmblck.narrowPeak"
		remove_backlist2 = "bedtools intersect -a NAME_bdgpeaks.bed -b "+Blacklist+" -v > NAME_bdgpeaks.rmblck.bed"
		

		annot_feature = "annot_gene_features.py -f NAME_peaks.rmblck.narrowPeak -o NAME_peaks.rmblck.narrowPeak.annot.tsv -g %s"%(str(self.args.genome).lower())
		pie_chart = "pie_plot.py -f NAME_peaks.rmblck.narrowPeak.annot.tsv --order Exon,Promoter,5UTR,3UTR,Intron,Intergenic --use_col -1 --header -o NAME_peaks.rmblck.narrowPeak.genomic_features.pie_chart"
		if self.args.subcmd == "atac_seq":
			frip = "module load python/2.7.13;/home/yli11/HemTools/bin/calculate_FRiP.py treatment_uid.bam NAME_peaks.rmblck.narrowPeak ATAC"
		else:
			frip = "module load python/2.7.13;/home/yli11/HemTools/bin/calculate_FRiP.py treatment_uid.bam NAME_peaks.rmblck.narrowPeak CHIP"
		commands.append(macs2_pealcall)
		commands.append(macs2_bdgcmp_qpois)
		commands.append(macs2_bdgcmp_FE)
		commands.append(macs2_bdgcmp_logLR)
		commands.append(macs2_bdgpeakcall)
		# if self.args.genome == "mm9" or self.args.genome == "mm10":
			# commands = [x+" -g mm" for x in commands]
		commands.append(remove_backlist1)
		commands.append(remove_backlist2)
		commands.append(annot_feature)
		commands.append(pie_chart)
		commands.append(frip)
		
		cut_run_markdup = "SEACR_1.3.sh treatment_uid.markdup.fragments.bdg control_uid.markdup.fragments.bdg norm relaxed NAME.SEACR.markdup"
		cut_run_rmdup = "SEACR_1.3.sh treatment_uid.rmdup.fragments.bdg control_uid.rmdup.fragments.bdg norm relaxed NAME.SEACR.rmdup"
		cut_run_rmdup_uq = "SEACR_1.3.sh treatment_uid.rmdup.uq.fragments.bdg control_uid.rmdup.uq.fragments.bdg norm relaxed NAME.SEACR.rmdup.uq"
		cut_run_commands = []
		
		cut_run_commands.append(cut_run_markdup)
		cut_run_commands.append(cut_run_rmdup)
		cut_run_commands.append(cut_run_rmdup_uq)
		for i in [0,1,2]:
			tmp = command_wigToBigWig.replace("INFILE",in_files[i])
			tmp = tmp.replace("OUTFILE",out_files[i])
			commands.append(tmp)

		myList1 = map(lambda x:x.replace("treatment_uid","${COL1}").replace("control_uid","${COL2}").replace("NAME","${COL3}"),commands)	
		myList_atac = map(lambda x:x.replace("treatment_uid","${COL3}.markdup.rmchrM").replace("control_uid","${COL2}.markdup.rmchrM").replace("NAME","${COL3}.markdup.rmchrM"),commands)	
		myList2 = map(lambda x:x.replace("treatment_uid","${COL1}.rmdup.uq.rmchrM").replace("control_uid","${COL2}.rmdup.uq.rmchrM").replace("NAME","${COL3}.rmdup.uq.rmchrM"),commands)	
		myList4 = map(lambda x:x.replace("treatment_uid","${COL1}.markdup.uq").replace("control_uid","${COL2}.markdup.uq").replace("NAME","${COL3}.markdup.uq"),commands)	
		myList3 = map(lambda x:x.replace("treatment_uid","${COL1}.rmdup").replace("control_uid","${COL2}.rmdup").replace("NAME","${COL3}.rmdup"),commands)	

		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "atac_seq":
			commands = myList_atac + map(lambda x:x.replace("treatment_uid","${COL3}.rmdup.uq.rmchrM").replace("NAME","${COL3}.rmdup.uq.rmchrM"),commands)	+ map(lambda x:x.replace("treatment_uid","${COL3}.rmdup").replace("NAME","${COL3}.rmdup"),commands)
			if self.args.single_end:
				commands += map(lambda x:x.replace("treatment_uid","${COL3}.markdup.uq").replace("NAME","${COL3}.markdup.uq"),commands)	
				commands = [x.replace("${COL3}","${COL2}") for x in commands]
		elif no_control_flag:
			commands = [x.replace("-c control_uid.bam","") for x in commands]
			
			commands = map(lambda x:x.replace("treatment_uid","${COL3}").replace("NAME","${COL3}"),commands) + map(lambda x:x.replace("treatment_uid","${COL3}.rmdup.uq.rmchrM").replace("NAME","${COL3}.rmdup.uq.rmchrM"),commands)	+ map(lambda x:x.replace("treatment_uid","${COL3}.rmdup").replace("NAME","${COL3}.rmdup"),commands) +  map(lambda x:x.replace("treatment_uid","${COL3}.markdup.uq").replace("NAME","${COL3}.markdup.uq"),commands)	
			if "single" in self.args.subcmd:
				commands = [x.replace("${COL3}","${COL2}") for x in commands]
		else:
			commands = myList1 + myList2+myList3+myList4
		if "cut_run" in self.args.subcmd:
			commands += [x.replace("treatment_uid","${COL1}").replace("control_uid","${COL2}").replace("NAME","${COL3}") for x in cut_run_commands]
		# for stranded
		try:
			if self.args.stranded:
				macs2_bdgpeakcall = "module load macs2/2.1.1;macs2 bdgpeakcall --no-trackline -i ${COL1}.rmdup.pos.RPGC.norm.bdg -c 5 -l 100 -g 75 -o ${COL1}.rmdup.pos.RPGC.norm.peak.bed"
				remove_backlist = "bedtools intersect -a ${COL1}.rmdup.pos.RPGC.norm.peak.bed -b "+Blacklist+" -v > ${COL1}.rmdup.pos.RPGC.norm.rmblck.peak.bed"
				commands.append(macs2_bdgpeakcall)
				commands.append(remove_backlist)
				macs2_bdgpeakcall = "module load macs2/2.1.1;macs2 bdgpeakcall --no-trackline -i ${COL1}.rmdup.neg.RPGC.norm.bdg -c 5 -l 100 -g 75 -o ${COL1}.rmdup.neg.RPGC.norm.peak.bed"
				remove_backlist = "bedtools intersect -a ${COL1}.rmdup.neg.RPGC.norm.peak.bed -b "+Blacklist+" -v > ${COL1}.rmdup.neg.RPGC.norm.rmblck.peak.bed"
				commands.append(macs2_bdgpeakcall)
				commands.append(remove_backlist)
				tmp = command_wigToBigWig.replace("INFILE","${COL1}.rmdup.pos.RPGC.norm.bdg")
				tmp = tmp.replace("OUTFILE","${COL1}.rmdup.pos.RPGC.norm.bw")
				commands.append(tmp)
				tmp = command_wigToBigWig.replace("INFILE","${COL1}.rmdup.neg.RPGC.norm.bdg")
				tmp = tmp.replace("OUTFILE","${COL1}.rmdup.neg.RPGC.norm.bw")
				commands.append(tmp)
		except:
			pass
		# commands = [x.replace(".uq.rmchrM","") for x in commands] # just for Li's data
		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".macs2.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		# self.parameter_dict['queue']="normal"
		self.parameter_dict['memory_request']=80000
		self.parameter_dict['job_id']="macs2"
		try:
			if self.args.design_matrix != None:
				self.parameter_dict['sample_list']=self.args.design_matrix
				self.parameter_dict['number_lines']=len(self.design_matrix)		
			else:
				self.parameter_dict['sample_list']=self.args.input
				self.parameter_dict['number_lines']=len(self.fastq_input_lines)			
		except:
			self.parameter_dict['sample_list']=self.args.input
			self.parameter_dict['number_lines']=len(self.fastq_input_lines)
		
		# depends on rmdup
		dependencies = ['module load samtools/1.10','module load bedtools/2.29.2','module load conda3','source activate /home/yli11/.conda/envs/py2','module load gcc/6.3.0','export PATH=$PATH:"/home/yli11/HemTools/bin"','module load R/3.5.1']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['rmdup']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".macs2.lsf"

		# submit job
		self.submit_array_job()

		# organize output
		try:
			myList = self.design_matrix
		except:
			myList = self.fastq_input_lines
		for i in myList:
			self.outputs_dict['peak_files'].append(i[-1]+"*.bed")
			self.outputs_dict['peak_files'].append(i[-1]+"*.xls")
			self.outputs_dict['peak_files'].append(i[-1]+"*narrowPeak*")
			self.outputs_dict['peak_files'].append(i[0]+"*peak.bed*")
			self.outputs_dict['bdg_files'].append(i[-1]+"*.bdg")
			self.outputs_dict['QC_files'].append(i[-1]+"*.FRiP.tsv")
			self.outputs_dict['QC_files'].append(i[-1]+"*png")
			self.outputs_dict['QC_files'].append(i[0]+"*png")
			self.outputs_dict['QC_files'].append(i[1]+"*png")
			self.outputs_dict['QC_files'].append(i[-1]+"*.tsv")
			self.outputs_dict['QC_files'].append(i[0]+"*.tsv")
			self.outputs_dict['QC_files'].append(i[1]+"*.tsv")
			if self.args.subcmd == "atac_seq":

				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_treat_pvalue.bw")	
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_treat_pvalue.bw")	
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_treat_pvalue.bw")							
				# self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_FE.bw")
				# self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_logLR.bw")
				# self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_treat_pvalue.bw")							
			# elif "single" in self.args.subcmd:
				# self.outputs_dict['bw_files'].append(i[-1]+"_FE.bw")
				# self.outputs_dict['bw_files'].append(i[-1]+"_logLR.bw")
				# self.outputs_dict['bw_files'].append(i[-1]+"_treat_pvalue.bw")	
			else:
				self.outputs_dict['bw_files'].append(i[-1]+"_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_treat_pvalue.bw")			
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_treat_pvalue.bw")			
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.uq_treat_pvalue.bw")			
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup_treat_pvalue.bw")
			self.outputs_dict['log_files'].append(i[-1]+"*.r")
			self.outputs_dict['log_files'].append("*.auc")
			self.outputs_dict['log_files'].append("*.auc.bed")

		pass

	def run_macs2_broadPeak(self):
		commands = []
		# atac-seq parameter is different
		# paired-end and single-end are different
		# paired-end data will call twice, one for raw bam (atac-seq will be raw.rmchrM.bam), one for uq.rmdup.bam
		# single-end data will always use raw bam

		macs2_pealcall = "macs2 callpeak -t treatment_uid.bam -c control_uid.bam -f BAMPE -B --keep-dup all -n NAME --broad"
		# local bias track and call peak
		# ref: https://github.com/taoliu/MACS/wiki/Build-Signal-Track#Run_MACS2_bdgcmp_to_generate_foldenrichment_and_logLR_track
		# ref: https://github.com/taoliu/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands#Step_2_Decide_the_fragment_length_d
		macs2_bdgcmp_qpois = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m qpois -o NAME_treat_pvalue.bdg"
		macs2_bdgcmp_FE = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m FE -o NAME_FE.bdg"
		macs2_bdgcmp_logLR = "macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg -m logLR -o NAME_logLR.bdg -p 0.01"
		macs2_bdgpeakcall = "macs2 bdgbroadcall -i NAME_treat_pvalue.bdg -c 1.301 -l 100 -g 75 -o NAME_bdgpeaks.bed"
		in_files = ['NAME_treat_pvalue.bdg','NAME_FE.bdg','NAME_logLR.bdg']
		out_files = ['NAME_treat_pvalue.bw','NAME_FE.bw','NAME_logLR.bw']
		
		# have to manually adjust , because the short queue is rhel6, and by default, the node is rhel7.
		if self.args.chrom_size != p_dir+"../hg19/hg19.chrom.sizes":
			chrom_size = self.args.chrom_size
		else:
			chrom_size = p_dir+"../hg19/hg19.chrom.sizes".replace("hg19",str(self.args.genome).lower())
		if self.args.short:
			command_wigToBigWig = p_dir+"../bin/wigToBigWig INFILE " + chrom_size + " OUTFILE"
		else:
			command_wigToBigWig = p_dir+"../bin/wigToBigWig INFILE " + chrom_size + " OUTFILE"
		if self.args.Blacklist != p_dir+'../hg19/hg19.blacklist.bed':
			Blacklist = self.args.Blacklist
		else:
			Blacklist = p_dir+"../hg19/hg19.blacklist.bed".replace("hg19",str(self.args.genome).lower())			
			
		remove_backlist1 = "bedtools intersect -a NAME_peaks.broadPeak -b "+Blacklist+" -v -wa > NAME_peaks.rmblck.broadPeak"
		remove_backlist2 = "bedtools intersect -a NAME_bdgpeaks.bed -b "+Blacklist+" -v -wa > NAME_bdgpeaks.rmblck.bed"
		
		commands.append(macs2_pealcall)
		commands.append(macs2_bdgcmp_qpois)
		commands.append(macs2_bdgcmp_FE)
		commands.append(macs2_bdgcmp_logLR)
		commands.append(macs2_bdgpeakcall)
		commands.append(remove_backlist1)
		commands.append(remove_backlist2)
		for i in [0,1,2]:
			tmp = command_wigToBigWig.replace("INFILE",in_files[i])
			tmp = tmp.replace("OUTFILE",out_files[i])
			commands.append(tmp)

		myList1 = map(lambda x:x.replace("treatment_uid","${COL1}").replace("control_uid","${COL2}").replace("NAME","${COL3}"),commands)	
		myList_atac = map(lambda x:x.replace("treatment_uid","${COL3}.markdup.rmchrM").replace("control_uid","${COL2}.markdup.rmchrM").replace("NAME","${COL3}.markdup.rmchrM"),commands)	
		myList2 = map(lambda x:x.replace("treatment_uid","${COL1}.rmdup.uq.rmchrM").replace("control_uid","${COL2}.rmdup.uq.rmchrM").replace("NAME","${COL3}.rmdup.uq.rmchrM"),commands)	

		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "atac_seq":
			commands = myList_atac + map(lambda x:x.replace("treatment_uid","${COL3}.rmdup.uq.rmchrM").replace("NAME","${COL3}.rmdup.uq.rmchrM"),commands)	
		elif "single" in self.args.subcmd:
			commands = myList1
		else:
			commands = myList1 + myList2

		# define LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".macs2.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		# self.parameter_dict['queue']="normal"
		self.parameter_dict['memory_request']=30000
		self.parameter_dict['job_id']="macs2"
		try:
			self.parameter_dict['sample_list']=self.args.design_matrix
			self.parameter_dict['number_lines']=len(self.design_matrix)			
		except:
			self.parameter_dict['sample_list']=self.args.input
			self.parameter_dict['number_lines']=len(self.fastq_input_lines)
		# depends on rmdup
		dependencies = ['module load macs2/2.1.1','module load bedtools/2.25.0']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['rmdup']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".macs2.lsf"

		# submit job
		self.submit_array_job()

		# organize output
		try:
			myList = self.design_matrix
		except:
			myList = self.fastq_input_lines
		for i in myList:
			self.outputs_dict['peak_files'].append(i[-1]+"*.bed")
			self.outputs_dict['peak_files'].append(i[-1]+"*.xls")
			self.outputs_dict['peak_files'].append(i[-1]+"*.broadPeak")
			self.outputs_dict['bdg_files'].append(i[-1]+"*.bdg")
			if self.args.subcmd == "atac_seq":
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".markdup.rmchrM_treat_pvalue.bw")	
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_treat_pvalue.bw")							
			elif "single" in self.args.subcmd:
				self.outputs_dict['bw_files'].append(i[-1]+"_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_treat_pvalue.bw")	
			else:
				self.outputs_dict['bw_files'].append(i[-1]+"_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+"_treat_pvalue.bw")			
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_FE.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_logLR.bw")
				self.outputs_dict['bw_files'].append(i[-1]+".rmdup.uq.rmchrM_treat_pvalue.bw")			
			self.outputs_dict['log_files'].append(i[-1]+"*.r")

		pass


	def custom_kallisto_count(self):
		"""specifically designed for the insulator project"""
		pass
	def run_hommer(self):

		pass

	def run_mEpigram(self):

		pass

	def send_email_command2(self):
		# this function should be applied for any API
		command = 'echo "Hi User_name,\n\nYour JOB_ID is finished. Please see the attachment for a summary of results. The results have been generated the following directory:\n RESULT_DIR \n\n If no attachment can be found, it might be caused by an error. In such case, please go to the result directory (above) and do HemTools report_bug. \n\n BigWiggle_URL" | mailx {{attachments}} -s "JOB_ID is finished" -- User_name@stjude.org'
		command = command.replace("JOB_ID",self.args.jid)
		# print self.attachments
		# real_attachments = []
		# for item in self.attachments:
			# if os.path.isfile(item):
				# real_attachments.append(real_attachments)
		# attachments_string = map(lambda x: '-a "'+x+'" ',real_attachments)
		attachments_string = map(lambda x: '-a "'+x+'" ',self.attachments)
		command = command.replace("{{attachments}}","".join(attachments_string))
		# command = command.replace("attachment_file",self.attachments)
		command = command.replace("User_name",username)
		command = command.replace("RESULT_DIR",os.getcwd()+"/"+self.args.jid+"/")
		if len(self.urlDict) > 0:
			outlines = []
			for k in self.urlDict:
				outlines.append(k+" tracks can be viewed at: "+self.urlDict[k])
			command = command.replace("BigWiggle_URL","\n\n".join(outlines))
		else:
			command = command.replace("BigWiggle_URL","")
		return 	command
		
	def send_email_command(self):
		# this function should be applied for any API
		outFile = str(uuid.uuid4()).split("-")[-1]+".bwDict"
		out = open(outFile,"wb")
		pickle.dump([self.attachments,self.args.jid,self.urlDict],out)			
		out.close()
		commands = p_dir+'send_email.py '+outFile
		
		self.outputs_dict['log_files'].append(outFile)	
		return 	commands


		

	def gather_email_attachments(self):

		pass

	def QC_report_command(self):
		commands = p_dir+'NGS.report.py '+self.args.jid + " "+self.args.subcmd+" "+",".join(self.fastq_input_list) + " " + ",".join(list(map(lambda x:x[-1],self.fastq_input_lines)))
		if str(self.args.subcmd).lower() == "rna_seq":
			commands = p_dir+'FASTQC.report.py '+self.args.jid + " "+self.args.subcmd+" "+",".join(self.fastq_input_list) + " " + ",".join(list(map(lambda x:x[-1],self.fastq_input_lines)))				
		# step 4
		# organize output

		self.outputs_dict['figures_tables'].append(self.args.jid+".report.html")
		self.attachments.append(self.args.jid+".report.html")			

		return commands

		

	def run_output_organization(self,report_flag=True,email_flag=True):
		# this function should be the most generalized function
		# add send email, html report command here
		# Should be the last step of NGS_pipeline
		# print "report_flag",report_flag
		# print "last step"
		commands = ['mkdir '+self.args.jid]
		if report_flag:
			# commands = [self.QC_report_command(),self.send_email_command(),'mkdir '+self.args.jid]
			commands += [self.QC_report_command()]
			# merge multiqc table
			# FRIP
			commands.append("mqc_tsv_merge.py FRiP.tsv FRIP_mqc.tsv")
			self.outputs_dict['QC_files'].append("FRIP_mqc.tsv")
			# chip-seq Qtag
			commands.append("mqc_tsv_merge.py spp.log.tsv SPP_mqc.tsv")
			self.outputs_dict['QC_files'].append("SPP_mqc.tsv")
			# lib complexity
			commands.append("mqc_tsv_merge.py lib.complexity.tsv Lib.complexity_mqc.tsv")
			self.outputs_dict['QC_files'].append("Lib.complexity_mqc.tsv")
			# atac-seq tss enrichment
			commands.append("mqc_tsv_merge.py TSS_enrichment.tsv TSS_enrichment_mqc.tsv")
			self.outputs_dict['QC_files'].append("TSS_enrichment_mqc.tsv")
		# if not report_flag:
			# if email_flag:
				# commands = [self.send_email_command(),'mkdir '+self.args.jid]
				# commands = ['mkdir '+self.args.jid]
			# else:
				# commands = ['mkdir '+self.args.jid]

		# print "attachments",self.attachments
		commands.append("module load conda3")
		commands.append("source activate /home/yli11/.conda/envs/multiQC/")
		commands.append("export LC_ALL=en_US.utf-8")
		commands.append("export LANG=en_US.utf-8")
		# commands.append("cd "+self.args.jid)
		commands.append("cp /home/yli11/HemTools/share/NGS_pipeline/multiqc_config.yaml .")
		commands.append("sed -i 's/{{jid}}/%s/g' multiqc_config.yaml"%(self.args.subcmd))
		commands.append("multiqc .")
		# self.outputs_dict['QC_files'].append("TSS_enrichment_mqc.tsv")
		self.attachments.append('multiqc_report.html')
		if email_flag:
			commands.append(self.send_email_command())
			
		for k in self.outputs_dict:
			if len(self.outputs_dict[k]) == 0:
				continue
			commands.append('mkdir '+self.args.jid+"/"+k)
			for v in self.outputs_dict[k]:
				commands.append('mv '+v+" "+self.args.jid+"/"+k)

		commands.append("mv multiqc_report.html %s"%(self.args.jid))
		commands.append("mv multiqc_config.yaml %s"%(self.args.jid))
		commands.append("mv multiqc_data/ %s"%(self.args.jid))
		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['job_script_file']=self.args.jid +".html.lsf"
		commands.append('mv '+self.parameter_dict['job_script_file']+" "+self.args.jid+"/log_files")
		commands.append('mv '+self.args.jid+".log "+self.args.jid+"/log_files")
		commands.append('HemTools my_dir -a '+",".join([self.args.jid,os.getcwd(),self.args.subcmd]))						
		self.parameter_dict['output_message']=self.args.jid+"/log_files/"+self.args.jid +".html.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['memory_request']=8000
		self.parameter_dict['job_id']="html"
		self.parameter_dict['queue']="standard"
		self.parameter_dict['number_lines']=1
		try:
			self.parameter_dict['sample_list']=self.args.input
		except:
			self.parameter_dict['sample_list']=self.args.design_matrix
		# depends on rmdup
		dependencies = ['module load python/2.7.13']
		if self.dep:
			tmp = []
			for v in self.submission_id_dict.values():
				tmp.append('ended('+v+')')
			dependencies=['"'+" && ".join(tmp)+'"']+dependencies
		self.parameter_dict['dependencies']='#BSUB -w '+"\n".join(dependencies)

		self.parameter_dict['commands']="\n".join(commands)
		# self.parameter_dict['job_script_file']=self.args.jid +".html.lsf"

		# step 3
		# submit job
		self.submit_array_job()
		
		self.logger.info("All jobs are submitted. You will be notified by email with the final report attached! Bye!")
		self.logger.info("User input command: %s"%(" ".join(sys.argv)))

		
		pass
		
	def create_CRISPR_bw_files(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		commands = []
		RRA_files = []
		norm_count_files = []
		for k in self.pairwise_comparisons:
			RRA_files.append(k+"_RRA_results.sgrna_summary.txt")
			norm_count_files.append(k+"_RRA_results.normalized.txt")
		
		RRA_bw_command = p_dir+"create_mageck_RRA_LFC_FDR_bw_files.py " +\
			",".join(RRA_files) + \
			" " + self.args.bed + \
			" " + self.args.genome
		commands.append(RRA_bw_command)
		count_bw_command = p_dir+"create_mageck_count_bw_file.py " +\
			",".join(norm_count_files) + \
			" " + self.args.bed + \
			" " + self.args.genome
		commands.append(count_bw_command)
		commands.append("chmod 777 -R /home/yli11/.conda/envs/mageck/lib/R/")
		commands.append("chmod 777 -R /research/rgs01/home/clusterHome/yli11/.conda/envs/py2/lib/python3.7/site-packages/")
		commands.append("chmod a+rx -R /research/rgs01/home/clusterHome/yli11/.conda/envs/py2/")
		
			

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".CRISPR_bw.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=4000
		self.parameter_dict['job_id']="CRISPR_bw"
		self.parameter_dict['sample_list']=self.args.design_matrix
		self.parameter_dict['number_lines']=1
		# depends on rmdup
		dependencies = ['module load bedops','module load python/2.7.13']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['mageck_RRA']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".CRISPR_bw.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output
		for i in self.fastq_input_list:
			# self.outputs_dict['bw_files'].append(i+".crispr.norm_count.bw")
			self.outputs_dict['bw_files'].append(i+".norm_count.bw")
		for k in self.pairwise_comparisons:
			self.outputs_dict['bw_files'].append(k+"_LFC.crispr.bw")
			self.outputs_dict['bw_files'].append(k+"_logFDR.crispr.bw")

		pass				  

						

		  
																
	def run_kallisto(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		number_cores=6
		# command = """kallisto quant -i {{genome_index}} --threads={{number_cores}} --bootstrap-samples=100 --output-dir=${COL3}_kallisto --pseudobam ${COL1} ${COL2} | samtools view -Sb - > ${COL3}.bam \n mv ${COL3}_kallisto/abundance.tsv ${COL3}.abundance.tsv \n mv ${COL3}_kallisto/run_info.json ${COL3}.run_info.json \n rm -rf ${COL3}_kallisto"""
		command = """kallisto quant -i {{genome_index}} --threads={{number_cores}} --bootstrap-samples=100 --output-dir=${COL3}_kallisto ${COL1} ${COL2} \n mv ${COL3}_kallisto/abundance.tsv ${COL3}.abundance.tsv \n mv ${COL3}_kallisto/run_info.json ${COL3}.run_info.json \n rm -rf ${COL3}_kallisto"""
		genome_index = p_dir+"../hg19/kallisto/hg19.idx".replace("hg19",str(self.args.genome).lower())
		command = multireplace(command, {'genome_index':genome_index,\
			'number_cores':number_cores,\
		})		
		
			

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".kallisto.message"
		self.parameter_dict['number_cores']=number_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=2000
		self.parameter_dict['job_id']="kallisto"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=str(len(self.fastq_input_lines))
		# depends on rmdup
		dependencies = ['module load kallisto/0.43.1','module load samtools/1.7']
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']=command
		self.parameter_dict['job_script_file']=self.args.jid +".kallisto.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output
		self.outputs_dict['kallisto_files'] = []
		for i in self.fastq_input_lines:
			self.outputs_dict['log_files'].append(i[-1]+".run_info.json")
			self.outputs_dict['kallisto_files'].append(i[-1]+".abundance.tsv")
			# self.outputs_dict['bam_files'].append(i[-1]+".bam")

		self.submission_id_dict['BWA'] = self.submission_id_dict['kallisto']
		pass				  


															 
																			   
																		
	def run_combine_kallisto(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		number_cores=1
		command = p_dir+"""combine_kallisto.py {{exp_files}} {{labels}} {{gene_info}} {{output_csv_file}} \n zip {{output_csv_file}}.zip {{output_csv_file}}"""
		gene_info = p_dir+"../hg19/kallisto/hg19_gene_info.tsv".replace("hg19",str(self.args.genome).lower())
		output_csv_file = self.args.jid + ".transcript.tpm.csv"
		command = multireplace(command, {'gene_info':gene_info,\
			'output_csv_file':output_csv_file,\
			'exp_files':",".join(self.outputs_dict['kallisto_files']),\
			'labels':",".join(map(lambda x:x[-1],self.fastq_input_lines))\
		})		
		
			

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".combine.message"
		self.parameter_dict['number_cores']=number_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=2000
		self.parameter_dict['job_id']="combine"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=1
		# depends on rmdup
		dependencies = ['module load python/2.7.13']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['kallisto']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']=command
		self.parameter_dict['job_script_file']=self.args.jid +".combine.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output

		self.outputs_dict['kallisto_files'].append(output_csv_file)
		self.outputs_dict['kallisto_files'].append(output_csv_file+".zip")
		self.attachments.append(output_csv_file+".zip")
		pass				  

													  
									 
														   
															  

		

	def run_STJtracks2(self):
		# dump bw files
		## if a user submit the same subcmd multiple times, same bw_types will overwrite each other.
		## add a random string	, 4/11/2019	
		random_string = "."+str(uuid.uuid4()).split("-")[-1]
		myURL = "https://ppr.stjude.org/?study=HemPipelines/UserName/UserJobID/{{tracks}}.json"
		myURL = myURL.replace("UserName",username)
		myURL = myURL.replace("UserJobID",self.args.jid)
				
		myDict = {}
		myDict['treat_pvalue.bw'] = []
		myDict['FE.bw'] = []
		myDict['logLR.bw'] = []
		myDict['all.bw'] = []
		myDict['rmdup.bw'] = []
		myDict['rmdup.uq.bw'] = []
		myDict['crispr.bw'] = []
		bw_files = self.outputs_dict['bw_files']
		
		for k in myDict:

						
			for n in bw_files:
				if k in n:
					myDict[k].append(n)
					continue
			if len(myDict[k]) == 0:
				continue
			tmp = myURL.replace("{{tracks}}",k+random_string)
			self.urlDict[k] = tmp					

		outFile = str(uuid.uuid4()).split("-")[-1]+".bwDict"
		out = open(outFile,"wb")
		pickle.dump([self.args.jid,myDict,random_string],out)			
		out.close()

		commands = p_dir+'STJtracks.py '+outFile


		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".tracks.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['memory_request']=4000
		self.parameter_dict['job_id']="tracks"
		self.parameter_dict['queue']="standard"
		self.parameter_dict['number_lines']=1
		try:
			self.parameter_dict['sample_list']=self.args.design_matrix		
		except:
			self.parameter_dict['sample_list']=self.args.input													 
		# depends on rmdup
		dependencies = ['module load python/2.7.13']
		# if self.dep:
			# tmp = []
			# for v in self.submission_id_dict.values():
				# tmp.append('#BSUB -w "done('+v+')"')
			# dependencies=tmp+dependencies
		if self.dep:
			tmp = []
			for v in self.submission_id_dict.values():
				tmp.append('ended('+v+')')
			dependencies=['"'+" && ".join(tmp)+'"']+dependencies
		self.parameter_dict['dependencies']='#BSUB -w '+"\n".join(dependencies)			
			
		# self.parameter_dict['dependencies']="\n".join(dependencies)

		self.parameter_dict['commands']=commands
		self.parameter_dict['job_script_file']=self.args.jid +".tracks.lsf"

		# step 3
		# submit job
		self.submit_array_job()		

		# step 4
		# organize output

		self.outputs_dict['log_files'].append(outFile)			
		self.outputs_dict['log_files'].append(self.args.jid+".url")			

		return commands
		
	def run_STJtracks(self):
		# dump bw files
		# upload_tracks.py asds 1047946_Hudep1_CTCFIP.all.bw all
		## if a user submit the same subcmd multiple times, same bw_types will overwrite each other.
		## add a random string	, 4/11/2019	
		commands=[]
		myURL = "https://ppr.stjude.org/?study=HemPipelines/UserName/random_string/{{tracks}}.json"
		myURL = myURL.replace("UserName",username)

				
		myDict = {}
		myDict['treat_pvalue.bw'] = []
		myDict['FE.bw'] = []
		myDict['logLR.bw'] = []
		myDict['all.bw'] = []
		myDict['rmdup.bw'] = []
		myDict['DRIPc-seq'] = []
		myDict['rmdup.uq.bw'] = []
		myDict['crispr.bw'] = []
		bw_files = self.outputs_dict['bw_files']
		
		for k in myDict:			
			for n in bw_files:
				if k in n:
					myDict[k].append(n)
					continue
				if k =="DRIPc-seq":
					if "rmdup.neg.RPGC.norm.bw" in n:
						myDict[k].append(n)
						continue
					if "rmdup.pos.RPGC.norm.bw" in n:
						myDict[k].append(n)
						continue
			if len(myDict[k]) == 0:
				continue
			tmp = myURL.replace("{{tracks}}",k)		
			random_string = str(uuid.uuid4()).split("-")[-1]
			tmp = tmp.replace("random_string",random_string)
			try:
				commands.append(p_dir+'../bin/upload_tracks.py '+random_string+" "+",".join(myDict[k])+" "+k+" "+self.args.genome)
			except:
				commands.append(p_dir+'../bin/upload_tracks.py '+random_string+" "+",".join(myDict[k])+" "+k)
			self.urlDict[k] = tmp					

		
		# upload_tracks.py asds 1047946_Hudep1_CTCFIP.all.bw all

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".tracks.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['memory_request']=2000
		self.parameter_dict['job_id']="tracks"
		self.parameter_dict['queue']="standard"
		self.parameter_dict['number_lines']=1
		try:
			self.parameter_dict['sample_list']=self.args.design_matrix		
		except:
			self.parameter_dict['sample_list']=self.args.input													 
		# depends on rmdup
		dependencies = ['module load python/2.7.13']

		if self.dep:
			tmp = []
			for v in self.submission_id_dict.values():
				tmp.append('ended('+v+')')
			dependencies=['"'+" && ".join(tmp)+'"']+dependencies
		self.parameter_dict['dependencies']='#BSUB -w '+"\n".join(dependencies)			
			
		# self.parameter_dict['dependencies']="\n".join(dependencies)

		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".tracks.lsf"

		# step 3
		# submit job
		self.submit_array_job()		


		return commands

	def run_DESEQ2(self):

		pass

	def parse_paired_fastq(self):
		# In paired-end input tsv, only len(cols) == 3 will be accepted
		# len(cols) == 2 won't raise Error, but it will be ignored
		my_list = []
		my_paired_files = []
		dos2unix(self.args.input)
		with open(self.args.input) as f:
			for line in f:
				line = line.strip().split()
				if len(line) == 3:
					my_paired_files.append([line[0],line[1],line[2]])
					my_list.append(line[0])
					my_list.append(line[1])
		self.fastq_input_list = my_list	
		self.fastq_input_lines = my_paired_files	

		pass

	def parse_single_fastq(self):
		# In single-end input tsv, only len(cols) == 2 will be accepted
		# len(cols) == 3 won't raise Error, but it will be ignored
		my_list = []
		my_single_files = []
		dos2unix(self.args.input)
		with open(self.args.input) as f:
			for line in f:
				line = line.strip().split()
				if len(line) == 2:
					my_single_files.append([line[0],line[1]])
					my_list.append(line[0])
		self.fastq_input_list = my_list	
		self.fastq_input_lines = my_single_files	
		pass

	def parse_treatment_control(self):
		my_list = []
		dos2unix(self.args.design_matrix)
		with open(self.args.design_matrix) as f:
			for line in f:
				line = line.strip().split()
				if len(line) == 3:
					my_list.append([line[0],line[1],line[2]])
		self.design_matrix = my_list
		self.pairwise_comparisons={}
		try:
			self.args.input
		except:
			for row in my_list:
				if not self.pairwise_comparisons.has_key(row[0]):
					self.pairwise_comparisons[row[0]] = {}
					self.pairwise_comparisons[row[0]]['control'] = []
					self.pairwise_comparisons[row[0]]['treatment'] = []
				self.pairwise_comparisons[row[0]][row[1]]=row[2]
				self.fastq_input_list+=row[2].split(",")				
		# try:			
			# for row in my_list:
				# if not self.pairwise_comparisons.has_key(row[0]):
					# self.pairwise_comparisons[row[0]] = {}
					# self.pairwise_comparisons[row[0]]['control'] = []
					# self.pairwise_comparisons[row[0]]['treatment'] = []
				# self.pairwise_comparisons[row[0]][row[1]]=row[2]
				# self.fastq_input_list+=row[2].split(",")							
		# except:
			# return 1

	def dry_run(self):
		# subcmd independent check-up
		# only check for file existence and sample ID matches
		self.logger.info("DRY RUN: check if all input files exist")
		file_exist_flag = True
		for fname in self.fastq_input_list:
			if not os.path.isfile(fname):
				
				# self.logger.info(fname+" NOT FOUND. Please check your input file: "+self.args.input)
				self.logger.info(fname+" NOT FOUND. Please check your input file: ")
				file_exist_flag = False
		try:
			self.design_matrix
			flag = True
		except:
			flag = False
		try:
			self.fastq_input_lines
		except:
			flag = False			   
		if flag:
			self.logger.info("DRY RUN: check if sample IDs match in the input fastq tsv and design matrix.")
			reference = list(map(lambda x:x[-1],self.fastq_input_lines))
			sample_id_match_flag = True
			for i in self.design_matrix:
				if not i[0] in reference:
					sample_id_match_flag = False
					self.logger.info(i[0]+" can't be found in "+self.args.input)
				if not i[1] in reference:
					sample_id_match_flag = False
					self.logger.info(i[1]+" can't be found in "+self.args.input)
			good_flag = file_exist_flag and sample_id_match_flag
		good_flag = file_exist_flag
		if good_flag:
				self.logger.info("LOOKS GOOD! Dry run PASSED! Analyses continue...")
		else:
			self.logger.info("Something is wrong! Please check the above error messages.")
			# self.logger.close()
			os.system("rm "+self.args.jid+".log")		
			exit()
		

	def run_example(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		commands = ""

		# step 1.5 adjust commands to different subcmd
		## DESIGN: if new API uses BWA, this might need to be adjusted
		if self.args.subcmd == "chip_seq_pair":
			commands = myList1 + myList2		

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".BWA.message"
		self.parameter_dict['number_cores']=Num_cores
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=4000
		self.parameter_dict['job_id']="BWA"
		self.parameter_dict['sample_list']=self.args.input
		self.parameter_dict['number_lines']=len(self.fastq_input_lines)
		# depends on rmdup
		dependencies = ['module load bedtools/2.25.0']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['BWA']+'[*])"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".BWA.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output
		for i in self.fastq_input_lines:
			self.outputs_dict['bam_files'].append(i[-1]+".bam")
			if self.args.subcmd == "chip_seq_pair":
				self.outputs_dict['bam_files'].append(i[-1]+".read1.bam")			


		pass

	def run_mageck_count(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		# tmp_files = glob.glob("*.fastq.gz")+glob.glob("*.fq.gz")+glob.glob("*.fastq")
		files = []
		df = pd.read_csv(self.args.design_matrix,sep="\t",header=None,usecols=[0,1,2])
		for f in df[2]:
			files += f.split(",")
		files = list(set(files))
		# for f in tmp_files:
			# if "undetermined" in f.lower():
				# continue
			# files.append(f)	
		self.fastq_input_list = files	
		if self.args.kallisto:
			# convert sgRNA library to sgRNA fasta sequence
			sgRNA_table_to_fasta(self.args.gRNA_library,self.args.jid)
			# run kallisto
			if self.args.PE_data:
				commands = ["module load python/2.7.13","/home/yli11/HemTools/bin/custom_kallisto_count.py -f {0}.fa -o {0}.kallisto.tsv --paired ".format(self.args.jid)+" ".join(self.fastq_input_list)]
			else:
				commands = ["module load python/2.7.13","/home/yli11/HemTools/bin/custom_kallisto_count.py -f {0}.fa -o {0}.kallisto.tsv ".format(self.args.jid)+" ".join(self.fastq_input_list)]
			commands.append("mageck count -k {0}.kallisto.tsv --pdf-report --output-prefix {0}_raw_counts".format(self.args.jid))
			
		else:
			try:
				commands = "mageck count --output-prefix "+self.args.jid+"_raw_counts -l "+self.args.gRNA_library+" --fastq "+" ".join(files)+" --sample-label "+",".join(files) + " --pdf-report --norm-method control --control-sgrna "+self.args.control_gRNAs
			except:
				commands = "mageck count --output-prefix "+self.args.jid+"_raw_counts -l "+self.args.gRNA_library+" --fastq "+" ".join(files)+" --sample-label "+",".join(files) + " --pdf-report "
			# outputs
			if self.args.trim_5:
				commands += " --trim-5 "+ self.args.trim_5
			commands = [commands]
		# all.count_normalized.txt  all_countsummary.R  all_countsummary.Rnw  all.countsummary.txt  all.count.txt  all.log
		number_of_files = len(files)
		if number_of_files < 5:
			rmd_to_html = "/home/yli11/HemTools/bin/rmd_to_html.R %s_raw_counts.count_report.Rmd %s"%(self.args.jid,7)
		else:
			rmd_to_html = "/home/yli11/HemTools/bin/rmd_to_html.R %s_raw_counts.count_report.Rmd %s"%(self.args.jid,number_of_files+2)
		PCA_command = "module load conda3;source activate /home/yli11/.conda/envs/py2;module load gcc/6.3.0;/home/yli11/HemTools/bin/plot_PCA.py -f %s_raw_counts.count_normalized.txt --index --transpose --remove_cols Gene --transpose --log2_transform --header --dbscan_label -o %s_total_count_PCA_plot"%(self.args.jid,self.args.jid)
		UMAP_command = "module load conda3;source activate /home/yli11/.conda/envs/py2;module load gcc/6.3.0;/home/yli11/HemTools/bin/plot_PCA.py -f %s_raw_counts.count_normalized.txt --index --transpose --remove_cols Gene --transpose --log2_transform --header --dbscan_label --UMAP -o %s_total_count_UMAP_plot"%(self.args.jid,self.args.jid)
		
		
		commands.append(rmd_to_html)
		commands.append(PCA_command)
		commands.append(UMAP_command)

		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".mageck_count.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=8000
		self.parameter_dict['job_id']="mageck_count"
		self.parameter_dict['sample_list']=self.args.design_matrix
		self.parameter_dict['number_lines']=1
		# depends on rmdup
		# dependencies = ['module load python/2.7.15-rhel7']
		dependencies = ['rm -r /home/yli11/.conda/envs/mageck/lib/R/etc/','cp -r /home/yli11/HemTools/share/misc/etc /home/yli11/.conda/envs/mageck/lib/R/etc','chmod 777 -R  /home/yli11/.conda/envs/mageck/lib/R/etc','module load conda3/202011','source activate /home/yli11/.conda/envs/mageck/','module load texlive/20190410']
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="\n".join(commands)
		self.parameter_dict['job_script_file']=self.args.jid +".mageck_count.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output



		self.outputs_dict['log_files'].append(self.args.jid+"_raw_counts_countsummary.R")
		self.outputs_dict['log_files'].append(self.args.jid+"_raw_counts_countsummary.Rnw")
		self.outputs_dict['log_files'].append(self.args.jid+"_raw_counts.log")
		self.outputs_dict['mageck_count_files'] = []

		self.outputs_dict['mageck_count_files'].append(self.args.jid+".kallisto.tsv")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+".fa")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts.countsummary.txt")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts_countsummary.aux")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts_countsummary.pdf")
		self.attachments.append(self.args.jid+"_raw_counts_countsummary.pdf")
		# self.attachments.append(self.args.jid+".report.html")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts_countsummary.toc")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts_countsummary.tex")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts_countsummary.log")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts.count.txt")
		self.outputs_dict['mageck_count_files'].append("%s_total_count_PCA_plot.html"%(self.args.jid))
		self.outputs_dict['mageck_count_files'].append("%s_total_count_UMAP_plot_minDist_0.01_metric_euclidean_nNeighbor_15.html"%(self.args.jid))
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts.count_normalized.txt")
		self.outputs_dict['mageck_count_files'].append(self.args.jid+"_raw_counts.count_report.*")
		self.attachments.append(self.args.jid+"_raw_counts.count_report.html")
		# self.attachments.append("%s_total_count_PCA_plot.html"%(self.args.jid))
		self.attachments.append("%s_total_count_UMAP_plot_minDist_0.01_metric_euclidean_nNeighbor_15.html"%(self.args.jid))
		for f in files:
			label = f.split(".")[0]
			self.outputs_dict['mageck_count_files'].append(label+"*.count")


		pass		
	def run_mageck_RRA(self,paired_flag):

		# step 1
		# Note: all the commands should be directly executable in bash
		commands=[]
		try:
			command = "mageck test --normcounts-to-file --count-table "+self.args.jid+"_raw_counts.count.txt"+" --norm-method control --output-prefix {{name}}_RRA_results --control-sgrna "+self.args.control_gRNAs+" -t {{treatment_col_names}} -c {{control_col_names}} "+' --additional-rra-parameters " --permutation 10000" '
		except:
			command = "mageck test --normcounts-to-file --count-table "+self.args.jid+"_raw_counts.count.txt"+" --norm-method median --output-prefix {{name}}_RRA_results "+" -t {{treatment_col_names}} -c {{control_col_names}} "+' --additional-rra-parameters " --permutation 10000" '
		if paired_flag:
			command += " --paired "
		## PCA UMAP plot
		PCA_command = "source activate /home/yli11/.conda/envs/py2 \t module load gcc/6.3.0 \t plot_PCA.py -f {{name}}_RRA_results.normalized.txt -o {{name}}_PCA  --index --transpose --remove_cols Gene --transpose --log2_transform --header --kmeans_label 2 --text"
		UMAP_command = "plot_PCA.py -f {{name}}_RRA_results.normalized.txt -o {{name}}  --index --transpose --remove_cols Gene --transpose --log2_transform --header --kmeans_label 2 --text --UMAP"			
		for k in self.pairwise_comparisons:
			t = self.pairwise_comparisons[k]['treatment']
			c = self.pairwise_comparisons[k]['control']
			name = k
			tmp = multireplace(command, {'name':name,\
				'treatment_col_names':t,\
				'control_col_names':c\
				})
			mageckflute_command = "source activate /home/yli11/.conda/envs/R4 \t /home/yli11/HemTools/bin/run_mageck_flute_RRA.R {{name}}_RRA_results.gene_summary.txt {{name}}_RRA_results.sgrna_summary.txt {{name}} hsa"
			tmp2 = multireplace(mageckflute_command, {'name':name})			
			commands.append(tmp+"\t"+tmp2+"\t"+multireplace(PCA_command, {'name':name})+"\t"+multireplace(UMAP_command, {'name':name}))
			# commands.append(tmp2)

		
		
		# prepare input list, the LSF job is to go through every line
		command_input = self.args.jid + ".RRA.input"
		write_file(command_input,"\n".join(commands))
		dos2unix(command_input)
		self.outputs_dict['log_files'].append(command_input)

		# outputs
		
		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".mageck_RRA.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=10000
		self.parameter_dict['job_id']="mageck_RRA"
		self.parameter_dict['sample_list']=command_input
		self.parameter_dict['number_lines']=len(commands)
		# depends on mageck_count
		# dependencies = ['module load python/2.7.15-rhel7']
		dependencies = ['module load conda3','source activate /home/yli11/.conda/envs/mageck/','module load texlive/20190410']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['mageck_count']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		# self.parameter_dict['commands']="${LINE}"
		self.parameter_dict['commands']="eval ${COL1}\n${COL2}\n${COL3}\n${COL4}\n${COL5}\n${COL6}"
		self.parameter_dict['job_script_file']=self.args.jid +".mageck_RRA.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output
		self.outputs_dict['mageck_RRA_files'] = []
		for k in self.pairwise_comparisons:
			name = k
			self.outputs_dict['log_files'].append(name+"_RRA_results.log")
			self.outputs_dict['log_files'].append(name+"_RRA_results_summary.Rnw")
			self.outputs_dict['log_files'].append(name+"_RRA_results.R")
			self.outputs_dict['log_files'].append(name+"_RRA_results*Rmd")
			self.outputs_dict['mageck_RRA_files'].append(name+"_RRA_results.gene_summary.txt")
			self.outputs_dict['mageck_RRA_files'].append(name+"_RRA_results.sgrna_summary.txt")
			self.outputs_dict['mageck_RRA_files'].append(name+"_RRA_results.normalized.txt")
			# self.outputs_dict['mageck_RRA_files'].append(name+"_Flute_Results")
			self.outputs_dict['mageck_RRA_files'].append("MAGeCKFlute_"+name)
			self.outputs_dict['mageck_RRA_files'].append("hsa*png")
			self.outputs_dict['mageck_RRA_files'].append("%s_PCA.html"%(name))
			self.outputs_dict['mageck_RRA_files'].append("%s_minDist_0.01_metric_euclidean_nNeighbor_15.html.html"%(name))



		pass		

	def run_mageck_MLE(self):

		# step 1
		# Note: all the commands should be directly executable in bash
		commands=[]
		try:
			command = "mageck mle -d {{design_matrix}} --count-table "+self.args.jid+"_raw_counts.count.txt"+" --norm-method control --output-prefix {{name}}_MLE_results --control-sgrna "+self.args.control_gRNAs	+ " -i {{selected_samples}}" 
		except:
			command = "mageck mle -d {{design_matrix}} --count-table "+self.args.jid+"_raw_counts.count.txt"+" --norm-method median --output-prefix {{name}}_MLE_results "+ " -i {{selected_samples}}" 
			
		# to_mageck_design_matrix(treatment_list,control_list,count_df,comparison_id):

		for k in self.pairwise_comparisons:
			t = self.pairwise_comparisons[k]['treatment']
			c = self.pairwise_comparisons[k]['control']
			print (k,t,c)
			design_matrix = ['1,0']*len(c.split(","))+['1,1']*len(t.split(","))
			design_matrix = ";".join(design_matrix)
			
			name = k
			tmp = multireplace(command, {'name':name,\
				'design_matrix':design_matrix,\
				'selected_samples':t+","+c\
				})
			
			commands.append(tmp)

		# prepare input list, the LSF job is to go through every line
		command_input = self.args.jid + ".MLE.input"
		write_file(command_input,"\n".join(commands))
		dos2unix(command_input)
		self.outputs_dict['log_files'].append(command_input)
		
		# step 2 
		# define 10 LSF parameters
		self.parameter_dict['output_message']=self.args.jid +".mageck_MLE.message"
		self.parameter_dict['number_cores']=1
		self.parameter_dict['queue']="standard"
		self.parameter_dict['memory_request']=20000
		self.parameter_dict['job_id']="mageck_MLE"
		self.parameter_dict['sample_list']=command_input
		self.parameter_dict['number_lines']=len(commands)
		# depends on mageck_count
		# dependencies = ['module load python/2.7.15-rhel7']
		dependencies = ['module load conda3','source activate /home/yli11/.conda/envs/mageck/','module load texlive/20190410']
		if self.dep:
			dependencies=['#BSUB -w "ended('+self.submission_id_dict['mageck_count']+')"']+dependencies
		self.parameter_dict['dependencies']="\n".join(dependencies)
		self.parameter_dict['commands']="${LINE}"
		self.parameter_dict['job_script_file']=self.args.jid +".mageck_MLE.lsf"

		# step 3
		# submit job
		self.submit_array_job()

		# step 4
		# organize output
		# test_MLE.gene_summary.txt  test_MLE.log  test_MLE.sgrna_summary.txt

		self.outputs_dict['mageck_MLE_files'] = []
		for k in self.pairwise_comparisons:
			name = k
			self.outputs_dict['log_files'].append(name+"_MLE_results.log")
			self.outputs_dict['log_files'].append(name+"_MLE_results.sgrna_summary.txt")
			self.outputs_dict['mageck_MLE_files'].append(name+"_MLE_results.gene_summary.txt")




		pass

	# def run_fastqc(self):

	# 	pass



def run_upload_tracks(jid,bw_files,bw_types,random_string,dir=False,genome="hg19"):
	# username = getpass.getuser()
	file = open(p_dir+"../misc/AT847CE", 'rb')
	password = pickle.load(file)
	file.close()
	paramiko.util.log_to_file("/tmp/"+jid)
	print "connecting to server"
	ssh_client =paramiko.SSHClient()
	ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	hostname = "10.220.17.3"
	ssh_client.connect(hostname=hostname,username="yli11",password=password)

	if dir:
		print "creating user's dir"
		stdin,stdout,stderr=ssh_client.exec_command("mkdir /research/rgs01/resgen/legacy/gb_customTracks/tp/HemPipelines/"+username)
		print "wait 5"
		time.sleep(5)
		print "wait 5 done"
		stdin,stdout,stderr=ssh_client.exec_command("mkdir /research/rgs01/resgen/legacy/gb_customTracks/tp/HemPipelines/"+username+"/"+jid)
	user_dir = "/research/rgs01/resgen/legacy/gb_customTracks/tp/HemPipelines/"+username+"/"+jid+"/"
	time.sleep(5)
	print "generating json file"
	tracks_template = '{"type":"bigwig","scale":{"auto": 1},"file": "{{relative_path}}","stackheight":20,"stackspace":1,"onerow":1,"name":"{{track_name}}"}'
	tmp_file = "/tmp/"+str(uuid.uuid4()).split("-")[-1]
	tmp_file_handle = open(tmp_file,"wb")
	tracks_json_list=[]
	for b in bw_files:
		tmp = tracks_template.replace("{{relative_path}}","HemPipelines/"+username+"/"+jid+"/"+b)
		tmp = tmp.replace("{{track_name}}",b.replace("."+bw_types,""))
		tracks_json_list.append(tmp)
	lines = open(p_dir+"../misc/template2_browser.json").readlines()
	lines = "".join(lines)
	try:
		if genome == "hg19_20copy":
			genome="hgvirus"
		if genome == "t2t":
			genome="hg38"
		if genome == "hg19_lenti":
			genome="hg19"
		if genome == "hg19_dctcf":
			genome="hg19"
		lines = lines.replace("{{genome_version}}",genome)
	except:
		lines = lines.replace("{{genome_version}}","hg19")
	lines = lines.replace("{{tracks_json_list}}",",\n".join(tracks_json_list))	
	print >>tmp_file_handle,lines
	tmp_file_handle.close()
	print "transfering file"
	print "wait 5" 
	time.sleep(5)
	print "wait 5 done"
	ftp_client=ssh_client.open_sftp()
	time.sleep(5)
	for b in bw_files:
		
		if os.path.isfile(b):
			print "trasferring",b
			ftp_client.put(b,user_dir+b)
			print "wait 5" 
			time.sleep(5)
			print "wait 5 done"
		else:
			print b,"not exist"
	ftp_client.put(tmp_file,user_dir+bw_types+random_string+".json")
	ftp_client.close()
	ssh_client.close()
	myURL = "https://ppr.stjude.org/?study=HemPipelines/UserName/UserJobID/{{tracks}}.json"
	myURL = myURL.replace("UserName",username)
	myURL = myURL.replace("UserJobID",jid)
	myURL = myURL.replace("{{tracks}}",bw_types+random_string)
	return myURL	


def fastq_pair(x,y):
	length_difference = len(x)-len(y)
	# print x
	# print y
	item_differences = 0
	paired_flag = 0
	if length_difference!= 0:
		return 0
	for i in range(len(x)):
		item_x = x[i]
		item_y = y[i]
		# print item_x,item_y
		if item_x!=item_y:
			item_differences+=1
			if item_x == "R1" and item_y == "R2":
				paired_flag = 1
			if item_y == "R1" and item_x == "R2":
				paired_flag = 1
	# print "item_differences",item_differences
	# print "paired_flag",paired_flag
	if item_differences == 1 and paired_flag == 1:
		return 1
	else:
		return 0
def define_fastq_label(x):
	label = []
	for i in re.split('_|-|\.',x):
		if i == "R1" or i == "R2":
			return "_".join(label)
		label.append(i)
	label = "_".join(label)
	# print label
	label = label.replace("_fastq","").replace("_gz","")
	return label

def prepare_paired_end_input():
	
	# observation:
	# 1. seperators are: - _ .
	# 2. the first element is often times unique for each sample
	# 	- 1631306_RFA001_R1.fastq.gz
	# 	- 110C-GATA1_S10_R1.fastq.gz
	# 	- 1047946_Hudep1_CTCFIP_R1.fastq.gz
	# 3. paired-end data, the only difference is R1 and R2
	# 	- all other items are the same 
	# 	- also, after spliting the seperators, the length of the array should be the same
	# logic: 
	# 1. split the string
	# 2. infer putative pairs
	# 	rules:
	# 		- length is the same 
	# 		- only difference is R1 and R2


	tmp_files = glob.glob("*.fastq.gz")+glob.glob("*.fastq")
	files = []
	for f in tmp_files:
		if "undetermined" in f.lower():
			continue
		files.append(f)
	files_array = map(lambda x:re.split('_|-|\.',x),files)
	myDict = {}
	
	for i in range(len(files_array)):
		myDict[files[i]] = []
		for j in range(len(files_array)):
			if fastq_pair(files_array[i],files_array[j]):
				myDict[files[i]].append(files[j])
	flag = False
	fname = "fastq.tsv"
	if os.path.isfile(fname):
		print fname,"exists!"
		fname = "fastq."+str(uuid.uuid4()).split("-")[-1]+".tsv"
		print "Will use new file name:",fname
	f = open(fname,"wb")
	used_files = []
	for k in myDict:
		if len(myDict[k]) != 1:
			print "FILE:",k,"didn't find a pair",myDict[k]
			flag = 	True
		else:

			f1 = k
			f2 = myDict[k][0]
			if f1 in used_files:
				continue
			if f2 in used_files:
				continue			
			label = define_fastq_label(f1)
			if len(label) == 0:
				label = f1[:5]
			## 7-9-2019, rna-seq variant call, input has to be order 
			mySeps = re.split('_|-|\.',f1)
			# print f1
			# print mySeps
			skip_flag = True
			for x in mySeps:
				if x == "R1":
					skip_flag = False
			if skip_flag:
				continue
			print >>f,"\t".join([f1,f2,label])
			used_files.append(f1)
			used_files.append(f2)
	unused_files = list(set(files) - set(used_files))	
	f.close()
	flag = False
	if len(unused_files) == 0:
		print "Input fastq files preparation complete! ALL GOOD!"
		print "Please check if you like the computer-generated labels in :",fname
		flag = True
	else:
		print "Input fastq files preparation complete! There are some unmatched files."
		for f in unused_files:
			print f
	return flag,fname

def prepare_single_end_input():
	
	tmp_files = glob.glob("*.fastq.gz")+glob.glob("*.fastq")
	files = []
	for f in tmp_files:
		if "undetermined" in f.lower():
			continue
		files.append(f)	
	fname = "fastq.tsv"
	if os.path.isfile(fname):
		print fname,"exists!"
		fname = "fastq."+str(uuid.uuid4()).split("-")[-1]+".tsv"
		print "Will use new file name:",fname
	f = open(fname,"wb")
	for fastq in files:
		label = define_fastq_label(fastq)
		if len(label) == 0:
			label = fastq[:5]		
		print >>f,"\t".join([fastq,label])
	print "Input fastq files preparation complete! ALL GOOD!"
	f.close()
	return True,fname

def find_control2(myList):
	for x in myList:
		if "igg" in x.lower():
			return x
		if 'input' in x.lower():
			return x
	return "no_control_found"
def prepare_design_matrix2(file):
	# any item matched to "igg" or "input" will be used as control, all other samples will be treatment
	items = map(lambda x:x.strip().split()[-1],open(file).readlines())
	control = find_control(items)
	if control == "no_control_found":
		print "No control sample found."
		return
	fname = "peakcall.tsv"
	if os.path.isfile(fname):
		print fname,"exists!"
		fname = "peakcall."+str(uuid.uuid4()).split("-")[-1]+".tsv"
		print "Will use new file name:",fname
	f = open(fname,"wb")


	for i in items:
		if i == control:
			continue	
		print >>f,"\t".join([i,control,i+".vs."+control])
	print "Input peakcall file preparation complete! File name:",fname	
	f.close()
	pass


def find_control(myList):
	return_list = []
	for x in myList:
		if "igg" in x.lower():
			return_list.append(x)
		if 'input' in x.lower():
			return_list.append(x)
	return return_list
def prepare_design_matrix(file):
	# 4-10, adjust for multiple controls
	# any item matched to "igg" or "input" will be used as control, all other samples will be treatment
	items = map(lambda x:x.strip().split()[-1],open(file).readlines())
	control = find_control(items)
	print len(control)
	if len(control) == 0:
		print "No control sample found."
		return
	fname = "peakcall.tsv"
	if os.path.isfile(fname):
		print fname,"exists!"
		fname = "peakcall."+str(uuid.uuid4()).split("-")[-1]+".tsv"
		print "Will use new file name:",fname
	f = open(fname,"wb")

	if len(control) > 1:
		print "Multiple control files found. Computer-generated design matrix could be inaccurate."
		treatment_list = list(set(items)-set(control))
		df = pd.DataFrame()
		for c in control:
			myList = map(lambda x:similar(c, x),treatment_list)
			df[c] = myList
		df.index = treatment_list
		df['Max'] = df.idxmax(axis=1)
		for i in df.index:
			c = df.at[i,'Max']
			print >>f,"\t".join([i,c,i+".vs."+c])
	else:
		for i in items:
			if i == control[0]:
				continue	
			print >>f,"\t".join([i,control[0],i+".vs."+control[0]])
	print "Input peakcall file preparation complete! File name:",fname	
	f.close()
	pass
def similar(a, b):
	return SequenceMatcher(None, a, b).ratio()
from itertools import tee, count
from collections import Counter # Counter counts the number of occurrences of each item

def uniquify(seq, suffs = ('_%s'%(x) for x in range(1, 1000))):
    """Make all the items unique by adding a suffix (1, 2, etc).

    `seq` is mutable sequence of strings.
    `suffs` is an optional alternative suffix iterable.
    """
    not_unique = [k for k,v in Counter(seq).items() if v>1] # so we have: ['name', 'zip']
    # suffix generator dict - e.g., {'name': <my_gen>, 'zip': <my_gen>}
    suff_gens = dict(zip(not_unique, tee(suffs, len(not_unique))))  
    for idx,s in enumerate(seq):
        try:
            suffix = str(next(suff_gens[s]))
        except KeyError:
            # s was unique
            continue
        else:
            seq[idx] += suffix
def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()
def sgRNA_table_to_fasta(input,jid):
	df = pd.read_csv(input)
	df.columns = [0,1,2]
	mylist = df[0].tolist()
	uniquify(mylist)
	df[2] = mylist
	df = df.set_index(2)
	write_fasta("%s.fa"%(jid),df[1].to_dict())


