#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2023.03.13
# @Author  : RichenLee 
# @Email   : lisheng@webmail.hzau.edu.cn



import logging
import argparse
import sys
import os
import gzip

__version__='1.0'
logging.basicConfig(
					 format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
					 datefmt='%a, %d %b %Y %H:%M:%S',
					 stream=sys.stderr,
					 filemode="w"
					 )
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

def get_logo():
	return (r'''
 ###########################      WELCOME     #########################################
                                                                                     
 多基因编辑拆分系统                                                                    
 MultipleGenesExtractor                                                               
 A program for cutting barcode in sequencing files                                
                                                                                     
                                                                                     
 Huazhong Agricultural University                                                    
 Version V1.0                                                                        
 LAST REVISED: 2023.03.13                                                            
 Usage python MultipleGenesExtractor.py -i R1.fq.gz,R2.fq.gz -b barcode.txt -o result 
    or python MultipleGenesExtractor.py -i R1.fq.gz -b barcode.txt                    
    or python GUI.py                                                                 
 ######################################################################################
''')

def get_header():
	term_width = 65
	logo = get_logo()
	description_str = logo+ "\n"
	description_str += ('[MultipleGenesExtractor version ' + __version__ + ']').center(term_width)+'\n'
	return description_str

def read_4lines(ff):
	ff1=ff.readline()
	if ff1:
		ff2=ff.readline()
		ff3=ff.readline()
		ff4=ff.readline()
		#print([ff1,ff2,ff3,ff4])
		return([ff1,ff2,ff3,ff4])
	else:
		return False

def get_sample(barcode_info_file):
	barcode_dic={}
	out_dic={}
	with open(barcode_info_file) as ff:
		for line in ff:
			sp=line.strip().split('\t')
			if len(sp)==3:
				sample=sp[0]
				F_barcode=sp[1]
				R_barcode=sp[2]
			elif len(sp)==2:
				sample=sp[0]
				F_barcode=sp[1]
				R_barcode=''				
			barcode_dic[F_barcode,R_barcode]=sample
			out_dic[sample]=[sample+'_R1.fq',sample+'_R2.fq']
	return barcode_dic,out_dic

def trim(content1,content2,barcode_dic,out_dic,out_path):
	write_content1=''
	write_content2=''
	read_sequenceR1 = content1[1].strip()
	read_sequenceR2 = content2[1].strip()
	#total_reads += 1
	for [F,R] in barcode_dic.keys():
		lst1=[read_sequenceR1.find(R),read_sequenceR2.find(F)]
		lst2=[read_sequenceR2.find(R),read_sequenceR1.find(F)]
		if (-1) not in lst1:
			content1[1]=read_sequenceR1[lst1[0]:]
			content2[1]=read_sequenceR2[lst1[1]:]
			write_content1=''.join(content1)
			write_content2=''.join(content2)
			with open(os.path.join(out_path,out_dic[barcode_dic[F,R]][0]),'at') as out:
				out.write(write_content1)
			with open(os.path.join(out_path,out_dic[barcode_dic[F,R]][1]),'at') as out:
				out.write(write_content2)			
			break
		elif (-1) not in lst2:
			content1[1]=read_sequenceR1[lst2[1]:]+'\n'
			content2[1]=read_sequenceR2[lst2[0]:]+'\n'
			write_content1=''.join(content1)
			write_content2=''.join(content2)
			with open(os.path.join(out_path,out_dic[barcode_dic[F,R]][0]),'at') as out:
				out.write(write_content1)
			with open(os.path.join(out_path,out_dic[barcode_dic[F,R]][1]),'at') as out:
				out.write(write_content2)
			break

def setcallback(x):
	with open('Result.txt','a+') as ff:
		ff.write(x)

def pairend_mode(barcode_dic,out_dic,fastqR1,fastqR2,out_path):
	for file in out_dic.values():
		out=gzip.open(os.path.join(out_path,file[0]),'wt')
		out.close()
		out=gzip.open(os.path.join(out_path,file[1]),'wt')
		out.close()

	infileR1=gzip.open(fastqR1, 'rt')
	infileR2=gzip.open(fastqR2, 'rt')
	total_reads=0
	while True:
		total_reads+=1
		content1=read_4lines(infileR1)
		content2=read_4lines(infileR2)			
		if not content1:
			break
		trim(content1,content2,barcode_dic,out_dic,out_path)
	infileR2.close()
	infileR1.close()

def singlend_mode(barcode_dic,out_dic,fastqR1,out_path):
	for file in out_dic.values():
		out=gzip.open(os.path.join(out_path,file[0]),'wt')
		out.close()
	temp_list=[]
	write_content=''
	n=0
	with gzip.open(fastqR1,'rt') as ff:
		for line in ff.readlines():
			n+=1
			temp_list.append(line)
			if n==4:
				read_sequence=temp_list[1].strip()
				for [F,R] in barcode_dic.keys():
					loc=read_sequence.find(F)
					if (-1) != loc:
						temp_list[1]=read_sequence[loc:]+'\n'
						write_content=''.join(temp_list)
						with open(os.path.join(out_path,out_dic[barcode_dic[F,R]][0]),'at') as out:
							out.write(write_content)
						break			
				n=0
				temp_list=[]

def main():
	print(get_header())
	parser = argparse.ArgumentParser(description='MultipleGenesExtractor Parameters')
	parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
	parser.add_argument('-i', '--input', type=str,  help='Sequencing fast file,pair-end files should be separated by commas', default='')
	parser.add_argument('-b', '--bar', type=str,  help='barcode file', default='')
	parser.add_argument('-o', '--output',  help='Output path', default='.')
	args = parser.parse_args()
	sys.stdout.flush() 
	if not args.input or not args.bar:
		parser.print_help()
		exit(1)
	barcode_file=args.bar
	info('Loading barcode file %s.'%barcode_file)
	try:
		barcode_dic,out_dic=get_sample(barcode_file)
	except:
		error('Wrong barcode file, please check!')
		exit(1)
	out_path=args.output
	if not os.path.exists(out_path):
		error('%s not exist,please check!'%out_path)
		exit(1)
	if len(args.input.split(','))==1:
		info('Input 1 file, Using sigle-end mode.')
		fastqR1=args.input
		info('Pocessing file %s.'%fastqR1)
		singlend_mode(barcode_dic,out_dic,fastqR1,out_path)
		info('Done!')
	elif len(args.input.split(','))==2:
		info('Input 2 files, Using pair-end mode.')
		fastqR1=args.input.split(',')[0]
		fastqR2=args.input.split(',')[1]
		info('Pocessing file %s,%s.'%(fastqR1,fastqR2))
		pairend_mode(barcode_dic,out_dic,fastqR1,fastqR2,out_path)

		info('Done!')


if __name__ == '__main__':
	main()