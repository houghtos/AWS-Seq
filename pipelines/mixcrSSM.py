#Runs [MiXCR](https://mixcr.readthedocs.io/en/master/) on paired fastq files.
#Set variables on lines 26, 73 for desired reference files. 
#These files will be copied from an S3 directory.
#This example shows mm10 (mouse) reference sequence.

import json
def mixcrSSM(file_list, s3_output):
	with open('config.json') as g:
		configVals = json.load(g)
		keyVar = configVals['awsCreds']['access_key']
		secretKeyVar = configVals['awsCreds']['aws_secret_key']
		regionVar = configVals['awsCreds']['region']
		bucket = configVals['S3bucket']

	if (len(file_list) % 2) == 0: 
		
		ssm_command_list = 	[
							"sudo mkdir /home/ec2-user/refs/",
							"sudo mkdir /home/ec2-user/fastq_RNA/",
							"sudo mkdir /home/ec2-user/output/",
							"sudo chmod 777 /home/ec2-user/refs/",
							"sudo chmod 777 /home/ec2-user/fastq_RNA/",
							"sudo chmod 777 /home/ec2-user/output/",
							"aws configure set aws_access_key_id {}".format(keyVar),
							"aws configure set aws_secret_access_key {}".format(secretKeyVar),
							"aws configure set default.region {}".format(regionVar),
							"aws s3 cp s3://yourbucket/Refs/mm10/ /home/ec2-user/refs/ --recursive", # <<< Set your S3 address for reference genome.  In this example, we are using mm10.
							]

		ssm_processing = str()
		ssm_mid_commands = str()
		iteration = 1
		
		for pair1 in range(0,len(file_list),2):
			pair2 = pair1 + 1

			aws_fastq1 = file_list[pair1]
			aws_fastq2 = file_list[pair2]

			fastq1 = aws_fastq1.split('/')[-1]
			fastq2 = aws_fastq2.split('/')[-1]


			fastq1_trimmed = fastq1.split('.')
			fastq1_trimmed.insert(1,'trimmed')
			fastq1_trimmed = ".".join(fastq1_trimmed)
			
			fastq1_trimmed_unpaired = fastq1_trimmed.split('.')
			fastq1_trimmed_unpaired.insert(2,'unpaired')
			fastq1_trimmed_unpaired = ".".join(fastq1_trimmed_unpaired)

			fastq2_trimmed = fastq2.split('.')
			fastq2_trimmed.insert(1,'trimmed')
			fastq2_trimmed = ".".join(fastq2_trimmed)

			fastq2_trimmed_unpaired = fastq2_trimmed.split('.')
			fastq2_trimmed_unpaired.insert(2,'unpaired')
			fastq2_trimmed_unpaired = ".".join(fastq2_trimmed_unpaired)

			iteration += 1

			fastq1_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq1, fastq1)
			fastq2_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq2, fastq2)
      
			trimmomatic = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 java -jar /home/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} ILLUMINACLIP:/home/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:6 MINLEN:28".format(fastq1,fastq2,fastq1_trimmed,fastq1_trimmed_unpaired,fastq2_trimmed,fastq2_trimmed_unpaired)
			
			mixcr = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true -r log.v3.txt /data/fastq_RNA/{} /data/fastq_RNA/{} data/output/T_Cell.v3.alignments.vdjca".format(fastq1_trimmed, fastq2_trimmed)
			mixcr_step1 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr assemblePartial data/output/T_Cell.v3.alignments.vdjca data/output/T_Cell.v3.alignments_rescued_1.vdjca'
			mixcr_step2 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr assemblePartial data/output/T_Cell.v3.alignments_rescued_1.vdjca data/output/T_Cell.v3.alignments_rescued_2.vdjca'
			mixcr_step3 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr extend data/output/T_Cell.v3.alignments_rescued_2.vdjca data/output/T_Cell.v3.alignments_rescued_2_extended.vdjca'
			mixcr_step4 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr assemble data/output/T_Cell.v3.alignments_rescued_2_extended.vdjca data/output/T_Cell.v3.clones.clns'
			mixcr_step5 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr exportClones -c TRB data/output/T_Cell.v3.clones.clns data/output/T_Cell.v3.clones.TRB.txt'
			mixcr_step6 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr exportClones -c TRA data/output/T_Cell.v3.clones.clns data/output/T_Cell.v3.clones.TRA.txt'
			mixcr_step7 = 'sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version3 /home/tools/mixcr/mixcr exportClones -c IGH data/output/T_Cell.v3.clones.clns data/output/T_Cell.v3.clones.IGH.txt'

			copy = "aws s3 cp /home/ec2-user/output/ " + "s3://" + bucket + "/{}/{}/ --recursive".format(s3_output,star_output_prefix[:-1])
			removeFastqFolder = "sudo rm -rf /home/ec2-user/fastq_RNA/"
			removeOutputFolder = "sudo rm -rf /home/ec2-user/output/"


			mkfastq = "sudo mkdir /home/ec2-user/fastq_RNA/"
			mkout = "sudo mkdir /home/ec2-user/output/"
			cmodfq = "sudo chmod 777 /home/ec2-user/fastq_RNA/"
			cmodout = "sudo chmod 777 /home/ec2-user/output/"


			ssm_command_list.append(fastq1_transfer)
			ssm_command_list.append(fastq2_transfer)
			ssm_command_list.append(trimmomatic)
			ssm_command_list.append(mixcr)
			ssm_command_list.append(mixcr_step1)
			ssm_command_list.append(mixcr_step2)
			ssm_command_list.append(mixcr_step3)
			ssm_command_list.append(mixcr_step4)
			ssm_command_list.append(mixcr_step5)
			ssm_command_list.append(mixcr_step6)
			ssm_command_list.append(mixcr_step7)			
			ssm_command_list.append(copy)
			ssm_command_list.append(removeFastqFolder)
			ssm_command_list.append(removeOutputFolder)
			ssm_command_list.append(mkfastq)
			ssm_command_list.append(mkout)
			ssm_command_list.append(cmodfq)
			ssm_command_list.append(cmodout)



	else:
		print("Non-even number of fastq pairs submitted.  Since this processes paired FASTQ files, there must be an even number of files submitted to an instance.  Please destroy infastrucutre and review your csv input.")
		pass

	
	ssm_command_list.append("sudo shutdown -h now")

	return(ssm_command_list)
