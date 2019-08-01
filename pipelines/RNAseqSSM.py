#Set variables on lines 25, 73, 74
#This example uses mm10 (mouse) reference sequence.

import json
def RNA_SSM(file_list, s3_output):
	with open('config.json') as g:
		configVals = json.load(g)
		keyVar = configVals['awsCreds']['access_key']
		secretKeyVar = configVals['awsCreds']['aws_secret_key']
		regionVar = configVals['awsCreds']['region']

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

			star_output_prefix = fastq1.split('.')[0]
			
			out_sam = star_output_prefix + ".Aligned.out.sam"
			aligned_out_sam = star_output_prefix + ".Aligned.out.sorted.sam"
			aligned_out_bam = star_output_prefix + ".Aligned.out.sorted.bam"
			gene_counts = star_output_prefix + ".gene.counts"

			iteration += 1

			fastq1_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq1, fastq1)
			fastq2_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq2, fastq2)
			trimmomatic = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 java -jar /home/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} ILLUMINACLIP:/home/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:6 MINLEN:28".format(fastq1,fastq2,fastq1_trimmed,fastq1_trimmed_unpaired,fastq2_trimmed,fastq2_trimmed_unpaired)
			star_align = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /home/tools/STAR-2.6.0a/source/STAR --runThreadN 4 --genomeDir /data/refs/ --readFilesCommand zcat --readFilesIn /data/fastq_RNA/{} /data/fastq_RNA/{} --outReadsUnmapped Fastx --outFileNamePrefix /data/output/{}".format(fastq1_trimmed,fastq2_trimmed,star_output_prefix)
			sam_sort = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /usr/local/bin/bin/samtools sort -n -@ 6 -m 4G /data/output/{} -o /data/output/{}".format(out_sam,aligned_out_bam)
			sam_view = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /usr/local/bin/bin/samtools view /data/output/{} -o /data/output/{}".format(aligned_out_bam,aligned_out_sam)
			htseq = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /usr/local/bin/htseq-count --stranded=no --order=name --idattr=gene_name --mode=intersection-nonempty /data/output/{} /data/refs/Mus_musculus.GRCm38.75.gtf > /home/ec2-user/output/{}".format(aligned_out_sam,gene_counts) #<<<<< The .GTF file here is set from the reference files (mm10 example here)
			copy = "aws s3 cp /home/ec2-user/output/ s3://yourbucket/{}/{}/ --recursive".format(s3_output,star_output_prefix[:-1]) #<<<<< Set your bucket name output here
			removeFastqFolder = "sudo rm -rf /home/ec2-user/fastq_RNA/"
			removeOutputFolder = "sudo rm -rf /home/ec2-user/output/"


			mkfastq = "sudo mkdir /home/ec2-user/fastq_RNA/"
			mkout = "sudo mkdir /home/ec2-user/output/"
			cmodfq = "sudo chmod 777 /home/ec2-user/fastq_RNA/"
			cmodout = "sudo chmod 777 /home/ec2-user/output/"


			ssm_command_list.append(fastq1_transfer)
			ssm_command_list.append(fastq2_transfer)
			ssm_command_list.append(trimmomatic)
			ssm_command_list.append(star_align)
			ssm_command_list.append(sam_sort)
			ssm_command_list.append(sam_view)
			ssm_command_list.append(htseq)
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
