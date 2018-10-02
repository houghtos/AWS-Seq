#Set variables on lines 14, 15, 15, 16, 66
#line 17 should include the hg19 reference files.


#Input a list of files.  Output a string of commands to be passed to instance ID. 
def RNA_SSM(file_list):
	if (len(file_list) % 2) == 0: 
		
		ssm_command_list = 	[
							"sudo mkdir /home/ec2-user/refs/",
							"sudo mkdir /home/ec2-user/fastq_RNA/",
							"sudo chmod 777 /home/ec2-user/refs/",
							"sudo chmod 777 /home/ec2-user/fastq_RNA/",
							"aws configure set aws_access_key_id YourKeyHere",
							"aws configure set aws_secret_access_key YourSecretKeyHere",
							"aws configure set default.region us-west-1",
							"aws s3 cp s3://YourBucket/hg19Refs/ /home/ec2-user/refs/ --recursive",
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
			star_output_prefix = star_output_prefix + "_ITER_" + str(iteration) + "."
			
			out_sam = star_output_prefix + "out.sam"
			aligned_out_sam = star_output_prefix + "Aligned.out.sorted.sam"
			aligned_out_bam = star_output_prefix + "Aligned.out.sorted.bam"
			gene_counts = star_output_prefix + "gene.counts"

			iteration += 1

			fastq1_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq1, fastq1)
			fastq2_transfer = "aws s3 cp {} /home/ec2-user/fastq_RNA/{}".format(aws_fastq2, fastq2)
			trimmomatic = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version2 java -jar /home/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} /data/fastq_RNA/{} ILLUMINACLIP:/home/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:6 MINLEN:28".format(fastq1,fastq2,fastq1_trimmed,fastq1_trimmed_unpaired,fastq2_trimmed,fastq2_trimmed_unpaired)
			star_align = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version2 /home/tools/STAR-2.6.0a/source/STAR --runThreadN 4 --genomeDir /data/refs/ --readFilesCommand zcat --readFilesIn /data/fastq_RNA/{} /data/fastq_RNA/{} --outReadsUnmapped Fastx --outFileNamePrefix {}".format(fastq1_trimmed,fastq2_trimmed,star_output_prefix)
			sam_sort = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version2 /usr/local/bin/bin/samtools sort -n {} > {}".format(out_sam,aligned_out_bam)
			sam_view = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version2 /usr/local/bin/bin/samtools view {} > {}".format(aligned_out_bam,aligned_out_sam)
			htseq = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version2 /usr/local/bin/htseq-count --stranded=no --order=name --idattr=gene_name --mode=intersection-nonempty {} /home/ec2-user/refs/Homo_sapiens.GRCh38.83.gtf > {}".format(aligned_out_sam,gene_counts)
			copy = "aws s3 cp /home/ec2-user/fastq_RNA/ s3://yourbucket/yourprefix/{}/ --recursive".format(star_output_prefix[:-1]) ############## <-  Set your output S3 bucket and prefix here
			remove_folder = "sudo rm -rf /home/ec2-user/fastq_RNA/"


			ssm_command_list.append(fastq1_transfer)
			ssm_command_list.append(fastq2_transfer)
			ssm_command_list.append(trimmomatic)
			ssm_command_list.append(star_align)
			ssm_command_list.append(sam_sort)
			ssm_command_list.append(sam_view)
			ssm_command_list.append(htseq)
			ssm_command_list.append(copy)
			ssm_command_list.append(remove_folder)
			



	else:
		pass

	
	ssm_command_list.append("sudo shutdown -h now")

	return(ssm_command_list)
