#Input a list of files.  Output a string of commands to be passed to instance ID. 
import json
def RNA_SSM(file_list, s3_output):
	with open('config.json') as g:
		configVals = json.load(g)
		keyVar = configVals['awsCreds']['access_key']
		secretKeyVar = configVals['awsCreds']['aws_secret_key']
		regionVar = configVals['awsCreds']['region']
		bucket = configVals['S3bucket']



  ssm_command_list = 	[
            "sudo mkdir /home/ec2-user/bam_HLA/",
            "sudo chmod 777 /home/ec2-user/bam_HLA/",
            "aws configure set aws_access_key_id {}".format(keyVar),
            "aws configure set aws_secret_access_key {}".format(secretKeyVar),
            "aws configure set default.region {}".format(regionVar)
            ]

  ssm_processing = str()
  ssm_mid_commands = str()
  iteration = 1

  for iterBam in range(0,len(file_list)):

    aws_bam = file_list[iterBam]

    bamFile = aws_bam.split('/')[-1] 

    folderOutput = 'HLA_DNA_' + bamFile 

    bamBaiFile = bamFile + '.bai'

    chr6Trim = bamFile
    chr6Trim = chr6Trim.split('.')
    chr6Trim.insert(1,'chr6')
    chr6Trim = ".".join(chr6Trim)

    bamFastq = bamFile.replace('.bam','.fastq')

    iteration += 1

    bam_transfer = "aws s3 cp {} /home/ec2-user/bam_HLA/{}".format(aws_bam, bamFile)
    index = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /usr/local/bin/bin/samtools index /data/bam_HLA/{}".format(bamFile)
    trim6 = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 /usr/local/bin/bin/samtools view -bS /data/bam_HLA/{} 6 > /home/ec2-user/bam_HLA/{}".format(bamFile,chr6Trim)
    bam2fastq = "sudo docker run -v /home/ec2-user/:/data houghtos/immunotools:version4 bedtools bamtofastq -i /data/bam_HLA/{} -fq /data/bam_HLA/{}".format(chr6Trim, bamFastq)
    callOptitypeDocker = "sudo docker pull fred2/optitype"
    runHLA = "sudo docker run -v /home/ec2-user/bam_HLA/:/data/ -t fred2/optitype -i {} -d -o /data/{}/".format(bamFastq, folderOutput)

    rmBAM = 'sudo rm /home/ec2-user/bam_HLA/{}'.format(bamFile)
    rmFASATQ = 'sudo rm /home/ec2-user/bam_HLA/{}'.format(bamFastq)
    rmBAI = 'sudo rm /home/ec2-user/bam_HLA/{}'.format(bamBaiFile)
    rmCHR6 = 'sudo rm /home/ec2-user/bam_HLA/{}'.format(chr6Trim)

    copyFolder = 'aws s3 cp /home/ec2-user/bam_HLA/ s3://{}/{}/{}/ --recursive'.format(bucket,s3_output,runHLA)

    removeFastqFolder = "sudo rm -rf /home/ec2-user/bam_HLA/"
    
    mkfastq = "sudo mkdir /home/ec2-user/bam_HLA/"
    cmodfq = "sudo chmod 777 /home/ec2-user/bam_HLA/"


    ssm_command_list.append(bam_transfer)
    ssm_command_list.append(index)
    ssm_command_list.append(trim6)
    ssm_command_list.append(bam2fastq)
    ssm_command_list.append(callOptitypeDocker)
    ssm_command_list.append(runHLA)
    ssm_command_list.append(rmBAM)
    ssm_command_list.append(rmFASATQ)
    ssm_command_list.append(rmBAI)
    ssm_command_list.append(rmCHR6)
    ssm_command_list.append(copyFolder)
    ssm_command_list.append(removeFastqFolder)
    ssm_command_list.append(mkfastq)
    ssm_command_list.append(cmodfq)
	
	ssm_command_list.append("sudo shutdown -h now")

	return(ssm_command_list)
