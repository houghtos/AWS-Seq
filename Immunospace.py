# To output SSM commands to AWS S3, uncomment lines 96-98 and fill in your relevant information.

#!/usr/bin/python
from MainTFSource import mainGen
from RNAseqSSM import RNA_SSM
from spotPriceModule import comparePrice 
import boto3
import os,sys
import re
import argparse
import csv

#Function to write the new Main.tf file in your terraform directory.
def writeTFmain(main):
	f= open("Main.tf","w+")
	f.write(main)
	f.close()
	return("**** Main TF file written *****")

#Regex method to obtain instance IDs.  Used to be compatible accross multiple OS systems.
#Can substitute w/ grep "spot_instance_id" terraform.tfstate
def getInstanceIDs():
	id_list = []
	file = open('terraform.tfstate', "r")
	for line in file:
		if re.search("spot_instance_id", line):
			line = line.replace(" ", "")
			line = line.replace('"', "")
			line = line.replace('spot_instance_id:', "")
			line = line.replace(',', "")
			line = line.replace('\n', "")
			id_list.append(line.replace(" ", ""))
	return(id_list)

#Shell commands to run Main.tf file from your Terraform installed directory.
def tfSTart():
	#print(os.getcwd())
	print("******  Initializing terraform ********")
	os.system("Terraform init")
	print("******  Running terraform plan ********")
	os.system("Terraform plan")
	print("******  Starting Terraform Infastructure Build ********")
	os.system("Terraform apply -auto-approve")
	return()

#Shell command to destory all Terraform infastructure (Spot EC2 instances) 
def destroy():
	os.system("Terraform destroy -auto-approve")
	return("Destroying all infastructure provisioned by terraform.")


#Checks what active instances are open that can be accessed by SSM and are properly configured.
def confirmSSM_Agent_Installed():
	os.system("aws ssm describe-instance-information --output text")
	return()

#Takes input of instance ID and returns command status for each one.  
def returnCommandStatus(instance_id):
	client = boto3.client('ssm')
	response = client.list_command_invocations(
	InstanceId=instance_id,
	MaxResults=50,
	)

	command_invocations = response["CommandInvocations"]
	for item in command_invocations:
		print("****** Commands invoked **********")
		print("CommandId: ",item["CommandId"])
		print("InstanceId ",item["InstanceId"])		
		print("InstanceName ",item["InstanceName"]) 	
		print("Comment ",item["Comment"])
		print("RequestedDateTime ",item["RequestedDateTime"])
		print("Status ",item["Status"])
		print("StatusDetails ",item["StatusDetails"])
		print()
	return()


#Requires input of instance ID <str> and SSM command <list>
#Should make input two lists: list of file names + instance ids.  Iterate over both.

def sendSSMCommand(instance_id, SSM_command):
	ssm_client = boto3.client('ssm')
	response = ssm_client.send_command(
			InstanceIds=[instance_id],
			DocumentName="AWS-RunShellScript",
			Comment='RNAseq aligning for paired fastq files', 
			Parameters={'commands': SSM_command }, 
			#OutputS3BucketName='MyBucket',  ########### <---- Specify your output AWS S3 bucket
			#OutputS3KeyPrefix='RNA_Files',            ########### <---- Specify your AWS S3 prefix
			#OutputS3Region='us-west-1',          ########### <---- Specify your AWS S3 region
			)
	return()


def process_argument():
	# Create argument parser
	parser = argparse.ArgumentParser(description=r"Script initiate and send comands to AWS EC2 spot instance.  Can also review historical spot pricing.")

	# Create sub parser for each function menus
	subparsers = parser.add_subparsers(dest='options', help='choose script action')

	#Call routine to describe EC2 spot prices (Amazon LINUX/UNIX) using CLI.
	price_parser = subparsers.add_parser('price', help='historical spot pricing for instance')
	price_parser.add_argument('-t', '--instance_type', help='EC2 instance type (e.g. t3.2xlarge)')
	price_parser.add_argument('-m', '--months_delta', type=int, help='Number of months prior to examine EC2 spot costs for')

	tf_pasers = subparsers.add_parser('start', help='Start terraform and pass commands to each instance ID')
	tf_pasers.add_argument('-n', '--instance_numbers', type=int, help='Number of EC2 instances to spin up')
	#tf_pasers.add_argument('-a', '--analysis_type', help='Analysis type to run... EG RNAseq paired end fastq')
	tf_pasers.add_argument('-t', '--instance_type', help='EC2 instance type (e.g. t3.2xlarge) to spin up')
	tf_pasers.add_argument('-b', '--bid', help='Maximum bid (e.g. 20 cents is 0.20) for spot instances')
	tf_pasers.add_argument('-s', '--space', type=int, help='Amount of space (GiB) to allocate to each instance')
	tf_pasers.add_argument('-f', '--rna_AWS_file', help='Path to the .csv containg the addresses of your RNA fastq.  Seee readme on github (houghtos) for more information')

	command_parser = subparsers.add_parser('command', help='Shows all commands invoked on all EC2 spot instances created by terraform')

	instance_id_Parser = subparsers.add_parser('instances', help='Display list of all instance IDs for active instances')

	dest = subparsers.add_parser('destroy', help='Destroys all Terraform infastructure')

	return (parser.parse_args())


if __name__ == "__main__":
	args = process_argument()
    
    # Get arguments on main program. Calculates prior months. Finally submits AWS CLI Describe command.
	if args.options == 'price':
		listPriceObject = comparePrice(args.instance_type, args.months_delta)
		listPriceObject.monthDelta()
		listPriceObject.describeSpotPrice()

	elif args.options == 'command':
		listReturn = getInstanceIDs()
		for instance_id in listReturn:
			returnCommandStatus(instance_id)

	elif args.options == 'instances':
		listReturn = getInstanceIDs()
		print(listReturn)

	elif args.options == 'start':

		f = open(args.rna_AWS_file, 'r') #rna_AWS_Files.csv kept in same directory as python file.
		cvs_reader = csv.reader(f)
		rna_list = []
		for row in cvs_reader:
			rna_list.append(row)
		f.close()

		writeTFmain(mainGen(args.instance_numbers, args.instance_type, args.bid, args.space))
		tfSTart()
		listReturn = getInstanceIDs()

		for instnace_ID in len(listReturn):
			RNA_Commands = RNA_SSM(rna_list[instnace_ID])
			sendSSMCommand(listReturn[instnace_ID],RNA_Commands) 

	elif args.options == 'destroy':
		destroy()
	
	else:
		print("Incorrect arguments given.  See argparse help")




