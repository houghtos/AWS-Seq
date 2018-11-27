#!/usr/bin/python

from RNAseqSSM import RNA_SSM
from spotPriceModule import comparePrice 
import boto3
import os,sys
import re
import argparse
import csv
import datetime
from writeConfig import runConfig
import json

#Function to write the new Main.tf file in your terraform directory.
def writeTFmain(instanceNumbers, AMI_type, spotBid, volumeSize):
	with open('config.json') as g:
		configVals = json.load(g)
		regionConfig = configVals['awsCreds']['region']
		vpcConfig = configVals['awsEC2']['secGroup']
		
		pemConfig = configVals['awsEC2']['pemAddress']
		pemNameConfig = pemConfig.split('/')
		pemNameConfig = pemNameConfig[-1]
		pemNameConfig = pemNameConfig.replace('.pem','')

	f=open("tfSource.txt", "r")
	if f.mode == 'r':
		contents =f.read()
		contents = re.sub('COUNT', str(instanceNumbers), contents)
		contents = re.sub('BID///', str(spotBid), contents)
		contents = re.sub('TYPE///', str(AMI_type), contents)
		contents = re.sub('SIZE///', str(volumeSize), contents)
		contents = re.sub('REGION///', str(regionConfig), contents)
		contents = re.sub('VPC///', str(vpcConfig), contents)
		contents = re.sub('PEM_ADDRESS///', str(pemConfig), contents)
		contents = re.sub('PEM_NAME///', str(pemNameConfig), contents)

		#print(contents)
		f.close()

	w= open("Main.tf","w+")
	w.write(contents)
	w.close()

	TFvars= open("variables.tf","w+")
	TFvars.write(r'variable "access_key" {}' +'\n')
	TFvars.write(r'variable "secret_key" {}' +'\n')
	TFvars.write('variable "region" {default = "%s"}'% regionConfig)
	TFvars.close()

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
	with open('config.json') as f:
		configVals = json.load(f)
		print("******  Initializing terraform ********")
		os.system("{} init".format(configVals['OS']))
		print("******  Running terraform plan ********")
		os.system("{} plan".format(configVals['OS']))
		print("******  Starting Terraform Infastructure Build ********")
		os.system("{} apply -auto-approve".format(configVals['OS']))
	return()

#Shell command to destory all Terraform infastructure (Spot EC2 instances) 
def destroy():
	with open('config.json') as f:
		configVals = json.load(f)
		os.system("{} destroy -auto-approve".format(configVals['OS']))
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
		command_ID = item["CommandId"]
		instance_ID = item["InstanceId"]
		comments =  item["Comment"]
		request_DT = item["RequestedDateTime"]
		request_status = item["Status"]
		status_details = item["StatusDetails"]	

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
def sendSSMCommand(instance_id, SSM_command, s3_output):
	with open('config.json') as f:
		configVals = json.load(f)


		nowDateTime = datetime.datetime.now()
		nowDateTime = str(nowDateTime).replace(' ', '_')
		ssmComment = 'Paired FASTQ Pipeline' + 	nowDateTime    

		s3_output = s3_output + '.Logs.' + nowDateTime
		
		ssm_client = boto3.client('ssm')
		response = ssm_client.send_command(
				InstanceIds=[instance_id],
				DocumentName="AWS-RunShellScript",
				Comment=ssmComment, 
				Parameters={'commands': SSM_command, "executionTimeout":["36000"] }, #7200 
				OutputS3BucketName=configVals['S3bucket'], 
				OutputS3KeyPrefix=s3_output,
				OutputS3Region= configVals['awsCreds']['region'] ,
				)
	return()


def process_argument():
	### Create argument parser
	parser = argparse.ArgumentParser(description=r"Script initiate and send comands to AWS EC2 spot instance.  Can also review historical spot pricing.")

	### Create sub parser for each function menus
	subparsers = parser.add_subparsers(dest='options', help='choose script action')

	### Call routine to describe EC2 spot prices (Amazon LINUX/UNIX) using CLI.
	price_parser = subparsers.add_parser('price', help='historical spot pricing for instance')
	price_parser.add_argument('-t', '--instance_type', help='EC2 instance type (e.g. t3.2xlarge)')
	price_parser.add_argument('-m', '--months_delta', type=int, help='Number of months prior to examine EC2 spot costs for')

	### Call routine to start EC2 spot instance and run pipeline on given samples.
	tf_pasers = subparsers.add_parser('start', help='Start terraform and pass commands to each instance ID')
	tf_pasers.add_argument('-n', '--instance_numbers', type=int, help='Number of EC2 instances to spin up')
	tf_pasers.add_argument('-t', '--instance_type', help='EC2 instance type (e.g. t3.2xlarge) to spin up')
	tf_pasers.add_argument('-b', '--bid', help='Maximum bid (e.g. 20 cents is 0.20) for spot instances')
	tf_pasers.add_argument('-s', '--space', type=int, help='Amount of space (GiB) to allocate to each instance')
	tf_pasers.add_argument('-f', '--rna_AWS_file', help='Path to the .csv containg the addresses of your RNA fastq.  Seee readme on github (houghtos) for more information')
	tf_pasers.add_argument('-c', '--s3_output', help='S3 output for processed files and the SSM run log')

	### List all commands currently invoked
	command_parser = subparsers.add_parser('command', help='Shows all commands invoked on all EC2 spot instances created by terraform')

	### List all active instances
	instance_id_Parser = subparsers.add_parser('instances', help='Display list of all instance IDs for active instances')

	### Terraform Destroy without any prompt.
	dest = subparsers.add_parser('destroy', help='Destroys all Terraform infastructure')

	### Invokes routine to configure credentials for immunospace.
	immunoConfig = subparsers.add_parser('configure', help='Destroys all Terraform infastructure')

	return (parser.parse_args())


if __name__ == "__main__":

	
	args = process_argument()
    
    # Get arguments on main program
	if args.options == 'price':
		listPriceObject = comparePrice(args.instance_type, args.months_delta)
		listPriceObject.monthDelta() #Calculates prior months
		listPriceObject.describeSpotPrice()

	elif args.options == 'command':
		listReturn = getInstanceIDs()
		for instance_id in listReturn:
			returnCommandStatus(instance_id)

	elif args.options == 'instances':
		listReturn = getInstanceIDs()
		print(listReturn)

	elif args.options == 'start':

		f = open(args.rna_AWS_file, 'r') #rna_AWS_Files.csv
		cvs_reader = csv.reader(f)
		rna_list = []
		for row in cvs_reader:
			rna_list.append(row)
		f.close()

		writeTFmain(args.instance_numbers, args.instance_type, args.bid, args.space)

		tfSTart()
		listReturn = getInstanceIDs()
		if len(rna_list) == len(listReturn): 
			for instance_ID in range(len(rna_list)):
				RNA_Commands = RNA_SSM(rna_list[instance_ID], args.s3_output)
				print("Sending SSM command.")
				print("ID: " + listReturn[instance_ID])
				print("S3 arg " + args.s3_output)
				print(RNA_Commands)
				sendSSMCommand(listReturn[instance_ID],RNA_Commands, args.s3_output) 
				
		else:
			stop = 0
			while stop == 0:
				mismatch = input('Number of instances does not match number of sample groups inputted.  Do you wish to continue? (y/n)?  If no, all terraform infastrucutre will be destroyed.' + '\n')
				if mismatch.lower() in ('y','yes'):
					for instance_ID in range(len(rna_list)):
						RNA_Commands = RNA_SSM(rna_list[instance_ID], args.s3_output)
					try:
						sendSSMCommand(listReturn[instance_ID],RNA_Commands, args.s3_output) 
					except:
						print("SSM Command failed")
						pass
					stop = 1
				elif mismatch.lower() in ('n','no'):
					print("Shutting down all Terraform infastructure")
					destroy()
					stop = 1
				else:
					print("Please select yes/no or y/n.  Captilization does not matter.")


	elif args.options == 'destroy':
		destroy()

	elif args.options == 'configure':
		runConfig()

	else:
		print("Incorrect arguments given.  See argparse help")










