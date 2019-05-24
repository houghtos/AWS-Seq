#!/usr/bin/python


#Command line prompt (similar to AWS/Google CLI) for writing and storing AWS credentials in config.json
#Opens and writes credentials to config.json. 
#Command line questions (similar to AWS/Google CLI) for configuring variables.
#coding utf-8

def runConfig():

	import json
	import io


	with open('config.json') as f:
		priorConfig = json.load(f)
		tempDat = priorConfig['awsCreds']['access_key'] 
		if tempDat == '':
			aws_key = input('AWS Access Key: ')
		else:
			tempDat = tempDat[-4:]
			aws_key = input('AWS Access Key: '+ '[****************' + tempDat + ']: ')
			if aws_key == '':
				aws_key = priorConfig['awsCreds']['access_key']
			else:
				pass

		tempDat = priorConfig['awsCreds']['aws_secret_key'] 
		if tempDat == '':
			aws_secret_key = input('AWS Secret Access Key: ')
		else:
			tempDat = tempDat[-4:]
			aws_secret_key = input('AWS Secret Access Key: ' + '[****************' + tempDat + ']: ')
			if aws_secret_key == '':
				aws_secret_key = priorConfig['awsCreds']['aws_secret_key']
			else:
				pass
		
		tempDat = priorConfig['awsCreds']['region'] 
		if tempDat == '':
			region = input('AWS region: ')
		else:
			region = input('AWS region [' + tempDat + ']: ')
			if region == '':
				region = priorConfig['awsCreds']['region'] 
			else:
				pass


		S3bucket = input('AWS S3 Bucket Output: ')
		if S3bucket == '':
			S3bucket = priorConfig['S3bucket']
		else:
			pass


		secGroup = input('Input security group: ')
		if secGroup == '':
			secGroup = priorConfig['awsEC2']['secGroup']
		else:
			pass

		pemAddress = input('Input Private Security Key (.pem) address: ')
		if pemAddress == '':
			pemAddress = priorConfig['awsEC2']['secGroup']
		else:
			pass

		osInput = input('Linux, Mac or Windows Operating System (accpetable answers are "linux", "mac", "unix", or "windows"): ' )
		if osInput == '':
			osInput = priorConfig['OS']
		else:
			if osInput.lower() in ('linux','mac', 'unix'):
				osType = 'terraform'
			else:
				osType = 'Terraform'

		userInput = input('Input username to associate with your Immunospace activity: ' )

	#Empty data structure for credentials.  
	#jsonInput variables will be used for writing to the config.json file.
	jsonInput = {
			'awsCreds': {
						'access_key': '',
						'aws_secret_key': '',
						'region': '',
							},
			'S3bucket': '',
			'awsEC2':{
				'secGroup':'',
				'pemAddress' : ''
				},
			'OS':'',
			"User": ''
			}

	#Embedding command line input credentials with jsonInput data structure.
	jsonInput['awsCreds']['access_key'] = aws_key
	jsonInput['awsCreds']['aws_secret_key'] = aws_secret_key
	jsonInput['awsCreds']['region'] = region 
	jsonInput['S3bucket'] = S3bucket
	jsonInput['awsEC2']['secGroup'] = secGroup
	jsonInput['awsEC2']['pemAddress'] = pemAddress
	jsonInput['OS'] = osType
	jsonInput['User'] = userInput

	#Implement try-catch block for checking unicode.
	#This ensures writing to a file is in in unicode :
	try:
		to_unicode = unicode
	except NameError:
		to_unicode = str
	#Write config.json file in unicode standard utf-8. 
	with io.open('config.json', 'w', encoding='utf8') as outfile:
		strOut = json.dumps(jsonInput,
						indent=4, sort_keys=True,
						separators=(',', ': '), ensure_ascii=False)
		outfile.write(to_unicode(strOut))
	return()
