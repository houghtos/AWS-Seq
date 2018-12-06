#!/usr/bin/python
#Command line prompt (similar to AWS/Google CLI) for writing and storing AWS credentials in config.json
#coding utf-8
#Open and write to config.json. 

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

		tempDat = priorConfig['awsCreds']['aws_secret_key'] 
		if tempDat == '':
			aws_secret_key = input('AWS Secret Access Key: ')
		else:
			tempDat = tempDat[-4:]
			aws_secret_key = input('AWS Secret Access Key: ' + '[****************' + tempDat + ']: ')
		
		tempDat = priorConfig['awsCreds']['region'] 
		if tempDat == '':
			region = input('AWS region: ')
		else:
			region = input('AWS region [' + tempDat + ']: ')


		tempDat = priorConfig['S3bucket']
		S3bucket = input('AWS S3 Bucket Output: ')


		secGroup = input('Input AWS security group: ')

		pemAddress = input('Input Private Security Key (.pem) directory address: ')

		osInput = input('Linux, Mac or Windows Operating System (accpetable answers are "linux", "mac", "unix", or "windows"): ' )
		if osInput.lower() in ('linux','mac', 'unix'):
			osType = './terraform'
		else:
			osType = 'Terraform'

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
			'OS':''
			}

	#Embedding command line input credentials with jsonInput data structure.
	jsonInput['awsCreds']['access_key'] = aws_key
	jsonInput['awsCreds']['aws_secret_key'] = aws_secret_key
	jsonInput['awsCreds']['region'] = region 
	jsonInput['S3bucket'] = S3bucket
	jsonInput['awsEC2']['secGroup'] = secGroup
	jsonInput['awsEC2']['pemAddress'] = pemAddress
	jsonInput['OS'] = osType

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
