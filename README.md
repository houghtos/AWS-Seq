# AWS-Seq
###### Using Terraform, AWS CLI, and boto3 (python) to automate EC2 spot instance creation and task managing.

###### Heavily recommend launching EC2 instances from same region as S3 buckets.  Copying between an EC2 spot instance and S3 from the same region does not incur download charges.

## Dependencies:

  - Terraform
  - Python 3
  - boto3 python library
  - AWS CLI
  - SSM Agent (See AWS docs)


## Functionality:

  1. List historic AWS prices.
  2. Launch custom AWS spot EC2 instance types and pass bash commands to process paired RNA fastq files.
  3. Get status and attributes of commands invoked.
  4. Get instance IDs of all active instances.
  5. Deostray all terraform infastructure.

## Setting up:

  ###### Not Required: 
  
  1. **SpotPriceModule.py**: Filters for what EC2 pricing can be manually modified in line 33
  
  2. **Immunospace.py**: Fill in lines 89-91 with AWS S3 information to output SSM job status (standard input/output)
  
  ###### Required:
  
  1.  **MainTFSource.py**: Set lines:
                      
                      - 25, 81   -- with security key (.pem) name and file address (see examples) 
                      - 32       -- with AWS region
                      - 63       -- with the name of your AWS security group (requires SSH inbound)
  
  2.  **terraform.tvars**: Set your AWS key/secret key
  
  3.  **variables.tf**: Set your AWS region
  
  4. **RNAseqSSM.py** Set lines:
                      
                      - 14       -- with your AWS key
                      - 15       -- with your secret access key
                      - 16       -- with your S3 region
                      - 17       -- with your AWS S3 list of relevant hg19 reference files (written as "s3://YourBucket/hg19Refs/")
                      - 67       -- Output S3 for all files processed (written as "s3://yourbucket/yourprefix/")
