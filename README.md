# AWS-Seq
###### Uses Terraform, AWS CLI, and boto3 (python) to automate EC2 spot instance creation and task managing.  Current setup is for processing paired fastq RNAseq files for gene counts. 

###### Heavily recommend launching EC2 instances from same region as S3 buckets.  Copying between an EC2 spot instance and S3 from the same region does not incur download charges.

###### Credit: David Redmond helping develop RNAseq fastq paired pipeline. 

## Dependencies:

  - [Terraform](https://www.terraform.io/intro/getting-started/install.html)
  - [Python 3](https://www.python.org/downloads/) and [pip](https://pip.pypa.io/en/stable/installing/)
  - [boto3 python library](https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html#installation)
  - [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/installing.html)
  - [SSM Agent](https://docs.aws.amazon.com/systems-manager/latest/userguide/ssm-agent.html)


## Functionality:

  1. List historic and current AWS spot pricing.  
  2. Configure credentials for use in provisioning instance(s), sending SSM commands, etc.
  3. Launch custom AWS spot EC2 instance types and pass bash commands to spot EC2 instance(s) to process paired RNA fastq files.
  4. Get status and attributes of commands being invoked.
  5. Get instance IDs of all active instances.
  6. Destroy all terraform infastructure.

## Setting up:

  ###### Not Required: 
  
  1. **SpotPriceModule.py**: Filters for what EC2 pricing can be manually modified in line 33
  
  2. **RNAseqSSM.py** Set lines:
                  
    - 16       -- with your AWS key
    - 17       -- with your secret access key
    - 18       -- with your S3 region
    - 19       -- with your AWS S3 list of relevant hg19 reference files (written as "s3://YourBucket/hg19Refs/")
    - 66       -- Output S3 for all files processed (written as "s3://yourbucket/yourprefix/")

## Usage:

All commands must be made from your terraform directory.

    1. $ python Immunospace.py price -t <instance type> -m <number months prior to query>
      - Instance type (e.g. t3.2xlarge)
      - Months prior to today to query pricing
     
    2. $ python Immunospace.py 
    
    3. $ python Immunospace.py start -n <number of instances> -t <instance type> -b <maximum spot bid> -s <number GiB for root storage> -f <S3 file list> -c <S3 cloud prefix>
      - Number of instance to launch (integer)
      - Instance type (e.g. t3.2xlarge)
      - Maximum spot bid in dollars (e.g. 0.25)
      - GiB assigned to root storage (integer).  Must account for size of all files or run will fail.
      - S3 file list is a CSV of S3 file addresses.  Each row contains fastq pairs.  The first row contains pairs to be submitted to first EC2 instance.  The row column pairs for second EC2 instance. etc. See next section for more information.
      -  S3 prefix (assuming bucket is already configured) to write files and SSM logs to.
      
    4. $ python Immunospace command
      -No options.  Looks for all recent commands invoked.
    
    5. $ python Immunospace instances 
      - No options.  Displays all active instances launched by terraform as list.
    
    6. $ python Immunospace destroy
     - No options.  Destroys all Terraform infastructure that currently exists.
     
## csv file input:

    See example in AWS-Seq repository: "rna_AWS_files.csv" 
    
    Required input for python "Immunospace start ..."  CSV formatted file of S3 addresses for all fastq.gz pairs.  Each pair should be next to eachother in order (if possible).  The first row of pairs will be submitted to first EC2 instance.  The second row of pairs submitted to the second EC2 instance. etc.  Ensure the number of columns and instances aligns to evenly spread processing or your needs. 
 

