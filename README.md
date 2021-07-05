
Note: This is no longer maintained and likely incompatible with current versions of Terraform and possibly AWS SSM. 


![Immunospace](https://github.com/houghtos/Immunospace/blob/master/Immunospace.jpg)

# Immunospace 1.2
###### Cost effective processing bulk-RNAseq data through command line.  Can be modified for additional analysis types (e.g. ATAC/ChIP-seq).   Uses Terraform and AWS SDK/Tool dependencies described below. Monitor instance progress regularly & check outputs.  Do not run workflows longer than 24 hours.      

###### Please ensure EC2 instances are from same region as S3 buckets.  Copying between an EC2 instance and S3 from the same region does not incur download charges.  

######  Downloading large files from S3 can be costly. Designed for uploading larger files (e.g. fastq or bam) to S3 and downloading resulting files (e.g. .tsv or .txt). 

##### Table of Contents  
1. [Dependencies](#Dependencies)
2. [Functionality](#Functionality)
3. [Setting up](#Setting-up)
4. [Usage](#Usage)
5. [Patch Notes](#Patch-Notes)
6. [Pipelines](#Pipelines)
7. [Known Issues](#Known-Issues)

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
  6. Destroy all terraform infastructure.  Note: please remember to destroy all terraform infastrucutre after a run completes or fails.  Even if a run is complete, there will still be some leftover infastructure that should be destroyed before the next run as it can induce errors.

## Setting up:

  1. Install [Terraform](https://www.terraform.io/intro/getting-started/install.html), [Python 3](https://www.python.org/downloads/) and [pip](https://pip.pypa.io/en/stable/installing/), , the and [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/installing.html).  
  2. Configure local machine's AWS CLI with access key, secret key, and region.  The
  3. Download all source files to the terraform installed directory.  
  4. Run the "Immunospace.py configure" command (see usage below).  This will set permissions each EC2 spot instance will inherit.  Please configure with AWS keys from an IAM user which ONLY has S3 read/write permission to where the input data (e.g. s3://bucket/AWS.pipeline.results/file_R1.fastq.gz) is stored. 
  5. Download the [full hg19 (human), mm10 (mouse)](http://hgdownload.cse.ucsc.edu/downloads.html), or other species reference genome data sets.  Upload this to an S3 bucket/prefix readable by the Immunospace configured IAM user from prior step (4).

  
  ###### Required Code Changes: 
  - **RNAseqSSM.py** Set lines:
                  
    - 25       -- With your S3 address for reference genome files.  In the RNAseqSSM.py example, we use mm10.
    - 73       -- Set the .gtf reference file used in htseq analysis.  In the RNAseqSSM.py example, we use Mus_musculus.GRCm38.75.gtf.  This would be one of the files copied from your reference genome files.
     
  
  
  ###### Optional Code Settings: 
  
  - **SpotPriceModule.py**: Filters for what EC2 pricing can be manually modified in line 33
  
  
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
      
    4. $ python Immunospace.py command
      -No options.  Looks for all recent commands invoked.
    
    5. $ python Immunospace.py instances 
      - No options.  Displays all active instances launched by terraform as list.
    
    6. $ python Immunospace.py destroy
     - No options.  Destroys all Terraform infastructure that currently exists.
     
## csv file input:

    See example in AWS-Seq repository: "rna_AWS_files.csv" 
    
    The number of rows in the CSV should be the same as number of instances you are starting.  Each row of paired fastq files will be processed by one instance via sending SSM command.
    
## Patch Notes

##### 1.2 Update Notes:
  - Patched Immunospace.py write main terraform function to include writing terraform.tfvars file based on config.json setup. 
  - Updated SSM pipelines to auto-copy to bucket designated in config.json
  - Now using Docker container Immunotools to version4 which takes up ~270MB of more space (~750MB total).   
    - Immunotoools version4 docker container now contains [MiXCR](https://mixcr.readthedocs.io/en/master/) and bedtools  
  - Added MiXCR pipeline that is swappable with the RNAseq default.  Same inputs as RNAseq pipeline. 
  - Immunospace configure now prompts for a "username". This is can be an arbitrary selected name that will make the EC2 IAM and key permission names unique reducing errors from multiple users from the same account using the software at once. 
  
##### 1.1 Update Notes:
  - Added configure functionality so user has to make far fewer modifications to code.
  - Set SSM timeout to 10 hours per instance (rather than 1 hour default).  This can be adjusted in the immunospace.py file under the boto3 SSM function.
  - Minor changes to RNAseq pipeline to run smoother.
  - Added option '-c' to "$ python immunospace.py start ..." for cloud prefix output of files and SSM run logs.

## Pipelines:

1. HLA Typing: Input SINGLE bam files in the csv.
2. RNAseq: Input paired FASTQ files (see example input)
3. MIXCR: Input paired FASTQ files (see example input)
4. MIXCR + RNAseq: Input paired FASTQ files.
 

## Known Issues:
  - "aws_iam_role.ec2_access: Error creating IAM Role xxxxx: EntityAlreadyExists: Role with name xxxx already exists."  
  
  Error induced from multiple users attempting to use same IAM role name at the same time.  Also caused by incomplete deletion of Terraform infastructure after use. 
  
  *Solutions:* delete all infastructure for all users and re-run.  This bug was addressed in version 1.2 and should be less prevelant.
  
  - "aws_spot_instance_request.test_spot.1: Error requesting spot instances: InvalidSubnetID.NotFound: No default subnet for availability zone: 'null'."  
  
  Error is random and appears to be induced by increased AWS usage (and decreased capacitiy for spot usage) in given region.    
  
   *Solutions:* Ensure you have a default VPC and subnets available for the region you are running your infastructure.  Check to ensure there is capacity to launch spot instances in your specified region's availability areas.  A simple work around to ensure the availability zone is not overcapacity is to add "availability_zone" to tfSource.py resource "aws_spot_instance_request" per below example:
   
    resource "aws_spot_instance_request" "test_spot" {
    count         = n
    ami           = ami_xxxx
    availability_zone = "us-east-1f"
    }

