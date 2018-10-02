# AWS-Seq
Using Terraform, AWS CLI, and boto3 (python) to automate EC2 spot instance creation and task managing.

Dependencies:

  - Terraform
  - Python 3
  - boto3 python library
  - AWS CLI
  - SSM Agent (See AWS docs)


Functionality:

  1. List historic AWS prices.
  2. Launch custom AWS spot EC2 instance types and pass bash commands to process paired RNA fastq files.
  3. Get status and attributes of commands invoked.
  4. Get instance IDs of all active instances.
  5. Deostray all terraform infastructure.
