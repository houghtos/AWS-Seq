#!/usr/bin/python
#Set lines 23, 30, 65

def mainGen(instanceNumbers, AMI_type, spot_bid, volume_size):
  string_1 = """
  resource "aws_iam_role" "ec2_access" {
    name               = "ec2_spot"
    assume_role_policy = "${file("assume-role-policy.json")}"
  }


  resource "aws_iam_role_policy_attachment" "SSM_Profile" {
    role = "${aws_iam_role.ec2_access.name}"
    policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2RoleforSSM"
  }

  resource "aws_iam_instance_profile" "SSM_Profile" {
    name  = "SSM_Profile"
    role = "${aws_iam_role.ec2_access.name}"
  }

  variable "aws_keypair" {
    default = "MyEC2Key"   ############## <- set the name of your key. 
  }

  provider "aws" {
    ## UNABLE TO ACCEPT KEYS WITH CHARACTERS LIKE '/+'  in most recent version.  Need fix.  Try config.
    access_key = "${var.access_key}"
    secret_key = "${var.secret_key}"
    region = "us-west-1"  ############ <-  Set region
  }

  data "aws_ami" "amazon_linux" {
    most_recent = true

    filter {
      name   = "name"
      values = ["amzn2-ami-hvm-*"]
    }

    filter {
      name   = "virtualization-type"
      values = ["hvm"]
    }

    owners = ["137112412989"] # Canonical
  }

"""

  string_2 =  ' resource "aws_spot_instance_request" "test_spot" {\n'
  string_3 =  ' count         =  {}\n'.format(instanceNumbers)
  string_4 =  ' ami           = "${data.aws_ami.amazon_linux.id}"\n'
  string_5 =  ' spot_price    = {}\n'.format(spot_bid)
  string_6 =  ' key_name =    "${var.aws_keypair}"\n'
  string_7 =  ' instance_type = "{}"\n'.format(AMI_type)
  string_8 =  ' spot_type     = "one-time"\n'
  string_9 =  """  

  iam_instance_profile   = "${aws_iam_instance_profile.SSM_Profile.name}"\n
  associate_public_ip_address = true
  wait_for_fulfillment = true
  
  key_name = "ec2_test_group"
  vpc_security_group_ids = ["sg-1a2b3c4d"]  ###### <--- specify your security group.  Does not require SSH.
  user_data       = "${file("userdata.sh")}"
    """
  string_10  = 'root_block_device {volume_size = "%d"}' % (volume_size)
  string_11  = 'tags {Name = "Immuno Spot"}}'
      

  return(string_1 + string_2 + string_3 + string_4 + string_5 + string_6 + string_7 + string_8 + string_9 + string_10 + string_11)

