#!/usr/bin/python

#Set lines 26, 33, 63, 80

def mainGen(instanceNumbers, AMI_type, spot_bid, volume_size):
  main_string = """

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
    default = "MyEC2Key"   ############## <- set the name of your key.  See line 80
  }

  provider "aws" {
    ## UNABLE TO ACCEPT KEYS WITH CHARACTERS LIKE '/+'  in most recent version.  Need fix.  Try config.
    access_key = "${var.access_key}"
    secret_key = "${var.secret_key}"
    region = "us-wes-1"  ############ <-  Set region
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

  resource "aws_spot_instance_request" "spot_seq" {
    count         =  %s
    ami           = "${data.aws_ami.amazon_linux.id}"
    spot_price    = %s 
    key_name =    "${var.aws_keypair}"
    instance_type = "%s"
    spot_type     = "one-time"
    iam_instance_profile   = "${aws_iam_instance_profile.SSM_Profile.name}"
    associate_public_ip_address = true
    wait_for_fulfillment = true
    
    key_name = "MyEC2Key"
    vpc_security_group_ids = ["sg-9a3c57d0"]  # Amazon Linux 2 AMI (HVM), SSD Volume Type 


    root_block_device {
      volume_size = "%s"
      }

    tags {
      Name = "Spot Seq"
    }
  }

  resource "null_resource" "configure-aws-instances" {
    count =  %s

    connection {
      type        = "ssh"
      user        = "ec2-user"
      private_key = "${file("/home/users/Sean/AWS/EC2_Private_Keys/MyEC2Key.pem")}" ############ <-- Specify private key local path
      host  = "${element(aws_spot_instance_request.spot_seq.*.public_ip, count.index)}"
    }

    provisioner "remote-exec"{
      inline = [
        "sudo yum update -y",
        "sudo yum install -y python-pip",        
        "sudo yum install java-1.8.0-openjdk-devel -y",
        "sudo yum install docker -y",
        "sudo sudo service docker start",
        "sudo docker pull houghtos/immunotools:version2",
        "sudo docker run houghtos/immunotools:version2",
        "sudo systemctl enable amazon-ssm-agent",
        "sudo systemctl start amazon-ssm-agent"
      ]
    }
  }

  """ %  (instanceNumbers, spot_bid, AMI_type, volume_size, instanceNumbers)

  return(main_string)
