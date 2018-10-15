#!/bin/bash -v
sudo yum update -y
sudo yum install -y python-pip        
sudo yum install java-1.8.0-openjdk-devel -y
sudo yum install docker -y
sudo sudo service docker start
sudo docker pull houghtos/immunotools:version2
sudo docker run houghtos/immunotools:version2
sudo systemctl enable amazon-ssm-agent
sudo systemctl start amazon-ssm-agent
