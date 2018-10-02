### Edit filters on line 33
#Requires AWS credentials to be embedded

#datetime is the only dependency in this module
from datetime import datetime as dt
import datetime
import os,sys

class comparePrice:

	def __init__(self, ami, delta):
		self.ami = ami
		self.delta = abs(delta) * -1		
		self.date = datetime.datetime.now()
		self.dateDelta = '' 

		#Inputs current datetime (date) and how many months you wish to change it by (delta)
	def monthDelta(self):
		m, y = (self.date.month+self.delta) % 12, self.date.year + ((self.date.month)+self.delta-1) // 12
		if not m: m = 12
		d = min(self.date.day, [31,
			29 if y%4==0 and not y%400==0 else 28,31,30,31,30,31,31,30,31,30,31][m-1])
		#return date.replace(day=d,month=m, year=y)
		self.dateDelta = str(self.date.replace(day=d,month=m, year=y))
		self.dateDelta = self.dateDelta.replace(" ", "T")

		self.date = str(self.date)
		self.date = self.date.replace(" ", "T")
		return (self.dateDelta, self.date)



	def describeSpotPrice(self):
		os.system('aws ec2 describe-spot-price-history --instance-types {} --product-description "Linux/UNIX (Amazon VPC)" --start-time {} --end-time {}'.format(self.ami,self.dateDelta,self.date))

