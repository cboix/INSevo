#!/bin/bash
awk 'BEGIN{split("",a); a[1]="2L";a[2]="2R";a[3]="3L";a[4]="3R";a[5]="4";a[6]="X";a[7]="Y";}
	{b="";
		for (i=1;i<=7;i++)
			{
			if($2 !~ a[i]) b=sprintf("%s--not-chr %s ",b,a[i])
			};
		print "vcftools --vcf /net/home/carlesba/db/DGRP/f.vcf",b,"--chr",$2,"--from-bp",$5,"--to-bp",$6,"--out",$1,"--recode";}' dmel_targets > batch
