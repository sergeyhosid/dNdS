##################################################
# calculation dN sites of codon matrix
##################################################
{
	if(FNR==1)
		num_file++;
	if(num_file==1)
	{
		if(substr($1,4,1)=="/")
			for(i=2;i<=NF;i++)
			{
				arr_ami[toupper($i)]=substr($1,5,1);
#				print toupper($i), arr_ami[toupper($i)];
			}
	}
	if(num_file==2)
	{
		arr_num_nuc[$1]=$2;
		arr_num_comb[$1]=$3;
		num_comb=1;
		pos=1;
		for(i=5;i<=NF;i++)
		{
			if($i==":")
			{
				num_comb++;
				pos=0;
			}
			else
				arr_comb[$1,num_comb,pos]=$i;
			pos++;
		}
	}
}
END{
	arr_nuc[0]="A";
	arr_nuc[1]="T";
	arr_nuc[2]="G";
	arr_nuc[3]="C";

#	generration codon
	for(i=0;i<=3;i++)
		for(j=0;j<=3;j++)
			for(k=0;k<=3;k++)
			{
				pos=i*4*4+j*4+k;
				arr_cod[pos]=arr_nuc[i] arr_nuc[j] arr_nuc[k];
				arr_cod_ind[arr_cod[pos],1]=pos;
				arr_cod_ind[arr_cod[pos],2]=i;
				arr_cod_ind[arr_cod[pos],3]=j;
				arr_cod_ind[arr_cod[pos],4]=k;
#				if(arr_ami[arr_cod[pos]])
#					print pos, arr_cod[pos], arr_cod_ind[arr_cod[pos]], arr_ami[arr_cod[pos]];
			}
			
# calculation dN dS sites of each codon
	for(i=0;i<=63;i++)
	{
# dN and dS sites of codon
		arr_cod_ind[arr_cod[i],5]=arr_cod_ind[arr_cod[i],6]=0;
		for(j=1;j<=3;j++)
		{
			for(k=0;k<=3;k++)
				if(k!=arr_cod_ind[arr_cod[i],j+1])
				{
					if(j==1)
						cod_new=arr_nuc[k] arr_nuc[arr_cod_ind[arr_cod[i],3]] arr_nuc[arr_cod_ind[arr_cod[i],4]]; 
					else 	if(j==2)
						cod_new=arr_nuc[arr_cod_ind[arr_cod[i],2]] arr_nuc[k] arr_nuc[arr_cod_ind[arr_cod[i],4]]; 
					else 	if(j==3)
						cod_new=arr_nuc[arr_cod_ind[arr_cod[i],2]] arr_nuc[arr_cod_ind[arr_cod[i],3]] arr_nuc[k]; 

					arr_cod_ind[arr_cod[i],5]+=(arr_ami[arr_cod[i]]!=arr_ami[cod_new]) ? 1 : 0;
					arr_cod_ind[arr_cod[i],6]+=(arr_ami[arr_cod[i]]==arr_ami[cod_new]) ? 1 : 0;
#					print i, arr_cod[i], cod_new, arr_ami[arr_cod[i]], arr_ami[cod_new], arr_cod_ind[arr_cod[i],5], arr_cod_ind[arr_cod[i],6] > "test_cod_log";
				}
		}
		print i, arr_cod[i], arr_ami[arr_cod[i]], arr_cod_ind[arr_cod[i],5], arr_cod_ind[arr_cod[i],6] > "test_cod_log";
	}
	
	calculation codone matrix
	for(i=0;i<=63;i++)
	{
		for(j=0;j<=63;j++) # 63
			if(arr_ami[arr_cod[i]]&&arr_ami[arr_cod[j]])
			{
				pos=0;
				kof=1;
				for(k=3;k>=1;k--)
				{
					pos+=(substr(arr_cod[i],k,1)!=substr(arr_cod[j],k,1)) ? kof : 0;
					kof*=4;
				}
				print "start", i,j, arr_cod[i], arr_cod[j], arr_ami[arr_cod[i]], arr_ami[arr_cod[j]], pos, arr_num_nuc[pos], arr_num_comb[pos] > "test_cod_count";
				sum_S=sum_N=sum_sS=sum_sN=num=0;
				
# change combination [1,2,6]				
				for(k1=1;k1<=arr_num_comb[pos];k1++)
				{
					sum_sS+=arr_cod_ind[arr_cod[i],5];
					sum_sN+=arr_cod_ind[arr_cod[i],6];
					num++;
# test combination [1,2,3 nucleotides]
					arr_tmp[1]=arr_nuc[arr_cod_ind[arr_cod[i],2]];
					arr_tmp[2]=arr_nuc[arr_cod_ind[arr_cod[i],3]];
					arr_tmp[3]=arr_nuc[arr_cod_ind[arr_cod[i],4]];
					cod_tmp=arr_tmp[1] arr_tmp[2] arr_tmp[3];
#					print "#", i,j, arr_cod[i], arr_cod_ind[arr_cod[i],1],arr_cod_ind[arr_cod[i],2],arr_cod_ind[arr_cod[i],3],arr_cod_ind[arr_cod[i],4], cod_tmp,"#";
					for(k2=1;k2<=arr_num_nuc[pos];k2++)
					{
							printf("%d ",arr_comb[pos,k1,k2]) > "test_cod_count";
							arr_tmp[arr_comb[pos,k1,k2]]=arr_nuc[arr_cod_ind[arr_cod[j],arr_comb[pos,k1,k2]+1]];
							cod_tmp_new=arr_tmp[1] arr_tmp[2] arr_tmp[3];
							printf("%s %s %d %d; ", cod_tmp, cod_tmp_new, arr_cod_ind[cod_tmp,5], arr_cod_ind[cod_tmp,6]) > "test_cod_count";
							sum_sS+=(k2<arr_num_nuc[pos]) ? arr_cod_ind[cod_tmp,5] : 0;
							sum_sN+=(k2<arr_num_nuc[pos]) ? arr_cod_ind[cod_tmp,6] : 0;
							num+=(k2<arr_num_nuc[pos]) ? 1 : 0;
							cod_tmp=cod_tmp_new;
					}
					printf("\n") > "test_cod_count";
				}
				printf("sum: %d %d %d\n",sum_sS,sum_sN,num) > "test_cod_count";
				print i,j, arr_cod[i], arr_cod[j], sum_sS, sum_sN, num, arr_ami[arr_cod[i]], arr_ami[arr_cod[j]], pos, arr_num_nuc[pos], arr_num_comb[pos];

			}
#		if(i==5)	
#			exit;
	}
}
